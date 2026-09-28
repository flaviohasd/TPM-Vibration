from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize_scalar

from .energy import precompute_energy_state
from .mass import compute_mass_properties


class FastDynamicOptimizer:
    """Precomputed base TVMS and PCHIP state for fast optimization of TPM relief amplitude."""

    def __init__(
        self,
        system,
        step: int = 1000,
        friction_model: str = "constant",
        fixed_mu: Optional[float] = 0.1211,
        cycles: int = 1,
        method: str = "LSODA",
    ) -> None:
        self.system = system
        self.step = step
        self.cycles = cycles
        self.method = method
        self.geom = system.geometry
        self.b_mm = system.b_mm
        self.density = system.density
        self.F = system.operation["F"]
        self.cycle_time = system.operation["cycle"] * cycles
        self.t_span = (0.0, self.cycle_time)

        # Precompute base energy state ONCE
        self.base_energy = precompute_energy_state(
            system, step=step, friction_model=friction_model, fixed_mu=fixed_mu
        )
        self.k_te = self.base_energy["k_te_base"]
        self.c = self.base_energy["c"]
        delta_norm = self.base_energy["delta_normalized"]
        k12 = self.base_energy["k12"]
        self.forcing_relief_rate = k12[:, 0] * delta_norm[:, 0] + k12[:, 1] * delta_norm[:, 1]

        # Automatic physical reference bound: Delta_ref = F / min(Kte) [um]
        min_kte = float(np.min(self.k_te))
        self.delta_ref_um = float((self.F / min_kte) * 1e6)

        # Build monotonic cubic Hermite (PCHIP) spline interpolators
        self.t_grid = np.linspace(0.0, self.cycle_time, step)
        self.dt = self.t_grid[1] - self.t_grid[0]
        self.n_intervals = step - 1

        self.pchip_k = PchipInterpolator(self.t_grid, self.k_te)
        self.pchip_c = PchipInterpolator(self.t_grid, self.c)
        self.pchip_rate = PchipInterpolator(self.t_grid, self.forcing_relief_rate)

        self.c_k = self.pchip_k.c
        self.c_c = self.pchip_c.c
        self.c_r = self.pchip_rate.c

        # Evaluation grid for steady-state RMS calculation (matching MATLAB signal.m)
        n_eval = max(1000, step)
        self.t_eval = np.linspace(0.0, self.cycle_time, n_eval)
        self.eval_half = n_eval // 2
        self.k_eval = self.pchip_k(self.t_eval)
        self.c_eval = self.pchip_c(self.t_eval)
        self.r_eval = self.pchip_rate(self.t_eval)

    def evaluate_rms_acceleration(self, deltamax_um: float) -> float:
        """Evaluate RMS acceleration for a candidate deltamax_um using fast ODE integration."""
        deltamax_um = float(deltamax_um)
        memod = compute_mass_properties(
            self.geom, self.b_mm, self.density, deltamax_um, correct_centroid_typo=False
        ).memod

        dt = self.dt
        n_intervals = self.n_intervals
        c_k = self.c_k
        c_c = self.c_c
        c_r = self.c_r
        F = self.F

        def rhs(t, y):
            i = int(t / dt)
            if i < 0:
                i = 0
            elif i >= n_intervals:
                i = n_intervals - 1
            dx = t - i * dt
            k_val = ((c_k[0, i] * dx + c_k[1, i]) * dx + c_k[2, i]) * dx + c_k[3, i]
            c_val = ((c_c[0, i] * dx + c_c[1, i]) * dx + c_c[2, i]) * dx + c_c[3, i]
            r_val = ((c_r[0, i] * dx + c_r[1, i]) * dx + c_r[2, i]) * dx + c_r[3, i]
            f_val = F + r_val * deltamax_um
            return [y[1], (f_val - k_val * y[0] - c_val * y[1]) / memod]

        sol = solve_ivp(
            rhs,
            self.t_span,
            [0.0, 0.0],
            method=self.method,
            t_eval=self.t_eval,
            rtol=1e-3,
            atol=1e-6,
        )
        x = sol.y[0]
        v = sol.y[1]

        # Compute dynamic acceleration a(t) = [F_dyn - K*x - C*v] / memod
        f_vals = F + self.r_eval * deltamax_um
        a = (f_vals - self.k_eval * x - self.c_eval * v) / memod

        # RMS on the second half of the signal (steady-state, identical to MATLAB signal.m)
        rms_a = float(np.sqrt(np.mean(a[self.eval_half:] ** 2)))
        return rms_a


def _objective(deltamax_um: float, system, step: int = 400) -> float:
    """Convenience objective function for backward compatibility and testing."""
    optimizer = FastDynamicOptimizer(system, step=step)
    return optimizer.evaluate_rms_acceleration(deltamax_um)


def optimize_deltamax(
    system,
    bounds: Optional[Tuple[float, float]] = None,
    step: int = 1000,
    opt_tol: float = 1e-4,
    n_points: int = 8,
    method: str = "LSODA",
) -> Dict[str, float]:
    """Optimize the TPM relief amplitude by minimizing steady-state RMS acceleration.

    Parameters
    ----------
    system : GearSystem
        Gear system with full geometry, material, and operational parameters.
    bounds : Tuple[float, float], optional
        Search interval (min_um, max_um). If None, automatically computed from
        static tooth deflection under load: [0.6 * Delta_ref, 1.4 * Delta_ref].
    step : int, default 1000
        Number of mesh discretization points for TVMS and spline interpolation.
    opt_tol : float, default 1e-4
        Convergence tolerance on deltamax_um (analogous to MATLAB TolX = 1e-4).
    n_points : int, default 8
        Retained for API compatibility.
    method : str, default 'LSODA'
        ODE integration method ('LSODA', 'BDF', 'Radau').
    """
    optimizer = FastDynamicOptimizer(system, step=step, method=method)

    if bounds is None:
        delta_ref = optimizer.delta_ref_um
        bounds = (0.6 * delta_ref, 1.4 * delta_ref)

    result = minimize_scalar(
        optimizer.evaluate_rms_acceleration,
        bounds=bounds,
        method="bounded",
        options={"xatol": opt_tol, "maxiter": 50},
    )

    return {
        "deltamax_um": float(result.x),
        "objective": float(result.fun),
        "success": bool(result.success),
        "delta_ref_um": optimizer.delta_ref_um,
        "bounds": bounds,
        "nfev": int(result.nfev),
    }
