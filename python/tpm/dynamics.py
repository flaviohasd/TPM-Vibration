from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

import numpy as np
from scipy.integrate import solve_ivp


@dataclass
class DynamicResponse:
    time: np.ndarray
    displacement: np.ndarray
    velocity: np.ndarray
    acceleration: np.ndarray
    dmf: np.ndarray | None = None


def simulate_dynamic(
    mesh_state,
    mass: float,
    force: float,
    time_span: tuple[float, float],
    initial_state: tuple[float, float] = (0.0, 0.0),
    num_points: int | None = None,
) -> DynamicResponse:
    """Solve a 1-DOF gear dynamic system using precomputed TVMS, damping, and forcing."""

    if isinstance(mesh_state, dict):
        k_te = np.asarray(mesh_state.get("k_te_base", mesh_state.get("k_te")), dtype=float)
        c = np.asarray(mesh_state.get("c", np.zeros_like(k_te)), dtype=float)
        forcing = mesh_state.get("forcing", None)
        tpm_relief = mesh_state.get("tpm_relief", None)
    else:
        k_te = getattr(mesh_state, "k_te_base", getattr(mesh_state, "k_te", None))
        k_te = np.asarray(k_te, dtype=float)
        c = np.asarray(getattr(mesh_state, "c", np.zeros_like(k_te)), dtype=float)
        forcing = getattr(mesh_state, "forcing", None)
        tpm_relief = getattr(mesh_state, "tpm_relief", None)

    if forcing is None:
        forcing = np.full(len(k_te), force, dtype=float)
    else:
        forcing = np.asarray(forcing, dtype=float)
        if forcing.ndim == 0:
            forcing = np.full(len(k_te), float(forcing), dtype=float)

    time_index = np.linspace(time_span[0], time_span[1], len(k_te))

    def rhs(t, y):
        idx = int(np.interp(t, time_index, np.arange(len(k_te))))
        idx = min(max(idx, 0), len(k_te) - 1)
        k = k_te[idx]
        damping = c[idx]
        excitation = forcing[idx]
        x, v = y
        return [v, (excitation - k * x - damping * v) / mass]

    grid_pts = 400 if num_points is None else int(num_points)
    time_grid = np.linspace(time_span[0], time_span[1], grid_pts)
    sol = solve_ivp(rhs, time_span, [initial_state[0], initial_state[1]], t_eval=time_grid, rtol=1e-5, atol=1e-8)
    x = sol.y[0]
    v = sol.y[1]

    # Analytical acceleration from equations of motion: a = (F_ext - k*x - c*v) / m
    k_eval = np.interp(time_grid, time_index, k_te)
    c_eval = np.interp(time_grid, time_index, c)
    f_eval = np.interp(time_grid, time_index, forcing)
    a = (f_eval - k_eval * x - c_eval * v) / mass

    # Dynamic meshing force: DMF = k*x + c*v - TPM_relief_force
    if tpm_relief is None:
        tpm_eval = f_eval - force
    else:
        tpm_eval = np.interp(time_grid, time_index, np.asarray(tpm_relief, dtype=float))
    dmf = k_eval * x + c_eval * v - tpm_eval

    return DynamicResponse(
        time=time_grid,
        displacement=x,
        velocity=v,
        acceleration=a,
        dmf=dmf,
    )
