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


def simulate_dynamic(
    mesh_state,
    mass: float,
    force: float,
    time_span: tuple[float, float],
    initial_state: tuple[float, float] = (0.0, 0.0),
) -> DynamicResponse:
    """Solve a simple 1-DOF dynamic system using precomputed stiffness and damping."""

    if isinstance(mesh_state, dict):
        k_te = np.asarray(mesh_state["k_te"], dtype=float)
        c = np.asarray(mesh_state.get("c", np.zeros_like(k_te)), dtype=float)
        forcing = mesh_state.get("forcing", None)
    else:
        k_te = np.asarray(mesh_state.k_te, dtype=float)
        c = np.asarray(mesh_state.c, dtype=float)
        forcing = getattr(mesh_state, "forcing", None)

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
        return [v, -k * x / mass - damping * v / mass + excitation / mass]

    time_grid = np.linspace(time_span[0], time_span[1], len(k_te) * 10)
    sol = solve_ivp(rhs, time_span, [initial_state[0], initial_state[1]], t_eval=time_grid, rtol=1e-6, atol=1e-9)
    x = sol.y[0]
    v = sol.y[1]
    a = np.gradient(v, time_grid, edge_order=2)

    return DynamicResponse(
        time=time_grid,
        displacement=x,
        velocity=v,
        acceleration=a,
    )
