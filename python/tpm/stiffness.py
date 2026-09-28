from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

import numpy as np


@dataclass
class MeshState:
    theta: np.ndarray
    k_te: np.ndarray
    k_tem: np.ndarray
    c: np.ndarray
    cm: np.ndarray
    vrel: np.ndarray
    transmission_error: np.ndarray


def precompute_mesh_state(
    geom: Dict[str, float],
    z1: int,
    z2: int,
    step: int,
    w: float,
    force: float,
    elasticity_modulus: float,
    poisson_ratio: float,
    deltamax_um: float,
) -> MeshState:
    """Create a simple periodic TVMS and damping profile for the initial Python port."""
    theta = np.linspace(0.0, 2.0 * np.pi, step)
    base_k = 2.0e9
    k_te = base_k * (1.0 + 0.15 * np.sin(theta))
    k_tem = k_te * (1.0 - 0.02 * deltamax_um / 100.0)
    c = 200.0 + 100.0 * np.sin(theta + 0.3)
    cm = c * 0.9
    vrel = 0.01 + 0.002 * np.sin(theta)
    transmission_error = 1e-6 * (1.0 + 0.2 * np.sin(theta))

    return MeshState(
        theta=theta,
        k_te=k_te,
        k_tem=k_tem,
        c=c,
        cm=cm,
        vrel=vrel,
        transmission_error=transmission_error,
    )
