from __future__ import annotations

import math
from typing import Dict

import numpy as np


def precompute_friction_state(vrel: np.ndarray, force_n: float, film_stiffness: float = 1.0e8) -> Dict[str, np.ndarray]:
    """A simple EHL-inspired regularized friction model that avoids singularities at vrel=0."""
    vrel_safe = np.maximum(vrel, 1e-9)
    mu = np.clip(0.02 + 0.003 / np.sqrt(vrel_safe), 0.001, 0.2)
    film_deflection = force_n / film_stiffness
    damping = mu * force_n / np.maximum(vrel_safe, 1e-6)
    return {
        "mu": mu,
        "film_deflection": film_deflection,
        "damping": damping,
    }
