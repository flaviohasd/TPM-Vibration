from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict

import numpy as np

from .geometry import compute_geometry as _compute_geometry
from .geometry import compute_mesh as _compute_mesh
from .mass import compute_mass_properties


@dataclass
class GearSystem:
    z1: int
    z2: int
    module: float
    alpha_deg: float
    b_mm: float
    elasticity_modulus: float
    poisson_ratio: float
    density: float
    power_w: float
    rpm: float
    lubricant_density: float
    temperature_c: float
    deltamax_um: float = 0.0

    def __post_init__(self) -> None:
        self.geometry = self._compute_geometry()
        self.mesh = self._compute_mesh()
        self.mass_properties = self._compute_mass_properties()
        self.operation = self._compute_operation()
        self.tvms = self._compute_tvms()

    def _compute_geometry(self) -> Dict[str, float]:
        geometry = _compute_geometry(self.z1, self.z2, self.module, self.alpha_deg)
        geometry["alpha"] = float(math.radians(self.alpha_deg))
        return geometry

    def _compute_mesh(self) -> Dict[str, float]:
        return _compute_mesh(self.geometry)

    def _compute_mass_properties(self) -> Dict[str, float]:
        mass_properties = compute_mass_properties(
            geom=self.geometry,
            b_mm=self.b_mm,
            density=self.density,
            deltamax_um=self.deltamax_um,
        )
        return {
            "me": mass_properties.me,
            "memod": mass_properties.memod,
            "mt1": mass_properties.mt1,
            "mt2": mass_properties.mt2,
            "J1": mass_properties.J1,
            "J2": mass_properties.J2,
            "mt1mod": mass_properties.mt1mod,
            "mt2mod": mass_properties.mt2mod,
            "J1mod": mass_properties.J1mod,
            "J2mod": mass_properties.J2mod,
        }

    def _compute_operation(self) -> Dict[str, float]:
        w = self.rpm * 2.0 * math.pi / 60.0
        db = self.geometry["db1"] / 1000.0
        dp1 = self.geometry["dp1"] / 1000.0
        dp2 = self.geometry["dp2"] / 1000.0
        torque = self.power_w * 60.0 / (2.0 * math.pi * self.rpm)
        force = torque / (db / 2.0)
        freq = w / (2.0 * math.pi)
        cycle = 2.0 * math.pi / w
        return {
            "w": w,
            "w1": w,
            "w2": -w * dp1 / dp2,
            "T": torque,
            "F": force,
            "freq": freq,
            "cycle": cycle,
        }

    def _compute_tvms(self, step: int = 400) -> Dict[str, np.ndarray]:
        from .energy import precompute_energy_state

        return precompute_energy_state(self, step=step)
