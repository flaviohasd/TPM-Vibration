"""TPM-Vibration: Tooth Profile Modification & Gear Mesh Vibration Optimization."""

from .faithful import GearSystem
from .optimization import FastDynamicOptimizer, optimize_deltamax
from .energy import precompute_energy_state
from .geometry import compute_geometry, compute_mesh
from .mass import compute_mass_properties
from .plots import generate_validation_plots

__all__ = [
    "GearSystem",
    "FastDynamicOptimizer",
    "optimize_deltamax",
    "precompute_energy_state",
    "compute_geometry",
    "compute_mesh",
    "compute_mass_properties",
    "generate_validation_plots",
]
