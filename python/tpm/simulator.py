from __future__ import annotations

from typing import Dict

from tpm.dynamics import simulate_dynamic
from tpm.geometry import compute_geometry
from tpm.mass import compute_mass_properties
from tpm.stiffness import precompute_mesh_state


def run_simulation(
    z1: int = 27,
    z2: int = 35,
    module_mm: float = 3.0,
    alpha_deg: float = 20.0,
    width_mm: float = 25.0,
    density: float = 7850.0,
    force_n: float = 8000.0,
    deltamax_um: float = 40.0,
) -> Dict[str, object]:
    geom = compute_geometry(z1=z1, z2=z2, module=module_mm, alpha_deg=alpha_deg)
    mass_props = compute_mass_properties(geom=geom, b_mm=width_mm, density=density, deltamax_um=deltamax_um)
    mesh_state = precompute_mesh_state(
        geom=geom,
        z1=z1,
        z2=z2,
        step=400,
        w=200.0,
        force=force_n,
        elasticity_modulus=206e9,
        poisson_ratio=0.3,
        deltamax_um=deltamax_um,
    )
    response = simulate_dynamic(mesh_state=mesh_state, mass=mass_props.me, force=force_n, time_span=(0.0, 0.01))

    return {
        "geometry": geom,
        "mass_properties": mass_props,
        "mesh_state": mesh_state,
        "response": response,
    }
