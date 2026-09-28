from math import isfinite

import numpy as np
from pytest import approx

from tpm.faithful import GearSystem
from tpm.optimization import optimize_deltamax


def test_gear_system_matches_matlab_like_geometry_and_mesh():
    system = GearSystem(
        z1=27,
        z2=35,
        module=3.0,
        alpha_deg=20.0,
        b_mm=25.0,
        elasticity_modulus=206e9,
        poisson_ratio=0.3,
        density=7850.0,
        power_w=80e3,
        rpm=2000.0,
        lubricant_density=870.0,
        temperature_c=60.0,
    )

    # Geometry checks
    assert system.geometry["dp1"] == approx(81.0, rel=1e-6)
    assert system.geometry["dp2"] == approx(105.0, rel=1e-6)
    assert system.geometry["db1"] == approx(81.0 * np.cos(np.radians(20.0)), rel=1e-6)
    assert system.geometry["rb1"] == approx(38.05755, rel=1e-4)
    assert system.geometry["rb2"] == approx(49.33238, rel=1e-4)

    # Mesh relations
    assert system.mesh["cr"] == approx(1.65815, rel=1e-4)
    assert system.mesh["PTH"] == approx(14.6852, rel=1e-4)
    assert system.mesh["double"] == approx(5.8288, rel=1e-4)
    assert system.mesh["single"] == approx(3.0276, rel=1e-4)
    assert system.mesh["thetad"] == approx(0.15316, rel=1e-4)
    assert system.mesh["thetas"] == approx(0.07955, rel=1e-4)

    # Mass & TVMS
    assert isfinite(system.mass_properties["me"])
    assert system.mass_properties["me"] > 0.0
    assert system.tvms["k_te"].shape[0] > 0
    assert np.all(system.tvms["k_te"] > 0.0)


def test_tvms_energy_method_physics():
    system = GearSystem(
        z1=27,
        z2=35,
        module=3.0,
        alpha_deg=20.0,
        b_mm=25.0,
        elasticity_modulus=206e9,
        poisson_ratio=0.3,
        density=7850.0,
        power_w=80e3,
        rpm=2000.0,
        lubricant_density=870.0,
        temperature_c=60.0,
    )

    k_te = system.tvms["k_te"]
    # Mean mesh stiffness is around 3.5e8 N/m
    assert np.mean(k_te) == approx(3.54e8, rel=0.05)
    assert np.min(k_te) > 2.0e8
    assert np.max(k_te) < 5.0e8


def test_dynamic_optimization_finds_paper_optimal_deltamax():
    system = GearSystem(
        z1=27,
        z2=35,
        module=3.0,
        alpha_deg=20.0,
        b_mm=25.0,
        elasticity_modulus=206e9,
        poisson_ratio=0.3,
        density=7850.0,
        power_w=80e3,
        rpm=2000.0,
        lubricant_density=870.0,
        temperature_c=60.0,
    )

    result = optimize_deltamax(system, bounds=None, step=400, opt_tol=1e-4)

    assert result["success"]
    assert result["delta_ref_um"] == approx(43.88, rel=0.02)
    assert result["deltamax_um"] == approx(43.5, rel=0.01)
