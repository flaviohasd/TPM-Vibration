import numpy as np

from tpm.faithful import GearSystem
from tpm.optimization import _objective, optimize_deltamax


def test_optimize_deltamax_returns_finite_value():
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

    result = optimize_deltamax(system, bounds=(30.0, 50.0), n_points=6)

    assert np.isfinite(result["deltamax_um"])
    assert result["objective"] >= 0.0


def test_objective_changes_with_relief_amplitude():
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

    low_objective = _objective(30.0, system, step=80)
    high_objective = _objective(50.0, system, step=80)

    assert np.isfinite(low_objective)
    assert np.isfinite(high_objective)
    assert not np.isclose(low_objective, high_objective)


def test_optimize_deltamax_automatic_bounds():
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
    assert 20.0 < result["delta_ref_um"] < 60.0
    assert 40.0 < result["deltamax_um"] < 46.0
