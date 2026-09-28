import numpy as np

from tpm.energy import _integral_compliance, precompute_energy_state


def test_energy_state_is_positive_and_vectorized():
    from tpm.faithful import GearSystem

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

    state = precompute_energy_state(system, step=200)

    assert state["k_te"].shape == (200,)
    assert state["k_tem"].shape == (200,)
    assert state["c"].shape == (200,)
    assert state["vrel"].shape == (200,)
    assert state["transmission_error"].shape == (200,)
    assert (state["k_te"] > 0).all()
    assert (state["vrel"] > 0).all()


def test_energy_state_depends_on_profile_relief_amplitude():
    from tpm.faithful import GearSystem

    base_system = GearSystem(
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
        deltamax_um=20.0,
    )
    modified_system = GearSystem(
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
        deltamax_um=60.0,
    )

    base_state = precompute_energy_state(base_system, step=200)
    modified_state = precompute_energy_state(modified_system, step=200)

    assert np.max(np.abs(modified_state["k_te"] - base_state["k_te"])) > 1e6
    assert np.max(np.abs(modified_state["transmission_error"] - base_state["transmission_error"])) > 1e-8


def test_integral_compliance_matches_matlab_formula_structure():
    value = _integral_compliance("ikbi", 0.2, 0.3, 0.4, 206e9, 25.0, 0.3)
    assert np.isfinite(value)
    assert value > 0.0
