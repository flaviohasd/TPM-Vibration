import numpy as np

from tpm.friction import precompute_friction_state


def test_friction_state_is_regularized_and_positive():
    vrel = np.array([0.0, 0.001, 0.01])
    state = precompute_friction_state(vrel=vrel, force_n=8000.0)

    assert state["mu"].shape == vrel.shape
    assert state["damping"].shape == vrel.shape
    assert np.all(state["mu"] > 0.0)
    assert np.all(state["damping"] > 0.0)
