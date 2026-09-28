import numpy as np

from tpm.geometry import compute_geometry
from tpm.stiffness import precompute_mesh_state


def test_precompute_mesh_state_has_expected_shapes():
    geom = compute_geometry(z1=27, z2=35, module=3.0, alpha_deg=20.0)
    state = precompute_mesh_state(
        geom=geom,
        z1=27,
        z2=35,
        step=200,
        w=200.0,
        force=8000.0,
        elasticity_modulus=206e9,
        poisson_ratio=0.3,
        deltamax_um=40.0,
    )

    assert state.theta.shape == (200,)
    assert state.k_te.shape == (200,)
    assert state.k_tem.shape == (200,)
    assert state.c.shape == (200,)
    assert state.cm.shape == (200,)
    assert np.all(state.k_te > 0)
    assert np.all(state.vrel > 0)
