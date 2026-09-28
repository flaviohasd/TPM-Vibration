import numpy as np

from tpm.geometry import compute_geometry


def test_compute_geometry_basic_values():
    geom = compute_geometry(z1=27, z2=35, module=3.0, alpha_deg=20.0)

    assert geom["dp1"] > 0
    assert geom["dp2"] > 0
    assert geom["db1"] > 0
    assert geom["db2"] > 0
    assert geom["rb1"] > 0
    assert geom["rb2"] > 0
    assert np.isfinite(geom["contact_ratio"])
