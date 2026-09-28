from tpm.geometry import compute_geometry, compute_mesh


def test_compute_mesh_uses_shared_geometry_data():
    geometry = compute_geometry(z1=27, z2=35, module=3.0, alpha_deg=20.0)
    mesh = compute_mesh(geometry)

    assert mesh["cr"] > 1.0
    assert mesh["PTH"] > 0.0
    assert mesh["thetad"] > 0.0
    assert mesh["thetas"] > 0.0
