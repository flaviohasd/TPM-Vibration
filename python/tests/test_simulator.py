from tpm.simulator import run_simulation


def test_run_simulation_returns_full_pipeline():
    result = run_simulation()

    assert "geometry" in result
    assert "mass_properties" in result
    assert "mesh_state" in result
    assert "response" in result
    assert result["response"].displacement.size > 0
