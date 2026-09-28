import numpy as np

from tpm.dynamics import simulate_dynamic
from tpm.stiffness import precompute_mesh_state


class DummyMeshState:
    def __init__(self):
        self.k_te = np.ones(20) * 2.0e9
        self.c = np.ones(20) * 200.0


def test_simulate_dynamic_returns_physical_response():
    state = DummyMeshState()
    response = simulate_dynamic(mesh_state=state, mass=10.0, force=8000.0, time_span=(0.0, 0.01))

    assert response.time.shape == (400,)
    assert response.displacement.shape == (400,)
    assert response.velocity.shape == (400,)
    assert response.acceleration.shape == (400,)
    assert np.isfinite(response.acceleration).all()
