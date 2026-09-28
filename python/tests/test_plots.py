from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from tpm.faithful import GearSystem
from tpm.plots import generate_validation_plots


def test_generate_validation_plots_writes_files(tmp_path: Path) -> None:
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

    outputs = generate_validation_plots(system, deltamax_um=40.0, output_dir=tmp_path)

    assert set(outputs) == {
        "contact_relations",
        "interpolation_results",
        "meshing_dynamics",
        "normalized_comparisons",
        "frequency_analysis",
    }
    for path in outputs.values():
        assert path.exists()


def test_generate_validation_plots_with_resolution_fast(tmp_path: Path) -> None:
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

    outputs = generate_validation_plots(system, deltamax_um=40.0, output_dir=tmp_path, resolution="fast")
    for path in outputs.values():
        assert path.exists()

