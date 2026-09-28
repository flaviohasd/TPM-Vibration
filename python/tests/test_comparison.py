import json
from pathlib import Path

from tpm.comparison import generate_comparison_report
from tpm.faithful import GearSystem


def test_generate_comparison_report_writes_files(tmp_path: Path) -> None:
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

    outputs = generate_comparison_report(system, deltamax_um=40.0, output_dir=tmp_path)

    assert {"summary_md", "metrics_json", "metrics_csv"}.issubset(outputs.keys())

    metrics = json.loads((tmp_path / "comparison_metrics.json").read_text())
    assert "python" in metrics
    assert "base" in metrics["python"]
    assert "modified" in metrics["python"]
    assert metrics["python"]["base"]["rms_acceleration"] > 0.0
    assert metrics["python"]["modified"]["rms_acceleration"] > 0.0


def test_generate_comparison_report_with_baseline(tmp_path: Path) -> None:
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

    baseline = {
        "base": {"rms_acceleration": 1.0, "peak_displacement": 1.0, "peak_transmission_error": 1.0},
        "modified": {"rms_acceleration": 1.0, "peak_displacement": 1.0, "peak_transmission_error": 1.0},
    }

    outputs = generate_comparison_report(system, deltamax_um=40.0, output_dir=tmp_path, baseline_metrics=baseline)

    summary = (tmp_path / "summary.md").read_text()
    assert "Baseline divergence" in summary
    assert "percent" in summary.lower()
    assert outputs["summary_md"].exists()
