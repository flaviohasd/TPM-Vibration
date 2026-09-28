from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Dict, Optional

import numpy as np

from .dynamics import simulate_dynamic
from .energy import precompute_energy_state
from .friction import precompute_friction_state


def _normalise_metrics(metrics: Dict[str, Dict[str, float]]) -> Dict[str, Dict[str, float]]:
    return {
        "base": {k: float(v) for k, v in metrics.get("base", {}).items()},
        "modified": {k: float(v) for k, v in metrics.get("modified", {}).items()},
    }


def _build_python_metrics(system, deltamax_um: float) -> Dict[str, Dict[str, float]]:
    base_system = type(system)(
        z1=system.z1,
        z2=system.z2,
        module=system.module,
        alpha_deg=system.alpha_deg,
        b_mm=system.b_mm,
        elasticity_modulus=system.elasticity_modulus,
        poisson_ratio=system.poisson_ratio,
        density=system.density,
        power_w=system.power_w,
        rpm=system.rpm,
        lubricant_density=system.lubricant_density,
        temperature_c=system.temperature_c,
        deltamax_um=0.0,
    )
    modified_system = type(system)(
        z1=system.z1,
        z2=system.z2,
        module=system.module,
        alpha_deg=system.alpha_deg,
        b_mm=system.b_mm,
        elasticity_modulus=system.elasticity_modulus,
        poisson_ratio=system.poisson_ratio,
        density=system.density,
        power_w=system.power_w,
        rpm=system.rpm,
        lubricant_density=system.lubricant_density,
        temperature_c=system.temperature_c,
        deltamax_um=deltamax_um,
    )

    for candidate in (base_system, modified_system):
        energy_state = precompute_energy_state(candidate, step=200)
        mesh_state = energy_state
        mass = candidate.mass_properties["memod"] if getattr(candidate, "deltamax_um", 0.0) > 0 else candidate.mass_properties["me"]
        cycle_time = candidate.operation["cycle"]
        response = simulate_dynamic(mesh_state=mesh_state, mass=mass, force=candidate.operation["F"], time_span=(0.0, cycle_time))
        candidate.__dict__["_response"] = response
        candidate.__dict__["_energy"] = energy_state

    base_response = base_system._response
    modified_response = modified_system._response
    base_energy = base_system._energy
    modified_energy = modified_system._energy

    return {
        "python": {
            "base": {
                "rms_acceleration": float(np.sqrt(np.mean(np.square(base_response.acceleration)))),
                "peak_displacement": float(np.max(np.abs(base_response.displacement))),
                "peak_transmission_error": float(np.max(np.abs(base_energy["transmission_error"]))),
            },
            "modified": {
                "rms_acceleration": float(np.sqrt(np.mean(np.square(modified_response.acceleration)))),
                "peak_displacement": float(np.max(np.abs(modified_response.displacement))),
                "peak_transmission_error": float(np.max(np.abs(modified_energy["transmission_error"]))),
            },
        }
    }


def _build_divergence_lines(metrics: Dict[str, Dict[str, Dict[str, float]]], baseline_metrics: Optional[Dict[str, Dict[str, float]]] = None) -> list[str]:
    lines: list[str] = []
    if not baseline_metrics:
        return lines

    baseline = _normalise_metrics(baseline_metrics)
    lines.append("")
    lines.append("## Baseline divergence")
    lines.append("")
    for case in ("base", "modified"):
        for metric in sorted(baseline[case].keys()):
            python_value = metrics["python"][case][metric]
            baseline_value = baseline[case][metric]
            if baseline_value == 0:
                delta = float("inf")
                pct = float("inf")
            else:
                delta = python_value - baseline_value
                pct = 100.0 * delta / baseline_value
            lines.append(
                f"- {case} {metric}: python={python_value:.6e}, baseline={baseline_value:.6e}, delta={delta:.6e}, percent={pct:.2f}%"
            )
    return lines


def generate_comparison_report(
    system,
    deltamax_um: float = 40.0,
    output_dir: str | Path | None = None,
    baseline_metrics: Optional[Dict[str, Dict[str, float]]] = None,
) -> Dict[str, Path]:
    """Create a structured Python-to-baseline comparison report for validation."""
    output_dir = Path(output_dir or "./comparison")
    output_dir.mkdir(parents=True, exist_ok=True)

    metrics = _build_python_metrics(system, deltamax_um=deltamax_um)

    metrics_path = output_dir / "comparison_metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2))

    csv_path = output_dir / "comparison_metrics.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["engine", "case", "metric", "value"])
        for engine, cases in metrics.items():
            for case, values in cases.items():
                for metric, value in values.items():
                    writer.writerow([engine, case, metric, value])

    summary_lines = [
        "# Comparison report",
        "",
        f"- Delta max evaluated: {deltamax_um} μm",
        "",
        "## Python metrics",
        "",
        f"- Base RMS acceleration: {metrics['python']['base']['rms_acceleration']:.6e}",
        f"- Modified RMS acceleration: {metrics['python']['modified']['rms_acceleration']:.6e}",
        f"- Base peak displacement: {metrics['python']['base']['peak_displacement']:.6e}",
        f"- Modified peak displacement: {metrics['python']['modified']['peak_displacement']:.6e}",
    ]
    summary_lines.extend(_build_divergence_lines(metrics, baseline_metrics=baseline_metrics))
    summary_path = output_dir / "summary.md"
    summary_path.write_text("\n".join(summary_lines) + "\n")

    return {
        "summary_md": summary_path,
        "metrics_json": metrics_path,
        "metrics_csv": csv_path,
    }
