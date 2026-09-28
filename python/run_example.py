from __future__ import annotations

from pathlib import Path

from tpm.faithful import GearSystem
from tpm.optimization import optimize_deltamax
from tpm.plots import generate_validation_plots


def main() -> None:
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

    result = optimize_deltamax(system, bounds=None)
    output_dir = Path("./plots")
    outputs = generate_validation_plots(system, deltamax_um=result["deltamax_um"], output_dir=output_dir)
    print("Reference deflection delta_ref_um:", result["delta_ref_um"])
    print("Automatic bounds:", result["bounds"])
    print("Optimized deltamax_um:", result["deltamax_um"])
    print("Objective (RMS a):", result["objective"])
    print("Success:", result["success"])
    print("Function evaluations:", result["nfev"])
    print("Plots written to:", output_dir)
    for name, path in outputs.items():
        print(f"- {name}: {path}")


if __name__ == "__main__":
    main()
