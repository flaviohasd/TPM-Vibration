from __future__ import annotations

from pathlib import Path
from typing import Dict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .dynamics import simulate_dynamic
from .energy import precompute_energy_state
from .friction import precompute_friction_state


def generate_validation_plots(system, deltamax_um: float = 40.0, output_dir: str | Path | None = None) -> Dict[str, Path]:
    """Generate MATLAB-style validation plots for the faithful Python port."""
    output_dir = Path(output_dir or "./plots")
    output_dir.mkdir(parents=True, exist_ok=True)

    base_system = system
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

    base_energy = precompute_energy_state(modified_system if False else base_system, step=1000)
    modified_energy = precompute_energy_state(modified_system, step=1000)

    base_friction = precompute_friction_state(base_energy["vrel"], force_n=base_system.operation["F"], film_stiffness=1.0e8)
    modified_friction = precompute_friction_state(modified_energy["vrel"], force_n=modified_system.operation["F"], film_stiffness=1.0e8)

    base_mesh = {
        "k_te": base_energy["k_te"],
        "c": base_energy["c"],
        "vrel": base_energy["vrel"],
        "transmission_error": base_energy["transmission_error"],
        "damping": base_friction["damping"],
    }
    modified_mesh = {
        "k_te": modified_energy["k_te"],
        "c": modified_energy["c"],
        "vrel": modified_energy["vrel"],
        "transmission_error": modified_energy["transmission_error"],
        "damping": modified_friction["damping"],
    }

    base_response = simulate_dynamic(mesh_state=base_mesh, mass=base_system.mass_properties["me"], force=base_system.operation["F"], time_span=(0.0, 0.01))
    modified_response = simulate_dynamic(mesh_state=modified_mesh, mass=modified_system.mass_properties["me"], force=modified_system.operation["F"], time_span=(0.0, 0.01))

    n_points = len(base_energy["k_te"])
    base_response_resampled = np.interp(np.linspace(0.0, 1.0, n_points), np.linspace(0.0, 1.0, len(base_response.time)), base_response.displacement)
    modified_response_resampled = np.interp(np.linspace(0.0, 1.0, n_points), np.linspace(0.0, 1.0, len(modified_response.time)), modified_response.displacement)
    base_velocity_resampled = np.interp(np.linspace(0.0, 1.0, n_points), np.linspace(0.0, 1.0, len(base_response.time)), base_response.velocity)
    modified_velocity_resampled = np.interp(np.linspace(0.0, 1.0, n_points), np.linspace(0.0, 1.0, len(modified_response.time)), modified_response.velocity)
    base_accel_resampled = np.interp(np.linspace(0.0, 1.0, n_points), np.linspace(0.0, 1.0, len(base_response.time)), base_response.acceleration)
    modified_accel_resampled = np.interp(np.linspace(0.0, 1.0, n_points), np.linspace(0.0, 1.0, len(modified_response.time)), modified_response.acceleration)

    theta = np.linspace(0.0, 2.0 * np.pi, n_points)
    plt.figure(figsize=(8, 6))
    plt.plot(theta, base_energy["vrel"], label="Δmax = 0 μm", linewidth=1.5, color="#1f77b4")
    plt.plot(theta, modified_energy["vrel"], label=f"Δmax = {deltamax_um:.4g} μm", linewidth=1.5, color="#d62728")
    plt.title("Relative velocity at the contact point")
    plt.xlabel("θ₁ (rad)")
    plt.ylabel("v_rel (m/s)")
    plt.xlim([0.0, 2.0 * np.pi])
    plt.legend(loc="upper right")
    contact_path = output_dir / "contact_relations.png"
    plt.savefig(contact_path, dpi=300)
    plt.close()

    # Interpolation results — use real time axis (one gear revolution)
    t_period = 2.0 * np.pi / base_system.operation["w"]
    t_interp = np.linspace(0.0, t_period, n_points)
    plt.figure(figsize=(8, 6))
    plt.plot(t_interp, base_energy["k_te"], label="Δmax = 0 μm", linewidth=1.5, color="#1f77b4")
    plt.plot(t_interp, modified_energy["k_te"], label=f"Δmax = {deltamax_um:.4g} μm", linewidth=1.5, color="#d62728")
    plt.title("Interpolation results")
    plt.xlabel("t (s)")
    plt.ylabel("k(t) (N/m)")
    plt.xlim([0.0, t_period])
    plt.legend(loc="upper right")
    interpolation_path = output_dir / "interpolation_results.png"
    plt.savefig(interpolation_path, dpi=300)
    plt.close()

    plt.figure(figsize=(10, 8))
    plt.subplot(2, 2, 1)
    plt.plot(base_response.time, base_response.displacement, label="Δmax = 0 μm", linewidth=0.8, color="#1f77b4")
    plt.plot(modified_response.time, modified_response.displacement, label=f"Δmax = {deltamax_um:.4g} μm", linewidth=0.8, color="#d62728")
    plt.title("Dynamic transmission error")
    plt.xlabel("t (s)")
    plt.ylabel("DTE (m)")
    plt.xlim([0.0, 0.01])
    plt.legend(loc="upper right")

    plt.subplot(2, 2, 2)
    plt.plot(base_response.time, base_response.acceleration, label="Δmax = 0 μm", linewidth=0.8, color="#1f77b4")
    plt.plot(modified_response.time, modified_response.acceleration, label=f"Δmax = {deltamax_um:.4g} μm", linewidth=0.8, color="#d62728")
    plt.title("Meshing acceleration")
    plt.xlabel("t (s)")
    plt.ylabel("a (m/s²)")
    plt.xlim([0.0, 0.01])
    plt.legend(loc="upper right")

    plt.subplot(2, 2, 3)
    plt.plot(base_response.time, base_response.velocity, label="Δmax = 0 μm", linewidth=0.8, color="#1f77b4")
    plt.plot(modified_response.time, modified_response.velocity, label=f"Δmax = {deltamax_um:.4g} μm", linewidth=0.8, color="#d62728")
    plt.title("Meshing velocity")
    plt.xlabel("t (s)")
    plt.ylabel("v (m/s)")
    plt.xlim([0.0, 0.01])
    plt.legend(loc="upper right")

    plt.subplot(2, 2, 4)
    base_force = base_energy["k_te"] * np.maximum(base_response_resampled, 0.0)
    modified_force = modified_energy["k_te"] * np.maximum(modified_response_resampled, 0.0)
    plt.plot(np.linspace(0.0, 1.0, n_points), base_force, label="Δmax = 0 μm", linewidth=0.8, color="#1f77b4")
    plt.plot(np.linspace(0.0, 1.0, n_points), modified_force, label=f"Δmax = {deltamax_um:.4g} μm", linewidth=0.8, color="#d62728")
    plt.title("Dynamic meshing force")
    plt.xlabel("Γ")
    plt.ylabel("DMF (N)")
    plt.legend(loc="upper right")
    plt.tight_layout()
    meshing_path = output_dir / "meshing_dynamics.png"
    plt.savefig(meshing_path, dpi=300)
    plt.close()

    tau = np.linspace(-0.5, 0.5, len(base_energy["k_te"]))
    plt.figure(figsize=(10, 8))
    plt.subplot(3, 2, 1)
    plt.plot(tau, base_energy["k_te"], label="Δmax = 0 μm", linewidth=1.0, color="#1f77b4")
    plt.plot(tau, modified_energy["k_te"], label=f"Δmax = {deltamax_um:.4g} μm", linewidth=1.0, color="#d62728")
    plt.title("k(t) along normalized coordinate")
    plt.xlabel("Γ")
    plt.ylabel("k(t) (N/m)")
    plt.xlim([-0.5, 0.5])
    plt.legend(loc="upper right")

    plt.subplot(3, 2, 2)
    plt.plot(tau, base_energy["transmission_error"], label="Δmax = 0 μm", linewidth=1.0, color="#1f77b4")
    plt.plot(tau, modified_energy["transmission_error"], label=f"Δmax = {deltamax_um:.4g} μm", linewidth=1.0, color="#d62728")
    plt.title("Transmission error")
    plt.xlabel("Γ")
    plt.ylabel("TE (m)")
    plt.xlim([-0.5, 0.5])
    plt.legend(loc="upper right")

    plt.subplot(3, 2, 3)
    plt.plot(tau, base_response_resampled, label="Δmax = 0 μm", linewidth=1.0, color="#1f77b4")
    plt.plot(tau, modified_response_resampled, label=f"Δmax = {deltamax_um:.4g} μm", linewidth=1.0, color="#d62728")
    plt.title("DTE along normalized coordinate")
    plt.xlabel("Γ")
    plt.ylabel("DTE (m)")
    plt.xlim([-0.5, 0.5])
    plt.legend(loc="upper right")

    plt.subplot(3, 2, 4)
    plt.plot(tau, base_accel_resampled, label="Δmax = 0 μm", linewidth=1.0, color="#1f77b4")
    plt.plot(tau, modified_accel_resampled, label=f"Δmax = {deltamax_um:.4g} μm", linewidth=1.0, color="#d62728")
    plt.title("Acceleration along normalized coordinate")
    plt.xlabel("Γ")
    plt.ylabel("a (m/s²)")
    plt.xlim([-0.5, 0.5])
    plt.legend(loc="upper right")

    plt.subplot(3, 2, 5)
    plt.plot(tau, base_energy["k_te"] * np.maximum(base_response_resampled, 0.0), label="Δmax = 0 μm", linewidth=1.0, color="#1f77b4")
    plt.plot(tau, modified_energy["k_te"] * np.maximum(modified_response_resampled, 0.0), label=f"Δmax = {deltamax_um:.4g} μm", linewidth=1.0, color="#d62728")
    plt.title("TVMF")
    plt.xlabel("Γ")
    plt.ylabel("TVMF (N)")
    plt.xlim([-0.5, 0.5])
    plt.legend(loc="upper right")

    plt.subplot(3, 2, 6)
    plt.plot(tau, base_energy["k_te"] * np.maximum(base_response_resampled, 0.0), label="Δmax = 0 μm", linewidth=1.0, color="#1f77b4")
    plt.plot(tau, modified_energy["k_te"] * np.maximum(modified_response_resampled, 0.0), label=f"Δmax = {deltamax_um:.4g} μm", linewidth=1.0, color="#d62728")
    plt.title("TVMF comparison")
    plt.xlabel("Γ")
    plt.ylabel("TVMF (N)")
    plt.xlim([-0.5, 0.5])
    plt.legend(loc="upper right")
    plt.tight_layout()
    normalized_path = output_dir / "normalized_comparisons.png"
    plt.savefig(normalized_path, dpi=300)
    plt.close()

    # Frequency analysis — FFT of the resampled displacement signal
    # dt is derived from the simulated time span divided by number of samples
    t_sim = float(base_response.time[-1] - base_response.time[0])
    dt_dyn = t_sim / len(base_response_resampled)
    disp_base_ac = base_response_resampled - np.mean(base_response_resampled)
    disp_mod_ac = modified_response_resampled - np.mean(modified_response_resampled)
    freqs = np.fft.rfftfreq(len(disp_base_ac), d=dt_dyn)
    spectrum_base = np.abs(np.fft.rfft(disp_base_ac)) / len(disp_base_ac)
    spectrum_modified = np.abs(np.fft.rfft(disp_mod_ac)) / len(disp_mod_ac)
    plt.figure(figsize=(8, 6))
    plt.plot(freqs, spectrum_base, label="Δmax = 0 μm", linewidth=1.0, color="#1f77b4")
    plt.plot(freqs, spectrum_modified, label=f"Δmax = {deltamax_um:.4g} μm", linewidth=1.0, color="#d62728")
    plt.title("Frequency analysis")
    plt.xlabel("f (Hz)")
    plt.ylabel("Magnitude")
    plt.xlim([0.0, min(500.0, freqs[-1])])
    plt.legend(loc="upper right")
    frequency_path = output_dir / "frequency_analysis.png"
    plt.savefig(frequency_path, dpi=300)
    plt.close()

    return {
        "contact_relations": contact_path,
        "interpolation_results": interpolation_path,
        "meshing_dynamics": meshing_path,
        "normalized_comparisons": normalized_path,
        "frequency_analysis": frequency_path,
    }
