from __future__ import annotations

from pathlib import Path
from typing import Dict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .dynamics import simulate_dynamic
from .energy import precompute_energy_state


def generate_validation_plots(
    system,
    deltamax_um: float = 40.0,
    output_dir: str | Path | None = None,
) -> Dict[str, Path]:
    """Generate high-resolution validation plots with both full rotation (time domain)
    and single tooth cycle (normalized coordinate Gamma in [-0.5, 0.5]) visualizations."""
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
        deltamax_um=deltamax_um,
    )

    step = 1000
    base_energy = precompute_energy_state(base_system, step=step)
    modified_energy = precompute_energy_state(modified_system, step=step)

    w = float(base_system.operation["w"])
    t_rot = 2.0 * np.pi / w
    force_n = float(base_system.operation["F"])
    z1 = int(base_system.z1)

    base_mesh = {
        "k_te": base_energy["k_te"],
        "c": base_energy["c"],
        "vrel": base_energy["vrel"],
        "transmission_error": base_energy["transmission_error"],
        "forcing": base_energy.get("forcing", np.full(step, force_n)),
    }
    modified_mesh = {
        "k_te": modified_energy["k_te"],
        "c": modified_energy["cm"],
        "vrel": modified_energy["vrel"],
        "transmission_error": modified_energy["transmission_error"],
        "forcing": modified_energy.get("forcing", np.full(step, force_n)),
    }

    # Simulate dynamic response across full gear revolution (T_rot)
    num_pts = 10000
    base_response = simulate_dynamic(
        mesh_state=base_mesh,
        mass=base_system.mass_properties["me"],
        force=force_n,
        time_span=(0.0, t_rot),
        num_points=num_pts,
    )
    modified_response = simulate_dynamic(
        mesh_state=modified_mesh,
        mass=modified_system.mass_properties["memod"],
        force=force_n,
        time_span=(0.0, t_rot),
        num_points=num_pts,
    )

    t_grid = np.linspace(0.0, t_rot, step)
    base_disp_resampled = np.interp(t_grid, base_response.time, base_response.displacement)
    mod_disp_resampled = np.interp(t_grid, modified_response.time, modified_response.displacement)
    base_vel_resampled = np.interp(t_grid, base_response.time, base_response.velocity)
    mod_vel_resampled = np.interp(t_grid, modified_response.time, modified_response.velocity)
    base_acc_resampled = np.interp(t_grid, base_response.time, base_response.acceleration)
    mod_acc_resampled = np.interp(t_grid, modified_response.time, modified_response.acceleration)

    base_dmf = base_energy["k_te"] * np.maximum(base_disp_resampled, 0.0)
    mod_dmf = modified_energy["k_te"] * np.maximum(mod_disp_resampled, 0.0)

    # Slice out one representative tooth meshing cycle in steady-state (penultimate tooth)
    pd = step / z1
    a = int(np.floor(pd * (z1 - 2 + 0.05)))
    b = int(np.floor(pd * (z1 - 0.4)))
    tau = np.linspace(-0.5, 0.5, b - a)

    # Single tooth contact kinematics (Tooth Engagement Cycle: 0 to 2*thetad + thetas)
    thetad = float(base_system.mesh["thetad"])
    thetas = float(base_system.mesh["thetas"])
    thetab1 = float(base_system.mesh["thetab1"])
    alpha10 = float(base_system.mesh["alpha10"])
    rb1 = float(base_system.geometry["rb1"])
    rp1 = float(base_system.geometry["rp1"])
    rp2 = float(base_system.geometry["rp2"])
    L = float(base_system.mesh["L"])
    alpha_rad = float(base_system.geometry["alpha"])

    theta_cycle = np.linspace(0.0, 2.0 * thetad + thetas, 200)
    phi1_cycle = theta_cycle + alpha10
    x1_cycle = rb1 * (thetab1 + phi1_cycle)
    x2_cycle = L - x1_cycle
    vrel_single_cycle = np.abs(w * x1_cycle / 1000.0 - (rp1 / rp2) * w * x2_cycle / 1000.0)

    # -------------------------------------------------------------
    # 1. CONTACT RELATIONS (Tooth Engagement Cycle & Contact Dynamics)
    # -------------------------------------------------------------
    plt.figure(figsize=(9, 8))

    plt.subplot(3, 1, 1)
    plt.plot(np.rad2deg(theta_cycle), vrel_single_cycle, color="#1f77b4", linewidth=1.5, label="Contact pair")
    plt.title("Relative velocity at the contact point (Single tooth cycle)")
    plt.xlabel("θ₁ (°)")
    plt.ylabel("v_rel (m/s)")
    plt.xlim([0.0, np.rad2deg(2.0 * thetad + thetas)])
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend(loc="upper right")

    plt.subplot(3, 1, 2)
    plt.plot(tau, base_energy["r1"][a:b], color="#2ca02c", linewidth=1.5, label="Load sharing ratio R₁")
    plt.title("Load sharing ratio along normalized coordinate")
    plt.xlabel("Γ")
    plt.ylabel("R₁")
    plt.xlim([-0.5, 0.5])
    plt.ylim([0.0, 1.1])
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend(loc="lower center")

    plt.subplot(3, 1, 3)
    plt.plot(tau, modified_energy["delta"][a:b, 0] * 1e6, color="#d62728", linewidth=1.5, label="Δ₁ (Pair 1)")
    plt.plot(tau, modified_energy["delta"][a:b, 1] * 1e6, color="#ff7f0e", linewidth=1.5, linestyle="--", label="Δ₂ (Pair 2)")
    plt.title(f"Tooth profile modification along normalized coordinate (Δmax = {deltamax_um:.4g} μm)")
    plt.xlabel("Γ")
    plt.ylabel("Δ (μm)")
    plt.xlim([-0.5, 0.5])
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend(loc="upper right")

    plt.tight_layout()
    contact_path = output_dir / "contact_relations.png"
    plt.savefig(contact_path, dpi=300)
    plt.close()

    # -------------------------------------------------------------
    # 2. FULL ROTATION TIME-DOMAIN CHARACTERISTICS (k(t), c(t), vrel(t), TE(t))
    # -------------------------------------------------------------
    plt.figure(figsize=(10, 8))

    plt.subplot(2, 2, 1)
    plt.plot(t_grid, base_energy["k_te"], label="Δmax = 0 μm", linewidth=1.0, color="#1f77b4")
    plt.plot(t_grid, modified_energy["k_te"], label=f"Δmax = {deltamax_um:.4g} μm", linewidth=1.0, color="#d62728")
    plt.title("Time-varying mesh stiffness (Full rotation)")
    plt.xlabel("t (s)")
    plt.ylabel("k(t) (N/m)")
    plt.xlim([0.0, t_rot])
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend(loc="upper right")

    plt.subplot(2, 2, 2)
    plt.plot(t_grid, base_energy["c"], label="Δmax = 0 μm", linewidth=1.0, color="#1f77b4")
    plt.plot(t_grid, modified_energy["cm"], label=f"Δmax = {deltamax_um:.4g} μm", linewidth=1.0, color="#d62728")
    plt.title("Time-varying damping coefficient (Full rotation)")
    plt.xlabel("t (s)")
    plt.ylabel("c(t) (N·s/m)")
    plt.xlim([0.0, t_rot])
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend(loc="upper right")

    plt.subplot(2, 2, 3)
    plt.plot(t_grid, base_energy["transmission_error"] * 1e6, label="Δmax = 0 μm", linewidth=1.0, color="#1f77b4")
    plt.plot(t_grid, modified_energy["transmission_error"] * 1e6, label=f"Δmax = {deltamax_um:.4g} μm", linewidth=1.0, color="#d62728")
    plt.title("Static transmission error (Full rotation)")
    plt.xlabel("t (s)")
    plt.ylabel("TE (μm)")
    plt.xlim([0.0, t_rot])
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend(loc="upper right")

    plt.subplot(2, 2, 4)
    plt.plot(t_grid, base_energy["vrel"], label="Δmax = 0 μm", linewidth=1.0, color="#1f77b4")
    plt.title("Contact relative velocity (Full rotation)")
    plt.xlabel("t (s)")
    plt.ylabel("v_rel (m/s)")
    plt.xlim([0.0, t_rot])
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend(loc="upper right")

    plt.tight_layout()
    interpolation_path = output_dir / "interpolation_results.png"
    plt.savefig(interpolation_path, dpi=300)
    plt.close()

    # -------------------------------------------------------------
    # 3. MESHING DYNAMICS (Full Rotation Time-Domain Response)
    # -------------------------------------------------------------
    plt.figure(figsize=(10, 8))

    plt.subplot(2, 2, 1)
    plt.plot(base_response.time, base_response.displacement * 1e6, label="Δmax = 0 μm", linewidth=0.8, color="#1f77b4")
    plt.plot(modified_response.time, modified_response.displacement * 1e6, label=f"Δmax = {deltamax_um:.4g} μm", linewidth=0.8, color="#d62728")
    plt.title("Dynamic transmission error")
    plt.xlabel("t (s)")
    plt.ylabel("DTE (μm)")
    plt.xlim([0.0, t_rot])
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend(loc="upper right")

    plt.subplot(2, 2, 2)
    plt.plot(base_response.time, base_response.acceleration, label="Δmax = 0 μm", linewidth=0.8, color="#1f77b4")
    plt.plot(modified_response.time, modified_response.acceleration, label=f"Δmax = {deltamax_um:.4g} μm", linewidth=0.8, color="#d62728")
    plt.title("Meshing acceleration")
    plt.xlabel("t (s)")
    plt.ylabel("a (m/s²)")
    plt.xlim([0.0, t_rot])
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend(loc="upper right")

    plt.subplot(2, 2, 3)
    plt.plot(base_response.time, base_response.velocity, label="Δmax = 0 μm", linewidth=0.8, color="#1f77b4")
    plt.plot(modified_response.time, modified_response.velocity, label=f"Δmax = {deltamax_um:.4g} μm", linewidth=0.8, color="#d62728")
    plt.title("Meshing velocity")
    plt.xlabel("t (s)")
    plt.ylabel("v (m/s)")
    plt.xlim([0.0, t_rot])
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend(loc="upper right")

    plt.subplot(2, 2, 4)
    plt.plot(t_grid, base_dmf, label="Δmax = 0 μm", linewidth=0.8, color="#1f77b4")
    plt.plot(t_grid, mod_dmf, label=f"Δmax = {deltamax_um:.4g} μm", linewidth=0.8, color="#d62728")
    plt.title("Dynamic meshing force")
    plt.xlabel("t (s)")
    plt.ylabel("DMF (N)")
    plt.xlim([0.0, t_rot])
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend(loc="upper right")

    plt.tight_layout()
    meshing_path = output_dir / "meshing_dynamics.png"
    plt.savefig(meshing_path, dpi=300)
    plt.close()

    # -------------------------------------------------------------
    # 4. SINGLE MESHING CYCLE COMPARISONS (Normalized Coordinate Γ in [-0.5, 0.5])
    # -------------------------------------------------------------
    plt.figure(figsize=(11, 9))

    plt.subplot(3, 2, 1)
    plt.plot(tau, base_energy["k_te"][a:b], "-.", label="Δmax = 0 μm", linewidth=1.4, color="#1f77b4")
    plt.plot(tau, modified_energy["k_te"][a:b], label=f"Δmax = {deltamax_um:.4g} μm", linewidth=1.4, color="#d62728")
    plt.title("k(t) along normalized coordinate")
    plt.xlabel("Γ")
    plt.ylabel("k(t) (N/m)")
    plt.xlim([-0.5, 0.5])
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend(loc="upper right")

    plt.subplot(3, 2, 2)
    plt.plot(tau, base_energy["transmission_error"][a:b] * 1e6, "-.", label="TE", linewidth=1.4, color="#1f77b4")
    plt.plot(tau, modified_energy["transmission_error"][a:b] * 1e6, label=f"TE_M ({deltamax_um:.4g} μm)", linewidth=1.4, color="#d62728")
    plt.title("Static transmission error along normalized coordinate")
    plt.xlabel("Γ")
    plt.ylabel("TE (μm)")
    plt.xlim([-0.5, 0.5])
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend(loc="upper right")

    plt.subplot(3, 2, 3)
    plt.plot(tau, base_disp_resampled[a:b] * 1e6, "-.", label="Δmax = 0 μm", linewidth=1.4, color="#1f77b4")
    plt.plot(tau, mod_disp_resampled[a:b] * 1e6, label=f"Δmax = {deltamax_um:.4g} μm", linewidth=1.4, color="#d62728")
    plt.title("DTE along normalized coordinate")
    plt.xlabel("Γ")
    plt.ylabel("DTE (μm)")
    plt.xlim([-0.5, 0.5])
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend(loc="upper right")

    plt.subplot(3, 2, 4)
    plt.plot(tau, base_acc_resampled[a:b], "-.", label="Δmax = 0 μm", linewidth=1.4, color="#1f77b4")
    plt.plot(tau, mod_acc_resampled[a:b], label=f"Δmax = {deltamax_um:.4g} μm", linewidth=1.4, color="#d62728")
    plt.title("Acceleration along normalized coordinate")
    plt.xlabel("Γ")
    plt.ylabel("a (m/s²)")
    plt.xlim([-0.5, 0.5])
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend(loc="upper right")

    # Subplot 5: TVMF * R1 with twin axis for R1
    ax5 = plt.subplot(3, 2, 5)
    r1_slice = base_energy["r1"][a:b]
    ax5.plot(tau, base_dmf[a:b] * r1_slice, "--", label="TVMF (Δmax = 0)", linewidth=1.4, color="#1f77b4")
    ax5.plot(tau, mod_dmf[a:b] * r1_slice, label=f"TVMF (Δmax = {deltamax_um:.4g})", linewidth=1.4, color="#d62728")
    ax5.set_xlabel("Γ")
    ax5.set_ylabel("Tooth Force (N)")
    ax5.set_xlim([-0.5, 0.5])
    ax5.grid(True, linestyle=":", alpha=0.6)

    ax5_twin = ax5.twinx()
    ax5_twin.plot(tau, r1_slice, "k:", linewidth=1.2, label="R₁")
    ax5_twin.set_ylabel("R₁", color="k")
    ax5_twin.set_ylim([0.0, 1.2])
    ax5.set_title("Tooth dynamic force with load sharing")

    lines_5 = ax5.get_lines() + ax5_twin.get_lines()
    labels_5 = [line.get_label() for line in lines_5]
    ax5.legend(lines_5, labels_5, loc="lower center")

    plt.subplot(3, 2, 6)
    plt.plot(tau, base_dmf[a:b], "-.", label="Δmax = 0 μm", linewidth=1.4, color="#1f77b4")
    plt.plot(tau, mod_dmf[a:b], label=f"Δmax = {deltamax_um:.4g} μm", linewidth=1.4, color="#d62728")
    plt.title("Total dynamic meshing force (TVMF)")
    plt.xlabel("Γ")
    plt.ylabel("TVMF (N)")
    plt.xlim([-0.5, 0.5])
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend(loc="upper right")

    plt.tight_layout()
    normalized_path = output_dir / "normalized_comparisons.png"
    plt.savefig(normalized_path, dpi=300)
    plt.close()

    # -------------------------------------------------------------
    # 5. FREQUENCY ANALYSIS (Steady-state FFT spectrum)
    # -------------------------------------------------------------
    steady_slice = slice(len(base_response.time) // 2, len(base_response.time))
    dt_steady = float(base_response.time[1] - base_response.time[0])
    acc_steady_base = base_response.acceleration[steady_slice]
    acc_steady_mod = modified_response.acceleration[steady_slice]

    n_steady = len(acc_steady_base)
    freqs = np.fft.rfftfreq(n_steady, d=dt_steady)
    spectrum_base = np.abs(np.fft.rfft(acc_steady_base - np.mean(acc_steady_base))) / n_steady
    spectrum_mod = np.abs(np.fft.rfft(acc_steady_mod - np.mean(acc_steady_mod))) / n_steady

    plt.figure(figsize=(9, 6))
    plt.plot(freqs, spectrum_base, label="Δmax = 0 μm", linewidth=1.2, color="#1f77b4")
    plt.plot(freqs, spectrum_mod, label=f"Δmax = {deltamax_um:.4g} μm", linewidth=1.2, color="#d62728")
    plt.title("Meshing acceleration frequency spectrum (Steady-state)")
    plt.xlabel("f (Hz)")
    plt.ylabel("Acceleration Magnitude (m/s²)")
    plt.xlim([0.0, 1500.0])
    mask_vis = freqs <= 1500.0
    y_max = max(float(np.max(spectrum_base[mask_vis])), float(np.max(spectrum_mod[mask_vis])))
    plt.ylim([0.0, max(y_max * 1.15, 10.0)])
    plt.grid(True, linestyle=":", alpha=0.6)
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
