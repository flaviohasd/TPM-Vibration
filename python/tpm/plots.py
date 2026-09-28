from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .dynamics import simulate_dynamic
from .energy import precompute_energy_state
from .friction import precompute_friction_state


def generate_validation_plots(
    system,
    deltamax_um: float = 40.0,
    output_dir: str | Path | None = None,
    step: Optional[int] = None,
    points_per_tooth: Optional[int] = None,
    resolution: str = "ultra",
) -> Dict[str, Path]:
    """Generate MATLAB-style validation plots matching the TCC and paper."""
    output_dir = Path(output_dir or "./plots")
    output_dir.mkdir(parents=True, exist_ok=True)

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

    z1 = int(base_system.z1)
    if points_per_tooth is not None:
        step = int(z1 * points_per_tooth)
    elif step is not None:
        step = int(step)
    else:
        preset_map = {"fast": 50, "standard": 100, "fine": 200, "ultra": 400, "publication": 400}
        step = int(z1 * preset_map.get(resolution, 400))

    w = base_system.operation["w"]
    t_period = 2.0 * np.pi / w
    F = base_system.operation["F"]
    me = base_system.mass_properties["me"]
    memod = modified_system.mass_properties["memod"]

    base_energy = precompute_energy_state(base_system, step=step)
    modified_energy = precompute_energy_state(modified_system, step=step)

    # Full rotation dynamic simulation using physical TVMS and excitation
    base_mesh = {
        "k_te_base": base_energy["k_te_base"],
        "c": base_energy["c"],
        "forcing": base_energy["forcing"],
    }
    modified_mesh = {
        "k_te_base": modified_energy["k_te_base"],
        "c": modified_energy["c"],
        "forcing": modified_energy["forcing"],
    }

    base_response = simulate_dynamic(
        mesh_state=base_mesh,
        mass=me,
        force=F,
        time_span=(0.0, t_period),
        num_points=step,
    )
    modified_response = simulate_dynamic(
        mesh_state=modified_mesh,
        mass=memod,
        force=F,
        time_span=(0.0, t_period),
        num_points=step,
    )

    t = base_response.time
    DTE_base = base_response.displacement * 1e6  # um
    DTE_mod = modified_response.displacement * 1e6  # um
    a_base = base_response.acceleration
    a_mod = modified_response.acceleration
    v_base = base_response.velocity
    v_mod = modified_response.velocity
    dmf_base = base_response.dmf
    dmf_mod = modified_response.dmf

    te_base_raw = base_energy.get("transmission_error_base", base_energy.get("te", base_energy.get("transmission_error")))
    te_mod_raw = modified_energy.get("transmission_error_modified", modified_energy.get("tem", modified_energy.get("transmission_error")))
    te_base_um = np.asarray(te_base_raw, dtype=float) * 1e6
    te_mod_um = np.asarray(te_mod_raw, dtype=float) * 1e6

    # ----------------------------------------------------
    # 1. Contact relations (Single tooth cycle)
    # ----------------------------------------------------
    thetad = base_system.mesh["thetad"]
    thetas = base_system.mesh["thetas"]
    th_contact_max_deg = np.rad2deg(2.0 * thetad + thetas)

    th_single_vec = np.linspace(0.0, 2.0 * thetad + thetas, 500)
    rb1 = base_system.geometry["rb1"]
    rp1 = base_system.geometry["rp1"]
    rp2 = base_system.geometry["rp2"]
    L = base_system.mesh["L"]
    alpha10 = base_system.mesh["alpha10"]
    thetab1 = base_system.mesh["thetab1"]
    phi1_single = th_single_vec + alpha10
    x1_single = rb1 * (thetab1 + phi1_single)
    x2_single = L - x1_single
    vrel_single = np.abs(w * x1_single / 1000.0 - (rp1 / rp2) * w * x2_single / 1000.0)

    tau_single = np.linspace(-0.5, 0.5, 500)
    r1_single = np.zeros(500)
    for idx_tau, g in enumerate(tau_single):
        if abs(g) <= 0.1:
            r1_single[idx_tau] = 1.0
        elif g < -0.1:
            r1_single[idx_tau] = 0.43 + 0.12 * (g + 0.5) / 0.4
        else:
            r1_single[idx_tau] = 0.42 + 0.12 * (g - 0.1) / 0.4

    delta1_um = np.zeros(500)
    delta2_um = np.zeros(500)
    for idx_tau, g in enumerate(tau_single):
        if g < -0.1:
            delta1_um[idx_tau] = deltamax_um * (-0.1 - g) / 0.4
            delta2_um[idx_tau] = deltamax_um * (g + 0.5) / 0.4
        elif g > 0.1:
            delta1_um[idx_tau] = deltamax_um * (0.5 - g) / 0.4
            delta2_um[idx_tau] = deltamax_um * (g - 0.1) / 0.4
        else:
            delta1_um[idx_tau] = 0.0
            delta2_um[idx_tau] = 0.0

    fig, axes = plt.subplots(3, 1, figsize=(9, 8))
    axes[0].plot(np.rad2deg(th_single_vec), vrel_single, label="Contact pair", linewidth=1.5, color="#1f77b4")
    axes[0].set_title("Relative velocity at the contact point (Single tooth cycle)")
    axes[0].set_xlabel("θ₁ (°)")
    axes[0].set_ylabel("v_rel (m/s)")
    axes[0].set_xlim([0.0, th_contact_max_deg])
    axes[0].grid(True, linestyle=":", alpha=0.6)
    axes[0].legend(loc="upper right")

    axes[1].plot(tau_single, r1_single, label="Load sharing ratio R₁", linewidth=1.5, color="#2ca02c")
    axes[1].set_title("Load sharing ratio along normalized coordinate")
    axes[1].set_xlabel("Γ")
    axes[1].set_ylabel("R₁")
    axes[1].set_xlim([-0.5, 0.5])
    axes[1].set_ylim([0.0, 1.1])
    axes[1].grid(True, linestyle=":", alpha=0.6)
    axes[1].legend(loc="lower center")

    axes[2].plot(tau_single, delta1_um, label="Δ₁ (Pair 1)", linewidth=1.5, color="#d62728")
    axes[2].plot(tau_single, delta2_um, label="Δ₂ (Pair 2)", linewidth=1.5, linestyle="--", color="#ff7f0e")
    axes[2].set_title(f"Tooth profile modification along normalized coordinate (Δmax = {deltamax_um:.2f} μm)")
    axes[2].set_xlabel("Γ")
    axes[2].set_ylabel("Δ (μm)")
    axes[2].set_xlim([-0.5, 0.5])
    axes[2].grid(True, linestyle=":", alpha=0.6)
    axes[2].legend(loc="upper right")

    plt.tight_layout()
    contact_path = output_dir / "contact_relations.png"
    plt.savefig(contact_path, dpi=300)
    plt.close()

    # ----------------------------------------------------
    # 2. Interpolation results (Full rotation t in [0, T_rot])
    # ----------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    axes[0, 0].plot(t, base_energy["k_te_base"], label="Δmax = 0 μm", linewidth=1.0, color="#1f77b4")
    axes[0, 0].plot(t, modified_energy["k_tem"], label=f"Δmax = {deltamax_um:.2f} μm", linewidth=1.0, color="#d62728")
    axes[0, 0].set_title("Time-varying mesh stiffness (Full rotation)")
    axes[0, 0].set_xlabel("t (s)")
    axes[0, 0].set_ylabel("k(t) (N/m)")
    axes[0, 0].set_xlim([0.0, t_period])
    axes[0, 0].grid(True, linestyle=":", alpha=0.6)
    axes[0, 0].legend(loc="upper right")

    axes[0, 1].plot(t, base_energy["c"], label="Δmax = 0 μm", linewidth=1.0, color="#1f77b4")
    axes[0, 1].plot(t, modified_energy["cm"], label=f"Δmax = {deltamax_um:.2f} μm", linewidth=1.0, color="#d62728")
    axes[0, 1].set_title("Time-varying damping coefficient (Full rotation)")
    axes[0, 1].set_xlabel("t (s)")
    axes[0, 1].set_ylabel("c(t) (N·s/m)")
    axes[0, 1].set_xlim([0.0, t_period])
    axes[0, 1].grid(True, linestyle=":", alpha=0.6)
    axes[0, 1].legend(loc="upper right")

    axes[1, 0].plot(t, te_base_um, label="Δmax = 0 μm", linewidth=1.0, color="#1f77b4")
    axes[1, 0].plot(t, te_mod_um, label=f"Δmax = {deltamax_um:.2f} μm", linewidth=1.0, color="#d62728")
    axes[1, 0].set_title("Static transmission error (Full rotation)")
    axes[1, 0].set_xlabel("t (s)")
    axes[1, 0].set_ylabel("TE (μm)")
    axes[1, 0].set_xlim([0.0, t_period])
    axes[1, 0].grid(True, linestyle=":", alpha=0.6)
    axes[1, 0].legend(loc="upper right")

    # Construct exact full-rotation relative velocity profile with zero sampling jitter
    th_contact = 2.0 * thetad + thetas
    t_cycle = th_contact / w
    n_pts_per_cycle = max(100, int(round(step / z1)))
    n_full_cycles = int(np.floor(t_period / t_cycle))
    remainder_time = t_period - n_full_cycles * t_cycle

    th_single_cycle = np.linspace(0.0, th_contact, n_pts_per_cycle)
    phi_single_cycle = th_single_cycle + alpha10
    x1_sc = rb1 * (thetab1 + phi_single_cycle)
    x2_sc = L - x1_sc
    v_single_cycle = np.abs(w * x1_sc / 1000.0 - (rp1 / rp2) * w * x2_sc / 1000.0)

    t_vrel_list = []
    vrel_list = []
    for c in range(n_full_cycles):
        t_c = c * t_cycle + np.linspace(0.0, t_cycle, n_pts_per_cycle, endpoint=False)
        t_vrel_list.append(t_c)
        vrel_list.append(v_single_cycle)

    if remainder_time > 0:
        n_rem = max(10, int(round(n_pts_per_cycle * remainder_time / t_cycle)))
        th_rem = np.linspace(0.0, remainder_time * w, n_rem)
        phi_rem = th_rem + alpha10
        x1_rem = rb1 * (thetab1 + phi_rem)
        x2_rem = L - x1_rem
        v_rem = np.abs(w * x1_rem / 1000.0 - (rp1 / rp2) * w * x2_rem / 1000.0)
        t_rem = n_full_cycles * t_cycle + th_rem / w
        t_vrel_list.append(t_rem)
        vrel_list.append(v_rem)

    t_vrel = np.concatenate(t_vrel_list)
    vrel_clean = np.concatenate(vrel_list)

    axes[1, 1].plot(t_vrel, vrel_clean, label="Δmax = 0 μm", linewidth=1.0, color="#1f77b4")
    axes[1, 1].set_title("Contact relative velocity (Full rotation)")
    axes[1, 1].set_xlabel("t (s)")
    axes[1, 1].set_ylabel("v_rel (m/s)")
    axes[1, 1].set_xlim([0.0, t_period])
    axes[1, 1].grid(True, linestyle=":", alpha=0.6)
    axes[1, 1].legend(loc="upper right")

    plt.tight_layout()
    interpolation_path = output_dir / "interpolation_results.png"
    plt.savefig(interpolation_path, dpi=300)
    plt.close()

    # ----------------------------------------------------
    # 3. Meshing dynamics (Full rotation t in [0, T_rot])
    # ----------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    axes[0, 0].plot(t, DTE_base, label="Δmax = 0 μm", linewidth=0.8, color="#1f77b4")
    axes[0, 0].plot(t, DTE_mod, label=f"Δmax = {deltamax_um:.2f} μm", linewidth=0.8, color="#d62728")
    axes[0, 0].set_title("Dynamic transmission error")
    axes[0, 0].set_xlabel("t (s)")
    axes[0, 0].set_ylabel("DTE (μm)")
    axes[0, 0].set_xlim([0.0, t_period])
    axes[0, 0].grid(True, linestyle=":", alpha=0.6)
    axes[0, 0].legend(loc="upper right")

    axes[0, 1].plot(t, a_base, label="Δmax = 0 μm", linewidth=0.8, color="#1f77b4")
    axes[0, 1].plot(t, a_mod, label=f"Δmax = {deltamax_um:.2f} μm", linewidth=0.8, color="#d62728")
    axes[0, 1].set_title("Meshing acceleration")
    axes[0, 1].set_xlabel("t (s)")
    axes[0, 1].set_ylabel("a (m/s²)")
    axes[0, 1].set_xlim([0.0, t_period])
    axes[0, 1].grid(True, linestyle=":", alpha=0.6)
    axes[0, 1].legend(loc="upper right")

    axes[1, 0].plot(t, v_base, label="Δmax = 0 μm", linewidth=0.8, color="#1f77b4")
    axes[1, 0].plot(t, v_mod, label=f"Δmax = {deltamax_um:.2f} μm", linewidth=0.8, color="#d62728")
    axes[1, 0].set_title("Meshing velocity")
    axes[1, 0].set_xlabel("t (s)")
    axes[1, 0].set_ylabel("v (m/s)")
    axes[1, 0].set_xlim([0.0, t_period])
    axes[1, 0].grid(True, linestyle=":", alpha=0.6)
    axes[1, 0].legend(loc="upper right")

    axes[1, 1].plot(t, dmf_base, label="Δmax = 0 μm", linewidth=0.8, color="#1f77b4")
    axes[1, 1].plot(t, dmf_mod, label=f"Δmax = {deltamax_um:.2f} μm", linewidth=0.8, color="#d62728")
    axes[1, 1].set_title("Dynamic meshing force")
    axes[1, 1].set_xlabel("t (s)")
    axes[1, 1].set_ylabel("DMF (N)")
    axes[1, 1].set_xlim([0.0, t_period])
    axes[1, 1].grid(True, linestyle=":", alpha=0.6)
    axes[1, 1].legend(loc="upper right")

    plt.tight_layout()
    meshing_path = output_dir / "meshing_dynamics.png"
    plt.savefig(meshing_path, dpi=300)
    plt.close()

    # ----------------------------------------------------
    # 4. Normalized coordinate comparisons (Steady-state tooth cycle)
    # ----------------------------------------------------
    z1 = base_system.z1
    pd = step / z1
    a_idx = int(np.floor(pd * (z1 - 2 + 0.05)))
    b_idx = int(np.floor(pd * (z1 - 0.4)))
    n_slice = b_idx - a_idx
    tau = np.linspace(-0.5, 0.5, n_slice)

    r1_slice = base_energy["r1"][a_idx:b_idx]

    fig, axes = plt.subplots(3, 2, figsize=(11, 9))

    axes[0, 0].plot(tau, base_energy["k_te_base"][a_idx:b_idx], "-.", label="Δmax = 0 μm", linewidth=1.5, color="#1f77b4")
    axes[0, 0].plot(tau, modified_energy["k_tem"][a_idx:b_idx], label=f"Δmax = {deltamax_um:.2f} μm", linewidth=1.5, color="#d62728")
    axes[0, 0].set_title("k(t) along normalized coordinate")
    axes[0, 0].set_xlabel("Γ")
    axes[0, 0].set_ylabel("k(t) (N/m)")
    axes[0, 0].set_xlim([-0.5, 0.5])
    axes[0, 0].grid(True, linestyle=":", alpha=0.6)
    axes[0, 0].legend(loc="upper right")

    axes[0, 1].plot(tau, te_base_um[a_idx:b_idx], "-.", label="TE", linewidth=1.5, color="#1f77b4")
    axes[0, 1].plot(tau, te_mod_um[a_idx:b_idx], label=f"TE_M ({deltamax_um:.2f} μm)", linewidth=1.5, color="#d62728")
    axes[0, 1].set_title("Static transmission error along normalized coordinate")
    axes[0, 1].set_xlabel("Γ")
    axes[0, 1].set_ylabel("TE (μm)")
    axes[0, 1].set_xlim([-0.5, 0.5])
    axes[0, 1].grid(True, linestyle=":", alpha=0.6)
    axes[0, 1].legend(loc="upper right")

    axes[1, 0].plot(tau, DTE_base[a_idx:b_idx], "-.", label="Δmax = 0 μm", linewidth=1.5, color="#1f77b4")
    axes[1, 0].plot(tau, DTE_mod[a_idx:b_idx], label=f"Δmax = {deltamax_um:.2f} μm", linewidth=1.5, color="#d62728")
    axes[1, 0].set_title("DTE along normalized coordinate")
    axes[1, 0].set_xlabel("Γ")
    axes[1, 0].set_ylabel("DTE (μm)")
    axes[1, 0].set_xlim([-0.5, 0.5])
    axes[1, 0].grid(True, linestyle=":", alpha=0.6)
    axes[1, 0].legend(loc="upper right")

    axes[1, 1].plot(tau, a_base[a_idx:b_idx], "-.", label="Δmax = 0 μm", linewidth=1.5, color="#1f77b4")
    axes[1, 1].plot(tau, a_mod[a_idx:b_idx], label=f"Δmax = {deltamax_um:.2f} μm", linewidth=1.5, color="#d62728")
    axes[1, 1].set_title("Acceleration along normalized coordinate")
    axes[1, 1].set_xlabel("Γ")
    axes[1, 1].set_ylabel("a (m/s²)")
    axes[1, 1].set_xlim([-0.5, 0.5])
    axes[1, 1].grid(True, linestyle=":", alpha=0.6)
    axes[1, 1].legend(loc="upper right")

    # Subplot 5: Tooth dynamic force with load sharing (dual y-axis)
    ax_left = axes[2, 0]
    ax_right = ax_left.twinx()
    f_tooth_base = dmf_base[a_idx:b_idx] * r1_slice
    f_tooth_mod = dmf_mod[a_idx:b_idx] * r1_slice
    line1 = ax_left.plot(tau, f_tooth_base, "--", label="TVMF (Δmax = 0)", linewidth=1.5, color="#1f77b4")
    line2 = ax_left.plot(tau, f_tooth_mod, label=f"TVMF (Δmax = {deltamax_um:.2f})", linewidth=1.5, color="#d62728")
    line3 = ax_right.plot(tau, r1_slice, ":", label="R₁", linewidth=1.5, color="black")
    ax_left.set_title("Tooth dynamic force with load sharing")
    ax_left.set_xlabel("Γ")
    ax_left.set_ylabel("Tooth Force (N)")
    ax_right.set_ylabel("R₁")
    ax_left.set_xlim([-0.5, 0.5])

    # Automatic dual-axis proportional scaling:
    # Since optimized TVMF ≈ F (nominal static transmitted force), the individual tooth force is
    # F_tooth(Γ) = TVMF * R_1(Γ) ≈ F * R_1(Γ).
    # Setting y_left = F * y_right ensures that the plateau (R_1 = 1.0) and the slope of R_1
    # align exactly with the modified tooth dynamic force for all operating conditions.
    peak_force = max(float(np.max(f_tooth_base)), float(np.max(f_tooth_mod)))
    r_max = max(1.2, float(np.ceil((peak_force / F) / 0.2) * 0.2))
    r_ticks = np.arange(0.0, r_max + 1e-5, 0.2)
    ax_right.set_ylim([0.0, r_max])
    ax_left.set_ylim([0.0, r_max * F])
    ax_right.set_yticks(r_ticks)
    ax_left.set_yticks(r_ticks * F)
    ax_left.yaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: f"{int(round(val)):d}"))

    ax_left.grid(True, linestyle=":", alpha=0.6)
    lines = line1 + line2 + line3
    labels = [l.get_label() for l in lines]
    ax_left.legend(lines, labels, loc="upper right", fontsize=8)

    # Subplot 6: Total dynamic meshing force (TVMF)
    axes[2, 1].plot(tau, dmf_base[a_idx:b_idx], "-.", label="Δmax = 0 μm", linewidth=1.5, color="#1f77b4")
    axes[2, 1].plot(tau, dmf_mod[a_idx:b_idx], label=f"Δmax = {deltamax_um:.2f} μm", linewidth=1.5, color="#d62728")
    axes[2, 1].set_title("Total dynamic meshing force (TVMF)")
    axes[2, 1].set_xlabel("Γ")
    axes[2, 1].set_ylabel("TVMF (N)")
    axes[2, 1].set_xlim([-0.5, 0.5])
    axes[2, 1].grid(True, linestyle=":", alpha=0.6)
    axes[2, 1].legend(loc="upper right")

    plt.tight_layout()
    normalized_path = output_dir / "normalized_comparisons.png"
    plt.savefig(normalized_path, dpi=300)
    plt.close()

    # ----------------------------------------------------
    # 5. Frequency analysis (Steady-state acceleration FFT)
    # ----------------------------------------------------
    rpm = base_system.rpm
    z1 = base_system.z1
    f_mesh = z1 * (rpm / 60.0)
    t_mesh = 1.0 / f_mesh

    # Use a steady-state window with an exact integer number of tooth meshing cycles
    # to eliminate spectral leakage (picket-fence effect) around f_mesh.
    n_teeth_fft = 12
    t_win = n_teeth_fft * t_mesh
    mask_ss = t >= (t[-1] - t_win)
    t_ss = t[mask_ss]
    a_base_ss = a_base[mask_ss]
    a_mod_ss = a_mod[mask_ss]

    n_fft_pts = 2048
    t_uniform = np.linspace(t_ss[0], t_ss[-1], n_fft_pts, endpoint=False)
    a_base_interp = np.interp(t_uniform, t_ss, a_base_ss)
    a_mod_interp = np.interp(t_uniform, t_ss, a_mod_ss)
    dt_uniform = t_uniform[1] - t_uniform[0]

    freqs = np.fft.rfftfreq(n_fft_pts, d=dt_uniform)
    spec_base = 2.0 * np.abs(np.fft.rfft(a_base_interp)) / n_fft_pts
    spec_mod = 2.0 * np.abs(np.fft.rfft(a_mod_interp)) / n_fft_pts

    plt.figure(figsize=(9, 6))
    plt.plot(freqs, spec_base, label="Δmax = 0 μm", linewidth=1.5, color="#1f77b4")
    plt.plot(freqs, spec_mod, label=f"Δmax = {deltamax_um:.2f} μm", linewidth=1.5, color="#d62728")
    plt.axvline(f_mesh, color="gray", linestyle="--", linewidth=1.0, alpha=0.7, label=f"$f_{{mesh}} = {f_mesh:.0f}$ Hz")

    plt.title("Meshing acceleration frequency spectrum (Steady-state)")
    plt.xlabel("f (Hz)")
    plt.ylabel("Acceleration Magnitude (m/s²)")
    plt.xlim([0.0, 1500.0])

    # Scale y-axis based on the visible frequency range [0, 1500 Hz]
    visible_mask = (freqs >= 0.0) & (freqs <= 1500.0)
    y_max = max(float(np.max(spec_base[visible_mask])), float(np.max(spec_mod[visible_mask]))) * 1.15
    plt.ylim([0.0, y_max])

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
