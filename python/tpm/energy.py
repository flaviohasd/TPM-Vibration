from __future__ import annotations

import math
from typing import Dict, Optional, Tuple, Union

import numpy as np
from scipy.integrate import quad


def _integral_compliance(kind: str, alpha1: float, alpha2: float, *args) -> float:
    """Numerically evaluate the MATLAB-style compliance integrals used in the stiffness model."""
    if len(args) == 4:
        _, E, b_mm, v = args
    elif len(args) == 3:
        E, b_mm, v = args
    else:
        raise ValueError(f"Invalid number of arguments: {len(args)}")

    b_m = b_mm / 1000.0
    if kind == "ikbi":
        def integrand(x: float) -> float:
            num = 3.0 * ((1.0 + math.cos(alpha1) * ((alpha2 - x) * math.sin(x) - math.cos(x))) ** 2) * (alpha2 - x) * math.cos(x)
            denom = 2.0 * E * b_m * (math.sin(x) + (alpha2 - x) * math.cos(x)) ** 3
            return num / denom
    elif kind == "iksi":
        def integrand(x: float) -> float:
            num = 1.2 * (1.0 + v) * (alpha2 - x) * math.cos(x) * (math.cos(alpha1) ** 2)
            denom = E * b_m * (math.sin(x) + (alpha2 - x) * math.cos(x))
            return num / denom
    elif kind == "ikai":
        def integrand(x: float) -> float:
            num = (alpha2 - x) * math.cos(x) * (math.sin(alpha1) ** 2)
            denom = 2.0 * E * b_m * (math.sin(x) + (alpha2 - x) * math.cos(x))
            return num / denom
    else:
        raise ValueError(f"Unsupported compliance kind: {kind}")

    return float(quad(integrand, -alpha1, alpha2, limit=100)[0])


def _calc_compliances_vectorized(
    alpha1_vec: np.ndarray,
    alpha2_scalar: float,
    E: float,
    b_mm: float,
    v: float,
    n_points: int = 25,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Machine-precision Gauss-Legendre quadrature for compliance integrals across all angles."""
    b_m = b_mm / 1000.0
    gl_nodes, gl_weights = np.polynomial.legendre.leggauss(n_points)
    gl_nodes = gl_nodes[:, None]
    gl_weights = gl_weights[:, None]

    alpha1_row = alpha1_vec[None, :]
    mid = (alpha2_scalar - alpha1_row) / 2.0
    half_len = (alpha2_scalar + alpha1_row) / 2.0
    x = mid + half_len * gl_nodes
    w_mat = half_len * gl_weights

    # Ikbi
    term1 = 1.0 + np.cos(alpha1_row) * ((alpha2_scalar - x) * np.sin(x) - np.cos(x))
    num_b = 3.0 * (term1**2) * (alpha2_scalar - x) * np.cos(x)
    denom_b = 2.0 * E * b_m * (np.sin(x) + (alpha2_scalar - x) * np.cos(x)) ** 3
    ikbi = np.sum(w_mat * (num_b / denom_b), axis=0)

    # Iksi
    num_s = 1.2 * (1.0 + v) * (alpha2_scalar - x) * np.cos(x) * (np.cos(alpha1_row) ** 2)
    denom_s = E * b_m * (np.sin(x) + (alpha2_scalar - x) * np.cos(x))
    iksi = np.sum(w_mat * (num_s / denom_s), axis=0)

    # Ikai
    num_a = (alpha2_scalar - x) * np.cos(x) * (np.sin(alpha1_row) ** 2)
    denom_a = 2.0 * E * b_m * (np.sin(x) + (alpha2_scalar - x) * np.cos(x))
    ikai = np.sum(w_mat * (num_a / denom_a), axis=0)

    return ikbi, iksi, ikai


def precompute_energy_state(
    system,
    step: Optional[int] = None,
    points_per_tooth: Optional[int] = None,
    resolution: Optional[str] = None,
    friction_model: str = "constant",
    fixed_mu: Optional[float] = 0.1211,
) -> Dict[str, np.ndarray]:
    """Precompute the exact MATLAB energy-based TVMS, damping, and normalized profile modifications."""
    z1 = int(system.z1)
    z2 = int(system.z2)

    if points_per_tooth is not None:
        step = int(z1 * points_per_tooth)
    elif resolution is not None:
        preset_map = {"fast": 50, "standard": 100, "fine": 200, "ultra": 400, "publication": 400}
        step = int(z1 * preset_map.get(resolution, 400))
    elif step is not None:
        step = int(step)
    else:
        step = 400
    alpha = float(system.geometry["alpha"])
    pb = float(system.geometry["pb"])
    rb1 = float(system.geometry["rb1"])
    rb2 = float(system.geometry["rb2"])
    rp1 = float(system.geometry["rp1"])
    rp2 = float(system.geometry["rp2"])
    di1 = float(system.geometry["di1"])
    di2 = float(system.geometry["di2"])

    thetab1 = float(system.mesh["thetab1"])
    thetab2 = float(system.mesh["thetab2"])
    alpha10 = float(system.mesh["alpha10"])
    alpha20 = float(system.mesh["alpha20"])
    thetad = float(system.mesh["thetad"])
    thetas = float(system.mesh["thetas"])
    L = float(system.mesh["L"])
    md1 = float(system.mesh["md1"])
    md2 = float(system.mesh["md2"])
    double_len = float(system.mesh["double"])
    single_len = float(system.mesh["single"])

    E = float(system.elasticity_modulus)
    b_mm = float(system.b_mm)
    v = float(system.poisson_ratio)
    w = float(system.operation["w"])
    force = float(system.operation["F"])
    me = float(system.mass_properties["me"])
    memod = float(system.mass_properties["memod"])
    relief_amplitude_um = float(getattr(system, "deltamax_um", 0.0))

    rf1 = di1 * 0.2 / 2.0
    rf2 = di2 * 0.2 / 2.0
    hf1 = di1 / (2.0 * rf1)
    hf2 = di2 / (2.0 * rf2)

    # Foundation stiffness coefficients (Sainsot / Cai model)
    X = np.array(
        [
            [-5.574e-5, -1.9986e-3, -2.3015e-4, 4.7702e-3, 0.0271, 6.8045],
            [60.111e-5, 28.1e-3, -83.431e-4, -9.9256e-3, 0.1624, 0.9086],
            [-50.952e-5, 185.5e-3, 0.0538e-4, 53.3e-3, 0.2895, 0.9236],
            [-6.2042e-5, 9.0889e-3, -4.0964e-4, 7.8297e-3, -0.1472, 0.6904],
        ]
    )
    XL1 = float(np.dot(X[0], [1.0 / thetab1**2, hf1**2, hf1 / thetab1, 1.0 / thetab1, hf1, 1.0]))
    XL2 = float(np.dot(X[0], [1.0 / thetab2**2, hf2**2, hf2 / thetab2, 1.0 / thetab2, hf2, 1.0]))
    XM1 = float(np.dot(X[1], [1.0 / thetab1**2, hf1**2, hf1 / thetab1, 1.0 / thetab1, hf1, 1.0]))
    XM2 = float(np.dot(X[1], [1.0 / thetab2**2, hf2**2, hf2 / thetab2, 1.0 / thetab2, hf2, 1.0]))
    XP1 = float(np.dot(X[2], [1.0 / thetab1**2, hf1**2, hf1 / thetab1, 1.0 / thetab1, hf1, 1.0]))
    XP2 = float(np.dot(X[2], [1.0 / thetab2**2, hf2**2, hf2 / thetab2, 1.0 / thetab2, hf2, 1.0]))
    XQ1 = float(np.dot(X[3], [1.0 / thetab1**2, hf1**2, hf1 / thetab1, 1.0 / thetab1, hf1, 1.0]))
    XQ2 = float(np.dot(X[3], [1.0 / thetab2**2, hf2**2, hf2 / thetab2, 1.0 / thetab2, hf2, 1.0]))

    # Hertzian contact stiffness
    kh = (math.pi * E * (b_mm / 1000.0)) / (4.0 * (1.0 - v**2))

    t_vec = np.linspace(0.0, (2.0 * math.pi) / w, step)
    theta1_vec = w * t_vec
    theta1vg = np.rad2deg(theta1_vec)

    theta1corr = np.floor(theta1_vec / (2.0 * np.pi / z1)) * (2.0 * np.pi / z1)
    tdouble_0 = theta1corr
    tdouble_1 = theta1corr + thetad
    tsingle_0 = theta1corr + thetad
    tsingle_1 = theta1corr + (thetad + thetas)

    theta = theta1_vec - np.floor(theta1_vec / (thetad + thetas)) * (thetad + thetas)
    thetadsd = theta1_vec - np.floor(theta1_vec / (2.0 * thetad + thetas)) * (2.0 * thetad + thetas)

    phi1 = theta + alpha10
    phi1dsd = thetadsd + alpha10
    phi1d = phi1 + 2.0 * np.pi / z1

    h1_ct = rb1 * ((thetab1 + phi1) * np.cos(phi1) - np.sin(phi1))
    h1d_ct = rb1 * ((thetab1 + phi1d) * np.cos(phi1d) - np.sin(phi1d))
    l1_ct = rb1 * ((thetab1 + phi1) * np.sin(phi1) + np.cos(phi1) - math.cos(thetab1))
    l1d_ct = rb1 * ((thetab1 + phi1d) * np.sin(phi1d) + np.cos(phi1d) - math.cos(thetab1))
    x1 = rb1 * (thetab1 + phi1)
    x1d = rb1 * (thetab1 + phi1d)

    x2 = L - x1
    x2d = L - x1d
    phi2 = alpha20 - (z1 / z2) * theta
    phi2dsd = alpha20 - (z1 / z2) * thetadsd
    phi2d = phi2 - 2.0 * np.pi / z2

    l2_ct = rb2 * ((thetab2 + phi2) * np.sin(phi2) + np.cos(phi2) - math.cos(thetab2))
    l2d_ct = rb2 * ((thetab2 + phi2d) * np.sin(phi2d) + np.cos(phi2d) - math.cos(thetab2))
    h2_ct = rb2 * ((thetab2 + phi2) * np.cos(phi2) - np.sin(phi2))
    h2d_ct = rb2 * ((thetab2 + phi2d) * np.cos(phi2d) - np.sin(phi2d))

    h1dsd = rb1 * ((thetab1 + phi1dsd) * np.cos(phi1dsd) - np.sin(phi1dsd))
    l1dsd = rb1 * ((thetab1 + phi1dsd) * np.sin(phi1dsd) + np.cos(phi1dsd) - math.cos(thetab1))
    x1dsd = rb1 * (thetab1 + phi1dsd)
    x2dsd = L - x1dsd

    alphax1 = phi1 + np.arctan(h1_ct / (rb1 * math.cos(thetab1) + l1_ct))
    alphax1dsd = phi1dsd + np.arctan(h1dsd / (rb1 * math.cos(thetab1) + l1dsd))

    rx1 = rb1 / np.cos(alphax1)
    rx1dsd = rb1 / np.cos(alphax1dsd)
    dx = rx1 * np.sin(alphax1) - md1
    alphax1d = np.arccos(np.sqrt((md1 + dx)**2 + rb1**2) * np.cos(alphax1) / np.sqrt((md1 + pb + dx)**2 + rb1**2))
    rx1d = (rx1 * np.sin(alphax1) + pb) / np.sin(alphax1d)

    rx2 = np.sqrt(rb2**2 + (L - np.sqrt(rx1**2 - rb1**2))**2)
    alphax2 = np.arccos(rb2 / rx2)
    rx2d = np.sqrt(rb2**2 + (L - np.sqrt(rx1d**2 - rb1**2))**2)

    Tau = (x1dsd / rb1) / math.tan(alpha) - 1.0

    dx1 = x1 - md1
    dx2 = x2d - md2

    uf1 = l1_ct - h1_ct * np.tan(phi1)
    Sf1 = thetab1 * di1
    uf2 = l2_ct - h2_ct * np.tan(phi2)
    Sf2 = thetab2 * di2
    uf1d = l1d_ct - h1d_ct * np.tan(phi1d)
    uf2d = l2d_ct - h2d_ct * np.tan(phi2d)

    Kf1 = (np.cos(phi1)**2 / (E * b_mm / 1000.0)) * (XL1 * (uf1 / Sf1)**2 + XM1 * (uf1 / Sf1) + XP1 * (1.0 + XQ1 * np.tan(phi1)**2))
    Kf1d = (np.cos(phi1d)**2 / (E * b_mm / 1000.0)) * (XL1 * (uf1d / Sf1)**2 + XM1 * (uf1d / Sf1) + XP1 * (1.0 + XQ1 * np.tan(phi1d)**2))
    Kf2 = (np.cos(phi2)**2 / (E * b_mm / 1000.0)) * (XL2 * (uf2 / Sf2)**2 + XM2 * (uf2 / Sf2) + XP2 * (1.0 + XQ2 * np.tan(phi2)**2))
    Kf2d = (np.cos(phi2d)**2 / (E * b_mm / 1000.0)) * (XL2 * (uf2d / Sf2)**2 + XM2 * (uf2d / Sf2) + XP2 * (1.0 + XQ2 * np.tan(phi2d)**2))

    Ikbi1, Iksi1, Ikai1 = _calc_compliances_vectorized(phi1, thetab1, E, b_mm, v)
    Ikbi2, Iksi2, Ikai2 = _calc_compliances_vectorized(phi2, thetab2, E, b_mm, v)
    Ikbi1d, Iksi1d, Ikai1d = _calc_compliances_vectorized(phi1d, thetab1, E, b_mm, v)
    Ikbi2d, Iksi2d, Ikai2d = _calc_compliances_vectorized(phi2d, thetab2, E, b_mm, v)

    K1 = Ikbi1 + Ikai1 + Iksi1 + Kf1
    K2 = Ikbi2 + Ikai2 + Iksi2 + Kf2
    K1d = Ikbi1d + Ikai1d + Iksi1d + Kf1d
    K2d = Ikbi2d + Ikai2d + Iksi2d + Kf2d

    Kte = np.zeros(step)
    Ktem = np.zeros(step)
    R1 = np.zeros(step)
    delta = np.zeros((step, 2))
    delta_normalized = np.zeros((step, 2))
    K12 = np.zeros((step, 2))
    TE = np.zeros(step)
    TEm = np.zeros(step)
    dTE = np.zeros(step)

    mem = 0
    for i in range(step):
        th1 = theta1_vec[i]
        if th1 >= tdouble_0[i] and th1 < tdouble_1[i]:
            kte1 = float(1.0 / (1.0 / kh + K1[i] + K2[i]))
            kte2 = float(1.0 / (1.0 / kh + K1d[i] + K2d[i]))
            Kte[i] = kte1 + kte2
            K12[i, 0] = kte1
            K12[i, 1] = kte2

            norm1 = (double_len - dx1[i]) / (1e6 * double_len)
            norm2 = (double_len - dx2[i]) / (1e6 * double_len)
            delta_normalized[i, 0] = norm1
            delta_normalized[i, 1] = norm2

            deltax1 = (relief_amplitude_um / 1000.0) * (double_len - dx1[i]) / double_len
            deltax2 = (relief_amplitude_um / 1000.0) * (double_len - dx2[i]) / double_len
            delta[i, 0] = deltax1 / 1000.0
            delta[i, 1] = deltax2 / 1000.0

            Ktem[i] = force * Kte[i] / (force + kte1 * deltax1 / 1000.0 + kte2 * deltax2 / 1000.0)
            R1[i] = kte1 / Kte[i]
        elif th1 >= tsingle_0[i] and th1 <= tsingle_1[i]:
            kte1 = float(1.0 / (1.0 / kh + K1[i] + K2[i]))
            kte2 = 0.0
            Kte[i] = kte1 + kte2
            K12[i, 0] = kte1
            K12[i, 1] = kte2

            delta_normalized[i, 0] = 0.0
            delta_normalized[i, 1] = 0.0
            delta[i, 0] = 0.0
            delta[i, 1] = 0.0

            Ktem[i] = Kte[i]
            R1[i] = kte1 / Kte[i]
            mem = i

        if th1 > (thetad + thetas) and th1 <= (2.0 * thetad + thetas):
            R1[i] = 1.0 - R1[i - mem]

        TE[i] = force / (K12[i, 0] + K12[i, 1])
        dTE1 = (K12[i, 0] * delta[i, 0]) / Kte[i]
        dTE2 = (K12[i, 1] * delta[i, 1]) / Kte[i]
        dTE[i] = dTE1 + dTE2
        TEm[i] = TE[i] + dTE[i]

    # Friction and damping calculation
    Fn = force * np.cos(alphax1dsd)
    vrel = np.abs(w * x1dsd / 1000.0 - (rp1 / rp2) * w * x2dsd / 1000.0)
    Ra = 0.63 / 1e6
    t0 = 60.0
    v0 = 38.5 + (t0 - 100.0) * (38.5 - 320.0) / (100.0 - 40.0)

    if friction_model == "constant":
        if fixed_mu is not None:
            mu = np.full(step, float(fixed_mu))
        else:
            # Calculate mean of valid EHL curve outside the pitch point singularity
            valid = vrel > 1e-3
            wt = force * np.sin(alphax1dsd)
            ehl = 0.12 * ((np.abs(wt) * Ra) / (v0 * np.maximum(vrel, 1e-4) * (rx1dsd / 1000.0))) ** 0.25
            mean_mu = float(np.mean(ehl[valid])) if np.any(valid) else 0.1211
            mu = np.full(step, mean_mu)
    elif friction_model == "ehl_variable":
        wt = force * np.sin(alphax1dsd)
        vrel_reg = np.sqrt(vrel**2 + 0.01**2)  # Smooth regularization avoiding division by zero
        mu = np.clip(0.12 * ((np.abs(wt) * Ra) / (v0 * vrel_reg * (rx1dsd / 1000.0))) ** 0.25, 0.01, 0.15)
    else:
        mu = np.full(step, 0.1211)

    Fa = Fn * mu
    B0 = Fn / Kte
    Wa = 4.0 * Fa * B0
    U = (Fn**2) / (2.0 * Kte)
    zeta = Wa / (4.0 * math.pi * U)
    C = 2.0 * zeta * np.sqrt(me * Kte)
    Cm = 2.0 * zeta * np.sqrt(memod * Ktem)

    forcing = force + K12[:, 0] * delta[:, 0] + K12[:, 1] * delta[:, 1]
    active_kte = Ktem if relief_amplitude_um > 0 else Kte

    return {
        "theta": theta1_vec,
        "k_te": active_kte,
        "k_te_base": Kte,
        "k_tem": Ktem,
        "c": C,
        "cm": Cm,
        "mu": mu,
        "vrel": vrel,
        "transmission_error": TEm,
        "delta": delta,
        "delta_normalized": delta_normalized,
        "k12": K12,
        "forcing": forcing,
        "r1": R1,
        "theta1vg": theta1vg,
        "vrelv": vrel,
        "te": TE,
        "dte": dTE,
        "kh": kh,
    }
