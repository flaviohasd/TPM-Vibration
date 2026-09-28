from __future__ import annotations

import math
from typing import Dict


def compute_geometry(z1: int, z2: int, module: float, alpha_deg: float) -> Dict[str, float]:
    """Compute basic spur gear geometry from standard parameters."""
    if z1 <= 0 or z2 <= 0:
        raise ValueError("Number of teeth must be positive")

    alpha = math.radians(alpha_deg)
    if math.asin(math.sqrt(2.0 / z1)) > alpha:
        raise ValueError("Pressure angle is below the critical angle for the selected teeth")

    e = 0.167 * module

    dp1 = module * z1
    de1 = module * (z1 + 2.0)
    db1 = dp1 * math.cos(alpha)
    di1 = module * (z1 - 2.334)

    rd1 = di1 / 2.0
    rde1 = di1 / 2.0 + e
    ra1 = de1 / 2.0
    rb1 = db1 / 2.0
    ri1 = di1 / 2.0
    rp1 = dp1 / 2.0

    rf1 = di1 * 0.2 / 2.0
    hf1 = di1 / (2.0 * rf1)
    sc1 = module * z1 * math.sin((math.pi / 2.0) / z1)
    Sa1 = de1 * (sc1 / dp1 + math.tan(alpha) - alpha - math.tan(math.acos(dp1 * math.cos(alpha) / de1)) + math.acos(dp1 * math.cos(alpha) / de1))
    Sd1 = 2.0 * rb1 * math.sin((math.pi + 4.0 * 0.0 * math.tan(alpha)) / (2.0 * z1) + math.tan(alpha) - alpha)
    hr1 = math.sqrt(rb1**2 - (Sd1 / 2.0) ** 2) - math.sqrt(rd1**2 - (Sd1 / 2.0) ** 2)
    h1 = math.sqrt(ra1**2 - (Sa1 / 2.0) ** 2) - math.sqrt(rd1**2 - (Sd1 / 2.0) ** 2)
    hi1 = (h1 * Sd1 - hr1 * Sa1) / (Sd1 - Sa1)

    dp2 = module * z2
    de2 = module * (z2 + 2.0)
    db2 = dp2 * math.cos(alpha)
    di2 = module * (z2 - 2.334)

    rd2 = di2 / 2.0
    rde2 = di2 / 2.0 + e
    ra2 = de2 / 2.0
    rb2 = db2 / 2.0
    ri2 = di2 / 2.0
    rp2 = dp2 / 2.0

    rf2 = di2 * 0.2 / 2.0
    hf2 = di2 / (2.0 * rf2)
    sc2 = module * z2 * math.sin((math.pi / 2.0) / z2)
    Sa2 = de2 * (sc2 / dp2 + math.tan(alpha) - alpha - math.tan(math.acos(dp2 * math.cos(alpha) / de2)) + math.acos(dp2 * math.cos(alpha) / de2))
    Sd2 = 2.0 * rb2 * math.sin((math.pi + 4.0 * 0.0 * math.tan(alpha)) / (2.0 * z2) + math.tan(alpha) - alpha)
    hr2 = math.sqrt(rde2**2 - (Sd2 / 2.0) ** 2) - math.sqrt(rd2**2 - (Sd2 / 2.0) ** 2)
    h2 = math.sqrt(ra2**2 - (Sa2 / 2.0) ** 2) - math.sqrt(rd2**2 - (Sd2 / 2.0) ** 2)
    hi2 = (h2 * Sd2 - hr2 * Sa2) / (Sd2 - Sa2)

    pb = math.pi * db1 / z1
    L = math.sqrt((dp1 / 2.0) ** 2 - rb1**2) + math.sqrt((dp2 / 2.0) ** 2 - rb2**2)
    contact_ratio = (
        math.sqrt(ra1**2 - rb1**2) + math.sqrt(ra2**2 - rb2**2) - (dp1 + dp2) * math.sin(alpha) / 2.0
    ) / pb

    return {
        "z1": float(z1),
        "z2": float(z2),
        "alpha": float(alpha),
        "dp1": float(dp1),
        "de1": float(de1),
        "db1": float(db1),
        "di1": float(di1),
        "rd1": float(rd1),
        "ra1": float(ra1),
        "rb1": float(rb1),
        "ri1": float(ri1),
        "rp1": float(rp1),
        "dp2": float(dp2),
        "de2": float(de2),
        "db2": float(db2),
        "di2": float(di2),
        "rd2": float(rd2),
        "ra2": float(ra2),
        "rb2": float(rb2),
        "ri2": float(ri2),
        "rp2": float(rp2),
        "rf1": float(rf1),
        "rf2": float(rf2),
        "hf1": float(hf1),
        "hf2": float(hf2),
        "Sa1": float(Sa1),
        "Sa2": float(Sa2),
        "Sd1": float(Sd1),
        "Sd2": float(Sd2),
        "hr1": float(hr1),
        "hr2": float(hr2),
        "h1": float(h1),
        "h2": float(h2),
        "hi1": float(hi1),
        "hi2": float(hi2),
        "pb": float(pb),
        "L": float(L),
        "contact_ratio": float(contact_ratio),
        "rde1": float(rde1),
        "rde2": float(rde2),
    }


def compute_mesh(geometry: Dict[str, float]) -> Dict[str, float]:
    """Compute mesh-related quantities from already computed geometry."""
    alpha = geometry["alpha"]
    dp1 = geometry["dp1"]
    dp2 = geometry["dp2"]
    db1 = geometry["db1"]
    rb1 = geometry["rb1"]
    ra1 = geometry["ra1"]
    rb2 = geometry["rb2"]
    ra2 = geometry["ra2"]
    pb = geometry["pb"]
    z1 = int(geometry["z1"])
    z2 = int(geometry["z2"])

    L = math.sqrt((dp1 / 2.0) ** 2 - rb1**2) + math.sqrt((dp2 / 2.0) ** 2 - rb2**2)
    cr = (
        math.sqrt(ra1**2 - rb1**2) + math.sqrt(ra2**2 - rb2**2) - (dp1 + dp2) * math.sin(alpha) / 2.0
    ) / pb

    md1 = math.sqrt(ra1**2 - rb1**2) - cr * pb
    md2 = math.sqrt(ra2**2 - rb2**2) - cr * pb

    thetad = (
        math.tan(math.acos(z1 * math.cos(alpha) / (z1 + 2.0)))
        - 2.0 * math.pi / z1
        - math.tan(
            math.acos(
                z1 * math.cos(alpha)
                / math.sqrt(
                    (z2 + 2.0) ** 2
                    + (z1 + z2) ** 2
                    - 2.0 * (z2 + 2.0) * (z1 + z2) * math.cos(math.acos(z2 * math.cos(alpha) / (z2 + 2.0)) - alpha)
                )
            )
        )
    )
    thetas = 2.0 * math.pi / z1 - thetad
    thetab1 = math.pi / (2.0 * z1) + math.tan(alpha) - alpha
    thetab2 = math.pi / (2.0 * z2) + math.tan(alpha) - alpha

    alpha10 = (
        -math.pi / (2.0 * z1)
        - math.tan(alpha)
        + alpha
        + math.tan(
            math.acos(
                z1 * math.cos(alpha)
                / math.sqrt(
                    (z2 + 2.0) ** 2
                    + (z1 + z2) ** 2
                    - 2.0 * (z2 + 2.0) * (z1 + z2) * math.cos(math.acos(z2 * math.cos(alpha) / (z2 + 2.0)) - alpha)
                )
            )
        )
    )
    alpha20 = (
        math.tan(math.acos(z2 * math.cos(alpha) / (z2 + 2.0)))
        - math.pi / (2.0 * z2)
        - math.tan(alpha)
        + alpha
    )

    double = (cr - 1.0) * pb
    single = pb - double

    return {
        "pb": float(pb),
        "L": float(L),
        "cr": float(cr),
        "PTH": float(cr * pb),
        "double": float(double),
        "single": float(single),
        "md1": float(md1),
        "md2": float(md2),
        "thetad": float(thetad),
        "thetas": float(thetas),
        "alpha10": float(alpha10),
        "alpha20": float(alpha20),
        "thetab1": float(thetab1),
        "thetab2": float(thetab2),
    }
