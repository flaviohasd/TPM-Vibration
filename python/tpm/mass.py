from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict


@dataclass
class MassProperties:
    me: float
    memod: float
    mt1: float = 0.0
    mt2: float = 0.0
    J1: float = 0.0
    J2: float = 0.0
    mt1mod: float = 0.0
    mt2mod: float = 0.0
    J1mod: float = 0.0
    J2mod: float = 0.0


def compute_mass_properties(
    geom: Dict[str, float],
    b_mm: float,
    density: float,
    deltamax_um: float,
    correct_centroid_typo: bool = False,
) -> MassProperties:
    """Compute equivalent mass and modified equivalent mass using the MATLAB centroid/inertia formulation."""
    z1 = int(geom["z1"])
    z2 = int(geom["z2"])
    rb1_mm = geom["rb1"]
    rb2_mm = geom["rb2"]
    ri1_mm = geom["ri1"]
    ri2_mm = geom["ri2"]
    rf1_mm = geom["rf1"]
    rf2_mm = geom["rf2"]
    hr1_mm = geom["hr1"]
    hr2_mm = geom["hr2"]
    hi1_mm = geom["hi1"]
    hi2_mm = geom["hi2"]
    h1_mm = geom["h1"]
    h2_mm = geom["h2"]
    Sa1_mm = geom["Sa1"]
    Sa2_mm = geom["Sa2"]
    Sd1_mm = geom["Sd1"]
    Sd2_mm = geom["Sd2"]

    ycg1 = ((Sd1_mm * hr1_mm**2) / 2.0 + (Sd1_mm * (hi1_mm - hr1_mm) / 2.0) * ((hi1_mm - hr1_mm) / 3.0 + hr1_mm) - (Sa1_mm * (hi1_mm - h1_mm) / 2.0) * ((hi1_mm - h1_mm) / 3.0 + h1_mm)) / (
        Sd1_mm * hr1_mm + (Sd1_mm * (hi1_mm - hr1_mm) / 2.0) - (Sa1_mm * (hi1_mm - h1_mm) / 2.0)
    )
    ycg2 = ((Sd2_mm * hr2_mm**2) / 2.0 + (Sd2_mm * (hi2_mm - hr2_mm) / 2.0) * ((hi2_mm - hr2_mm) / 3.0 + hr2_mm) - (Sa2_mm * (hi2_mm - h2_mm) / 2.0) * ((hi2_mm - h2_mm) / 3.0 + h2_mm)) / (
        Sd2_mm * hr2_mm + (Sd2_mm * (hi2_mm - hr2_mm) / 2.0) - (Sa2_mm * (hi2_mm - h2_mm) / 2.0)
    )
    ycg1mod = (
        (Sd1_mm * hr1_mm**2) / 2.0
        + (Sd1_mm * (hi1_mm - hr1_mm) / 2.0) * ((hi1_mm - hr1_mm) / 3.0 + hr1_mm)
        - (Sa1_mm * (hi1_mm - h1_mm) / 2.0) * ((hi1_mm - h1_mm) / 3.0 + h1_mm)
        - (deltamax_um / 1000.0) * (h1_mm - hr1_mm) * (2.0 * (h1_mm - hr1_mm) / 3.0 + hr1_mm)
    ) / (
        Sd1_mm * hr1_mm
        + (Sd1_mm * (hi1_mm - hr1_mm) / 2.0)
        - (Sa1_mm * (hi1_mm - h1_mm) / 2.0)
        - (deltamax_um / 1000.0) * (h1_mm - hr1_mm)
    )
    h2_term = (h2_mm - hr2_mm) if correct_centroid_typo else (h1_mm - hr1_mm)
    ycg2mod = (
        (Sd2_mm * hr2_mm**2) / 2.0
        + (Sd2_mm * (hi2_mm - hr2_mm) / 2.0) * ((hi2_mm - hr2_mm) / 3.0 + hr2_mm)
        - (Sa2_mm * (hi2_mm - h2_mm) / 2.0) * ((hi2_mm - h2_mm) / 3.0 + h2_mm)
        - (deltamax_um / 1000.0) * (h2_mm - hr2_mm) * (2.0 * (h2_mm - hr2_mm) / 3.0 + hr2_mm)
    ) / (
        Sd2_mm * hr2_mm
        + (Sd2_mm * (hi2_mm - hr2_mm) / 2.0)
        - (Sa2_mm * (hi2_mm - h2_mm) / 2.0)
        - (deltamax_um / 1000.0) * h2_term
    )

    A1 = Sd1_mm * hr1_mm + Sd1_mm * (hi1_mm - hr1_mm) / 2.0 - Sa1_mm * (hi1_mm - h1_mm) / 2.0
    A2 = Sd2_mm * hr2_mm + Sd2_mm * (hi2_mm - hr2_mm) / 2.0 - Sa2_mm * (hi2_mm - h2_mm) / 2.0
    Ar1 = (deltamax_um / 1000.0) * (hi1_mm - hr1_mm)
    Ar2 = (deltamax_um / 1000.0) * (hi2_mm - hr2_mm)

    V1 = A1 * b_mm / 1e9
    V2 = A2 * b_mm / 1e9
    Vr1 = Ar1 * b_mm / 1e9
    Vr2 = Ar2 * b_mm / 1e9

    mD1 = V1 * density
    mD2 = V2 * density
    mDr1 = Vr1 * density
    mDr2 = Vr2 * density

    mDt1 = mD1 * z1
    mDt2 = mD2 * z2
    mDt1mod = (mD1 - mDr1) * z1
    mDt2mod = (mD2 - mDr2) * z2

    m1c = density * math.pi * (b_mm / 1000.0) * (ri1_mm / 1000.0) ** 2
    m2c = density * math.pi * (b_mm / 1000.0) * (ri2_mm / 1000.0) ** 2
    m1cr = density * math.pi * (b_mm / 1000.0) * (rf1_mm / 1000.0) ** 2
    m2cr = density * math.pi * (b_mm / 1000.0) * (rf2_mm / 1000.0) ** 2

    mt1 = mDt1 + m1c - m1cr
    mt2 = mDt2 + m2c - m2cr
    mt1mod = mDt1mod + m1c - m1cr
    mt2mod = mDt2mod + m2c - m2cr

    Id1 = mDt1 * (ycg1 / 1000.0 + ri1_mm / 1000.0) ** 2
    Id2 = mDt2 * (ycg2 / 1000.0 + ri2_mm / 1000.0) ** 2
    Id1mod = mDt1mod * (ycg1mod / 1000.0 + ri1_mm / 1000.0) ** 2
    Id2mod = mDt2mod * (ycg2mod / 1000.0 + ri2_mm / 1000.0) ** 2

    Ic1 = m1c * (ri1_mm / 1000.0) ** 2 / 2.0 - m1cr * (rf1_mm / 1000.0) ** 2 / 2.0
    Ic2 = m2c * (ri2_mm / 1000.0) ** 2 / 2.0 - m2cr * (rf2_mm / 1000.0) ** 2 / 2.0

    J1 = Id1 + Ic1
    J2 = Id2 + Ic2
    J1mod = Id1mod + Ic1
    J2mod = Id2mod + Ic2

    me = J1 * J2 / (J2 * (rb1_mm / 1000.0) ** 2 + J1 * (rb2_mm / 1000.0) ** 2)
    memod = J1mod * J2mod / (J2mod * (rb1_mm / 1000.0) ** 2 + J1mod * (rb2_mm / 1000.0) ** 2)

    return MassProperties(
        me=float(me),
        memod=float(memod),
        mt1=float(mt1),
        mt2=float(mt2),
        J1=float(J1),
        J2=float(J2),
        mt1mod=float(mt1mod),
        mt2mod=float(mt2mod),
        J1mod=float(J1mod),
        J2mod=float(J2mod),
    )
