import math

import pytest

from tpm.geometry import compute_geometry
from tpm.mass import compute_mass_properties


def test_mass_properties_are_positive():
    geom = compute_geometry(z1=27, z2=35, module=3.0, alpha_deg=20.0)
    props = compute_mass_properties(geom=geom, b_mm=25.0, density=7850.0, deltamax_um=40.0, correct_centroid_typo=False)
    props_corrected = compute_mass_properties(geom=geom, b_mm=25.0, density=7850.0, deltamax_um=40.0, correct_centroid_typo=True)

    assert props.me > 0.0
    assert props.memod > 0.0
    assert props_corrected.memod > 0.0


def test_mass_properties_match_matlab_geometry_model():
    geom = compute_geometry(z1=27, z2=35, module=3.0, alpha_deg=20.0)
    props_legacy = compute_mass_properties(geom=geom, b_mm=25.0, density=7850.0, deltamax_um=40.0, correct_centroid_typo=False)
    props_corrected = compute_mass_properties(geom=geom, b_mm=25.0, density=7850.0, deltamax_um=40.0, correct_centroid_typo=True)

    z1 = int(geom["z1"])
    z2 = int(geom["z2"])
    p = 7850.0
    b = 25.0
    deltamax = 40.0
    rb1 = geom["rb1"] / 1000.0
    rb2 = geom["rb2"] / 1000.0
    ri1_mm = geom["ri1"]
    ri2_mm = geom["ri2"]
    rf1 = geom["rf1"]
    rf2 = geom["rf2"]
    hr1 = geom["hr1"]
    hr2 = geom["hr2"]
    hi1 = geom["hi1"]
    hi2 = geom["hi2"]
    h1 = geom["h1"]
    h2 = geom["h2"]
    Sa1 = geom["Sa1"]
    Sa2 = geom["Sa2"]
    Sd1 = geom["Sd1"]
    Sd2 = geom["Sd2"]

    ycg1 = ((Sd1 * hr1**2) / 2.0 + (Sd1 * (hi1 - hr1) / 2.0) * ((hi1 - hr1) / 3.0 + hr1) - (Sa1 * (hi1 - h1) / 2.0) * ((hi1 - h1) / 3.0 + h1)) / (Sd1 * hr1 + (Sd1 * (hi1 - hr1) / 2.0) - (Sa1 * (hi1 - h1) / 2.0))
    ycg2 = ((Sd2 * hr2**2) / 2.0 + (Sd2 * (hi2 - hr2) / 2.0) * ((hi2 - hr2) / 3.0 + hr2) - (Sa2 * (hi2 - h2) / 2.0) * ((hi2 - h2) / 3.0 + h2)) / (Sd2 * hr2 + (Sd2 * (hi2 - hr2) / 2.0) - (Sa2 * (hi2 - h2) / 2.0))
    ycg1mod = ((Sd1 * hr1**2) / 2.0 + (Sd1 * (hi1 - hr1) / 2.0) * ((hi1 - hr1) / 3.0 + hr1) - (Sa1 * (hi1 - h1) / 2.0) * ((hi1 - h1) / 3.0 + h1) - (deltamax / 1000.0) * (h1 - hr1) * (2.0 * (h1 - hr1) / 3.0 + hr1)) / (Sd1 * hr1 + (Sd1 * (hi1 - hr1) / 2.0) - (Sa1 * (hi1 - h1) / 2.0) - (deltamax / 1000.0) * (h1 - hr1))
    ycg2mod = ((Sd2 * hr2**2) / 2.0 + (Sd2 * (hi2 - hr2) / 2.0) * ((hi2 - hr2) / 3.0 + hr2) - (Sa2 * (hi2 - h2) / 2.0) * ((hi2 - h2) / 3.0 + h2) - (deltamax / 1000.0) * (h2 - hr2) * (2.0 * (h2 - hr2) / 3.0 + hr2)) / (Sd2 * hr2 + (Sd2 * (hi2 - hr2) / 2.0) - (Sa2 * (hi2 - h2) / 2.0) - (deltamax / 1000.0) * (h1 - hr1))

    A1 = Sd1 * hr1 + Sd1 * (hi1 - hr1) / 2.0 - Sa1 * (hi1 - h1) / 2.0
    A2 = Sd2 * hr2 + Sd2 * (hi2 - hr2) / 2.0 - Sa2 * (hi2 - h2) / 2.0
    Ar1 = (deltamax / 1000.0) * (hi1 - hr1)
    Ar2 = (deltamax / 1000.0) * (hi2 - hr2)

    V1 = A1 * b / 1e9
    V2 = A2 * b / 1e9
    Vr1 = Ar1 * b / 1e9
    Vr2 = Ar2 * b / 1e9

    mD1 = V1 * p
    mD2 = V2 * p
    mDr1 = Vr1 * p
    mDr2 = Vr2 * p

    mDt1 = mD1 * z1
    mDt2 = mD2 * z2
    mDt1mod = (mD1 - mDr1) * z1
    mDt2mod = (mD2 - mDr2) * z2

    m1c = p * math.pi * (b / 1000.0) * (ri1_mm / 1000.0) ** 2
    m2c = p * math.pi * (b / 1000.0) * (ri2_mm / 1000.0) ** 2
    m1cr = p * math.pi * (b / 1000.0) * (rf1 / 1000.0) ** 2
    m2cr = p * math.pi * (b / 1000.0) * (rf2 / 1000.0) ** 2

    mt1 = mDt1 + m1c - m1cr
    mt2 = mDt2 + m2c - m2cr
    mt1mod = mDt1mod + m1c - m1cr
    mt2mod = mDt2mod + m2c - m2cr

    Id1 = mDt1 * (ycg1 / 1000.0 + ri1_mm / 1000.0) ** 2
    Id2 = mDt2 * (ycg2 / 1000.0 + ri2_mm / 1000.0) ** 2
    Id1mod = mDt1mod * (ycg1mod / 1000.0 + ri1_mm / 1000.0) ** 2
    Id2mod = mDt2mod * (ycg2mod / 1000.0 + ri2_mm / 1000.0) ** 2

    Ic1 = m1c * (ri1_mm / 1000.0) ** 2 / 2.0 - m1cr * (rf1 / 1000.0) ** 2 / 2.0
    Ic2 = m2c * (ri2_mm / 1000.0) ** 2 / 2.0 - m2cr * (rf2 / 1000.0) ** 2 / 2.0

    J1 = Id1 + Ic1
    J2 = Id2 + Ic2
    J1mod = Id1mod + Ic1
    J2mod = Id2mod + Ic2

    expected_me = J1 * J2 / (J2 * (rb1**2) + J1 * (rb2**2))
    expected_memod = J1mod * J2mod / (J2mod * (rb1**2) + J1mod * (rb2**2))

    assert props_legacy.me == pytest.approx(expected_me, rel=1e-6)
    assert props_legacy.memod == pytest.approx(expected_memod, rel=1e-6)
    assert props_corrected.memod == pytest.approx(expected_memod, rel=1e-4)
