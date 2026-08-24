"""Tests for the Hurley et al. (2000) main-sequence module (ZAMS-to-TMS
only -- see src/realta/stellar/main_sequence.py's module docstring for
scope). See docs/provenance.md and
docs/science/rlof-ce-classifier-proposal.md for the transcription
discipline these pin: a real coefficient bug (a40's gamma exponent and
a41's alpha exponent, both off by a factor of 10) was caught during
this session precisely because ms_radius() collapsed to near-zero for
M=5/20 Msun instead of growing smoothly -- these tests exist so a
regression of that kind fails loudly again.
"""

from itertools import pairwise

import pytest

from realta.stellar import main_sequence as ms
from realta.stellar.zams import zams_luminosity, zams_radius

Z_SOLAR = 0.02


def test_zams_luminosity_radius_solar_calibration():
    """1 Msun, Z=0.02 ZAMS should be close to (but below) 1 Lsun/Rsun --
    the Sun has evolved somewhat off the ZAMS over 4.6 Gyr, so a modest
    undershoot is physically expected, not a bug.
    """
    l_zams = zams_luminosity(1.0, Z_SOLAR)
    r_zams = zams_radius(1.0, Z_SOLAR)
    assert 0.5 < l_zams < 1.0
    assert 0.7 < r_zams < 1.0


def test_ms_lifetime_solar_calibration():
    """1 Msun, Z=0.02 MS lifetime should be of order 10 Gyr (textbook
    solar MS lifetime), not e.g. 10x or 0.1x off.
    """
    lifetime = ms.t_ms(1.0, Z_SOLAR)
    assert 8_000.0 < lifetime < 13_000.0


@pytest.mark.parametrize("mass", [0.6, 1.0, 3.0, 5.0, 10.0, 20.0, 50.0])
def test_ms_radius_and_luminosity_monotonically_increase(mass):
    """Radius and luminosity must increase smoothly and monotonically
    across the MS for every mass in Realta's realistic range (up to and
    beyond mcut=8 Msun) -- this is the direct regression test for the
    session's caught bug, where R_MS collapsed to near-planet-size for
    M=5/20 Msun for most of the track before snapping back only in the
    last ~1 per cent of t_MS.
    """
    lifetime = ms.t_ms(mass, Z_SOLAR)
    fractions = [1e-6, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
    radii = [ms.ms_radius(mass, Z_SOLAR, lifetime * f) for f in fractions]
    luminosities = [ms.ms_luminosity(mass, Z_SOLAR, lifetime * f) for f in fractions]

    for r_prev, r_next in pairwise(radii):
        assert (
            r_next > r_prev * 0.99
        ), f"radius must not collapse across the MS for M={mass}: {radii}"
    for l_prev, l_next in pairwise(luminosities):
        assert (
            l_next > l_prev * 0.99
        ), f"luminosity must not collapse across the MS for M={mass}: {luminosities}"
    # Sanity floor: never smaller than a tenth of the ZAMS radius --
    # catches the specific class of bug this session found (a
    # near-total collapse), without being so strict it rejects a
    # legitimate small early-MS contraction.
    r_zams = zams_radius(mass, Z_SOLAR)
    assert all(r > 0.1 * r_zams for r in radii)


@pytest.mark.parametrize("mass", [0.6, 1.0, 3.0, 5.0, 10.0, 20.0, 50.0])
def test_ms_radius_and_luminosity_converge_to_tms_endpoint(mass):
    """At t -> t_MS, R_MS/L_MS must converge to R_TMS/L_TMS exactly (by
    construction of eqs. 12-13's perturbation decomposition).
    """
    lifetime = ms.t_ms(mass, Z_SOLAR)
    t = lifetime * (1.0 - 1e-9)
    assert ms.ms_radius(mass, Z_SOLAR, t) == pytest.approx(
        ms.r_tms(mass, Z_SOLAR), rel=1e-4
    )
    assert ms.ms_luminosity(mass, Z_SOLAR, t) == pytest.approx(
        ms.l_tms(mass, Z_SOLAR), rel=1e-4
    )


def test_phase_reports_hg_past_ms_lifetime_and_raises_past_bgb():
    """phase() must report k=2 (HG) once a star has left the MS but
    not yet reached the GB, and raise -- not guess -- once it reaches
    the GB, which is out of scope for this module (see its docstring
    and tests/test_hertzsprung_gap.py for the HG-specific coverage).
    """
    lifetime = ms.t_ms(5.0, Z_SOLAR)
    bgb = ms.t_bgb(5.0, Z_SOLAR)
    assert ms.phase(5.0, Z_SOLAR, lifetime * 0.5) == 1
    assert ms.phase(0.5, Z_SOLAR, 1.0) == 0
    assert ms.phase(5.0, Z_SOLAR, lifetime * 1.001) == 2
    with pytest.raises(ValueError, match="not implemented"):
        ms.phase(5.0, Z_SOLAR, bgb * 1.001)


def test_pinned_values_z_solar():
    """Exact self-consistency regression pin (this implementation
    against itself, not an independent reference run -- see
    docs/provenance.md for why no independent numeric reference exists
    yet). A change to any coefficient or equation will move these.
    """
    assert ms.t_ms(1.0, Z_SOLAR) == pytest.approx(11003.121, rel=1e-4)
    assert ms.r_tms(1.0, Z_SOLAR) == pytest.approx(1.623978, rel=1e-4)
    assert ms.l_tms(1.0, Z_SOLAR) == pytest.approx(2.261420, rel=1e-4)

    assert ms.t_ms(5.0, Z_SOLAR) == pytest.approx(104.0162, rel=1e-4)
    assert ms.r_tms(5.0, Z_SOLAR) == pytest.approx(6.141333, rel=1e-4)

    assert ms.t_ms(20.0, Z_SOLAR) == pytest.approx(8.675435, rel=1e-4)
    assert ms.r_tms(20.0, Z_SOLAR) == pytest.approx(16.13165, rel=1e-4)


def test_low_mass_degenerate_radius_floor():
    """Eq. (24): very-low-mass MS radius has a degenerate floor,
    0.0258*(1+X)^(5/3)*M^(-1/3), X = 0.76 - 3Z.
    """
    mass = 0.1
    x_hydrogen = 0.76 - 3.0 * Z_SOLAR
    floor = 0.0258 * (1.0 + x_hydrogen) ** (5.0 / 3.0) * mass ** (-1.0 / 3.0)
    r = ms.ms_radius(mass, Z_SOLAR, 1.0)
    assert r >= floor
