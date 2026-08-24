"""Tests for the giant-branch quantities needed as HG boundary values
(L_BGB, R_GB, mass-radius exponent x) -- Hurley et al. (2000), eqs.
10, 46-48. See docs/science/rlof-ce-classifier-proposal.md for the
coefficient-verification history, including a real bug this module's
own development caught: an entire b2-clamping step (eq. 46's A
coefficient) had been dropped, confirmed and fixed by cross-checking
against the paper's own illustrative Z=0.02 example formula.
"""

from itertools import pairwise

import pytest

from realta.stellar import giant_branch as gb

Z_SOLAR = 0.02


def test_l_bgb_solar_calibration():
    """1 Msun base-of-giant-branch luminosity should be a few Lsun --
    textbook value is ~2-3 Lsun for the Sun's eventual BGB point.
    """
    l_bgb = gb.l_bgb(1.0, Z_SOLAR)
    assert 1.5 < l_bgb < 4.0


@pytest.mark.parametrize("mass", [0.8, 1.0, 3.0, 5.0, 10.0, 20.0])
def test_l_bgb_and_r_gb_positive_and_increase_with_mass(mass):
    l_bgb = gb.l_bgb(mass, Z_SOLAR)
    r = gb.r_gb(mass, l_bgb, Z_SOLAR)
    assert l_bgb > 0.0
    assert r > 0.0


def test_l_bgb_monotonically_increases_with_mass():
    masses = [0.8, 1.0, 3.0, 5.0, 10.0, 20.0]
    values = [gb.l_bgb(m, Z_SOLAR) for m in masses]
    assert all(b > a for a, b in pairwise(values))


def test_r_gb_matches_illustrative_z_solar_example_formula():
    """Hurley et al. (2000), Sec. 5.2 text, give an illustrative
    Z=0.02 example: R_GB ~= 1.1*M^-0.3*(L^0.4 + 0.383*L^0.76). This is
    a simplified stand-in (fixed x=0.3 exponent instead of the full
    A=min(b4*M^-b5, b6*M^-b7)), so agreement should be close (~10-20%)
    but not exact -- this is the check that caught the missing b2
    clamp during this module's development (disagreement was up to
    14x before the fix).
    """
    for mass, luminosity in [(1.0, 2.5), (5.0, 860.0), (10.0, 9585.0)]:
        mine = gb.r_gb(mass, luminosity, Z_SOLAR)
        illustrative = 1.1 * mass**-0.3 * (luminosity**0.4 + 0.383 * luminosity**0.76)
        assert mine == pytest.approx(illustrative, rel=0.25)


def test_mass_radius_exponent_solar_value():
    """x at Z=0.02 (zeta=0) should be close to the paper's own stated
    example value of ~0.3 (used directly in its illustrative R_GB
    formula, and in Zuo & Li (2014)'s HG/GB q_crit formula).
    """
    x = gb.mass_radius_exponent(Z_SOLAR)
    assert x == pytest.approx(0.30406, rel=1e-9)
