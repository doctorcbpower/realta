"""Tests for compact-remnant radius (Hurley et al. 2000, eq. 91, and
Sec. 6.3's core-radius definition for HG/GB stars) -- used as R_c1 in
the eventual CE binding-energy formula. See remnant.py's module
docstring for why this one wasn't put through a paste-verification
round the way the dense appendix tables were (clean text, plus a
strong physical sanity check).
"""

from itertools import pairwise

import pytest

from realta.stellar import giant_branch as gb
from realta.stellar import main_sequence as ms
from realta.stellar import remnant

Z_SOLAR = 0.02


def test_white_dwarf_radius_matches_known_real_value():
    """A 0.6 Msun white dwarf (Sirius B territory) has a well-known
    real radius of ~8000-9000 km (~0.0115-0.013 Rsun). This is the
    check that gave confidence in the formula before implementation.
    """
    r = remnant.white_dwarf_radius(0.6)
    assert 0.010 < r < 0.015


def test_white_dwarf_radius_decreases_with_mass():
    """The defining (and easy-to-get-backwards) feature of the
    degenerate mass-radius relation: more massive white dwarfs are
    SMALLER, not larger.
    """
    masses = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
    radii = [remnant.white_dwarf_radius(m) for m in masses]
    for r_prev, r_next in pairwise(radii):
        assert r_next < r_prev


def test_white_dwarf_radius_shrinks_toward_chandrasekhar_mass():
    """Approaching M_ch, the radius must shrink toward zero (the
    Chandrasekhar-limit collapse), not blow up or go negative.
    """
    r_near_limit = remnant.white_dwarf_radius(1.43)
    r_moderate = remnant.white_dwarf_radius(1.0)
    assert 0.0 < r_near_limit < r_moderate


def test_white_dwarf_radius_floored_at_neutron_star_radius():
    """eq. 91's max(R_NS, ...) floor -- must never return a radius
    smaller than the fixed 10 km neutron-star radius.
    """
    for mass in [0.1, 0.5, 1.0, 1.43]:
        assert remnant.white_dwarf_radius(mass) >= remnant.R_NEUTRON_STAR


def test_core_radius_matches_white_dwarf_radius_for_intermediate_mass_donor():
    """For donor_mass >= M_HeF (the only regime this codebase's
    core-mass tracking reaches), core_radius must equal
    white_dwarf_radius(core_mass) exactly -- no other physics applied.
    """
    donor_mass, core_mass = 5.0, 0.86
    assert donor_mass >= ms.m_hef(Z_SOLAR)
    assert remnant.core_radius(core_mass, donor_mass, Z_SOLAR) == pytest.approx(
        remnant.white_dwarf_radius(core_mass)
    )


def test_core_radius_raises_below_m_hef():
    low_mass = ms.m_hef(Z_SOLAR) - 0.5
    with pytest.raises(ValueError, match="M_HeF"):
        remnant.core_radius(0.3, low_mass, Z_SOLAR)


def test_core_radius_end_to_end_for_hg_donor():
    """Full chain: core_mass_hg -> core_radius, for a real HG donor
    within the supported mass range, must produce a finite, positive,
    physically reasonable (sub-white-dwarf-scale) radius.
    """
    mass = 5.0
    assert mass < gb.CORE_MASS_BGB_MAX_MASS
    t_ms_val = ms.t_ms(mass, Z_SOLAR)
    t_bgb_val = ms.t_bgb(mass, Z_SOLAR)
    mid_hg = (t_ms_val + t_bgb_val) / 2.0

    core_mass = ms.core_mass_hg(mass, Z_SOLAR, mid_hg)
    r_c = remnant.core_radius(core_mass, mass, Z_SOLAR)
    assert 0.0 < r_c < 0.1  # white-dwarf-scale, not stellar-scale
