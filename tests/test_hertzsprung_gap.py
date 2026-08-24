"""Tests for the Hertzsprung-Gap radius/luminosity extension to
main_sequence.py (Hurley et al. 2000, Sec. 5.1.2, eqs. 25-30), and the
M_FGB scope boundary (stars with no GB phase, needing L_HeI/R_HeI,
which this module does not implement).
"""

from itertools import pairwise

import pytest

from realta.stellar import main_sequence as ms

Z_SOLAR = 0.02


@pytest.mark.parametrize("mass", [0.8, 1.0, 3.0, 5.0, 10.0])
def test_hg_radius_grows_smoothly_from_tms_to_ehg(mass):
    """R_HG(t) must increase monotonically across the HG (the defining
    physical feature of the Hertzsprung Gap -- rapid radius expansion
    at roughly constant or declining L), converging exactly to
    R_TMS/R_EHG at the endpoints.
    """
    t_ms_val = ms.t_ms(mass, Z_SOLAR)
    t_bgb_val = ms.t_bgb(mass, Z_SOLAR)
    fractions = [1e-9, 0.2, 0.4, 0.6, 0.8, 1.0 - 1e-9]
    radii = [
        ms.hg_radius(mass, Z_SOLAR, t_ms_val + f * (t_bgb_val - t_ms_val))
        for f in fractions
    ]
    for r_prev, r_next in pairwise(radii):
        assert r_next > r_prev * 0.999

    assert radii[0] == pytest.approx(ms.r_tms(mass, Z_SOLAR), rel=1e-4)
    assert radii[-1] == pytest.approx(ms.r_ehg(mass, Z_SOLAR), rel=1e-3)


@pytest.mark.parametrize("mass", [0.8, 1.0, 3.0, 5.0, 10.0])
def test_hg_luminosity_converges_to_endpoints(mass):
    t_ms_val = ms.t_ms(mass, Z_SOLAR)
    t_bgb_val = ms.t_bgb(mass, Z_SOLAR)
    l_start = ms.hg_luminosity(mass, Z_SOLAR, t_ms_val + 1e-9)
    l_end = ms.hg_luminosity(mass, Z_SOLAR, t_bgb_val - 1e-9)
    assert l_start == pytest.approx(ms.l_tms(mass, Z_SOLAR), rel=1e-4)
    assert l_end == pytest.approx(ms.l_ehg(mass, Z_SOLAR), rel=1e-3)


def test_phase_returns_two_during_hg():
    mass = 5.0
    t_ms_val = ms.t_ms(mass, Z_SOLAR)
    t_bgb_val = ms.t_bgb(mass, Z_SOLAR)
    mid_hg = (t_ms_val + t_bgb_val) / 2.0
    assert t_ms_val < mid_hg < t_bgb_val
    assert ms.phase(mass, Z_SOLAR, mid_hg) == 2


def test_phase_raises_past_bgb():
    mass = 5.0
    t_bgb_val = ms.t_bgb(mass, Z_SOLAR)
    with pytest.raises(ValueError, match="not implemented"):
        ms.phase(mass, Z_SOLAR, t_bgb_val * 1.001)


def test_m_fgb_solar_sanity_value():
    """M_FGB at Z=0.02 should be of order 10-15 Msun (the mass above
    which a star ignites helium in the HG rather than on the GB).
    """
    assert 10.0 < ms.m_fgb(Z_SOLAR) < 16.0


def test_l_ehg_and_r_ehg_raise_above_m_fgb():
    """Stars with M >= M_FGB skip the GB entirely and need L_HeI/R_HeI
    (helium-ignition quantities, part of the CHeB machinery this
    module does not implement) -- must raise, not silently return a
    GB-based value that doesn't apply to these stars.
    """
    mass_above_fgb = ms.m_fgb(Z_SOLAR) + 1.0
    with pytest.raises(ValueError, match="M_FGB"):
        ms.l_ehg(mass_above_fgb, Z_SOLAR)
    with pytest.raises(ValueError, match="M_FGB"):
        ms.r_ehg(mass_above_fgb, Z_SOLAR)


def test_phase_still_reports_hg_above_m_fgb_but_radius_raises():
    """phase() identifies the evolutionary type without needing
    M_FGB, but hg_radius()/hg_luminosity() must still raise for those
    stars since they need L_HeI/R_HeI -- see module docstring.
    """
    mass_above_fgb = ms.m_fgb(Z_SOLAR) + 1.0
    t_ms_val = ms.t_ms(mass_above_fgb, Z_SOLAR)
    t_bgb_val = ms.t_bgb(mass_above_fgb, Z_SOLAR)
    mid_hg = (t_ms_val + t_bgb_val) / 2.0
    assert ms.phase(mass_above_fgb, Z_SOLAR, mid_hg) == 2
    with pytest.raises(ValueError, match="M_FGB"):
        ms.hg_radius(mass_above_fgb, Z_SOLAR, mid_hg)
