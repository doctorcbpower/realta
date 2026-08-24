"""Tests for HG core-mass tracking (Hurley et al. 2000, eqs. 2, 28-30,
44, 66) -- the CE-prerequisite piece added on top of the HG radius/
luminosity module. See docs/science/rlof-ce-classifier-proposal.md's
"Core mass: implementation note" for why `core_mass_bgb` uses the
paper's own stated large-mass asymptotic approximation
(0.098*M^1.35), rather than the full eq. 44 (which needs the entire
GB core-mass-luminosity relation, eqs. 31-43 -- out of scope), and why
it's capped at `CORE_MASS_BGB_MAX_MASS` (~7.3 Msun, derived from eq.
44's own c1/c2 constants): above that mass the approximation itself
goes super-Chandrasekhar, which was found and fixed during
implementation (a real, corrected scope-boundary mistake, not
speculation -- the original framing was "valid for M >= M_HeF" with
no upper bound).
"""

from itertools import pairwise

import pytest

from realta.stellar import giant_branch as gb
from realta.stellar import main_sequence as ms

Z_SOLAR = 0.02


def test_core_mass_bgb_positive_and_increases_with_mass():
    masses = [1.0, 3.0, 5.0, 7.0]
    values = [gb.core_mass_bgb(m) for m in masses]
    assert all(v > 0.0 for v in values)
    assert all(b > a for a, b in pairwise(values))


def test_core_mass_bgb_stays_sub_chandrasekhar_within_supported_range():
    for mass in [1.0, 3.0, 5.0, 7.0]:
        assert gb.core_mass_bgb(mass) < gb._M_CHANDRASEKHAR


def test_core_mass_bgb_raises_above_max_mass():
    """This is the real bug this session caught: without this raise,
    core_mass_bgb(10.0) silently returned a super-Chandrasekhar core
    mass (2.17 Msun), which broke white_dwarf_radius() (NaN under the
    sqrt) rather than being flagged as out of scope.
    """
    assert gb.CORE_MASS_BGB_MAX_MASS == pytest.approx(7.317, abs=0.01)
    with pytest.raises(ValueError, match="super-Chandrasekhar"):
        gb.core_mass_bgb(10.0)
    with pytest.raises(ValueError, match="super-Chandrasekhar"):
        gb.core_mass_bgb(gb.CORE_MASS_BGB_MAX_MASS)  # boundary itself excluded


def test_m_hef_solar_sanity_value():
    """M_HeF at Z=0.02 should be close to ~2 Msun (the classic
    threshold above which helium ignition is non-degenerate).
    """
    assert ms.m_hef(Z_SOLAR) == pytest.approx(1.995, rel=1e-3)


def test_core_mass_ehg_uses_core_mass_bgb_in_intermediate_mass_range():
    mass = 5.0
    assert ms.m_hef(Z_SOLAR) <= mass < ms.m_fgb(Z_SOLAR)
    assert ms.core_mass_ehg(mass, Z_SOLAR) == pytest.approx(
        gb.core_mass_bgb(mass), rel=1e-12
    )


def test_core_mass_ehg_raises_below_m_hef():
    low_mass = ms.m_hef(Z_SOLAR) - 0.5
    with pytest.raises(ValueError, match="M_HeF"):
        ms.core_mass_ehg(low_mass, Z_SOLAR)


def test_core_mass_ehg_raises_above_m_fgb():
    high_mass = ms.m_fgb(Z_SOLAR) + 1.0
    with pytest.raises(ValueError, match="M_FGB"):
        ms.core_mass_ehg(high_mass, Z_SOLAR)


@pytest.mark.parametrize("mass", [3.0, 5.0, 7.0])
def test_core_mass_hg_grows_monotonically_and_converges_to_endpoints(mass):
    t_ms_val = ms.t_ms(mass, Z_SOLAR)
    t_bgb_val = ms.t_bgb(mass, Z_SOLAR)
    fractions = [1e-9, 0.3, 0.6, 1.0 - 1e-9]
    core_masses = [
        ms.core_mass_hg(mass, Z_SOLAR, t_ms_val + f * (t_bgb_val - t_ms_val))
        for f in fractions
    ]
    for m_prev, m_next in pairwise(core_masses):
        assert m_next >= m_prev  # monotonically non-decreasing

    m_c_ehg = ms.core_mass_ehg(mass, Z_SOLAR)
    rho = ms._rho_coefficient(mass)
    assert core_masses[0] == pytest.approx(rho * m_c_ehg, rel=1e-4)
    assert core_masses[-1] == pytest.approx(m_c_ehg, rel=1e-4)


def test_core_mass_hg_starts_below_end_of_hg_value():
    """rho < 1 always (eq. 29's form guarantees this for M > 0), so
    the core mass at the start of HG must be strictly less than at the
    end -- confirms the growth direction, not just monotonicity.
    """
    mass = 5.0
    rho = ms._rho_coefficient(mass)
    assert 0.0 < rho < 1.0
