"""Tests for stellar/cak_wind.py -- the CAK radiatively-driven wind
physics feeding the wind-capture accretion module
(binaries/wind_capture.py). Source: El Mellah & Casse (2017,
arXiv:1609.01532), verified directly against the paper before
implementing.
"""

from itertools import pairwise

import pytest

from realta.stellar import cak_wind as cw

# Vela X-1 donor (HD 77598), Friend & Castor (1982, ApJ 261, 293)
# Table 1: T=26000K, R/Rsun=35, L/Lsun=5e5, M_evol/Msun=24,
# Gamma_*=0.50, v_inf observed 700-1700 km/s, Mdot observed
# 0.6-2e-6 Msun/yr.
VELA_X1_MASS = 24.0
VELA_X1_RADIUS = 35.0
VELA_X1_LUMINOSITY = 5.0e5
VELA_X1_ALPHA = 0.5  # within Shimada et al. (1994)'s alpha in [0.47,0.52] for T=20-30kK


def test_escape_velocity_matches_known_solar_value():
    """Sanity/calibration check: v_esc(1 Msun, 1 Rsun) should match
    the well-known solar escape velocity (~617.5-618 km/s)."""
    v_esc = (2.0 * cw.G_MSUN_RSUN_KM2_S2 * 1.0 / 1.0) ** 0.5
    assert v_esc == pytest.approx(617.8, abs=0.5)


def test_eddington_luminosity_and_factor_match_vela_x1():
    """Sanity/calibration check against Friend & Castor's own Table 1
    Gamma_*=0.50 for Vela X-1's donor."""
    gamma = cw.eddington_factor(VELA_X1_MASS, VELA_X1_LUMINOSITY)
    assert gamma == pytest.approx(0.50, rel=0.15)


def test_wind_terminal_velocity_matches_vela_x1_observed_range():
    """Sanity/calibration check: eq. (7)'s calibrated terminal
    velocity should fall within Vela X-1's own observed range
    (700-1700 km/s, Friend & Castor Table 1) -- confirms the eq. (7)
    vs. eq. (5) choice made while implementing this (see
    wind_terminal_velocity's own docstring for why eq. 5 alone
    undershoots this range for the same parameters)."""
    v_inf = cw.wind_terminal_velocity(VELA_X1_MASS, VELA_X1_RADIUS, VELA_X1_ALPHA)
    assert 700.0 <= v_inf <= 1700.0


def test_wind_terminal_velocity_eq5_alone_undershoots_observed_range():
    """Regression guard for the actual finding that justified using
    eq. (7) over eq. (5): the basic point-source CAK formula (eq. 5)
    gives a terminal velocity BELOW Vela X-1's observed range for the
    same parameters -- confirms this was a real, checked reason to
    prefer eq. (7), not an arbitrary choice."""
    v_esc = (2.0 * cw.G_MSUN_RSUN_KM2_S2 * VELA_X1_MASS / VELA_X1_RADIUS) ** 0.5
    v_inf_eq5 = v_esc * (VELA_X1_ALPHA / (1.0 - VELA_X1_ALPHA)) ** 0.5
    assert v_inf_eq5 < 700.0


def test_wind_velocity_increases_monotonically_and_approaches_terminal():
    v_inf = cw.wind_terminal_velocity(VELA_X1_MASS, VELA_X1_RADIUS, VELA_X1_ALPHA)
    radii = [VELA_X1_RADIUS * f for f in (1.01, 2.0, 5.0, 20.0, 1000.0)]
    velocities = [cw.wind_velocity(r, VELA_X1_RADIUS, v_inf) for r in radii]
    assert all(v2 > v1 for v1, v2 in pairwise(velocities))
    assert velocities[-1] == pytest.approx(v_inf, rel=1e-3)


def test_wind_velocity_zero_at_and_below_stellar_radius():
    v_inf = cw.wind_terminal_velocity(VELA_X1_MASS, VELA_X1_RADIUS, VELA_X1_ALPHA)
    assert cw.wind_velocity(VELA_X1_RADIUS, VELA_X1_RADIUS, v_inf) == 0.0
    assert cw.wind_velocity(VELA_X1_RADIUS * 0.5, VELA_X1_RADIUS, v_inf) == 0.0


def test_wind_mass_loss_rate_matches_vela_x1_order_of_magnitude():
    """Sanity/calibration check: order-of-magnitude agreement with
    Vela X-1's observed mass-loss rate (0.6-2e-6 Msun/yr), using
    Q~900 (their Sec. 3.3(ii) fiducial value for OB supergiants) --
    not an exact match (Q is not independently calibrated per star
    here, matching both papers' own acknowledged sensitivity to this
    parameter)."""
    gamma = cw.eddington_factor(VELA_X1_MASS, VELA_X1_LUMINOSITY)
    mdot = cw.wind_mass_loss_rate(
        VELA_X1_MASS, VELA_X1_LUMINOSITY, VELA_X1_ALPHA, 900.0, gamma
    )
    assert 1e-7 < mdot < 1e-4  # within ~2 orders of magnitude of observed


def test_wind_mass_loss_rate_increases_with_q_force():
    """Direction/sensitivity check: a larger Q force multiplier (more
    line-driving strength) should increase the mass-loss rate."""
    gamma = cw.eddington_factor(VELA_X1_MASS, VELA_X1_LUMINOSITY)
    mdot_low_q = cw.wind_mass_loss_rate(
        VELA_X1_MASS, VELA_X1_LUMINOSITY, VELA_X1_ALPHA, 500.0, gamma
    )
    mdot_high_q = cw.wind_mass_loss_rate(
        VELA_X1_MASS, VELA_X1_LUMINOSITY, VELA_X1_ALPHA, 2000.0, gamma
    )
    assert mdot_high_q > mdot_low_q


def test_wind_mass_loss_rate_increases_with_gamma():
    """Direction check: a higher Eddington factor (more luminous
    relative to Eddington) should drive a stronger wind."""
    mdot_low_gamma = cw.wind_mass_loss_rate(
        VELA_X1_MASS, VELA_X1_LUMINOSITY, VELA_X1_ALPHA, 900.0, 0.2
    )
    mdot_high_gamma = cw.wind_mass_loss_rate(
        VELA_X1_MASS, VELA_X1_LUMINOSITY, VELA_X1_ALPHA, 900.0, 0.6
    )
    assert mdot_high_gamma > mdot_low_gamma


@pytest.mark.parametrize("bad_alpha", [0.0, 1.0, -0.1, 1.5])
def test_alpha_out_of_range_rejected(bad_alpha):
    with pytest.raises(ValueError):
        cw.wind_terminal_velocity(VELA_X1_MASS, VELA_X1_RADIUS, bad_alpha)


@pytest.mark.parametrize("bad_gamma", [0.0, 1.0, -0.1, 1.5])
def test_gamma_out_of_range_rejected(bad_gamma):
    with pytest.raises(ValueError):
        cw.wind_mass_loss_rate(
            VELA_X1_MASS, VELA_X1_LUMINOSITY, VELA_X1_ALPHA, 900.0, bad_gamma
        )
