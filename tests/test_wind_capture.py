"""Tests for binaries/wind_capture.py -- converting a CAK wind
(stellar/cak_wind.py) into a compact-object accretion rate and
circularization-radius estimate. Source: El Mellah & Casse (2017,
arXiv:1609.01532), verified directly against the paper before
implementing.
"""

import pytest

from realta.binaries import wind_capture as wc

M_COMPACT = 1.4  # Msun, neutron star
V_WIND = 1280.0  # km/s, illustrative (Vela-X-1-like terminal velocity)
V_ORBITAL = 200.0  # km/s, illustrative
SEPARATION = 60.0  # Rsun, illustrative


def test_escape_velocity_matches_known_solar_value():
    assert wc.escape_velocity(1.0, 1.0) == pytest.approx(617.8, abs=0.5)


def test_relative_wind_velocity_is_quadrature_sum():
    v_rel = wc.relative_wind_velocity(V_WIND, V_ORBITAL)
    assert v_rel == pytest.approx((V_WIND**2 + V_ORBITAL**2) ** 0.5)
    assert v_rel > V_WIND  # orbital motion always adds some relative speed


def test_relative_wind_velocity_reduces_to_wind_speed_when_orbital_is_zero():
    assert wc.relative_wind_velocity(V_WIND, 0.0) == V_WIND


def test_accretion_radius_smaller_for_faster_relative_wind():
    """Direction check: a faster wind should be harder to capture
    (smaller accretion radius) -- R_acc ~ 1/v_rel^2."""
    r_acc_slow = wc.accretion_radius(M_COMPACT, 500.0)
    r_acc_fast = wc.accretion_radius(M_COMPACT, 1500.0)
    assert r_acc_fast < r_acc_slow


def test_accretion_radius_larger_for_more_massive_compact_object():
    r_acc_ns = wc.accretion_radius(1.4, V_WIND)
    r_acc_bh = wc.accretion_radius(10.0, V_WIND)
    assert r_acc_bh > r_acc_ns


def test_bhl_accretion_fraction_smaller_for_wider_orbit():
    """Direction check: a wider separation captures a smaller fraction
    of the wind (beta ~ (R_acc/a)^2)."""
    v_rel = wc.relative_wind_velocity(V_WIND, V_ORBITAL)
    r_acc = wc.accretion_radius(M_COMPACT, v_rel)
    beta_close = wc.bhl_accretion_fraction(r_acc, 20.0)
    beta_wide = wc.bhl_accretion_fraction(r_acc, 200.0)
    assert beta_close > beta_wide


def test_bhl_accretion_fraction_is_clamped_to_one():
    """An unphysically small separation (R_acc comparable to or
    exceeding it) must not give beta > 1."""
    huge_r_acc = 1000.0
    tiny_separation = 1.0
    beta = wc.bhl_accretion_fraction(huge_r_acc, tiny_separation)
    assert beta == 1.0


def test_wind_capture_rate_is_simple_product():
    assert wc.wind_capture_rate(1e-6, 0.1) == pytest.approx(1e-7)
    assert wc.wind_capture_rate(1e-6, 0.0) == 0.0


def test_bhl_accretion_rate_simple_matches_manual_composition():
    mdot_wind = 1.6e-5
    expected_v_rel = wc.relative_wind_velocity(V_WIND, V_ORBITAL)
    expected_r_acc = wc.accretion_radius(M_COMPACT, expected_v_rel)
    expected_beta = wc.bhl_accretion_fraction(expected_r_acc, SEPARATION)
    expected = wc.wind_capture_rate(mdot_wind, expected_beta)

    actual = wc.bhl_accretion_rate_simple(
        M_COMPACT, mdot_wind, V_WIND, V_ORBITAL, SEPARATION
    )
    assert actual == pytest.approx(expected)


def test_circularization_radius_fraction_matches_exact_pinned_value():
    """Precise pin for the exact (v_orbital/v_rel)^4 exponent and 1/4
    prefactor -- the order-of-magnitude test below is too loose to
    catch a wrong exponent (e.g. ^2 instead of ^4 still lands inside
    that test's 1e-5 to 1e-1 tolerance for these illustrative
    numbers). Pinned against an independently-computed literal value,
    not a re-derivation of the formula under test, so a change to
    either the exponent or the prefactor is actually caught -- I
    verified sensitivity directly: changed the exponent from ^4 to
    ^2, confirmed neither this nor the direction/order-of-magnitude
    tests originally caught it, added this pin, reverted."""
    v_rel = wc.relative_wind_velocity(V_WIND, V_ORBITAL)
    frac = wc.circularization_radius_fraction(V_ORBITAL, v_rel)
    assert frac == pytest.approx(0.00014199369139068394, rel=1e-9)


def test_circularization_radius_fraction_within_order_of_magnitude_of_paper():
    """Sanity check against El Mellah & Casse's own stated order of
    magnitude (Sec. 4.4: R_circ/R_acc ~ 1e-3 to 1e-2) -- see
    circularization_radius_fraction's own docstring for why this is a
    lower-confidence estimate than the rest of this module (not their
    own closed-form result -- they give none) and the real bug this
    caught (a first-draft formula off by four orders of magnitude)."""
    frac = wc.circularization_radius_fraction(
        V_ORBITAL, wc.relative_wind_velocity(V_WIND, V_ORBITAL)
    )
    assert 1e-5 < frac < 1e-1  # within ~1-2 orders of magnitude, not exact


def test_circularization_radius_fraction_increases_with_orbital_velocity():
    """Direction check: faster orbital motion (relative to the wind)
    should shear more angular momentum into the flow, raising
    R_circ/R_acc."""
    v_rel = wc.relative_wind_velocity(V_WIND, V_ORBITAL)
    frac_slow_orbit = wc.circularization_radius_fraction(50.0, v_rel)
    frac_fast_orbit = wc.circularization_radius_fraction(300.0, v_rel)
    assert frac_fast_orbit > frac_slow_orbit


def test_circularization_radius_is_fraction_times_accretion_radius():
    v_rel = wc.relative_wind_velocity(V_WIND, V_ORBITAL)
    r_acc = wc.accretion_radius(M_COMPACT, v_rel)
    r_circ = wc.circularization_radius(V_ORBITAL, v_rel, r_acc)
    expected = wc.circularization_radius_fraction(V_ORBITAL, v_rel) * r_acc
    assert r_circ == pytest.approx(expected)
    assert r_circ < r_acc  # circularization radius must be well inside R_acc


def test_circularization_radius_fraction_zero_orbital_velocity_gives_zero():
    """No orbital shear at all (v_orbital=0) must give zero
    circularization -- purely radial infall has no captured angular
    momentum in this model."""
    v_rel = wc.relative_wind_velocity(V_WIND, 0.0)
    assert wc.circularization_radius_fraction(0.0, v_rel) == 0.0


@pytest.mark.parametrize("bad_v_rel", [0.0, -1.0])
def test_accretion_radius_rejects_non_positive_v_rel(bad_v_rel):
    with pytest.raises(ValueError):
        wc.accretion_radius(M_COMPACT, bad_v_rel)


def test_bhl_accretion_fraction_rejects_non_positive_separation():
    with pytest.raises(ValueError):
        wc.bhl_accretion_fraction(1.0, 0.0)
