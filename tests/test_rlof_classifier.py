"""Tests for the RLOF outcome classifier (Stage 2 of
docs/science/rlof-ce-classifier-proposal.md). MS-donor scope only --
see src/realta/binaries/interaction.py's module docstring.
"""

import pytest

from realta.binaries.interaction import (
    Q_CRIT_MS,
    RLOFOutcome,
    apply_stable_mass_transfer,
    classify_rlof,
    find_rlof_onset,
    merge_stellar_masses,
    roche_lobe_radius,
)
from realta.stellar import main_sequence as ms

Z_SOLAR = 0.02


def test_roche_lobe_radius_equal_mass_sanity_value():
    """q1=1 (equal masses) should give R_L1/a ~ 0.38 -- the well-known
    literature value for the Eggleton (1983) fit at q=1.
    """
    assert roche_lobe_radius(1.0, 1.0) == pytest.approx(0.3789, abs=1e-3)


def test_roche_lobe_radius_scales_with_separation():
    assert roche_lobe_radius(10.0, 1.0) == pytest.approx(
        10.0 * roche_lobe_radius(1.0, 1.0)
    )


def test_roche_lobe_radius_rejects_non_positive_ratio():
    with pytest.raises(ValueError, match="mass_ratio"):
        roche_lobe_radius(1.0, 0.0)


def test_classify_detached_when_donor_radius_below_roche_lobe():
    """Wide separation: donor's Roche lobe is far larger than its
    actual radius -- no interaction.
    """
    donor_mass, companion_mass = 5.0, 3.0
    age = ms.t_ms(donor_mass, Z_SOLAR) * 0.9
    outcome = classify_rlof(donor_mass, companion_mass, 100.0, Z_SOLAR, age)
    assert outcome == RLOFOutcome.DETACHED


def test_classify_immediate_merger_for_high_mass_ratio_dynamical_donor():
    """q1 = 5/3 = 1.667 > Q_CRIT_MS -- a donor filling its Roche lobe
    with this mass ratio must merge dynamically (HTP02 Sec. 2.6.4).
    """
    donor_mass, companion_mass = 5.0, 3.0
    assert donor_mass / companion_mass > Q_CRIT_MS
    age = ms.t_ms(donor_mass, Z_SOLAR) * 0.9
    outcome = classify_rlof(donor_mass, companion_mass, 10.0, Z_SOLAR, age)
    assert outcome == RLOFOutcome.IMMEDIATE_MERGER


def test_classify_stable_mass_transfer_for_low_mass_ratio_donor():
    """q1 = 2/10 = 0.2 < Q_CRIT_MS -- a donor filling its Roche lobe
    with this mass ratio undergoes stable mass transfer.
    """
    donor_mass, companion_mass = 2.0, 10.0
    assert donor_mass / companion_mass < Q_CRIT_MS
    age = ms.t_ms(donor_mass, Z_SOLAR) * 0.9
    outcome = classify_rlof(donor_mass, companion_mass, 8.0, Z_SOLAR, age)
    assert outcome == RLOFOutcome.STABLE_MASS_TRANSFER


def test_classify_phase_not_modelled_past_ms_lifetime():
    """A donor past t_MS is HG/GB+ -- out of this module's scope; must
    return PHASE_NOT_MODELLED rather than guessing an outcome.
    """
    donor_mass, companion_mass = 5.0, 3.0
    age = ms.t_ms(donor_mass, Z_SOLAR) * 1.01
    outcome = classify_rlof(donor_mass, companion_mass, 10.0, Z_SOLAR, age)
    assert outcome == RLOFOutcome.PHASE_NOT_MODELLED


def test_classify_sensitivity_to_q_crit_override():
    """Overriding q_crit_ms must actually change the classification --
    confirms the threshold is live, not a dead parameter.
    """
    donor_mass, companion_mass = 5.0, 3.0
    age = ms.t_ms(donor_mass, Z_SOLAR) * 0.9
    default_outcome = classify_rlof(donor_mass, companion_mass, 10.0, Z_SOLAR, age)
    raised_outcome = classify_rlof(
        donor_mass, companion_mass, 10.0, Z_SOLAR, age, q_crit_ms=10.0
    )
    assert default_outcome == RLOFOutcome.IMMEDIATE_MERGER
    assert raised_outcome == RLOFOutcome.STABLE_MASS_TRANSFER


def test_merge_stellar_masses_is_conservative():
    assert merge_stellar_masses(5.0, 3.0) == pytest.approx(8.0)


def test_find_rlof_onset_wide_binary_never_interacts():
    """a=200 stays detached through both the MS and the (much more
    radius-expanded) HG -- a=100 does NOT (the donor's HG expansion
    eventually overflows it, confirmed once HG search was added), so
    "wide enough" now has to account for both phases, not just the MS.
    """
    t_rlof, outcome, _ = find_rlof_onset(5.0, 3.0, 200.0, Z_SOLAR)
    assert t_rlof == float("inf")
    assert outcome == RLOFOutcome.DETACHED


def test_find_rlof_onset_crosses_during_hg_not_ms():
    """a=100 is wide enough to stay detached through the whole MS, but
    the donor's radius expansion during HG eventually overflows it --
    confirms the HG search extension actually does something, not just
    reproduces the MS-only result.
    """
    t_rlof, outcome, donor_is_star1 = find_rlof_onset(5.0, 3.0, 100.0, Z_SOLAR)
    donor_mass = 5.0
    assert ms.t_ms(donor_mass, Z_SOLAR) < t_rlof < ms.t_bgb(donor_mass, Z_SOLAR)
    assert donor_is_star1
    assert outcome != RLOFOutcome.DETACHED


def test_find_rlof_onset_finds_crossing_and_correct_donor():
    """m1=5, m2=3, a=10: only m1 (the more massive, faster-expanding
    star) crosses its Roche lobe within the MS -- confirmed by direct
    radius/Roche-lobe evaluation at the returned crossing time.
    """
    t_rlof, outcome, donor_is_star1 = find_rlof_onset(5.0, 3.0, 10.0, Z_SOLAR)
    assert 0.0 < t_rlof < ms.t_ms(5.0, Z_SOLAR)
    assert donor_is_star1
    assert outcome == RLOFOutcome.IMMEDIATE_MERGER  # q1 = 5/3 > Q_CRIT_MS

    r_donor_at_cross = ms.ms_radius(5.0, Z_SOLAR, t_rlof)
    r_l1 = roche_lobe_radius(10.0, 5.0 / 3.0)
    assert r_donor_at_cross == pytest.approx(r_l1, rel=1e-6)


def test_find_rlof_onset_identifies_star2_as_donor_when_more_massive():
    """m1=2 (light), m2=10 (heavy): star 2 is the one that expands
    fast enough to fill its own Roche lobe first, despite being
    labelled 'star 2' -- the function must not assume m1 is always the
    donor.
    """
    t_rlof, outcome, donor_is_star1 = find_rlof_onset(2.0, 10.0, 8.0, Z_SOLAR)
    assert t_rlof < float("inf")
    assert not donor_is_star1
    assert outcome == RLOFOutcome.IMMEDIATE_MERGER  # q1 = 10/2 > Q_CRIT_MS


def test_find_rlof_onset_born_overflowing_returns_zero():
    """Separation small enough that the Roche lobe is already smaller
    than the ZAMS radius -- RLOF from t=0.
    """
    t_rlof, outcome, _ = find_rlof_onset(5.0, 3.0, 3.0, Z_SOLAR)
    assert t_rlof == 0.0
    assert outcome != RLOFOutcome.DETACHED


def test_find_rlof_onset_favours_the_more_massive_star_as_donor():
    """Emergent property, not an assumption: because R_L1/a increases
    monotonically with the donor's own mass ratio (see the roche_lobe
    sanity values above -- 0.207 at q1=0.1 vs. 0.578 at q1=10), the
    lighter star's Roche lobe is proportionally smaller and easier to
    fill, but in practice its slower radius growth means the *heavier*
    star still tends to reach its (larger) Roche lobe first across a
    wide range of mass/separation combinations -- confirmed by a
    parameter sweep during this module's development, which did not
    find a stable-mass-transfer case via automatic donor selection.
    This is a real, checked consequence of q_crit_ms=0.695 < 1, not
    something to special-case; classify_rlof() itself still returns
    STABLE_MASS_TRANSFER correctly when a low-q1 donor is specified
    directly (see the tests above). PHASE_NOT_MODELLED is also skipped
    here, the same way DETACHED is -- the (2.0, 10.0, 20.0) case
    crosses right at the MS/HG boundary with a 10 Msun donor, above
    CORE_MASS_BGB_MAX_MASS, so hg_q_crit() correctly can't classify it;
    that gap is covered directly by
    tests/test_hg_ce_classifier.py::test_classify_rlof_hg_donor_beyond_core_mass_range_not_modelled.
    """
    for m1, m2, a in [(2.0, 10.0, 20.0), (3.0, 4.0, 8.0), (1.5, 20.0, 8.0)]:
        _t, outcome, donor_is_star1 = find_rlof_onset(m1, m2, a, Z_SOLAR)
        if outcome in (RLOFOutcome.DETACHED, RLOFOutcome.PHASE_NOT_MODELLED):
            continue
        donor_mass = m1 if donor_is_star1 else m2
        companion_mass = m2 if donor_is_star1 else m1
        assert donor_mass >= companion_mass
        assert outcome == RLOFOutcome.IMMEDIATE_MERGER


def test_apply_stable_mass_transfer_conserves_mass_and_reaches_detachment():
    """New Roche-lobe radius, evaluated at the new (widened) separation
    and mass ratio, must exactly equal the donor's (unchanged) radius
    -- the self-consistent instantaneous-detachment point this
    function is defined to reach.
    """
    donor_mass, companion_mass, a = 2.0, 10.0, 8.0
    age = ms.t_ms(donor_mass, Z_SOLAR) * 0.9
    r_donor = ms.ms_radius(donor_mass, Z_SOLAR, age)

    new_donor, new_companion, new_a = apply_stable_mass_transfer(
        donor_mass, companion_mass, a, r_donor
    )

    assert new_donor + new_companion == pytest.approx(donor_mass + companion_mass)
    assert new_donor < new_companion  # direction: donor gets lighter, not heavier
    assert new_a > a  # orbit widens

    r_l1_after = roche_lobe_radius(new_a, new_donor / new_companion)
    assert r_l1_after == pytest.approx(r_donor, rel=1e-6)


def test_apply_stable_mass_transfer_rejects_wrong_direction():
    """Only reachable when donor_mass < companion_mass (the only regime
    classify_rlof() ever labels STABLE_MASS_TRANSFER) -- must reject
    the opposite ordering rather than silently doing something wrong.
    """
    with pytest.raises(ValueError, match="donor_mass < companion_mass"):
        apply_stable_mass_transfer(10.0, 2.0, 8.0, 3.0)
