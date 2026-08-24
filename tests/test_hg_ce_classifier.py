"""Tests for the HG-donor extension of classify_rlof() -- HTP02's own
GB q_crit formula (eqs. 56-57), reused for HG donors per Zuo & Li
(2014)'s eq. 1, and the resulting COMMON_ENVELOPE outcome (as opposed
to IMMEDIATE_MERGER for MS donors). See
docs/science/rlof-ce-classifier-proposal.md's "Literature findings"
note. The CE outcome itself (survive vs. merge) is not resolved by
this classifier -- that needs the not-yet-implemented energy-balance
solve.
"""

import pytest

from realta.binaries.interaction import (
    RLOFOutcome,
    classify_rlof,
    find_rlof_onset,
    hg_q_crit,
)
from realta.stellar import giant_branch as gb
from realta.stellar import main_sequence as ms

Z_SOLAR = 0.02


def _mid_hg_age(mass: float, z: float = Z_SOLAR) -> float:
    t_ms_val = ms.t_ms(mass, z)
    t_bgb_val = ms.t_bgb(mass, z)
    return (t_ms_val + t_bgb_val) / 2.0


def test_hg_q_crit_lower_than_htp02_crude_hg_default():
    """HTP02 itself uses a crude fixed q_crit=4 for HG donors, calling
    it "rather approximate." The refined, core-mass-dependent formula
    should give a materially different (in practice much lower) value
    for a typical HG donor -- confirms the refinement is actually
    doing something, not silently reducing to the same number.
    """
    donor_mass = 5.0
    age = _mid_hg_age(donor_mass)
    q_crit = hg_q_crit(donor_mass, Z_SOLAR, age)
    assert 0.0 < q_crit < 4.0


def test_hg_q_crit_matches_direct_formula_evaluation():
    donor_mass = 5.0
    age = _mid_hg_age(donor_mass)
    core_mass = ms.core_mass_hg(donor_mass, Z_SOLAR, age)
    x = gb.mass_radius_exponent(Z_SOLAR)
    expected = (1.67 - x + 2.0 * (core_mass / donor_mass) ** 5) / 2.13
    assert hg_q_crit(donor_mass, Z_SOLAR, age) == pytest.approx(expected)


def test_classify_rlof_hg_donor_unstable_gives_common_envelope_not_merger():
    """The key behavioural difference from MS donors: a dynamically
    unstable HG donor forms a common envelope, not an immediate
    merger -- HG donors are in HTP02's CE-eligible donor list, unlike
    MS donors.
    """
    donor_mass, companion_mass = 5.0, 3.0
    age = _mid_hg_age(donor_mass)
    q_crit = hg_q_crit(donor_mass, Z_SOLAR, age)
    assert donor_mass / companion_mass > q_crit  # confirms instability

    outcome = classify_rlof(donor_mass, companion_mass, 5.0, Z_SOLAR, age)
    assert outcome == RLOFOutcome.COMMON_ENVELOPE


def test_classify_rlof_hg_donor_stable_mass_transfer():
    donor_mass, companion_mass = 5.0, 20.0
    age = _mid_hg_age(donor_mass)
    q_crit = hg_q_crit(donor_mass, Z_SOLAR, age)
    assert donor_mass / companion_mass < q_crit  # confirms stability

    outcome = classify_rlof(donor_mass, companion_mass, 5.0, Z_SOLAR, age)
    assert outcome == RLOFOutcome.STABLE_MASS_TRANSFER


def test_classify_rlof_hg_donor_detached_when_wide():
    donor_mass, companion_mass = 5.0, 3.0
    age = _mid_hg_age(donor_mass)
    outcome = classify_rlof(donor_mass, companion_mass, 200.0, Z_SOLAR, age)
    assert outcome == RLOFOutcome.DETACHED


def test_classify_rlof_hg_donor_beyond_core_mass_range_not_modelled():
    """A donor still on the HG (k=2) but with mass >=
    CORE_MASS_BGB_MAX_MASS: phase() reports k=2, but hg_radius()/
    hg_q_crit() raise (core mass out of the supported range) -- must
    be caught and reported as PHASE_NOT_MODELLED, not crash or guess.
    """
    donor_mass = 10.0
    assert donor_mass >= gb.CORE_MASS_BGB_MAX_MASS
    age = _mid_hg_age(donor_mass)
    assert ms.phase(donor_mass, Z_SOLAR, age) == 2  # still identified as HG

    outcome = classify_rlof(donor_mass, 3.0, 20.0, Z_SOLAR, age)
    assert outcome == RLOFOutcome.PHASE_NOT_MODELLED


def test_find_rlof_onset_reaches_common_envelope_via_hg_search():
    """End-to-end: a binary wide enough to stay detached through the
    whole MS, but whose donor's HG expansion overflows it, must be
    found by the root-finder AND correctly classified as
    COMMON_ENVELOPE (not IMMEDIATE_MERGER) -- confirms find_rlof_onset
    actually reaches classify_rlof's new HG branch, not just that the
    branch exists in isolation.
    """
    donor_mass, companion_mass, a = 5.0, 3.0, 100.0
    t_rlof, outcome, donor_is_star1 = find_rlof_onset(
        donor_mass, companion_mass, a, Z_SOLAR
    )

    t_ms_val = ms.t_ms(donor_mass, Z_SOLAR)
    t_bgb_val = ms.t_bgb(donor_mass, Z_SOLAR)
    assert t_ms_val < t_rlof < t_bgb_val  # crosses during HG, not MS
    assert donor_is_star1
    assert outcome == RLOFOutcome.COMMON_ENVELOPE

    # Cross-check against a direct classify_rlof() call at the found time.
    assert classify_rlof(donor_mass, companion_mass, a, Z_SOLAR, t_rlof) in (
        RLOFOutcome.COMMON_ENVELOPE,
        RLOFOutcome.DETACHED,  # floating-point boundary tolerance
    )


def test_find_rlof_onset_reaches_stable_mass_transfer_via_hg_search():
    donor_mass, companion_mass, a = 5.0, 20.0, 100.0
    t_rlof, outcome, donor_is_star1 = find_rlof_onset(
        donor_mass, companion_mass, a, Z_SOLAR
    )

    t_ms_val = ms.t_ms(donor_mass, Z_SOLAR)
    t_bgb_val = ms.t_bgb(donor_mass, Z_SOLAR)
    assert t_ms_val < t_rlof < t_bgb_val
    assert donor_is_star1
    assert outcome == RLOFOutcome.STABLE_MASS_TRANSFER
