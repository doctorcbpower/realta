"""Tests for B3's MS-mass-gainer rejuvenation (docs/science/
paper1-detailed-work-breakdown.md) -- Tout, Aarseth, Pols & Eggleton
(1997, MNRAS 291, 732) eq. (41), verified directly against the paper
before implementing (see binaries/interaction.py::rejuvenate_ms_gainer's
own docstring for the transcription).
"""

import numpy as np
import pytest

from realta.binaries.interaction import RLOFOutcome, rejuvenate_ms_gainer
from realta.binaries.population import BinaryPopulation
from realta.config import SimulationConfig
from realta.stellar import main_sequence as ms

Z_SOLAR = 0.02


def test_radiative_core_reduces_to_simple_fractional_age_preservation():
    """For 0.3 < M < 1.3 Msun, eq. 41's mass-ratio factor collapses to
    1, reducing to the same simple fractional-age-preservation formula
    already verified against HTP02/Paper I Sec. 7.1 (t' = (t'_MS/t_MS)*t)."""
    mass_before, mass_after = 1.0, 1.05
    t_ms_before = ms.t_ms(mass_before, Z_SOLAR)
    age_before = 0.5 * t_ms_before

    frac = rejuvenate_ms_gainer(mass_before, mass_after, age_before, Z_SOLAR)

    t_ms_after = ms.t_ms(mass_after, Z_SOLAR)
    age_after_simple = (t_ms_after / t_ms_before) * age_before
    expected = 1.0 - age_after_simple / t_ms_after
    assert frac == pytest.approx(expected, rel=1e-9)


def test_convective_core_gives_extra_rejuvenation_beyond_simple_case():
    """For M > 1.3 Msun, eq. 41's surviving M/M' factor (< 1 for a
    gainer) makes the star appear younger STILL than fractional-age
    preservation alone -- the "mixes in unburnt fuel" effect Tout et
    al. describe. Direction check: convective-case remaining fraction
    must exceed what the simple (radiative-only) formula would give
    for the same masses/age."""
    mass_before, mass_after = 5.0, 6.0
    t_ms_before = ms.t_ms(mass_before, Z_SOLAR)
    age_before = 0.5 * t_ms_before

    frac = rejuvenate_ms_gainer(mass_before, mass_after, age_before, Z_SOLAR)

    t_ms_after = ms.t_ms(mass_after, Z_SOLAR)
    age_after_simple = (t_ms_after / t_ms_before) * age_before
    simple_only_frac = 1.0 - age_after_simple / t_ms_after
    assert frac > simple_only_frac


def test_no_mass_change_matches_ordinary_unrejuvenated_aging():
    """If mass_before == mass_after (no actual gain), the formula must
    reduce to plain, un-rejuvenated fractional aging -- a sanity check
    that eq. 41 is a genuine no-op when there is nothing to rejuvenate."""
    mass = 5.0
    t_ms_val = ms.t_ms(mass, Z_SOLAR)
    age_before = 0.3 * t_ms_val

    frac = rejuvenate_ms_gainer(mass, mass, age_before, Z_SOLAR)
    expected = 1.0 - age_before / t_ms_val
    assert frac == pytest.approx(expected, rel=1e-9)


def test_remaining_fraction_is_clamped_to_positive_floor():
    """An old age already close to/exceeding the new (shorter) MS
    lifetime must not give a negative or zero remaining fraction."""
    mass_before, mass_after = 5.0, 5.001  # negligible mass gain
    t_ms_before = ms.t_ms(mass_before, Z_SOLAR)
    age_before = t_ms_before * 0.999999  # essentially at end of life

    frac = rejuvenate_ms_gainer(mass_before, mass_after, age_before, Z_SOLAR)
    assert 0.0 < frac <= 1.0


def test_rejuvenation_direction_more_mass_gained_is_more_rejuvenated():
    """Sensitivity/direction check: for a fixed convective-core donor
    mass and age, a LARGER mass gain should rejuvenate the star more
    (larger remaining fraction) than a smaller one."""
    mass_before = 5.0
    age_before = 0.5 * ms.t_ms(mass_before, Z_SOLAR)

    frac_small_gain = rejuvenate_ms_gainer(mass_before, 5.5, age_before, Z_SOLAR)
    frac_large_gain = rejuvenate_ms_gainer(mass_before, 8.0, age_before, Z_SOLAR)
    assert frac_large_gain > frac_small_gain


def _hand_constructed_stable_mt_population(
    config, donor_mass, companion_mass, rlof_time, separation_rsun=8.0
):
    pop = BinaryPopulation(config)
    pop.m1 = np.array([donor_mass])
    pop.m2 = np.array([companion_mass])
    pop.a = np.array([separation_rsun / BinaryPopulation.RSUN_PER_AU])
    pop.period = np.array([100.0])
    pop.turnoff_time = np.array([1.0e6])
    pop.t2_lifetime = np.array([1.0e6])
    pop.nturn = np.zeros(1, dtype=np.int8)
    pop.is_survived = np.ones(1, dtype=bool)
    pop.lum_xray = np.zeros(1)
    pop.did_merge = np.zeros(1, dtype=bool)
    pop.merge_time = np.array([np.nan])
    pop.rlof_time = np.array([rlof_time])
    pop.rlof_outcome = np.array([RLOFOutcome.STABLE_MASS_TRANSFER], dtype=object)
    pop.rlof_donor_is_star1 = np.array([True])
    pop.rlof_processed = np.zeros(1, dtype=bool)
    return pop


def test_evolve_rejuvenates_ms_companion_instead_of_full_reset():
    """Wiring test: when the companion is genuinely MS-phase at
    rlof_time, evolve() must use rejuvenate_ms_gainer's remaining
    fraction (against LifetimeTable's own lifetime at the new mass),
    not the old full-reset simplification.

    Hand-constructs the scenario directly (bypassing find_rlof_onset's
    own donor selection) because a "companion still on its own MS at a
    non-trivial donor-discovered rlof_time" combination is essentially
    unreachable via find_rlof_onset in practice for STABLE_MASS_TRANSFER
    (donor lighter than companion): by the time a lighter donor's own
    much longer MS+HG evolution reaches RLOF, the heavier companion has
    almost always already died by its own (shorter) lifetime -- the
    same structural finding already documented in docs/provenance.md
    for why STABLE_MASS_TRANSFER rarely reaches evolve() at all. This
    test verifies the rejuvenation code path itself is correct;
    whether it is commonly reached in a real population is a separate,
    already-documented question.
    """
    config = SimulationConfig(
        ntot=10, imetal=2, use_rlof_classifier=True, fsur=1.0, iseed=1
    )
    # Masses/separation confirmed (via find_rlof_onset) to be within
    # the "born overflowing" zone for this pair -- donor=5 crosses at
    # t=0.0 exactly for a<=7 Rsun; rlof_time=0.05 here is a small,
    # physically-motivated nudge past that (the system is already
    # overflowing shortly after formation), not an arbitrary pick --
    # a genuinely later, natural crossing with the companion still
    # MS-phase does not exist for this outcome (see this test's own
    # docstring).
    donor_mass, companion_mass, rlof_time = 5.0, 7.246376811594203, 0.05
    # Confirm the companion really is MS-phase at this age first.
    assert ms.phase(companion_mass, 0.008, rlof_time) in (0, 1)

    pop = _hand_constructed_stable_mt_population(
        config, donor_mass, companion_mass, rlof_time, separation_rsun=5.0
    )
    pop.evolve(tnow=rlof_time, dt=1.0)

    companion_new_mass = pop.m2[0]
    assert companion_new_mass > companion_mass  # confirms it gained mass

    expected_fraction = rejuvenate_ms_gainer(
        companion_mass, companion_new_mass, rlof_time, 0.008
    )
    expected_t2 = rlof_time + expected_fraction * pop.lifetime_table.get_lifetime(
        companion_new_mass
    )
    assert pop.t2_lifetime[0] == pytest.approx(expected_t2, rel=1e-9)

    full_reset_would_be = rlof_time + pop.lifetime_table.get_lifetime(
        companion_new_mass
    )
    assert pop.t2_lifetime[0] != pytest.approx(full_reset_would_be)


def test_evolve_falls_back_to_full_reset_for_non_ms_companion():
    """When the companion is NOT MS-phase at rlof_time (e.g. already
    past its own t_BGB -- classify_rlof() places no phase constraint
    on the companion for STABLE_MASS_TRANSFER, so this is a real case),
    evolve() must fall back to the pre-existing full-reset
    simplification rather than misapplying the MS-only rejuvenation
    formula.
    """
    config = SimulationConfig(
        ntot=10, imetal=2, use_rlof_classifier=True, fsur=1.0, iseed=1
    )
    donor_mass, companion_mass, rlof_time = 2.0, 10.0, 875.0
    with pytest.raises(ValueError):
        ms.phase(companion_mass, 0.008, rlof_time)  # confirms it's past t_BGB

    pop = _hand_constructed_stable_mt_population(
        config, donor_mass, companion_mass, rlof_time
    )
    pop.evolve(tnow=rlof_time, dt=1.0)

    companion_new_mass = pop.m2[0]
    full_reset_expected = rlof_time + pop.lifetime_table.get_lifetime(
        companion_new_mass
    )
    assert pop.t2_lifetime[0] == pytest.approx(full_reset_expected, rel=1e-9)


def test_evolve_sensitivity_rejuvenation_actually_differs_from_full_reset():
    """Direct sensitivity check requested by this session's discipline:
    confirms the rejuvenated companion's new lifetime is measurably
    different from (and, for a genuinely young companion, later than)
    what the old full-reset code would have given -- i.e. this is a
    real behavioural change, not a no-op wearing new code."""
    config = SimulationConfig(
        ntot=10, imetal=2, use_rlof_classifier=True, fsur=1.0, iseed=1
    )
    donor_mass, companion_mass, rlof_time = 5.0, 7.246376811594203, 0.05

    pop = _hand_constructed_stable_mt_population(
        config, donor_mass, companion_mass, rlof_time, separation_rsun=5.0
    )
    pop.evolve(tnow=rlof_time, dt=1.0)

    companion_new_mass = pop.m2[0]
    full_reset_value = rlof_time + pop.lifetime_table.get_lifetime(companion_new_mass)
    assert pop.t2_lifetime[0] < full_reset_value
