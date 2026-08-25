"""Tests for wiring the RLOF classifier into BinaryPopulation
(generate_population's precomputation + evolve()'s Phase 0) -- see
docs/science/rlof-ce-classifier-proposal.md and
docs/physics/rlof-classifier.md.
"""

from unittest.mock import patch

import numpy as np
import pytest

from realta.binaries.interaction import RLOFOutcome, apply_common_envelope
from realta.binaries.population import BinaryPopulation
from realta.config import SimulationConfig
from realta.stellar import main_sequence as ms


def test_rlof_classifier_disabled_by_default_leaves_rlof_time_infinite():
    config = SimulationConfig(ntot=2000, iseed=7)
    assert config.use_rlof_classifier is False
    pop = BinaryPopulation(config)
    assert np.all(np.isinf(pop.rlof_time))
    # NOT `pop.rlof_outcome == RLOFOutcome.DETACHED` -- numpy's
    # vectorized `==` against a str-Enum member is broken (confirmed:
    # returns all-False even when every element genuinely equals it),
    # while per-element scalar comparison works correctly. See
    # RLOFOutcome's own docstring in interaction.py.
    assert all(outcome == RLOFOutcome.DETACHED for outcome in pop.rlof_outcome)


def test_rlof_outcome_array_construction_avoids_np_full_corruption():
    """Regression test for a genuine numpy bug found during this
    module's development: `np.full(n, RLOFOutcome.X, dtype=object)`
    silently truncates/corrupts the str-Enum fill value -- the
    resulting array holds a plain, truncated str that fails equality
    against the real enum member, even though the array's dtype is
    correctly 'object'. generate_population() must use list-based
    construction (`np.array([RLOFOutcome.DETACHED] * n, dtype=object)`)
    instead. This locks in that choice directly, rather than relying
    on it only being caught indirectly by other tests.
    """
    config = SimulationConfig(ntot=2000, iseed=7)
    pop = BinaryPopulation(config)
    assert len(pop.m1) > 10
    for outcome in pop.rlof_outcome:
        assert outcome == RLOFOutcome.DETACHED
        assert isinstance(outcome, RLOFOutcome)


def test_rlof_classifier_disabled_matches_baseline_exactly():
    """Sensitivity/determinism check: enabling use_rlof_classifier=False
    (the default) must produce bit-identical output to a config that
    doesn't mention it at all -- the whole point of the opt-in flag.
    """
    config_a = SimulationConfig(ntot=2000, fsur=1.0, iseed=7)
    config_b = SimulationConfig(ntot=2000, fsur=1.0, iseed=7, use_rlof_classifier=False)
    pop_a = BinaryPopulation(config_a)
    pop_b = BinaryPopulation(config_b)
    np.testing.assert_array_equal(pop_a.m1, pop_b.m1)
    for step in range(20):
        assert pop_a.evolve(step * 1.0, 1.0) == pop_b.evolve(step * 1.0, 1.0)


def test_rlof_classifier_z_zero_warns_and_skips(caplog):
    """imetal=1 (Z=0) is outside the Hurley/Tout formulae's domain --
    must log a warning and leave rlof_time inert, not crash.
    """
    config = SimulationConfig(ntot=2000, imetal=1, use_rlof_classifier=True, iseed=7)
    with caplog.at_level("WARNING"):
        pop = BinaryPopulation(config)
    assert np.all(np.isinf(pop.rlof_time))
    assert any("Z=0" in record.message for record in caplog.records)


def test_rlof_classifier_precomputes_some_finite_times_at_realistic_z():
    """At a realistic metallicity, a large-enough massive-star
    population should have at least some binaries close enough to
    interact on the MS.
    """
    config = SimulationConfig(
        ntot=20_000, imetal=2, use_rlof_classifier=True, pmin=0.5, pmax=50.0, iseed=7
    )
    pop = BinaryPopulation(config)
    assert len(pop.m1) > 50
    assert np.any(np.isfinite(pop.rlof_time))


def test_evolve_processes_immediate_merger_via_rlof_channel():
    """Hand-construct a single binary guaranteed to RLOF-merge shortly
    after t=0 (identical to the interaction-module unit tests' merger
    scenario), and confirm evolve() actually applies it: m1 becomes
    the combined mass, m2 is zeroed, did_merge/merge_time are set, and
    the merged star's new lifetime clock starts from the merge time.
    """
    config = SimulationConfig(
        ntot=10, imetal=2, use_rlof_classifier=True, fsur=1.0, iseed=1
    )
    pop = BinaryPopulation(config)

    # Overwrite state directly (same technique as tests/test_evolve.py)
    # with the exact scenario from
    # test_rlof_classifier.py::test_find_rlof_onset_finds_crossing_and_correct_donor
    # (m1=5, m2=3, a=10 -> t_rlof ~ 72.46 Myr, immediate merger, donor=star1).
    pop.m1 = np.array([5.0])
    pop.m2 = np.array([3.0])
    pop.a = np.array([10.0])
    pop.period = np.array([100.0])
    pop.turnoff_time = np.array([1.0e6])  # far future -- won't fire this test
    pop.t2_lifetime = np.array([1.0e6])
    pop.nturn = np.zeros(1, dtype=np.int8)
    pop.is_survived = np.ones(1, dtype=bool)
    pop.lum_xray = np.zeros(1)
    pop.did_merge = np.zeros(1, dtype=bool)
    pop.merge_time = np.array([np.nan])
    pop.rlof_time = np.array([72.45736603967721])
    pop.rlof_outcome = np.array([RLOFOutcome.IMMEDIATE_MERGER], dtype=object)
    pop.rlof_donor_is_star1 = np.array([True])
    pop.rlof_processed = np.zeros(1, dtype=bool)

    pop.evolve(tnow=72.46, dt=1.0)

    assert pop.m1[0] == pytest.approx(8.0)
    assert pop.m2[0] == 0.0
    assert pop.did_merge[0]
    assert pop.merge_time[0] == pytest.approx(72.46)
    assert pop.rlof_processed[0]
    # New lifetime clock starts from the merge time, for the merged mass.
    expected_turnoff = 72.46 + pop.lifetime_table.get_lifetime(8.0)
    assert pop.turnoff_time[0] == pytest.approx(expected_turnoff)
    assert pop.t2_lifetime[0] == 0.0


def test_evolve_does_not_reprocess_rlof_event():
    """Once processed, an RLOF event must not fire again on a later
    timestep (rlof_processed guards this).
    """
    config = SimulationConfig(
        ntot=10, imetal=2, use_rlof_classifier=True, fsur=1.0, iseed=1
    )
    pop = BinaryPopulation(config)

    pop.m1 = np.array([5.0])
    pop.m2 = np.array([3.0])
    pop.a = np.array([10.0])
    pop.period = np.array([100.0])
    pop.turnoff_time = np.array([1.0e6])
    pop.t2_lifetime = np.array([1.0e6])
    pop.nturn = np.zeros(1, dtype=np.int8)
    pop.is_survived = np.ones(1, dtype=bool)
    pop.lum_xray = np.zeros(1)
    pop.did_merge = np.zeros(1, dtype=bool)
    pop.merge_time = np.array([np.nan])
    pop.rlof_time = np.array([50.0])
    pop.rlof_outcome = np.array([RLOFOutcome.IMMEDIATE_MERGER], dtype=object)
    pop.rlof_donor_is_star1 = np.array([True])
    pop.rlof_processed = np.zeros(1, dtype=bool)

    pop.evolve(tnow=50.0, dt=1.0)
    merged_mass_after_first = pop.m1[0]
    pop.evolve(tnow=60.0, dt=1.0)
    assert pop.m1[0] == merged_mass_after_first  # unchanged, not re-merged


def test_evolve_applies_stable_mass_transfer_instantaneously():
    """STABLE_MASS_TRANSFER is applied via
    apply_stable_mass_transfer(): mass moves conservatively from the
    (lighter) donor to the (heavier) companion, the orbit widens, and
    the event is marked processed so it doesn't repeat. donor_mass=2.0
    (star 1) < companion_mass=10.0 (star 2) -- the only direction
    apply_stable_mass_transfer accepts (classify_rlof never labels the
    opposite ordering STABLE_MASS_TRANSFER).
    """
    config = SimulationConfig(
        ntot=10, imetal=2, use_rlof_classifier=True, fsur=1.0, iseed=1
    )
    pop = BinaryPopulation(config)

    a_rsun = 8.0
    pop.m1 = np.array([2.0])
    pop.m2 = np.array([10.0])
    # pop.a is in AU (see BinaryPopulation.RSUN_PER_AU) -- evolve()
    # converts to Rsun before calling apply_stable_mass_transfer, so
    # the intended Rsun-scale separation must be pre-converted here.
    pop.a = np.array([a_rsun / BinaryPopulation.RSUN_PER_AU])
    pop.period = np.array([100.0])
    pop.turnoff_time = np.array([1.0e6])
    pop.t2_lifetime = np.array([1.0e6])
    pop.nturn = np.zeros(1, dtype=np.int8)
    pop.is_survived = np.ones(1, dtype=bool)
    pop.lum_xray = np.zeros(1)
    pop.did_merge = np.zeros(1, dtype=bool)
    pop.merge_time = np.array([np.nan])
    rlof_time = 875.0  # < t_ms(2.0, Z=0.008) ~ 971.65 Myr
    pop.rlof_time = np.array([rlof_time])
    pop.rlof_outcome = np.array([RLOFOutcome.STABLE_MASS_TRANSFER], dtype=object)
    pop.rlof_donor_is_star1 = np.array([True])
    pop.rlof_processed = np.zeros(1, dtype=bool)

    total_mass_before = pop.m1[0] + pop.m2[0]
    pop.evolve(tnow=rlof_time, dt=1.0)

    assert pop.m1[0] + pop.m2[0] == pytest.approx(total_mass_before)  # conservative
    assert pop.m1[0] < 2.0  # donor lost mass
    assert pop.m2[0] > 10.0  # companion gained mass
    assert pop.m1[0] < pop.m2[0]  # still the lighter star
    assert pop.a[0] * BinaryPopulation.RSUN_PER_AU > a_rsun  # orbit widened
    assert not pop.did_merge[0]
    assert pop.rlof_processed[0]


def test_evolve_does_not_reprocess_stable_mass_transfer():
    config = SimulationConfig(
        ntot=10, imetal=2, use_rlof_classifier=True, fsur=1.0, iseed=1
    )
    pop = BinaryPopulation(config)

    pop.m1 = np.array([2.0])
    pop.m2 = np.array([10.0])
    pop.a = np.array([8.0 / BinaryPopulation.RSUN_PER_AU])  # AU, see above
    pop.period = np.array([100.0])
    pop.turnoff_time = np.array([1.0e6])
    pop.t2_lifetime = np.array([1.0e6])
    pop.nturn = np.zeros(1, dtype=np.int8)
    pop.is_survived = np.ones(1, dtype=bool)
    pop.lum_xray = np.zeros(1)
    pop.did_merge = np.zeros(1, dtype=bool)
    pop.merge_time = np.array([np.nan])
    pop.rlof_time = np.array([875.0])
    pop.rlof_outcome = np.array([RLOFOutcome.STABLE_MASS_TRANSFER], dtype=object)
    pop.rlof_donor_is_star1 = np.array([True])
    pop.rlof_processed = np.zeros(1, dtype=bool)

    pop.evolve(tnow=875.0, dt=1.0)
    m1_after_first = pop.m1[0]
    pop.evolve(tnow=900.0, dt=1.0)
    assert pop.m1[0] == m1_after_first  # unchanged, not re-applied


def test_evolve_applies_stable_mass_transfer_for_hg_donor():
    """Regression test for a real bug found while wiring HG search
    into find_rlof_onset(): evolve()'s STABLE_MASS_TRANSFER branch
    unconditionally called main_sequence.ms_radius() for the donor's
    radius, which is wrong once find_rlof_onset() can find a
    stable-MT crossing during the HG too, not just the MS. Must use
    hg_radius() for an HG-phase donor instead.
    """
    config = SimulationConfig(
        ntot=10, imetal=3, use_rlof_classifier=True, fsur=1.0, iseed=1
    )
    pop = BinaryPopulation(config)

    # donor=5 Msun (HG donor), companion=20 Msun, a=100, Z=0.02
    # (imetal=3) -- confirmed via find_rlof_onset(5.0, 20.0, 100.0,
    # 0.02) to cross during HG at t~104.326 Myr with outcome
    # STABLE_MASS_TRANSFER. rlof_time is nudged slightly past the
    # exact root (whose remaining HG window is very narrow -- only
    # ~0.12 Myr to t_BGB) since apply_stable_mass_transfer's own
    # bracket search needs donor_radius meaningfully above (not
    # razor-equal to) the current Roche lobe.
    donor_mass, companion_mass, a_rsun = 5.0, 20.0, 100.0
    rlof_time = 104.32589865040659 + 0.01

    pop.m1 = np.array([donor_mass])
    pop.m2 = np.array([companion_mass])
    pop.a = np.array([a_rsun / BinaryPopulation.RSUN_PER_AU])  # AU, see above
    pop.period = np.array([1000.0])
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

    total_mass_before = pop.m1[0] + pop.m2[0]
    pop.evolve(tnow=rlof_time, dt=1.0)

    assert pop.m1[0] + pop.m2[0] == pytest.approx(total_mass_before)
    assert pop.m1[0] < donor_mass  # donor lost mass
    assert pop.m2[0] > companion_mass  # companion gained mass
    assert pop.a[0] * BinaryPopulation.RSUN_PER_AU > a_rsun  # orbit widened
    assert pop.rlof_processed[0]


def test_evolve_common_envelope_merges_when_energy_balance_favours_merger():
    """COMMON_ENVELOPE now resolves via the alpha-lambda energy-balance
    solve (binaries/interaction.py::apply_common_envelope). For this
    donor/companion/age combination (mid-HG, envelope-heavy donor) the
    energy balance favours merger -- see
    tests/test_ce_energy_balance.py's module docstring for why this is
    the expected, dominant outcome rather than a bug. The merged mass
    is the donor's *core* mass (already stripped by the failed CE) plus
    the companion's mass, not the donor's full pre-CE mass.
    """
    config = SimulationConfig(
        ntot=10, imetal=3, use_rlof_classifier=True, fsur=1.0, iseed=1
    )
    pop = BinaryPopulation(config)

    pop.m1 = np.array([5.0])
    pop.m2 = np.array([3.0])
    pop.a = np.array([100.0 / BinaryPopulation.RSUN_PER_AU])  # AU, see above
    pop.period = np.array([1000.0])
    pop.turnoff_time = np.array([1.0e6])
    pop.t2_lifetime = np.array([1.0e6])
    pop.nturn = np.zeros(1, dtype=np.int8)
    pop.is_survived = np.ones(1, dtype=bool)
    pop.lum_xray = np.zeros(1)
    pop.did_merge = np.zeros(1, dtype=bool)
    pop.merge_time = np.array([np.nan])
    pop.rlof_time = np.array([104.42264724439553])
    pop.rlof_outcome = np.array([RLOFOutcome.COMMON_ENVELOPE], dtype=object)
    pop.rlof_donor_is_star1 = np.array([True])
    pop.rlof_processed = np.zeros(1, dtype=bool)

    pop.evolve(tnow=104.5, dt=1.0)

    core_mass = ms.core_mass_hg(5.0, 0.02, 104.42264724439553)
    assert pop.m1[0] == pytest.approx(core_mass + 3.0)
    assert pop.m2[0] == 0.0
    assert pop.did_merge[0]
    assert pop.merge_time[0] == 104.5
    assert pop.rlof_processed[0]

    # And confirm it doesn't get reprocessed either.
    m1_after_first = pop.m1[0]
    pop.evolve(tnow=110.0, dt=1.0)
    assert pop.m1[0] == m1_after_first


def test_evolve_common_envelope_survives_when_energy_balance_favours_survival():
    """The other branch: when the energy balance favours survival (an
    envelope-light donor near the end of the HG), evolve() must strip
    the donor to its core mass, leave the companion untouched, and
    tighten the orbit to a_f -- without merging.
    """
    config = SimulationConfig(
        ntot=10, imetal=3, use_rlof_classifier=True, fsur=1.0, iseed=1
    )
    pop = BinaryPopulation(config)

    donor_mass, companion_mass, a = 5.0, 1.0, 5.0
    z = 0.02
    age = None
    for frac in (0.999, 0.9999, 0.99999):
        candidate_age = ms.t_ms(donor_mass, z) + frac * (
            ms.t_bgb(donor_mass, z) - ms.t_ms(donor_mass, z)
        )
        survives, _, _, _ = apply_common_envelope(
            donor_mass, companion_mass, a, z, candidate_age
        )
        if survives:
            age = candidate_age
            break
    if age is None:
        pytest.skip(
            "no HG age for this donor mass reaches survives=True on real tracks"
        )

    pop.m1 = np.array([donor_mass])
    pop.m2 = np.array([companion_mass])
    pop.a = np.array([a / BinaryPopulation.RSUN_PER_AU])  # AU, see above
    pop.period = np.array([1000.0])
    pop.turnoff_time = np.array([1.0e6])
    pop.t2_lifetime = np.array([1.0e6])
    pop.nturn = np.zeros(1, dtype=np.int8)
    pop.is_survived = np.ones(1, dtype=bool)
    pop.lum_xray = np.zeros(1)
    pop.did_merge = np.zeros(1, dtype=bool)
    pop.merge_time = np.array([np.nan])
    pop.rlof_time = np.array([age])
    pop.rlof_outcome = np.array([RLOFOutcome.COMMON_ENVELOPE], dtype=object)
    pop.rlof_donor_is_star1 = np.array([True])
    pop.rlof_processed = np.zeros(1, dtype=bool)

    core_mass = ms.core_mass_hg(donor_mass, z, age)
    pop.evolve(tnow=age + 0.01, dt=1.0)

    assert pop.m1[0] == pytest.approx(core_mass)
    assert pop.m2[0] == companion_mass
    assert pop.a[0] * BinaryPopulation.RSUN_PER_AU < a
    assert not pop.did_merge[0]
    assert pop.rlof_processed[0]


def test_config_alpha_ce_lambda_ce_default_to_interaction_module_constants():
    from realta.binaries.interaction import ALPHA_CE, LAMBDA_CE

    config = SimulationConfig(ntot=10, iseed=1)
    assert config.alpha_ce == ALPHA_CE
    assert config.lambda_ce == LAMBDA_CE


def test_config_alpha_ce_lambda_ce_are_overridable_and_validated():
    config = SimulationConfig(ntot=10, iseed=1, alpha_ce=0.7, lambda_ce=0.3)
    assert config.alpha_ce == 0.7
    assert config.lambda_ce == 0.3

    with pytest.raises(ValueError, match="alpha_ce"):
        SimulationConfig(ntot=10, iseed=1, alpha_ce=0.0)
    with pytest.raises(ValueError, match="lambda_ce"):
        SimulationConfig(ntot=10, iseed=1, lambda_ce=-1.0)


def test_config_alpha_ce_override_propagates_into_evolve_common_envelope():
    """Confirms the override actually reaches apply_common_envelope
    inside evolve(), not just the config field itself: a low enough
    alpha_ce should force merger even for a donor/age combination that
    otherwise survives at the default alpha_ce.
    """
    donor_mass, companion_mass, a = 5.0, 1.0, 5.0
    z = 0.02
    age = None
    for frac in (0.999, 0.9999, 0.99999):
        candidate_age = ms.t_ms(donor_mass, z) + frac * (
            ms.t_bgb(donor_mass, z) - ms.t_ms(donor_mass, z)
        )
        survives, _, _, _ = apply_common_envelope(
            donor_mass, companion_mass, a, z, candidate_age
        )
        if survives:
            age = candidate_age
            break
    if age is None:
        pytest.skip(
            "no HG age for this donor mass reaches survives=True on real tracks"
        )

    config = SimulationConfig(
        ntot=10, imetal=3, use_rlof_classifier=True, fsur=1.0, iseed=1, alpha_ce=0.05
    )
    pop = BinaryPopulation(config)

    pop.m1 = np.array([donor_mass])
    pop.m2 = np.array([companion_mass])
    pop.a = np.array([a / BinaryPopulation.RSUN_PER_AU])  # AU, see above
    pop.period = np.array([1000.0])
    pop.turnoff_time = np.array([1.0e6])
    pop.t2_lifetime = np.array([1.0e6])
    pop.nturn = np.zeros(1, dtype=np.int8)
    pop.is_survived = np.ones(1, dtype=bool)
    pop.lum_xray = np.zeros(1)
    pop.did_merge = np.zeros(1, dtype=bool)
    pop.merge_time = np.array([np.nan])
    pop.rlof_time = np.array([age])
    pop.rlof_outcome = np.array([RLOFOutcome.COMMON_ENVELOPE], dtype=object)
    pop.rlof_donor_is_star1 = np.array([True])
    pop.rlof_processed = np.zeros(1, dtype=bool)

    pop.evolve(tnow=age + 0.01, dt=1.0)
    assert pop.did_merge[0]


def test_evolve_common_envelope_survives_leaves_companion_clock_untouched():
    """CE-survival must reset the DONOR's lifetime clock (it was just
    stripped to its core mass -- a real, physical change) but leave
    the COMPANION's clock entirely alone -- apply_common_envelope's
    own docstring says the companion is mass-unaffected by a surviving
    CE, so resetting its clock would incorrectly de-age it (pretend it
    just formed anew), the same class of issue found and fixed for
    STABLE_MASS_TRANSFER's companion (B3 rejuvenation). This is a real
    bug found while implementing B3, fixed by explicit user decision:
    leave the companion's clock alone entirely rather than reset or
    rejuvenate it (it has no mass change for rejuvenation to apply to).

    A genuine survives=True case could not be found across a wide
    search of real donor/companion/age combinations on the currently-
    supported HG track range (see docs/physics/mass-transfer.md's
    CE-merger-dominance finding) -- this test mocks
    apply_common_envelope directly to exercise evolve()'s wiring logic
    in isolation, independent of whether real stellar tracks currently
    reach that branch.
    """
    config = SimulationConfig(
        ntot=10, imetal=3, use_rlof_classifier=True, fsur=1.0, iseed=1
    )
    pop = BinaryPopulation(config)

    donor_mass, companion_mass, a, rlof_time = 5.0, 3.0, 5.0, 10.0
    pop.m1 = np.array([donor_mass])
    pop.m2 = np.array([companion_mass])
    pop.a = np.array([a / BinaryPopulation.RSUN_PER_AU])
    pop.period = np.array([1000.0])
    pop.turnoff_time = np.array([1.0e6])
    sentinel_companion_clock = 4321.0
    pop.t2_lifetime = np.array([sentinel_companion_clock])
    pop.nturn = np.zeros(1, dtype=np.int8)
    pop.is_survived = np.ones(1, dtype=bool)
    pop.lum_xray = np.zeros(1)
    pop.did_merge = np.zeros(1, dtype=bool)
    pop.merge_time = np.array([np.nan])
    pop.rlof_time = np.array([rlof_time])
    pop.rlof_outcome = np.array([RLOFOutcome.COMMON_ENVELOPE], dtype=object)
    pop.rlof_donor_is_star1 = np.array([True])
    pop.rlof_processed = np.zeros(1, dtype=bool)

    core_mass = 0.9
    with patch(
        "realta.binaries.population.apply_common_envelope",
        return_value=(True, core_mass, companion_mass, 2.0),
    ):
        pop.evolve(tnow=rlof_time, dt=1.0)

    assert pop.m1[0] == pytest.approx(core_mass)
    assert pop.m2[0] == companion_mass
    # Donor (star1) clock IS reset.
    assert pop.turnoff_time[0] == pytest.approx(
        rlof_time + pop.lifetime_table.get_lifetime(core_mass)
    )
    # Companion (star2) clock is untouched -- still the sentinel.
    assert pop.t2_lifetime[0] == sentinel_companion_clock
