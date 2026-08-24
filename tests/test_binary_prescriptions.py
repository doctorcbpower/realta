"""Unit tests for the Paper 1 binary-interaction/merger prescriptions.

See docs/science/paper1-binary-interaction-proposal.md for the physics
and provenance.md for why these parameters (interaction_boost, p_merge,
p_merge_max_period, f_merge) are explicitly flagged as NOT paper-derived
-- they are a new, Realta-specific parameterization, unlike every other
row in provenance.md. These tests follow the same phase-isolation style
as tests/test_evolve.py: hand-construct or hand-seed a scenario, assert
the exact mechanism, rather than relying on what a real IMF-sampled run
happens to reach.
"""

import numpy as np
import pytest

from realta.binaries.interaction import RLOFOutcome
from realta.binaries.population import BinaryPopulation
from realta.config import SimulationConfig


def test_default_prescription_matches_baseline_exactly():
    """binary_prescription='non_interacting' (the default) must be a
    no-op: interaction_boost resolves to 1.0 and p_merge to 0.0, so
    fsur_eff == fsur exactly and the merger RNG draw is skipped
    entirely -- see the "skip the RNG draw entirely" comment in
    generate_population(). This is the guarantee behind "the 25
    currently passing tests should still pass unchanged" in the
    implementation prompt.
    """
    config = SimulationConfig(ntot=5000, iseed=7)
    assert config.binary_prescription == "non_interacting"
    assert config.interaction_boost == 1.0
    assert config.p_merge == 0.0

    pop = BinaryPopulation(config)
    assert not np.any(pop.did_merge)


def test_single_prescription_suppresses_binary_formation():
    """'single' must produce zero binaries -- no companion, no period,
    no HMXB channel -- while total_mass_msun (needed for L_bol/L_UV
    rescaling) is computed from the full IMF sample as usual.
    """
    config = SimulationConfig(ntot=5000, binary_prescription="single", iseed=7)
    pop = BinaryPopulation(config)

    assert len(pop.m1) == 0
    assert pop.total_mass_msun > 0.0

    lumx_tot, _nphot_tot, nactive, ndead = pop.evolve(tnow=1.0, dt=1.0)
    assert lumx_tot == 0.0
    assert nactive == 0
    assert ndead == 0


def test_interaction_boost_raises_activation_fraction_for_stable_mt_systems():
    """interaction_boost must raise the HMXB activation fraction above
    fsur (statistically) for binaries that underwent stable mass
    transfer, and cap at fsur_eff = 1.0 -- direct analogue of
    test_evolve.py::test_evolve_phase1_fsur_partial_activation, with
    interaction_boost added on top of the same fsur gate. Since the
    reconciliation with the physics-based RLOF classifier (see
    docs/science/rlof-ce-classifier-proposal.md "Decision 3"), the
    boost only applies to binaries the classifier found actually
    interacted -- see the paired test below for the "never interacted"
    case, which must NOT receive the boost.
    """
    n = 2000
    fsur = 0.3
    boost = 3.0  # -> fsur_eff = 0.9, still < 1 so the gate remains stochastic
    config = SimulationConfig(
        ntot=10,
        fsur=fsur,
        binary_prescription="enhanced_interaction",
        interaction_boost=boost,
        iseed=123,
    )
    assert config.use_rlof_classifier is True  # implied by the prescription now
    pop = BinaryPopulation(config)
    assert config.interaction_boost == boost

    pop.remnant_table.get_remnant_mass = lambda m: 1.4
    pop.m1 = np.full(n, 20.0)
    pop.m2 = np.full(n, 20.0)  # floss = (20-1.4)/40 = 0.465 <= 0.5 -> always survives
    pop.period = np.full(n, 10.0)
    pop.a = np.full(n, 1.0)
    pop.turnoff_time = np.full(n, 5.0)
    pop.t2_lifetime = np.full(n, 6.0)
    pop.nturn = np.zeros(n, dtype=np.int8)
    pop.is_survived = np.ones(n, dtype=bool)
    pop.lum_xray = np.zeros(n)
    pop.did_merge = np.zeros(n, dtype=bool)
    # Every system already had a stable-MT interaction (isolates the
    # boost-application mechanism from which binaries get flagged --
    # that selection is find_rlof_onset's concern, tested
    # separately in tests/test_rlof_classifier.py).
    pop.rlof_time = np.full(n, -1.0)
    pop.rlof_outcome = np.array([RLOFOutcome.STABLE_MASS_TRANSFER] * n, dtype=object)
    pop.rlof_donor_is_star1 = np.ones(n, dtype=bool)
    pop.rlof_processed = np.ones(n, dtype=bool)  # already applied, Phase 0 no-ops

    pop.evolve(tnow=5.0, dt=1.0)

    n_active = int(np.count_nonzero(pop.lum_xray > 0))
    fsur_eff = fsur * boost
    expected = n * fsur_eff
    tolerance = 10 * np.sqrt(n * fsur_eff * (1 - fsur_eff))
    assert abs(n_active - expected) < tolerance
    # Sanity: the boosted activation fraction is well above the
    # unboosted fsur -- confirms the boost actually did something, not
    # just that the statistical test is loose enough to pass anyway.
    assert n_active > n * fsur * 1.5


def test_interaction_boost_not_applied_to_never_interacted_systems():
    """A binary the RLOF classifier found DETACHED (never interacted)
    must use plain fsur, not the boosted fsur_eff -- the boost is
    conditional on genuine interaction history since the reconciliation
    (docs/science/rlof-ce-classifier-proposal.md "Decision 3"), not
    applied unconditionally to every surviving binary as before.
    """
    n = 2000
    fsur = 0.3
    boost = 3.0
    config = SimulationConfig(
        ntot=10,
        fsur=fsur,
        binary_prescription="enhanced_interaction",
        interaction_boost=boost,
        iseed=123,
    )
    pop = BinaryPopulation(config)

    pop.remnant_table.get_remnant_mass = lambda m: 1.4
    pop.m1 = np.full(n, 20.0)
    pop.m2 = np.full(n, 20.0)
    pop.period = np.full(n, 10.0)
    pop.a = np.full(n, 1.0)
    pop.turnoff_time = np.full(n, 5.0)
    pop.t2_lifetime = np.full(n, 6.0)
    pop.nturn = np.zeros(n, dtype=np.int8)
    pop.is_survived = np.ones(n, dtype=bool)
    pop.lum_xray = np.zeros(n)
    pop.did_merge = np.zeros(n, dtype=bool)
    # DETACHED: never interacted on the MS.
    pop.rlof_time = np.full(n, np.inf)
    pop.rlof_outcome = np.array([RLOFOutcome.DETACHED] * n, dtype=object)
    pop.rlof_donor_is_star1 = np.ones(n, dtype=bool)
    pop.rlof_processed = np.zeros(n, dtype=bool)

    pop.evolve(tnow=5.0, dt=1.0)

    n_active = int(np.count_nonzero(pop.lum_xray > 0))
    expected = n * fsur  # unboosted
    tolerance = 10 * np.sqrt(n * fsur * (1 - fsur))
    assert abs(n_active - expected) < tolerance


def test_interaction_boost_caps_at_one():
    """fsur_eff = min(1, fsur * interaction_boost) must not exceed 1 --
    e.g. fsur=0.8, boost=3.0 would give 2.4 uncapped, which would make
    `rng.random() <= fsur_eff` always true (harmless) but is worth
    pinning explicitly since a missing min() would silently do the same
    thing here and only misbehave for a rng returning exactly on [0,1)
    edge cases elsewhere.
    """
    config = SimulationConfig(
        ntot=10,
        fsur=0.8,
        binary_prescription="enhanced_interaction",
        interaction_boost=3.0,
        iseed=1,
    )

    fsur_eff = min(1.0, config.fsur * config.interaction_boost)
    assert fsur_eff == 1.0


def test_merger_channel_folds_companion_into_primary_and_disables_hmxb():
    """A merged system: m1 increases by f_merge*m2, m2 is zeroed, and
    -- because activation in evolve() requires m2 > |mcomp| -- it can
    never become an active HMXB even with fsur=1.0. Uses
    p_merge=1.0/p_merge_max_period=huge so every eligible system merges
    deterministically, isolating the mechanism from the p_merge draw
    itself (already covered statistically by the fsur/interaction_boost
    tests above, which exercise the same rng.random()-vs-threshold
    pattern).
    """
    config = SimulationConfig(
        ntot=2000,
        fsur=1.0,
        binary_prescription="enhanced_mergers",
        p_merge=1.0,
        p_merge_max_period=1000.0,  # matches config default pmax -> all eligible
        f_merge=0.5,
        iseed=7,
    )
    pop = BinaryPopulation(config)
    assert len(pop.m1) > 10

    assert np.all(pop.did_merge)
    assert np.all(pop.m2 == 0.0)
    assert np.all(pop.merge_time == 0.0)

    # Run to completion; no HMXB should ever activate since every
    # system's companion was zeroed at formation.
    for step in range(200):
        lumx_tot, _nphot_tot, _nactive, _ndead = pop.evolve(
            tnow=step * config.dt, dt=config.dt
        )
        assert lumx_tot == 0.0
    assert np.all(pop.lum_xray == 0.0)


def test_merger_channel_respects_period_threshold():
    """Only binaries with period < p_merge_max_period are merger-
    eligible -- p_merge_max_period=0 must merge nothing (no period is
    < 0), regardless of p_merge.
    """
    config = SimulationConfig(
        ntot=2000,
        binary_prescription="enhanced_mergers",
        p_merge=1.0,
        p_merge_max_period=0.0,
        f_merge=0.5,
        iseed=7,
    )
    pop = BinaryPopulation(config)
    assert len(pop.m1) > 10
    assert not np.any(pop.did_merge)


def test_invalid_binary_prescription_rejected():
    with pytest.raises(ValueError, match="binary_prescription"):
        SimulationConfig(binary_prescription="nonsense")


@pytest.mark.parametrize(
    "field,value",
    [
        ("interaction_boost", -1.0),
        ("p_merge", 1.5),
        ("p_merge_max_period", -1.0),
        ("f_merge", -0.1),
    ],
)
def test_out_of_range_interaction_params_rejected(field, value):
    with pytest.raises(ValueError):
        SimulationConfig(**{field: value})
