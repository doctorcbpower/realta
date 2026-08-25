"""Tests for A1 (docs/science/paper1-detailed-work-breakdown.md):
independently configurable binary_fraction, mass_ratio_distribution,
and period_distribution in BinaryPopulation.generate_population().

Config-level defaults/validation are covered in
tests/test_binary_prescriptions.py; this file covers the actual
sampling behaviour at population scale.
"""

import numpy as np
import pytest

from realta.binaries.interaction import RLOFOutcome
from realta.binaries.population import BinaryPopulation
from realta.config import SimulationConfig

BASE_KWARGS = {"ntot": 20_000, "iseed": 42, "mmin": 0.1, "mmax": 100.0, "mcut": 8.0}


def test_binary_fraction_default_pairs_every_massive_star():
    """binary_fraction=1.0 (default) must still pair every M >= mcut
    star -- the pre-existing baseline, unchanged."""
    config = SimulationConfig(**BASE_KWARGS)
    pop = BinaryPopulation(config)
    assert np.all(pop.m2 > 0.0)
    assert np.all(pop.period > 0.0)
    assert np.all(pop.a > 0.0)


def test_binary_fraction_zero_gives_no_companions():
    config = SimulationConfig(**BASE_KWARGS, binary_fraction=0.0)
    pop = BinaryPopulation(config)
    assert np.all(pop.m2 == 0.0)
    assert np.all(pop.period == 0.0)
    assert np.all(pop.a == 0.0)
    # m1 (and therefore SN/lifetime tracking) is still populated --
    # this generalizes, rather than replaces, the "single"
    # binary_prescription's own array-emptying shortcut.
    assert len(pop.m1) > 0


def test_binary_fraction_intermediate_value_matches_expected_fraction():
    """Statistical check: the observed companion fraction should sit
    close to the configured binary_fraction for a large-enough sample.
    """
    config = SimulationConfig(**BASE_KWARGS, binary_fraction=0.5)
    pop = BinaryPopulation(config)
    observed_fraction = np.mean(pop.m2 > 0.0)
    assert 0.4 < observed_fraction < 0.6


def test_binary_fraction_one_is_rng_neutral_vs_baseline():
    """Sensitivity check for the "skip the RNG draw entirely at
    binary_fraction=1.0" claim: explicitly passing binary_fraction=1.0
    must produce bit-identical output to not passing it at all --
    confirms the new Bernoulli draw truly isn't consuming RNG state in
    the default case, not just that the *fraction* happens to come out
    right.
    """
    config_default = SimulationConfig(**BASE_KWARGS)
    config_explicit = SimulationConfig(**BASE_KWARGS, binary_fraction=1.0)
    pop_default = BinaryPopulation(config_default)
    pop_explicit = BinaryPopulation(config_explicit)

    assert np.array_equal(pop_default.m1, pop_explicit.m1)
    assert np.array_equal(pop_default.m2, pop_explicit.m2)
    assert np.array_equal(pop_default.period, pop_explicit.period)
    assert np.array_equal(pop_default.a, pop_explicit.a)


def test_no_companion_stars_excluded_from_rlof_classifier():
    """Regression guard for the divide-by-zero/inf-q1 risk: a star with
    no companion (m2=0) must never reach find_rlof_onset -- its
    rlof_time must stay at the untouched np.inf default, and its
    rlof_outcome must stay DETACHED (the untouched default), not
    crash or silently produce garbage.
    """
    config = SimulationConfig(
        **BASE_KWARGS,
        binary_fraction=0.5,
        use_rlof_classifier=True,
        imetal=3,
    )
    pop = BinaryPopulation(config)
    no_companion = pop.m2 == 0.0
    assert np.any(no_companion)  # confirms the scenario actually arises
    assert np.all(np.isinf(pop.rlof_time[no_companion]))
    assert all(
        outcome == RLOFOutcome.DETACHED for outcome in pop.rlof_outcome[no_companion]
    )


def test_mass_ratio_distribution_flat_q_gives_uniform_mass_ratio():
    config = SimulationConfig(**BASE_KWARGS, mass_ratio_distribution="flat_q")
    pop = BinaryPopulation(config)
    q = pop.m2 / pop.m1
    assert np.all(q >= 0.0) and np.all(q <= 1.0)
    # Uniform(0,1) mean ~0.5 -- statistical check, large-N sample.
    assert 0.45 < q.mean() < 0.55


def test_mass_ratio_distribution_default_matches_baseline_shape():
    """The default "uniform" distribution (m2 ~ Uniform(mcomp, m1),
    clipped to m1) concentrates more probability near m1 than flat_q's
    uniform-in-q shape does, for a fixed mcomp floor -- confirms the
    two distributions are genuinely different shapes, not the same
    thing under two names."""
    config = SimulationConfig(**BASE_KWARGS)
    pop = BinaryPopulation(config)
    q = pop.m2 / pop.m1
    assert q.mean() > 0.5


def test_period_distribution_log_normal_stays_within_bounds():
    config = SimulationConfig(**BASE_KWARGS, period_distribution="log_normal")
    pop = BinaryPopulation(config)
    assert np.all(pop.period >= config.pmin)
    assert np.all(pop.period <= config.pmax)


def test_period_distribution_log_normal_centres_near_geometric_mean():
    """Sanity/calibration check for the pmin/pmax-derived mu: the mean
    of log10(period) should sit close to the midpoint of
    log10(pmin)/log10(pmax) -- confirms the derivation formula is
    actually producing the intended centring, not just "some"
    distribution."""
    config = SimulationConfig(**BASE_KWARGS, period_distribution="log_normal")
    pop = BinaryPopulation(config)
    expected_mu = 0.5 * (np.log10(config.pmin) + np.log10(config.pmax))
    assert np.log10(pop.period).mean() == pytest.approx(expected_mu, abs=0.05)


def test_period_distribution_log_normal_is_narrower_than_log_uniform():
    """Direction check: a truncated normal centred in log-period space
    should have a smaller spread than the full-range log-uniform
    default -- confirms log_normal isn't just relabelled log_uniform.
    """
    config_uniform = SimulationConfig(**BASE_KWARGS)
    config_normal = SimulationConfig(**BASE_KWARGS, period_distribution="log_normal")
    pop_uniform = BinaryPopulation(config_uniform)
    pop_normal = BinaryPopulation(config_normal)
    assert np.log10(pop_normal.period).std() < np.log10(pop_uniform.period).std()


def test_binary_sampling_pinned_regression_non_default():
    """Numeric regression pin combining non-default binary_fraction,
    mass_ratio_distribution, and period_distribution -- see
    docs/physics/binary-sampling.md. Values captured by an actual run,
    cross-checked for run-to-run determinism before being pinned, not
    hand-derived analytically, following tests/test_regression.py's
    own discipline.
    """
    config = SimulationConfig(
        **BASE_KWARGS,
        binary_fraction=0.5,
        mass_ratio_distribution="flat_q",
        period_distribution="log_normal",
    )
    pop = BinaryPopulation(config)

    assert pop.total_mass_msun == pytest.approx(17891.377637950565, rel=1e-9)
    assert len(pop.m1) == 251
    assert np.sum(pop.m2 > 0.0) == 122
    assert pop.m1[0] == pytest.approx(96.11973553566165, rel=1e-9)
    assert pop.m2[0] == pytest.approx(52.73633665578549, rel=1e-9)
    assert pop.period[0] == pytest.approx(3.497299979643711, rel=1e-9)
