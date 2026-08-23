"""Numeric regression test pinning the Power et al. (2009) baseline output.

Closes the "no numeric regression test" gap flagged in
docs/provenance.md's "Known gaps" section: `test_reference_cluster_run`
below only checked that a run didn't crash or produce NaNs, so nothing
protected the ported physics (SN survival, HMXB activation, X-ray
luminosity, ionising photon rate) from silent drift under refactoring.

This test pins exact per-timestep trajectory values plus population-
generation summary stats, for a fixed config and iseed, against values
captured by an actual run (see the values below -- generated and cross-
checked for run-to-run determinism before being pinned; not hand-derived
analytically). A run is fully deterministic for a fixed `iseed` (see
docs/provenance.md Section 5, "Random number generation") -- if this
test fails, either a real bug was introduced, or a deliberate change to
the ported physics was made and these pinned values need updating
*with an explicit note on what changed and why*, per the brief's
"do not silently change scientific behaviour" principle. Do not update
the pinned values without first identifying which specific change moved
them and confirming it was the change you intended to make.

Two configs are pinned (see REGRESSION_CASES below): one with fsur=1.0
(every surviving binary becomes active -- exercises HMXB-activation and
X-ray-luminosity code paths unconditionally) and one with fsur=0.5
(exercises the *stochastic rejection* branch of the activation gate,
which fsur=1.0 can never reach since `rng.random() <= 1.0` is always
true). See docs/provenance.md Section 2 for why fsur<1 coverage matters
here specifically -- tests/test_evolve.py's
test_evolve_phase1_fsur_partial_activation additionally unit-tests that
gate in isolation, with a hand-controlled scenario rather than a full
population run.
"""

import numpy as np
import pytest

from realta import ClusterSimulation, SimulationConfig

# Population-generation summary stats and per-timestep
# (time, lumx_tot, nphot_tot, nactive, ndead) trajectory, pinned at
# every 5th step (dt=1.0 Myr) of a 50-step run, for each config below.
# ntot/tmax chosen to run in a few seconds while still producing a
# nontrivial number of active HMXBs (regression-relevant, not a full
# production run). Values captured by an actual run each, cross-checked
# for run-to-run determinism before being pinned -- not hand-derived
# analytically.
REGRESSION_CASES = {
    "fsur=1.0": {
        "config": dict(
            ntot=20_000,
            imf_type=2,  # Kroupa
            tmax=50.0,
            dt=1.0,
            fsur=1.0,
            lxmin=33.0,
            lxmax=38.0,
            xray_distribution="weibull",
            iseed=42,
        ),
        "total_mass_msun": 15836.446390529192,
        "n_massive": 225,
        "trajectory": [
            (0.0, 0.0, 0.0, 0, 0),
            (5.0, 2.8352774824704705e39, 5.461419996548588e49, 12, 5),
            (10.0, 1.2896603831968012e40, 2.4841931871198718e50, 7, 30),
            (20.0, 1.533413364849877e40, 2.9537195091285727e50, 9, 56),
            (30.0, 1.8777611171450152e40, 3.617015328242748e50, 4, 78),
            (40.0, 1.8631054345297103e40, 3.5887849914968177e50, 3, 96),
            (50.0, 1.5607903733997897e40, 3.00645415075169e50, 0, 110),
        ],
    },
    "fsur=0.5": {
        # Same population (ntot/iseed/imf_type unchanged, so mass/counts
        # are identical to fsur=1.0 above -- floss/survival/lifetime are
        # all deterministic and unaffected by fsur; only which surviving
        # binaries actually get activated as HMXBs changes) -- confirms
        # the activation gate itself is what fsur is changing here, not
        # something upstream.
        "config": dict(
            ntot=20_000,
            imf_type=2,
            tmax=50.0,
            dt=1.0,
            fsur=0.5,
            lxmin=33.0,
            lxmax=38.0,
            xray_distribution="weibull",
            iseed=42,
        ),
        "total_mass_msun": 15836.446390529192,
        "n_massive": 225,
        "trajectory": [
            (0.0, 0.0, 0.0, 0, 0),
            (5.0, 5.065385723578021e37, 9.75713983975777e47, 12, 5),
            (10.0, 4.7929800527375766e39, 9.232421611260102e49, 7, 30),
            (20.0, 6.743327946233837e39, 1.2989256366102418e50, 9, 56),
            (30.0, 6.65390586443703e39, 1.2817008129844234e50, 4, 78),
            (40.0, 6.899024011734441e39, 1.3289164086163744e50, 3, 96),
            (50.0, 4.872783620101247e39, 9.386142296903928e49, 0, 110),
        ],
    },
}


def test_reference_cluster_run():
    """Basic sanity: a run completes and produces finite output.

    Kept separate from the numeric-pinning test below so a genuine
    crash/NaN is reported distinctly from a physics-value drift.
    """
    config = SimulationConfig(ntot=1000, tmax=10.0)
    cluster = ClusterSimulation(config)
    results = cluster.run(output_dir="tests/output_tmp")

    assert len(results) > 0
    assert results[0]["time"] == 0.0
    assert not np.isnan(results[-1]["lumx_tot"])


def test_run_is_deterministic_for_fixed_seed():
    """A fixed iseed must reproduce bit-identical output across runs.

    This is the actual claim docs/provenance.md Section 5 makes (single
    seeded np.random.Generator threaded through population generation,
    SN survival, HMXB activation, and XRayLuminosity.get_lumx) -- distinct
    from test_reference_cluster_run_pinned_trajectory, which pins *what*
    the values are; this pins that they don't vary run-to-run at all.
    Uses a small ntot for speed since only reproducibility, not physics
    coverage, is being tested here.
    """
    config_kwargs = dict(ntot=500, tmax=5.0, dt=1.0, fsur=1.0, iseed=7)

    def run():
        sim = ClusterSimulation(SimulationConfig(**config_kwargs))
        results = sim.run(output_dir="tests/output_tmp")
        return sim.population.total_mass_msun, results

    mass_a, results_a = run()
    mass_b, results_b = run()

    assert mass_a == mass_b
    assert results_a == results_b


@pytest.mark.parametrize(
    "case_name", list(REGRESSION_CASES.keys()), ids=list(REGRESSION_CASES.keys())
)
def test_reference_cluster_run_pinned_trajectory(case_name):
    """Pin exact simulation output against the values captured above.

    See module docstring for why this exists (including why fsur=0.5 is
    pinned alongside fsur=1.0) and what a failure means.
    """
    case = REGRESSION_CASES[case_name]
    config = SimulationConfig(**case["config"])
    sim = ClusterSimulation(config)
    results = sim.run(output_dir="tests/output_tmp")

    assert sim.population.total_mass_msun == pytest.approx(
        case["total_mass_msun"], rel=1e-9
    )
    assert len(sim.population.m1) == case["n_massive"]

    dt = case["config"]["dt"]
    for time, lumx_tot, nphot_tot, nactive, ndead in case["trajectory"]:
        idx = round(time / dt)
        r = results[idx]
        assert r["time"] == pytest.approx(time)
        assert r["lumx_tot"] == pytest.approx(lumx_tot, rel=1e-9)
        assert r["nphot_tot"] == pytest.approx(nphot_tot, rel=1e-9)
        assert r["nactive"] == nactive
        assert r["ndead"] == ndead
