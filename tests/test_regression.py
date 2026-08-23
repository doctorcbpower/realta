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
        "config": {
            "ntot": 20_000,
            "imf_type": 2,  # Kroupa
            "tmax": 50.0,
            "dt": 1.0,
            "fsur": 1.0,
            "lxmin": 33.0,
            "lxmax": 38.0,
            "xray_distribution": "weibull",
            "iseed": 42,
        },
        "total_mass_msun": 17891.377637950565,
        "n_massive": 251,
        "trajectory": [
            (0.0, 0.0, 0.0, 0, 0),
            (5.0, 1.0792994993518963e40, 2.0789880018689555e50, 13, 6),
            (10.0, 1.9367829502108998e40, 3.73070544193776e50, 8, 30),
            (20.0, 1.9838986565644173e40, 3.821461518696216e50, 9, 67),
            (30.0, 2.0917421183153941e40, 4.029193722032664e50, 7, 91),
            (40.0, 2.071274531037113e40, 3.989768272095641e50, 0, 117),
            (50.0, 1.9281491935918848e40, 3.7140747695129444e50, 0, 132),
        ],
    },
    "fsur=0.5": {
        # Same population (ntot/iseed/imf_type unchanged, so mass/counts
        # are identical to fsur=1.0 above -- floss/survival/lifetime are
        # all deterministic and unaffected by fsur; only which surviving
        # binaries actually get activated as HMXBs changes) -- confirms
        # the activation gate itself is what fsur is changing here, not
        # something upstream.
        "config": {
            "ntot": 20_000,
            "imf_type": 2,
            "tmax": 50.0,
            "dt": 1.0,
            "fsur": 0.5,
            "lxmin": 33.0,
            "lxmax": 38.0,
            "xray_distribution": "weibull",
            "iseed": 42,
        },
        "total_mass_msun": 17891.377637950565,
        "n_massive": 251,
        "trajectory": [
            (0.0, 0.0, 0.0, 0, 0),
            (5.0, 1.2746482371385436e40, 2.455276216847799e50, 13, 6),
            (10.0, 1.1261812399102986e40, 2.1692934047586806e50, 8, 30),
            (20.0, 9.946778743922909e39, 1.915986589290412e50, 9, 67),
            (30.0, 7.986502664821549e39, 1.5383907086984483e50, 7, 91),
            (40.0, 8.747347011508785e39, 1.6849474586091476e50, 0, 117),
            (50.0, 7.806512213812698e39, 1.503720259177895e50, 0, 132),
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
    config_kwargs = {"ntot": 500, "tmax": 5.0, "dt": 1.0, "fsur": 1.0, "iseed": 7}

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
