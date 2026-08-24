"""Regression coverage for the Paper 1 RLOF-classifier pipeline end to
end -- closes a real coverage gap found this session: every existing
unit test for `binaries/interaction.py` uses hand-picked, internally
self-consistent Rsun-scale separations, so none of them could have
caught the AU/Rsun unit mismatch (`self.a` is AU, but
`classify_rlof`/`find_rlof_onset`/`apply_common_envelope` expect Rsun)
that made Figure 2 completely degenerate until it was found by
manually running `scripts/run_paper1_experiment.py` (see
docs/provenance.md Section 6's "Residual, understood degeneracy" note
and Section 10's units-bug row). Nothing in the existing suite would
have failed if that conversion were silently removed again.

Two layers, matching `tests/test_regression.py`'s existing pattern for
the plain `fsur`-only baseline:

1. `test_standard_interaction_pinned_trajectory` -- a numeric
   regression pin (population summary stats, RLOF-outcome-distribution
   counts, and a per-timestep lumx_tot/nphot_tot/nactive/ndead
   trajectory) for a `standard_interaction` (`use_rlof_classifier`
   opt-in) run at a smaller-but-still-representative scale than the
   full Paper 1 config, for test speed. If this test starts failing
   with lumx_tot/nphot_tot suddenly collapsing toward the values a
   near-total-`IMMEDIATE_MERGER` population would give, that is exactly
   the AU/Rsun-regression symptom to look for first.
2. `test_run_paper1_experiment_produces_nondegenerate_output` -- an
   actual smoke test of `scripts/run_paper1_experiment.py`'s own
   functions (loaded directly from the script file, since `scripts/`
   is not an importable package), using a small config, confirming the
   RLOF-classifier prescriptions produce non-zero L_X/Q_H (not the
   degenerate all-zero output the units bug produced) and that both
   figure files actually get written.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

from realta import ClusterSimulation, SimulationConfig
from realta.binaries.interaction import RLOFOutcome

SCRIPT_PATH = Path(__file__).parent.parent / "scripts" / "run_paper1_experiment.py"


def _load_run_paper1_experiment_module():
    spec = importlib.util.spec_from_file_location("run_paper1_experiment", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


# Population summary stats, RLOF-outcome-distribution counts, and
# per-timestep (time, lumx_tot, nphot_tot, nactive, ndead) trajectory,
# pinned at every 10th step (dt=1.0 Myr) of a 50-step run -- captured
# by an actual run (cross-checked for run-to-run determinism, per the
# same discipline as tests/test_regression.py) after the AU/Rsun fix,
# not hand-derived analytically. mcut/pmin/pmax/mcomp mirror
# configs/paper1_basic_experiment.yml's own values (the corrected
# pmin=1.0, not the global pmin=0.1 default) so this test exercises the
# same physical regime; ntot is smaller than the full Paper 1 config
# for test speed.
PINNED_CONFIG = {
    "ntot": 20_000,
    "imf_type": 2,  # Kroupa
    "tmax": 50.0,
    "dt": 1.0,
    "mcut": 8.0,
    "pmin": 1.0,
    "pmax": 1000.0,
    "mcomp": 0.5,
    "fsur": 0.5,
    "imetal": 2,
    "binary_prescription": "standard_interaction",
    "iseed": 42,
}
PINNED_TOTAL_MASS_MSUN = 17891.377637950565
PINNED_N_MASSIVE = 251
PINNED_OUTCOME_COUNTS = {
    RLOFOutcome.IMMEDIATE_MERGER: 49,
    RLOFOutcome.DETACHED: 102,
    RLOFOutcome.STABLE_MASS_TRANSFER: 10,
    RLOFOutcome.PHASE_NOT_MODELLED: 90,
}
PINNED_TRAJECTORY = [
    (0.0, 0.0, 0.0, 0, 0),
    (5.0, 1.3921928623100108e37, 2.6816951724411323e47, 9, 6),
    (10.0, 7.811108920381723e39, 1.5046056943894237e50, 11, 37),
    (20.0, 8.559642571565061e39, 1.6487911110175045e50, 13, 80),
    (30.0, 9.34231218345191e39, 1.7995519270391257e50, 11, 109),
    (40.0, 1.037294342621823e40, 1.9980760613826616e50, 1, 140),
    (50.0, 9.840851312015995e39, 1.89558244195823e50, 0, 154),
]


def test_standard_interaction_pinned_trajectory():
    """Numeric regression pin for the RLOF-classifier pipeline (not
    just the plain fsur-only baseline tests/test_regression.py covers).
    See module docstring for what a failure here would mean.
    """
    config = SimulationConfig(**PINNED_CONFIG)
    sim = ClusterSimulation(config)
    results = sim.run(output_dir="tests/output_tmp")

    assert sim.population.total_mass_msun == pytest.approx(
        PINNED_TOTAL_MASS_MSUN, rel=1e-9
    )
    assert len(sim.population.m1) == PINNED_N_MASSIVE

    # Per-element comparison, not vectorized `==` -- see
    # RLOFOutcome's own docstring for why a vectorized comparison
    # against an object array of str-Enum members is unreliable.
    outcome_counts = {}
    for outcome in sim.population.rlof_outcome:
        outcome_counts[outcome] = outcome_counts.get(outcome, 0) + 1
    assert outcome_counts == PINNED_OUTCOME_COUNTS

    dt = PINNED_CONFIG["dt"]
    for time, lumx_tot, nphot_tot, nactive, ndead in PINNED_TRAJECTORY:
        idx = round(time / dt)
        r = results[idx]
        assert r["time"] == pytest.approx(time)
        assert r["lumx_tot"] == pytest.approx(lumx_tot, rel=1e-9)
        assert r["nphot_tot"] == pytest.approx(nphot_tot, rel=1e-9)
        assert r["nactive"] == nactive
        assert r["ndead"] == ndead


def test_run_paper1_experiment_produces_nondegenerate_output(tmp_path):
    """Smoke test of the actual script's own functions (not a
    reimplementation) at small scale: confirms an RLOF-classifier
    prescription produces non-zero L_X/Q_H at late times (the
    regression guard for the AU/Rsun-mismatch class of bug, which made
    every RLOF-classifier prescription produce all-zero L_X/Q_H via
    near-total IMMEDIATE_MERGER) and that both figure files are
    actually written.
    """
    module = _load_run_paper1_experiment_module()

    base = {
        "ntot": 5_000,
        "imf_type": 2,
        "tmax": 20.0,
        "dt": 1.0,
        "mcut": 8.0,
        "pmin": 1.0,
        "pmax": 1000.0,
        "mcomp": 0.5,
        "fsur": 0.5,
        "imetal": 2,
        "iseed": 7,
    }

    results_non_interacting = module.run_variant(base, "non_interacting", tmp_path)
    results_standard = module.run_variant(base, "standard_interaction", tmp_path)

    # non_interacting doesn't use the RLOF classifier at all -- a
    # sanity check that this population produces HMXBs, so a later
    # all-zero result for standard_interaction can't be blamed on the
    # underlying population being too small/short to form any.
    assert np.max(results_non_interacting["l_x"]) > 0.0

    # The actual regression guard: before the AU/Rsun fix, this was
    # zero for every RLOF-classifier-enabled prescription across the
    # whole run.
    assert np.max(results_standard["l_x"]) > 0.0
    assert np.max(results_standard["q_h"]) > 0.0

    figures_dir = tmp_path / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    fig1_path = figures_dir / "figure1_population_evolution.png"
    fig2_path = figures_dir / "figure2_xray_uv_evolution.png"
    module.make_figure1([results_non_interacting, results_standard], fig1_path)
    module.make_figure2([results_non_interacting, results_standard], fig2_path)

    assert fig1_path.exists() and fig1_path.stat().st_size > 0
    assert fig2_path.exists() and fig2_path.stat().st_size > 0
