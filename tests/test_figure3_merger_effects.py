"""Smoke test for scripts/figure3_merger_effects.py (C3, docs/science/
paper1-detailed-work-breakdown.md). Not a numeric regression pin --
this is a comparison figure built from already-tested underlying
quantities (lumx_tot via A2/A3, merge_time via the existing merger
bookkeeping); this test guards against the script's own aggregation
logic (cumulative merger counting, config parsing) breaking silently.
"""

import importlib.util
import sys
from pathlib import Path

import numpy as np

SCRIPT_PATH = Path(__file__).parent.parent / "scripts" / "figure3_merger_effects.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("figure3_merger_effects", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_load_experiment_config_parses_the_real_config():
    module = _load_module()
    base, prescriptions = module.load_experiment_config(
        "configs/figure3_merger_effects.yml"
    )
    assert base["mcut"] == 8.0
    assert prescriptions == [
        "non_interacting",
        "standard_interaction",
        "enhanced_mergers",
    ]


def test_no_mergers_prescription_gives_zero_cumulative_mergers(tmp_path):
    module = _load_module()
    base = {
        "ntot": 5_000,
        "mmin": 0.1,
        "mmax": 100.0,
        "mcut": 8.0,
        "tmax": 20.0,
        "dt": 2.0,
        "pmin": 1.0,
        "pmax": 1000.0,
        "mcomp": 0.5,
        "fsur": 0.5,
        "imetal": 2,
        "iseed": 1,
    }
    r = module.run_variant(base, "non_interacting", tmp_path)
    assert np.all(r["cumulative_mergers"] == 0)


def test_enhanced_mergers_produces_real_mergers_and_cumulative_count_is_monotonic(
    tmp_path,
):
    module = _load_module()
    base = {
        "ntot": 20_000,
        "mmin": 0.1,
        "mmax": 100.0,
        "mcut": 8.0,
        "tmax": 30.0,
        "dt": 1.0,
        "pmin": 1.0,
        "pmax": 1000.0,
        "mcomp": 0.5,
        "fsur": 0.5,
        "imetal": 2,
        "iseed": 42,
    }
    r = module.run_variant(base, "enhanced_mergers", tmp_path)
    assert r["cumulative_mergers"][-1] > 0
    # A cumulative count must never decrease.
    assert np.all(np.diff(r["cumulative_mergers"]) >= 0)
    assert np.all(np.diff(r["ndead"]) >= 0)


def test_make_figure3_produces_a_file(tmp_path):
    module = _load_module()
    base = {
        "ntot": 5_000,
        "mmin": 0.1,
        "mmax": 100.0,
        "mcut": 8.0,
        "tmax": 10.0,
        "dt": 2.0,
        "pmin": 1.0,
        "pmax": 1000.0,
        "mcomp": 0.5,
        "fsur": 0.5,
        "imetal": 2,
        "iseed": 1,
    }
    all_results = [
        module.run_variant(base, p, tmp_path)
        for p in ["non_interacting", "enhanced_mergers"]
    ]
    out_path = tmp_path / "figure3_test.png"
    module.make_figure3(all_results, out_path)
    assert out_path.exists() and out_path.stat().st_size > 0
