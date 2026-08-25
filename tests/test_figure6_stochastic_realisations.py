"""Smoke test for scripts/figure6_stochastic_realisations.py (Figure
6, docs/science/research-programme.md). Not a numeric regression pin
-- this is a lightweight first look at stochasticity, not a physics
path with a citable exact value; this test guards against the script's
own sweep/config-parsing/plotting logic breaking silently.
"""

import importlib.util
import sys
from pathlib import Path

import numpy as np

SCRIPT_PATH = (
    Path(__file__).parent.parent / "scripts" / "figure6_stochastic_realisations.py"
)


def _load_module():
    spec = importlib.util.spec_from_file_location(
        "figure6_stochastic_realisations", SCRIPT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_load_sweep_config_parses_the_real_config():
    module = _load_module()
    base, ntot_values, n_realizations, selected_ages, base_iseed = (
        module.load_sweep_config("configs/figure6_stochastic_realisations.yml")
    )
    assert base["mcut"] == 8.0
    assert len(ntot_values) >= 2
    assert n_realizations > 0
    assert len(selected_ages) >= 1
    assert isinstance(base_iseed, int)


def test_load_sweep_config_rejects_forced_fields():
    module = _load_module()
    import tempfile

    import yaml

    for forced_field in ("ntot", "iseed"):
        bad_config = {
            "base": {"mcut": 8.0, forced_field: 100},
            "ntot_values": [1000],
            "n_realizations": 3,
            "selected_ages": [10.0],
            "base_iseed": 1,
        }
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yml", delete=False) as f:
            yaml.dump(bad_config, f)
            path = f.name
        try:
            try:
                module.load_sweep_config(path)
                raised = False
            except ValueError:
                raised = True
            assert raised
        finally:
            import os

            os.remove(path)


def test_build_distributions_runs_and_produces_a_figure(tmp_path):
    module = _load_module()
    base = {
        "mmin": 0.1,
        "mmax": 100.0,
        "mcut": 8.0,
        "tmax": 20.0,
        "dt": 4.0,
        "pmin": 1.0,
        "pmax": 1000.0,
        "mcomp": 0.5,
        "fsur": 0.5,
        "binary_prescription": "non_interacting",
        "imetal": 2,
    }
    ntot_values = [500, 3000]
    n_realizations = 6
    selected_ages = [10.0]

    sweep = module.build_distributions(
        base, ntot_values, n_realizations, selected_ages, 1, tmp_path
    )
    assert set(sweep.keys()) == {500, 3000}
    # Larger ntot must give a genuinely larger realized cluster mass.
    assert sweep[3000]["mean_mass"] > sweep[500]["mean_mass"]

    for ntot in ntot_values:
        age_entry = sweep[ntot]["ages"][10.0]
        n_finite = age_entry["log10_ratios"].size
        n_quiescent = round(age_entry["quiescent_fraction"] * n_realizations)
        assert n_finite + n_quiescent == n_realizations

    fig_path = tmp_path / "figure6_test.png"
    module.make_figure6(sweep, selected_ages, fig_path)
    assert fig_path.exists() and fig_path.stat().st_size > 0


def test_smaller_cluster_has_a_higher_or_equal_quiescent_fraction(tmp_path):
    """Direction/sanity check: a smaller cluster has fewer massive
    stars, so it should be at least as likely (never strictly less
    likely) to have zero active HMXBs at a given age -- the actual
    physical point of this figure. Uses fsur=1.0 and a wide period
    range so activation itself isn't the bottleneck, isolating pure
    finite-N sampling noise (a cluster with zero massive companions
    at all vs. one with several)."""
    module = _load_module()
    base = {
        "mmin": 0.1,
        "mmax": 100.0,
        "mcut": 8.0,
        "tmax": 15.0,
        "dt": 5.0,
        "pmin": 1.0,
        "pmax": 1000.0,
        "mcomp": 0.5,
        "fsur": 1.0,
        "binary_prescription": "non_interacting",
        "imetal": 2,
    }
    ntot_values = [200, 20000]
    n_realizations = 15
    selected_ages = [10.0]

    sweep = module.build_distributions(
        base, ntot_values, n_realizations, selected_ages, 5, tmp_path
    )
    small_q = sweep[200]["ages"][10.0]["quiescent_fraction"]
    large_q = sweep[20000]["ages"][10.0]["quiescent_fraction"]
    assert small_q >= large_q
    assert small_q > 0.0
    assert np.isclose(large_q, 0.0)
