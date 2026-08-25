"""Smoke test for scripts/figure4_imf_binary_grid.py (C2, docs/science/
paper1-detailed-work-breakdown.md) -- guards against the script
silently breaking due to an API change elsewhere. Not a numeric
regression pin: this is a one-off degeneracy-illustration figure, not
a physics path whose exact output should be pinned.
"""

import importlib.util
import sys
from pathlib import Path

import numpy as np

SCRIPT_PATH = Path(__file__).parent.parent / "scripts" / "figure4_imf_binary_grid.py"


def _load_module():
    spec = importlib.util.spec_from_file_location(
        "figure4_imf_binary_grid", SCRIPT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_load_grid_config_parses_the_real_config():
    module = _load_module()
    base, alpha_imf_values, f_bin_values, selected_ages = module.load_grid_config(
        "configs/figure4_imf_binary_grid.yml"
    )
    assert base["mcut"] == 8.0
    assert len(alpha_imf_values) > 1
    assert len(f_bin_values) > 1
    assert len(selected_ages) >= 1


def test_load_grid_config_rejects_forced_fields():
    module = _load_module()
    import tempfile

    import yaml

    bad_config = {
        "base": {"ntot": 100, "imf_slope": 2.0},
        "alpha_imf_values": [2.0],
        "f_bin_values": [1.0],
        "selected_ages": [10.0],
    }
    with tempfile.NamedTemporaryFile(mode="w", suffix=".yml", delete=False) as f:
        yaml.dump(bad_config, f)
        path = f.name
    try:
        try:
            module.load_grid_config(path)
            raised = False
        except ValueError:
            raised = True
        assert raised
    finally:
        import os

        os.remove(path)


def test_build_grid_runs_and_produces_a_figure(tmp_path):
    module = _load_module()
    base = {
        "ntot": 3_000,
        "mmin": 0.1,
        "mmax": 100.0,
        "mcut": 8.0,
        "tmax": 20.0,
        "dt": 2.0,
        "pmin": 1.0,
        "pmax": 1000.0,
        "mcomp": 0.5,
        "fsur": 1.0,
        "imetal": 2,
        "iseed": 7,
    }
    alpha_imf_values = [2.0, 2.5]
    f_bin_values = [0.0, 1.0]
    selected_ages = [10.0]

    grid = module.build_grid(
        base, alpha_imf_values, f_bin_values, selected_ages, tmp_path
    )
    assert grid.shape == (1, 2, 2)
    # f_bin=0.0 (no binaries at all) must give L_X/L_UV == 0 or NaN
    # (no HMXB channel possible), regardless of alpha_imf.
    assert np.all(np.nan_to_num(grid[0, 0, :]) == 0.0)

    fig_path = tmp_path / "figure4_test.png"
    module.make_figure4(grid, alpha_imf_values, f_bin_values, selected_ages, fig_path)
    assert fig_path.exists() and fig_path.stat().st_size > 0
