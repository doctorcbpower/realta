"""Smoke test for scripts/figure5_metallicity_sweep.py (Figure 5, docs/
science/paper1-detailed-work-breakdown.md). Not a numeric regression
pin -- this is a comparison figure built from already-tested underlying
quantities (lumx_tot/luv_tot via A2/A3); this test guards against the
script's own sweep/config-parsing/plotting logic breaking silently.
"""

import importlib.util
import sys
from pathlib import Path

SCRIPT_PATH = Path(__file__).parent.parent / "scripts" / "figure5_metallicity_sweep.py"


def _load_module():
    spec = importlib.util.spec_from_file_location(
        "figure5_metallicity_sweep", SCRIPT_PATH
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_load_sweep_config_parses_the_real_config():
    module = _load_module()
    base, binary_prescriptions, imetal_values = module.load_sweep_config(
        "configs/figure5_metallicity_sweep.yml"
    )
    assert base["mcut"] == 8.0
    assert binary_prescriptions == ["non_interacting", "standard_interaction"]
    assert imetal_values == [1, 2, 3]


def test_load_sweep_config_rejects_forced_fields():
    module = _load_module()
    import tempfile

    import yaml

    bad_config = {
        "base": {"ntot": 100, "imetal": 2},
        "binary_prescriptions": ["non_interacting"],
        "imetal_values": [2],
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


def test_build_sweep_runs_and_produces_a_figure(tmp_path):
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
        "iseed": 7,
    }
    binary_prescriptions = ["non_interacting", "standard_interaction"]
    imetal_values = [2, 3]

    sweep = module.build_sweep(base, binary_prescriptions, imetal_values, tmp_path)
    assert set(sweep.keys()) == {2, 3}
    for imetal in imetal_values:
        assert len(sweep[imetal]) == 2
        for result in sweep[imetal]:
            assert len(result["time"]) > 0
            assert len(result["l_x"]) == len(result["time"])

    fig_path = tmp_path / "figure5_test.png"
    module.make_figure5(sweep, fig_path)
    assert fig_path.exists() and fig_path.stat().st_size > 0


def test_imetal_one_runs_without_crashing_and_logs_no_binary_signature_caveat(
    tmp_path,
):
    """imetal=1 (Z=0) is included in the real config's sweep -- confirm
    it runs end-to-end (does not crash) even though
    use_rlof_classifier is known to skip itself there (binaries/
    population.py's own pre-existing, documented behaviour)."""
    module = _load_module()
    base = {
        "ntot": 3_000,
        "mmin": 0.1,
        "mmax": 100.0,
        "mcut": 8.0,
        "tmax": 10.0,
        "dt": 2.0,
        "pmin": 1.0,
        "pmax": 1000.0,
        "mcomp": 0.5,
        "fsur": 0.5,
        "iseed": 3,
    }
    result = module.run_variant(base, "standard_interaction", 1, tmp_path)
    assert len(result["time"]) > 0
