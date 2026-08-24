"""Smoke test for scripts/xu2025_smc_crosscheck.py (B2, docs/science/
paper1-detailed-work-breakdown.md) -- guards against the script
silently breaking due to an API change elsewhere (e.g. RLOFOutcome
values, apply_common_envelope's signature, BinaryPopulation.RSUN_PER_AU).
Not a numeric regression pin: the script's whole point is a one-off
literature comparison (see its own module docstring and
docs/provenance.md Section 10 for the actual result and discussion),
not a physics path whose exact output fractions should be pinned here.
"""

import importlib.util
import sys
from pathlib import Path

SCRIPT_PATH = Path(__file__).parent.parent / "scripts" / "xu2025_smc_crosscheck.py"


def test_xu2025_crosscheck_runs_without_crashing(capsys):
    spec = importlib.util.spec_from_file_location("xu2025_smc_crosscheck", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)

    # Small ntot for test speed -- exercises the script's real logic
    # (population generation, q-window filtering, COMMON_ENVELOPE
    # resolution via apply_common_envelope, fraction reporting), not a
    # re-derived duplicate.
    module.main(ntot=5_000)

    out = capsys.readouterr().out
    assert "post-mass-transfer" in out
    assert "merger" in out
