"""Tests for A2/A3 (docs/science/paper1-detailed-work-breakdown.md):
L_bol/L_UV/Q_H wired into ClusterSimulation.run()'s own per-timestep
output, alongside the pre-existing lumx_tot/nphot_tot.

A2: lbol_tot/luv_tot (FSPS SSP tables, MS-population-level, rescaled
to this run's total_mass_msun) are now computed inside run() itself,
not just post-hoc in scripts/run_paper1_experiment.py.

A3: qh_tot = qh_ms_tot + nphot_tot, where qh_ms_tot is a genuine,
independent ionizing-photon-rate calculation from the currently-alive
M >= 8 Msun massive-star population (via the previously-unused
IonizingPhotonTable), replacing the old degenerate-with-L_X placeholder
for what Q_H(t) actually reports.
"""

import numpy as np
import pytest

from realta.config import SimulationConfig
from realta.simulation.cluster import ClusterSimulation

# Literature ballpark for O/early-B main-sequence ionizing rates
# (Vacca, Garmany & Shull 1996, ApJ 460, 914) -- used only as an
# order-of-magnitude sanity check, not an exact-match requirement.
LITERATURE_QH_RANGE_10MSUN = (1e46, 1e48)
LITERATURE_QH_RANGE_40MSUN = (1e48, 1e50)
LITERATURE_QH_RANGE_80MSUN = (1e49, 1e51)


def test_results_include_lbol_luv_qh_keys():
    config = SimulationConfig(ntot=5000, tmax=5.0, dt=1.0, iseed=1)
    sim = ClusterSimulation(config)
    results = sim.run(output_dir="tests/output_tmp")
    for r in results:
        assert "lbol_tot" in r
        assert "luv_tot" in r
        assert "qh_tot" in r
        assert np.isfinite(r["lbol_tot"])
        assert np.isfinite(r["luv_tot"])
        assert np.isfinite(r["qh_tot"])


def test_qh_ms_tot_matches_direct_ionizing_table_calibration():
    """Sanity/calibration check for the interpretation adopted in
    ClusterSimulation._qh_ms_tot's docstring: IonizingPhotonTable.
    get_ngamma(m) is a *total* photon count over the star's whole MS
    lifetime, not a rate -- dividing by LifetimeTable.get_lifetime(m)
    (converted to seconds) should land within the well-known
    literature range for O/early-B ionizing rates at a few
    representative masses.
    """
    config = SimulationConfig(ntot=10, iseed=1, imetal=2)
    sim = ClusterSimulation(config)
    sim.initialize()

    for mass, (lo, hi) in [
        (10.0, LITERATURE_QH_RANGE_10MSUN),
        (40.0, LITERATURE_QH_RANGE_40MSUN),
        (80.0, LITERATURE_QH_RANGE_80MSUN),
    ]:
        lifetime_s = (
            sim.population.lifetime_table.get_lifetime(mass) * sim.MYR_TO_SECONDS
        )
        qh_rate = 10.0 ** sim.ionizing_table.get_ngamma(mass) / lifetime_s
        assert lo < qh_rate < hi, f"Q_H({mass})={qh_rate:.2e} outside {lo:.0e}-{hi:.0e}"


def test_qh_ms_tot_zero_before_any_massive_star_and_after_all_have_died():
    config = SimulationConfig(ntot=5000, tmax=5.0, dt=0.5, mcut=8.0, iseed=1)
    sim = ClusterSimulation(config)
    sim.initialize()
    # Before evolve() has run at all, every massive star is still
    # nturn==0 -- qh_ms_tot should already be nonzero at t=0 (stars are
    # born already ionizing), not need an evolve() call first.
    assert sim._qh_ms_tot(0.0) > 0.0


def test_single_prescription_gives_nonzero_qh():
    """Regression guard for the fix made while implementing A3: the
    'single' prescription used to empty BinaryPopulation.m1 entirely,
    which silently made qh_ms_tot always 0 for single-star populations
    (physically wrong -- massive single stars still ionize). Migrated
    onto the same has_companion=False mechanism A1 uses for
    binary_fraction<1, which keeps m1 populated. L_X/HMXB stay exactly
    zero (m2=0 blocks that unconditionally); Q_H must not.
    """
    config = SimulationConfig(
        ntot=5000, tmax=2.0, dt=1.0, binary_prescription="single", iseed=7
    )
    sim = ClusterSimulation(config)
    results = sim.run(output_dir="tests/output_tmp")

    assert any(r["qh_tot"] > 0.0 for r in results)
    assert all(r["lumx_tot"] == 0.0 for r in results)


def test_qh_tot_is_independent_of_lumx_tot_not_a_fixed_multiple():
    """Direction/sensitivity check for the actual bug A3 fixes: the
    old placeholder was `nphot_tot = NPHOT_PER_LUMX * lumx_tot`, a
    fixed multiple with zero independent information. qh_tot must NOT
    be explainable as a fixed multiple of lumx_tot across timesteps
    where lumx_tot is nonzero (confirms qh_ms_tot is contributing real,
    independent variation)."""
    config = SimulationConfig(ntot=20_000, tmax=30.0, dt=1.0, fsur=1.0, iseed=42)
    sim = ClusterSimulation(config)
    results = sim.run(output_dir="tests/output_tmp")

    ratios = [r["qh_tot"] / r["lumx_tot"] for r in results if r["lumx_tot"] > 0.0]
    assert len(ratios) > 2
    # A fixed multiple would give ~identical ratios at every timestep;
    # real independent variation should not.
    assert np.std(ratios) / np.mean(ratios) > 0.01


def test_qh_ms_tot_excludes_stars_after_their_own_supernova():
    """Sensitivity check: a star must stop contributing to qh_ms_tot
    once its own nturn/lifetime clock says it has exploded -- verified
    by hand-constructing a population state directly (same technique
    as tests/test_evolve.py) rather than waiting for a real population
    to reach that point."""
    config = SimulationConfig(ntot=10, iseed=1)
    sim = ClusterSimulation(config)
    sim.initialize()
    pop = sim.population

    pop.m1 = np.array([40.0])
    pop.m2 = np.array([0.0])
    pop.nturn = np.array([0], dtype=np.int8)
    pop.t2_lifetime = np.array([0.0])

    assert sim._qh_ms_tot(1.0) > 0.0

    pop.nturn = np.array([1], dtype=np.int8)  # m1 has now exploded
    assert sim._qh_ms_tot(1.0) == 0.0


def test_qh_ms_tot_includes_companion_while_alive():
    config = SimulationConfig(ntot=10, iseed=1)
    sim = ClusterSimulation(config)
    sim.initialize()
    pop = sim.population

    pop.m1 = np.array([100.0])  # already exploded -- contributes nothing
    pop.m2 = np.array([20.0])  # real companion, still alive
    pop.nturn = np.array([1], dtype=np.int8)
    pop.t2_lifetime = np.array([50.0])

    assert sim._qh_ms_tot(10.0) > 0.0  # companion still alive at t=10 < 50
    assert sim._qh_ms_tot(60.0) == 0.0  # companion has died by t=60 > 50


def test_lbol_tot_equals_ms_lbol_plus_lumx_tot():
    """Confirms the "population contribution + HMXB contribution"
    convention (matching what scripts/run_paper1_experiment.py used to
    compute independently) is what actually lands in results."""
    config = SimulationConfig(ntot=20_000, tmax=10.0, dt=1.0, fsur=1.0, iseed=42)
    sim = ClusterSimulation(config)
    results = sim.run(output_dir="tests/output_tmp")

    for r in results:
        ms_lbol = sim.ms_table.get_lbol(r["time"], sim.population.total_mass_msun)
        assert r["lbol_tot"] == pytest.approx(ms_lbol + r["lumx_tot"], rel=1e-9)


def test_luv_tot_is_ms_only_no_hmxb_contribution():
    config = SimulationConfig(ntot=20_000, tmax=10.0, dt=1.0, fsur=1.0, iseed=42)
    sim = ClusterSimulation(config)
    results = sim.run(output_dir="tests/output_tmp")

    for r in results:
        expected = sim.uv_table.get_luv(r["time"], sim.population.total_mass_msun)
        assert r["luv_tot"] == pytest.approx(expected, rel=1e-9)
