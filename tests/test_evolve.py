"""Unit tests for individual BinaryPopulation.evolve() phases, isolated
from the IMF-sampling/lifetime-table stochasticity that drives which
systems reach each phase in a real run.

These complement, not replace, tests/test_regression.py's pinned
end-to-end trajectory: that test proves the *combined* pipeline produces
known output for one config+seed, but (per docs/provenance.md's status
notes) doesn't isolate which phase is responsible if it ever fails, and
can't reach branches a given config doesn't happen to exercise (notably
fsur < 1's stochastic rejection, closed by
test_evolve_phase1_fsur_partial_activation below). Each test here
constructs a BinaryPopulation via the normal config path (so
xray_calc/remnant_table/etc. are wired up correctly) and then overwrites
its internal state arrays directly to set up an exact, deterministic
scenario for one evolve() phase -- monkeypatching
`remnant_table.get_remnant_mass` where needed so the SN-mass-loss
arithmetic is exactly known rather than depending on the bundled data
table (which is covered separately, transitively, by the pinned
regression test).
"""

import numpy as np

from realta.binaries.population import BinaryPopulation
from realta.config import SimulationConfig


def test_evolve_phase1_survival_criterion():
    """floss <= 0.5 must survive; floss > 0.5 must be disrupted.

    Power et al. (2009) Sec. 2.1's deterministic sudden-mass-loss
    criterion -- see binaries/population.py::evolve, Phase 1, and
    docs/provenance.md Section 2. Two systems, identical primary mass
    and remnant mass (so deltam is identical), differing only in
    companion mass m2 -- chosen so one lands just inside floss <= 0.5
    and the other just outside it:
        deltam = 20.0 - 1.4 = 18.6
        system 0: m2=20.0 -> floss = 18.6 / 40.0 = 0.465  (survives)
        system 1: m2=5.0  -> floss = 18.6 / 25.0 = 0.744  (disrupted)
    """
    config = SimulationConfig(ntot=10, fsur=1.0, iseed=1)
    pop = BinaryPopulation(config)

    pop.remnant_table.get_remnant_mass = lambda m: 1.4
    pop.m1 = np.array([20.0, 20.0])
    pop.m2 = np.array([20.0, 5.0])
    pop.period = np.array([10.0, 10.0])
    pop.a = np.array([1.0, 1.0])
    pop.turnoff_time = np.array([5.0, 5.0])
    pop.t2_lifetime = np.array([6.0, 6.0])
    pop.nturn = np.zeros(2, dtype=np.int8)
    pop.is_survived = np.ones(2, dtype=bool)
    pop.lum_xray = np.zeros(2)

    lumx_tot, nphot_tot, nactive, ndead = pop.evolve(tnow=5.0, dt=1.0)

    # Both undergo the primary SN and are relabelled to the remnant mass,
    # regardless of survival outcome.
    assert pop.m1[0] == 1.4
    assert pop.m1[1] == 1.4
    assert pop.nturn[0] == 1
    assert pop.nturn[1] == 1
    assert pop.turnoff_time[0] == pop.t2_lifetime[0]

    # Only the survival outcome differs.
    assert pop.is_survived[0]
    assert not pop.is_survived[1]

    # fsur=1.0 and m2 > mcomp_abs -> the survivor is activated as an
    # HMXB (persistent L_X drawn); the disrupted system is never
    # considered for activation at all.
    assert pop.lum_xray[0] > 0.0
    assert pop.lum_xray[1] == 0.0

    # nactive counts primary-SN *events* this step (formation-rate, not
    # a census of currently-active HMXBs -- both systems count, even
    # the disrupted one) -- see the Phase 3 comment in population.py.
    assert nactive == 2
    assert ndead == 0
    assert lumx_tot == pop.lum_xray[0] * config.lunit


def test_evolve_phase2_death():
    """Secondary SN must mark the system permanently dead.

    binaries/population.py::evolve, Phase 2. A single active HMXB
    (nturn=1, mid-lifetime) reaching its secondary's turnoff should:
    get a remnant mass for m2, drop turnoff_time to 0 (no further
    transitions), be marked nturn=2 (dead), have is_survived cleared,
    and have its persistent X-ray luminosity zeroed.
    """
    config = SimulationConfig(ntot=10, fsur=1.0, iseed=1)
    pop = BinaryPopulation(config)

    pop.remnant_table.get_remnant_mass = lambda m: 1.2
    pop.m1 = np.array([1.4])
    pop.m2 = np.array([5.0])
    pop.period = np.array([10.0])
    pop.a = np.array([1.0])
    pop.turnoff_time = np.array([6.0])
    pop.t2_lifetime = np.array([6.0])
    pop.nturn = np.array([1], dtype=np.int8)
    pop.is_survived = np.array([True])
    pop.lum_xray = np.array([1.0e5])

    lumx_tot, nphot_tot, nactive, ndead = pop.evolve(tnow=6.0, dt=1.0)

    assert pop.m2[0] == 1.2
    assert pop.turnoff_time[0] == 0.0
    assert pop.nturn[0] == 2
    assert not pop.is_survived[0]
    assert pop.lum_xray[0] == 0.0
    assert nactive == 0  # no primary SN this step -- system was already past that
    assert ndead == 1
    assert lumx_tot == 0.0


def test_evolve_phase1_fsur_partial_activation():
    """fsur < 1 must sometimes reject HMXB activation for a surviving binary.

    Statistical check, not an exact-count assertion -- the activation
    gate (`self.np_rng.random() <= self.config.fsur`) is a real draw
    from the same seeded stream also used for the X-ray luminosity draw
    itself (main.f: `ran3(iseed).le.fbin`). This closes the gap noted in
    docs/provenance.md Section 2: previously no automated test exercised
    fsur < 1 at all (the pinned regression config uses fsur=1.0
    specifically so activation is unconditional). Uses a fixed iseed so
    the result is reproducible, not flaky; the tolerance (10 binomial
    std devs) is wide enough that normal RNG variation won't trip it,
    while a real fsur bug (e.g. the gate being ignored, or inverted)
    fails it clearly.
    """
    n = 2000
    fsur = 0.3
    config = SimulationConfig(ntot=10, fsur=fsur, iseed=123)
    pop = BinaryPopulation(config)

    # Identical, guaranteed-survivor HMXB progenitors -- isolates the
    # fsur gate from the floss/survival logic already covered by
    # test_evolve_phase1_survival_criterion above.
    pop.remnant_table.get_remnant_mass = lambda m: 1.4
    pop.m1 = np.full(n, 20.0)
    pop.m2 = np.full(n, 20.0)  # floss = (20-1.4)/40 = 0.465 <= 0.5 -> always survives
    pop.period = np.full(n, 10.0)
    pop.a = np.full(n, 1.0)
    pop.turnoff_time = np.full(n, 5.0)
    pop.t2_lifetime = np.full(n, 6.0)
    pop.nturn = np.zeros(n, dtype=np.int8)
    pop.is_survived = np.ones(n, dtype=bool)
    pop.lum_xray = np.zeros(n)

    pop.evolve(tnow=5.0, dt=1.0)

    assert np.all(pop.is_survived)  # floss criterion unaffected by fsur

    n_active = int(np.count_nonzero(pop.lum_xray > 0))
    expected = n * fsur
    tolerance = 10 * np.sqrt(n * fsur * (1 - fsur))  # ~10 sigma binomial
    assert abs(n_active - expected) < tolerance

    # Sanity bounds: fsur=0.3 must reject a real fraction (not activate
    # everyone) and activate a real fraction (not reject everyone).
    assert 0 < n_active < n
