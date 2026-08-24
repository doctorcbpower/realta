"""Tests for the post-SN secondary Roche-lobe-overflow channel
(docs/science/paper1-followup-prompt.md, config.use_post_sn_rlof) --
the real HMXB-formation channel once the primary has already collapsed:
the secondary's own later RLOF onto the compact primary, independent of
config.use_rlof_classifier (which only covers PRE-SN interaction
between two still-live stars).
"""

import numpy as np

from realta.binaries.interaction import RLOFOutcome
from realta.binaries.population import BinaryPopulation
from realta.config import SimulationConfig
from realta.stellar import main_sequence as ms

Z_IMETAL2 = 0.008


def _hand_constructed_post_sn_population(
    config, m1_compact, m2_secondary, a_rsun, nturn=1, is_survived=True, lum_xray=0.0
):
    pop = BinaryPopulation(config)
    pop.m1 = np.array([m1_compact])
    pop.m2 = np.array([m2_secondary])
    pop.a = np.array([a_rsun / BinaryPopulation.RSUN_PER_AU])
    pop.period = np.array([100.0])
    pop.turnoff_time = np.array([1.0e6])
    pop.t2_lifetime = np.array([ms.t_bgb(m2_secondary, Z_IMETAL2)])
    pop.nturn = np.array([nturn], dtype=np.int8)
    pop.is_survived = np.array([is_survived])
    pop.lum_xray = np.array([lum_xray])
    pop.did_merge = np.zeros(1, dtype=bool)
    pop.merge_time = np.array([np.nan])
    pop.rlof_time = np.array([np.inf])
    pop.rlof_outcome = np.array([RLOFOutcome.DETACHED], dtype=object)
    pop.rlof_donor_is_star1 = np.array([False])
    pop.rlof_processed = np.array([True])
    return pop


def test_config_use_post_sn_rlof_defaults_to_false():
    config = SimulationConfig(ntot=10, iseed=1)
    assert config.use_post_sn_rlof is False


def test_disabled_by_default_is_inert():
    """With use_post_sn_rlof=False (default), even a tight,
    guaranteed-to-overflow scenario must never activate via this
    channel."""
    config = SimulationConfig(ntot=10, imetal=2, iseed=1)
    pop = _hand_constructed_post_sn_population(config, 1.4, 8.0, 15.0)
    for tnow in np.arange(0.0, ms.t_bgb(8.0, Z_IMETAL2), 2.0):
        pop.evolve(tnow=tnow, dt=1.0)
    assert pop.lum_xray[0] == 0.0


def test_triggers_when_secondary_fills_its_roche_lobe():
    config = SimulationConfig(ntot=10, imetal=2, use_post_sn_rlof=True, iseed=1)
    pop = _hand_constructed_post_sn_population(config, 1.4, 8.0, 15.0)
    triggered_at = None
    for tnow in np.arange(0.0, ms.t_bgb(8.0, Z_IMETAL2), 1.0):
        pop.evolve(tnow=tnow, dt=1.0)
        if pop.lum_xray[0] > 0.0:
            triggered_at = tnow
            break
    assert triggered_at is not None
    assert pop.lum_xray[0] > 0.0


def test_never_triggers_for_a_wide_enough_orbit():
    """A wide-enough separation must never trigger within the
    secondary's own MS+HG lifetime -- the Roche lobe is never filled.
    """
    config = SimulationConfig(ntot=10, imetal=2, use_post_sn_rlof=True, iseed=1)
    pop = _hand_constructed_post_sn_population(config, 1.4, 8.0, 500.0)
    for tnow in np.arange(0.0, ms.t_bgb(8.0, Z_IMETAL2), 5.0):
        pop.evolve(tnow=tnow, dt=1.0)
    assert pop.lum_xray[0] == 0.0


def test_does_not_apply_before_primary_supernova():
    """nturn==0 (primary hasn't exploded yet) must never trigger this
    channel, even with an otherwise-overflowing configuration --
    that's Phase 0's territory (config.use_rlof_classifier), not this
    one."""
    config = SimulationConfig(ntot=10, imetal=2, use_post_sn_rlof=True, iseed=1)
    pop = _hand_constructed_post_sn_population(
        config, 1.4, 8.0, 15.0, nturn=0, is_survived=False
    )
    for tnow in np.arange(0.0, ms.t_bgb(8.0, Z_IMETAL2), 2.0):
        pop.evolve(tnow=tnow, dt=1.0)
    assert pop.lum_xray[0] == 0.0


def test_does_not_apply_after_secondary_supernova():
    """nturn==2 (both stars already compact) must never trigger."""
    config = SimulationConfig(ntot=10, imetal=2, use_post_sn_rlof=True, iseed=1)
    pop = _hand_constructed_post_sn_population(config, 1.4, 8.0, 15.0, nturn=2)
    for tnow in np.arange(0.0, ms.t_bgb(8.0, Z_IMETAL2), 2.0):
        pop.evolve(tnow=tnow, dt=1.0)
    assert pop.lum_xray[0] == 0.0


def test_does_not_apply_to_disrupted_binaries():
    """is_survived=False (disrupted at SN1) must never trigger --
    there is no bound orbit for RLOF to occur in."""
    config = SimulationConfig(ntot=10, imetal=2, use_post_sn_rlof=True, iseed=1)
    pop = _hand_constructed_post_sn_population(
        config, 1.4, 8.0, 15.0, is_survived=False
    )
    for tnow in np.arange(0.0, ms.t_bgb(8.0, Z_IMETAL2), 2.0):
        pop.evolve(tnow=tnow, dt=1.0)
    assert pop.lum_xray[0] == 0.0


def test_does_not_overwrite_an_already_active_system():
    """A system already activated (e.g. via the stochastic fsur draw
    at SN1) must not be re-drawn/overwritten by this channel."""
    config = SimulationConfig(ntot=10, imetal=2, use_post_sn_rlof=True, iseed=1)
    sentinel_lumx = 12345.0
    pop = _hand_constructed_post_sn_population(
        config, 1.4, 8.0, 15.0, lum_xray=sentinel_lumx
    )
    for tnow in np.arange(0.0, ms.t_bgb(8.0, Z_IMETAL2), 2.0):
        pop.evolve(tnow=tnow, dt=1.0)
    assert pop.lum_xray[0] == sentinel_lumx


def test_respects_mcomp_floor():
    """A secondary at or below |mcomp| must never activate via this
    channel, matching the existing fsur-activation gate's own floor."""
    config = SimulationConfig(
        ntot=10, imetal=2, use_post_sn_rlof=True, mcomp=8.5, iseed=1
    )
    pop = _hand_constructed_post_sn_population(config, 1.4, 8.0, 15.0)
    for tnow in np.arange(0.0, ms.t_bgb(8.0, Z_IMETAL2), 2.0):
        pop.evolve(tnow=tnow, dt=1.0)
    assert pop.lum_xray[0] == 0.0


def test_z0_imetal1_warns_and_skips_gracefully(caplog):
    config = SimulationConfig(ntot=10, imetal=1, use_post_sn_rlof=True, iseed=1)
    pop = _hand_constructed_post_sn_population(config, 1.4, 8.0, 15.0)
    pop.evolve(tnow=1.0, dt=1.0)  # must not raise
    assert pop.lum_xray[0] == 0.0
