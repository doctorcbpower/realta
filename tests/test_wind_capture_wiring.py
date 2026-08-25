"""Wiring tests for use_wind_capture (docs/science/paper1-followup-
prompt.md): BinaryPopulation.evolve()'s Phase 1.5 branch that converts
a still-detached secondary donor's CAK wind into a deterministic
accretion luminosity, using stellar/cak_wind.py and
binaries/wind_capture.py -- both already unit-tested on their own
(tests/test_cak_wind.py, tests/test_wind_capture.py). This file only
tests the *wiring*: that evolve() calls them correctly, with the right
gating, hand-off to use_post_sn_rlof, and Eddington capping.
"""

from __future__ import annotations

from unittest.mock import patch

import numpy as np
import pytest

from realta.binaries.interaction import RLOFOutcome
from realta.binaries.population import BinaryPopulation
from realta.config import SimulationConfig

M_COMPACT = 1.4  # Msun, neutron star
M_DONOR = 24.0  # Msun, Vela-X-1-like OB donor (Friend & Castor 1982 Table 1)
Z_IMETAL2 = 0.008
TNOW = 5.0  # Myr, donor still comfortably on the MS at this mass/Z


def _detached_hmxb_population(config, a_rsun, m2=M_DONOR, m1=M_COMPACT):
    pop = BinaryPopulation(config)
    pop.m1 = np.array([m1])
    pop.m2 = np.array([m2])
    pop.a = np.array([a_rsun / BinaryPopulation.RSUN_PER_AU])
    pop.period = np.array([8.96])
    pop.turnoff_time = np.array([1.0e6])
    pop.t2_lifetime = np.array([1.0e6])
    pop.nturn = np.array([1], dtype=np.int8)
    pop.is_survived = np.ones(1, dtype=bool)
    pop.lum_xray = np.zeros(1)
    pop.did_merge = np.zeros(1, dtype=bool)
    pop.merge_time = np.array([np.nan])
    pop.rlof_time = np.array([np.inf])
    pop.rlof_outcome = np.array([RLOFOutcome.DETACHED], dtype=object)
    pop.rlof_donor_is_star1 = np.array([False])
    pop.rlof_processed = np.ones(1, dtype=bool)
    return pop


def test_config_default_is_disabled():
    config = SimulationConfig(ntot=10)
    assert config.use_wind_capture is False


def test_disabled_is_inert_even_with_a_qualifying_detached_donor():
    config = SimulationConfig(
        ntot=10, imetal=2, use_wind_capture=False, fsur=0.0, iseed=1
    )
    pop = _detached_hmxb_population(config, a_rsun=60.0)
    pop.evolve(tnow=TNOW, dt=1.0)
    assert pop.lum_xray[0] == 0.0


def test_enabled_produces_a_positive_sub_eddington_luminosity():
    config = SimulationConfig(
        ntot=10, imetal=2, use_wind_capture=True, fsur=0.0, iseed=1
    )
    pop = _detached_hmxb_population(config, a_rsun=60.0)
    pop.evolve(tnow=TNOW, dt=1.0)

    assert pop.lum_xray[0] > 0.0
    ledd = pop.xray_calc.eddington_luminosity(M_COMPACT)
    assert pop.lum_xray[0] <= ledd


def test_capture_luminosity_increases_for_a_closer_orbit():
    """Direction/sensitivity check: R_acc/a increases as the orbit
    tightens (both because R_acc/a itself grows and because the wind
    hasn't fully accelerated yet at smaller r), so a closer binary
    should capture more of the wind."""
    config = SimulationConfig(
        ntot=10, imetal=2, use_wind_capture=True, fsur=0.0, iseed=1
    )
    lx_wide = _detached_hmxb_population(config, a_rsun=150.0)
    lx_wide.evolve(tnow=TNOW, dt=1.0)
    lx_close = _detached_hmxb_population(config, a_rsun=30.0)
    lx_close.evolve(tnow=TNOW, dt=1.0)

    assert lx_close.lum_xray[0] > lx_wide.lum_xray[0]


def test_post_sn_rlof_takes_over_once_donor_fills_its_roche_lobe():
    """Wiring/hand-off check: with both channels enabled, a donor that
    already fills its Roche lobe must go through the certain-activation
    use_post_sn_rlof path (xray_calc.get_lumx's draw), not the
    wind-capture formula -- driven by shrinking the separation until
    the donor (radius fixed by its own mass/age) overflows."""
    config = SimulationConfig(
        ntot=10,
        imetal=2,
        use_wind_capture=True,
        use_post_sn_rlof=True,
        fsur=0.0,
        iseed=1,
    )
    # A tight enough orbit that the donor's own MS radius exceeds its
    # Roche lobe -- confirmed by shrinking `a` until this triggers.
    pop = _detached_hmxb_population(config, a_rsun=12.0)
    with patch.object(
        pop.xray_calc, "get_lumx", wraps=pop.xray_calc.get_lumx
    ) as mock_get_lumx:
        pop.evolve(tnow=TNOW, dt=1.0)
    # Confirms the RLOF branch (the only caller of get_lumx in Phase
    # 1.5) fired, not the wind-capture one -- a direct check of which
    # code path ran, rather than trying to bound the stochastic RLOF
    # draw's magnitude (xray_calc.get_lumx's weibull draw is not
    # confined to [lxmin, lxmax] itself -- see xray/luminosity.py).
    mock_get_lumx.assert_called_once()
    assert pop.lum_xray[0] > 0.0


def test_wind_capture_applies_only_while_donor_is_below_its_roche_lobe():
    """Companion check to the hand-off test above: with
    use_post_sn_rlof left OFF, the same tight-orbit (Roche-lobe-
    filling) configuration must NOT activate at all via wind-capture
    -- donor_radius < r_l2 is a hard gate on that branch, not merely a
    preference."""
    config = SimulationConfig(
        ntot=10,
        imetal=2,
        use_wind_capture=True,
        use_post_sn_rlof=False,
        fsur=0.0,
        iseed=1,
    )
    pop = _detached_hmxb_population(config, a_rsun=12.0)
    pop.evolve(tnow=TNOW, dt=1.0)
    assert pop.lum_xray[0] == 0.0


def test_zero_metallicity_skips_gracefully_with_a_warning(caplog):
    config = SimulationConfig(
        ntot=10, imetal=1, use_wind_capture=True, fsur=0.0, iseed=1
    )
    pop = _detached_hmxb_population(config, a_rsun=60.0)
    with caplog.at_level("WARNING"):
        pop.evolve(tnow=TNOW, dt=1.0)
    assert pop.lum_xray[0] == 0.0
    assert any("Z=0" in record.message for record in caplog.records)


def test_does_not_double_activate_once_lum_xray_is_already_set():
    """Same "activate once, hold fixed" convention as the rest of this
    module (fsur, use_post_sn_rlof) -- a subsequent evolve() call must
    not overwrite an already-set lum_xray."""
    config = SimulationConfig(
        ntot=10, imetal=2, use_wind_capture=True, fsur=0.0, iseed=1
    )
    pop = _detached_hmxb_population(config, a_rsun=60.0)
    pop.evolve(tnow=TNOW, dt=1.0)
    first_value = pop.lum_xray[0]
    assert first_value > 0.0

    pop.evolve(tnow=TNOW + 1.0, dt=1.0)
    assert pop.lum_xray[0] == pytest.approx(first_value)


def test_mcomp_floor_prevents_activation_for_a_negligible_donor():
    config = SimulationConfig(
        ntot=10, imetal=2, use_wind_capture=True, fsur=0.0, mcomp=5.0, iseed=1
    )
    pop = _detached_hmxb_population(config, a_rsun=60.0, m2=1.0)
    pop.evolve(tnow=TNOW, dt=1.0)
    assert pop.lum_xray[0] == 0.0


def test_sensitivity_wind_cak_alpha_changes_the_result():
    """Direct sensitivity check requested by this project's testing
    discipline: confirms config.wind_cak_alpha is actually wired
    through to the CAK wind calculation, not silently ignored."""
    config_low = SimulationConfig(
        ntot=10,
        imetal=2,
        use_wind_capture=True,
        fsur=0.0,
        wind_cak_alpha=0.45,
        iseed=1,
    )
    config_high = SimulationConfig(
        ntot=10,
        imetal=2,
        use_wind_capture=True,
        fsur=0.0,
        wind_cak_alpha=0.65,
        iseed=1,
    )
    pop_low = _detached_hmxb_population(config_low, a_rsun=60.0)
    pop_low.evolve(tnow=TNOW, dt=1.0)
    pop_high = _detached_hmxb_population(config_high, a_rsun=60.0)
    pop_high.evolve(tnow=TNOW, dt=1.0)

    assert pop_low.lum_xray[0] != pytest.approx(pop_high.lum_xray[0])
