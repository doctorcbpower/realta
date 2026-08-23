"""Unit tests for XRayLuminosity, covering the "uniform" distribution
path directly -- see docs/provenance.md Section 3: this branch was
previously only reachable indirectly (config.xray_distribution="uniform"
is not exercised by any pinned regression config, which use "weibull"
throughout), and was not itself tested in isolation.
"""

import numpy as np
import pytest

from realta.xray.luminosity import XRayLuminosity


def test_uniform_distribution_bounds_and_log_uniform_shape():
    """distribution="uniform": flat log-uniform draw on [lxmin, lxmax].

    get_lumx.f's dead debug branch (only reachable at iseed==-1, never
    hit by a real reference run -- see the class docstring); Realta
    exposes it as an explicit opt-in rather than reproducing it as
    reachable-but-dead code.
    """
    xray = XRayLuminosity(lxmin=1.0, lxmax=100.0, lunit=1.0e33, distribution="uniform")
    rng = np.random.default_rng(11)
    draws = np.array(
        [xray.get_lumx(20.0, 5.0, 10.0, 1.0, rng=rng) for _ in range(5000)]
    )

    assert np.all(draws >= xray.lxmin)
    assert np.all(draws <= xray.lxmax)

    log_draws = np.log10(draws)
    expected_mean = (np.log10(xray.lxmin) + np.log10(xray.lxmax)) / 2.0
    n = len(draws)
    sigma_uniform = (np.log10(xray.lxmax) - np.log10(xray.lxmin)) / np.sqrt(12.0)
    tolerance = 5.0 * sigma_uniform / np.sqrt(n)
    assert abs(log_draws.mean() - expected_mean) < tolerance


def test_uniform_distribution_is_not_eddington_limited_unlike_weibull():
    """The "uniform" branch has no Eddington-limit rejection step.

    This is a real, documented asymmetry in get_lumx() -- "weibull"
    rejection-samples below the Eddington luminosity, but "uniform"
    does not (see xray/luminosity.py: get_lumx()'s uniform branch
    returns directly, with no `if lumx <= ledd` check). Confirmed here
    by choosing an artificially low primary mass so its Eddington limit
    falls inside [lxmin, lxmax] -- if a future change accidentally made
    "uniform" Eddington-limited too (or "weibull" stopped being so),
    this test would catch it.
    """
    massp = 1.0e-4  # synthetic value chosen only to put ledd inside [lxmin, lxmax]
    xray_uniform = XRayLuminosity(
        lxmin=1.0, lxmax=1.0e4, lunit=1.0e33, distribution="uniform"
    )
    xray_weibull = XRayLuminosity(
        lxmin=1.0, lxmax=1.0e4, lunit=1.0e33, distribution="weibull"
    )
    ledd = xray_uniform.eddington_luminosity(massp)
    assert 1.0 < ledd < 1.0e4  # sanity: ledd actually falls inside the draw range

    rng = np.random.default_rng(11)
    uniform_draws = np.array(
        [xray_uniform.get_lumx(massp, 5.0, 10.0, 1.0, rng=rng) for _ in range(2000)]
    )
    rng = np.random.default_rng(11)
    weibull_draws = np.array(
        [xray_weibull.get_lumx(massp, 5.0, 10.0, 1.0, rng=rng) for _ in range(2000)]
    )

    assert np.any(uniform_draws > ledd)  # uncapped -- some draws exceed ledd
    assert np.all(weibull_draws <= ledd)  # Eddington-limited by rejection sampling


def test_lxmin_equals_lxmax_returns_fixed_value_for_both_distributions():
    """lxmin == lxmax is a degenerate case handled before either branch."""
    for distribution in ("uniform", "weibull"):
        xray = XRayLuminosity(
            lxmin=5.0, lxmax=5.0, lunit=1.0e33, distribution=distribution
        )
        rng = np.random.default_rng(1)
        assert xray.get_lumx(20.0, 5.0, 10.0, 1.0, rng=rng) == 5.0


def test_unknown_distribution_rejected():
    with pytest.raises(ValueError):
        XRayLuminosity(lxmin=1.0, lxmax=100.0, distribution="not_a_real_distribution")
