import numpy as np
import pytest

from realta.imf.chabrier import ChabrierIMF
from realta.imf.factory import get_imf
from realta.imf.kroupa import KroupaIMF
from realta.imf.salpeter import SalpeterIMF


def test_kroupa_cdf():
    imf = KroupaIMF()
    # CDF should evaluate to 0 at lower bound and 1 at upper bound
    assert np.isclose(imf.cdf(0.01, mmin=0.01, mmax=100.0), 0.0)
    assert np.isclose(imf.cdf(100.0, mmin=0.01, mmax=100.0), 1.0)


def test_kroupa_sampling():
    imf = KroupaIMF()
    rng = np.random.default_rng(42)
    samples = imf.sample(1000, mmin=0.01, mmax=100.0, rng=rng)

    assert len(samples) == 1000
    assert np.all(samples >= 0.01)
    assert np.all(samples <= 100.0)


def test_salpeter_cdf():
    # config.imf_type=1 -- see imf/salpeter.py for the Salpeter (1955)
    # provenance and reference-Fortran (salpeter.f) match.
    imf = SalpeterIMF()
    assert np.isclose(imf.cdf(0.01, mmin=0.01, mmax=100.0), 0.0)
    assert np.isclose(imf.cdf(100.0, mmin=0.01, mmax=100.0), 1.0)
    # CDF must be monotonically nondecreasing
    masses = np.linspace(0.01, 100.0, 50)
    cdf_vals = [imf.cdf(m, mmin=0.01, mmax=100.0) for m in masses]
    assert np.all(np.diff(cdf_vals) >= 0)


def test_salpeter_sampling():
    imf = SalpeterIMF()
    rng = np.random.default_rng(42)
    samples = imf.sample(1000, mmin=0.01, mmax=100.0, rng=rng)

    assert len(samples) == 1000
    assert np.all(samples >= 0.01)
    assert np.all(samples <= 100.0)


def test_get_imf_salpeter_slope_override_reaches_alpha():
    """config.imf_slope (A4, docs/science/paper1-detailed-work-
    breakdown.md) overrides SalpeterIMF's own alpha via
    imf/factory.py::get_imf's `slope` parameter -- confirms it's
    actually wired through, not a dead parameter.
    """
    imf = get_imf(imf_type=1, slope=1.8)
    assert isinstance(imf, SalpeterIMF)
    assert imf.alpha == 1.8


def test_get_imf_salpeter_no_slope_matches_default():
    imf_default = get_imf(imf_type=1)
    imf_explicit_none = get_imf(imf_type=1, slope=None)
    assert imf_default.alpha == imf_explicit_none.alpha == 2.35


def test_get_imf_slope_ignored_for_non_salpeter():
    """slope only applies to imf_type=1 (Salpeter) -- Kroupa/Chabrier
    have no single power-law slope to override, per get_imf's own
    docstring; passing slope for them must be silently ignored, not
    raise or otherwise error.
    """
    kroupa = get_imf(imf_type=2, slope=1.8)
    assert isinstance(kroupa, KroupaIMF)
    chabrier = get_imf(imf_type=3, slope=1.8)
    assert isinstance(chabrier, ChabrierIMF)


def test_salpeter_shallower_slope_favours_higher_mass_sampling():
    """Sanity/direction check: a shallower slope (smaller alpha) means
    relatively more high-mass stars -- the sampled mean mass should
    increase as alpha decreases, holding mmin/mmax/seed fixed."""
    rng_steep = np.random.default_rng(7)
    rng_shallow = np.random.default_rng(7)
    steep = get_imf(imf_type=1, slope=3.0).sample(2000, 0.5, 100.0, rng_steep)
    shallow = get_imf(imf_type=1, slope=1.5).sample(2000, 0.5, 100.0, rng_shallow)
    assert shallow.mean() > steep.mean()


def test_get_imf_slope_singularity_raises_uncaught():
    """The imf_slope=1.0 rejection lives in SimulationConfig.__post_init__
    (SalpeterIMF's CDF denominator, mmax**0 - mmin**0, is identically
    0.0 there) -- get_imf() itself does not validate slope, so calling
    it directly with slope=1.0 constructs a SalpeterIMF whose cdf()
    raises ZeroDivisionError (0.0/0.0 in plain Python floats), a
    documented pre-existing edge case now guarded against at the
    config boundary instead (see
    tests/test_binary_prescriptions.py::test_out_of_range_interaction_params_rejected)."""
    imf = get_imf(imf_type=1, slope=1.0)
    with pytest.raises(ZeroDivisionError):
        imf.cdf(10.0, mmin=0.5, mmax=100.0)


def test_chabrier_cdf():
    # config.imf_type=3 -- see imf/chabrier.py for the Chabrier (2003)
    # provenance and reference-Fortran (log_normal_IMF.f) match.
    imf = ChabrierIMF()
    assert np.isclose(imf.cdf(0.01, mmin=0.01, mmax=100.0), 0.0)
    assert np.isclose(imf.cdf(100.0, mmin=0.01, mmax=100.0), 1.0)
    masses = np.linspace(0.01, 100.0, 50)
    cdf_vals = [imf.cdf(m, mmin=0.01, mmax=100.0) for m in masses]
    assert np.all(np.diff(cdf_vals) >= 0)


def test_chabrier_sampling():
    imf = ChabrierIMF()
    rng = np.random.default_rng(42)
    samples = imf.sample(1000, mmin=0.01, mmax=100.0, rng=rng)

    assert len(samples) == 1000
    assert np.all(samples >= 0.01)
    assert np.all(samples <= 100.0)


def test_imf_slope_pinned_population_regression():
    """Numeric regression pin for a non-default imf_slope (A4), at
    population-generation scale -- see docs/physics/binary-sampling.md.
    Values captured by an actual run, cross-checked for run-to-run
    determinism before being pinned, not hand-derived analytically,
    following tests/test_regression.py's own discipline.
    """
    from realta.binaries.population import BinaryPopulation
    from realta.config import SimulationConfig

    config = SimulationConfig(
        ntot=20_000,
        imf_type=1,
        imf_slope=1.8,
        iseed=42,
        mmin=0.1,
        mmax=100.0,
        mcut=8.0,
    )
    pop = BinaryPopulation(config)

    assert pop.total_mass_msun == pytest.approx(24970.24511970111, rel=1e-9)
    assert len(pop.m1) == 532
