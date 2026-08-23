import numpy as np

from realta.imf.chabrier import ChabrierIMF
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
