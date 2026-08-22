import numpy as np
from realta.imf.kroupa import KroupaIMF


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
