"""Unit tests for BinaryPopulation.generate_population()'s orbital
parameters (period, semi-major axis) -- closes the gap noted in
docs/provenance.md: these never feed `floss`/survival, so they are the
one part of population generation the pinned integration test
(tests/test_regression.py) genuinely cannot see (a bug here would not
move `lumx_tot`/`nphot_tot`/`nactive`/`ndead` at all).
"""

import numpy as np
import pytest

from realta.binaries.population import BinaryPopulation
from realta.config import SimulationConfig


def test_orbital_period_bounds_and_log_uniform_distribution():
    """Period must be log-flat between pmin and pmax (make_stars.f).

    Bounds are checked exactly; the log-uniform shape is checked via the
    sample mean of log(period), which should sit near the midpoint of
    [log(pmin), log(pmax)] within a wide (5 sigma) tolerance from the
    CLT for a uniform distribution -- a real distribution bug (e.g. a
    linear instead of log-uniform draw, or a shape skewed toward one
    bound) would move the mean far outside this, while ordinary sampling
    noise won't.
    """
    config = SimulationConfig(ntot=10_000, pmin=1.0, pmax=1000.0, iseed=7)
    pop = BinaryPopulation(config)
    assert len(pop.period) > 50  # sanity: population actually formed binaries

    assert np.all(pop.period >= config.pmin)
    assert np.all(pop.period <= config.pmax)

    log_p = np.log(pop.period)
    expected_mean = (np.log(config.pmin) + np.log(config.pmax)) / 2.0
    n = len(pop.period)
    sigma_uniform = (np.log(config.pmax) - np.log(config.pmin)) / np.sqrt(12.0)
    tolerance = 5.0 * sigma_uniform / np.sqrt(n)
    assert abs(log_p.mean() - expected_mean) < tolerance


def test_semi_major_axis_matches_reference_formula():
    """a = AFAC * m1^(1/3) * (1 + m2/m1)^(1/3) * period^(2/3).

    Reference: main.f, `a(n)=afac*mass(1,n)**(1./3.)*(1.+mass(2,n)
    /mass(1,n))**(1./3.)` then `a(n)=a(n)*period(n)**(2./3.)` --
    AFAC=0.0193852859 is the reference's own literal hardcoded constant
    (in AU), not independently re-derived from Kepler's third law here:
    a from-scratch physical derivation (G=4*pi^2 AU^3/(Msun*yr^2), using
    PFAC=365.229126 as the days-per-year conversion) gives a
    ~1% *different* prefactor (~0.019571), which is the reference's own
    precision/convention, not a bug in either constant -- matching the
    reference's literal value is what matters here, not re-deriving
    physics independently.
    """
    config = SimulationConfig(ntot=5000, iseed=7)
    pop = BinaryPopulation(config)
    assert len(pop.m1) > 10

    assert pop.AFAC == pytest.approx(0.0193852859)

    q = pop.m2 / pop.m1
    expected_a = (
        pop.AFAC
        * (pop.m1 ** (1.0 / 3.0))
        * ((1.0 + q) ** (1.0 / 3.0))
        * (pop.period ** (2.0 / 3.0))
    )
    np.testing.assert_allclose(pop.a, expected_a, rtol=1e-12)
