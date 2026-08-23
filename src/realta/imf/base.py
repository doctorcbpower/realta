from abc import ABC, abstractmethod

import numpy as np


class IMF(ABC):
    """Abstract base class for Initial Mass Functions (IMFs).

    An IMF describes the relative number of stars formed per unit mass
    at birth. Realta uses the CDF form (`cdf`) exclusively, sampled via
    inverse-transform (`sample`) to draw a population of stellar masses.

    config.imf_type selects the concrete IMF (see
    realta.imf.factory.get_imf and SimulationConfig.imf_type):
    1=Salpeter, 2=Kroupa, 3=Chabrier.
    """

    @abstractmethod
    def cdf(self, m: float, mmin: float, mmax: float) -> float:
        """Cumulative distribution function P(<m), normalized on [mmin, mmax].

        m, mmin, mmax: stellar mass in Msun. Returns a value in [0, 1].
        """

    def sample(
        self, n: int, mmin: float, mmax: float, rng: np.random.Generator
    ) -> np.ndarray:
        """Sample n masses from the IMF using inverse transform sampling.

        For each of n uniform draws u ~ U(0,1), bisects `cdf(m, mmin,
        mmax)` on [mmin, mmax] for the mass m with P(<m) = u, to a
        relative tolerance of 1e-10 in m or 100 iterations, whichever
        comes first.

        Returns an array of n masses in Msun.
        """
        u = np.array([rng.random() for _ in range(n)])
        masses = []

        for ui in u:
            low, high = mmin, mmax
            for _ in range(100):
                mid = (low + high) / 2
                cdf_mid = self.cdf(mid, mmin, mmax)
                if cdf_mid < ui:
                    low = mid
                else:
                    high = mid
                if abs(high - low) < 1e-10:
                    break
            masses.append((low + high) / 2)

        return np.array(masses)
