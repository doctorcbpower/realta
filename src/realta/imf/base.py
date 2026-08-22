from abc import ABC, abstractmethod

import numpy as np

from realta.random import RandomGenerator


class IMF(ABC):
    """Abstract base class for Initial Mass Functions."""

    @abstractmethod
    def cdf(self, m: float, mmin: float, mmax: float) -> float:
        """Cumulative distribution function P(<m)."""

    def sample(
        self, n: int, mmin: float, mmax: float, rng: RandomGenerator
    ) -> np.ndarray:
        """Sample n masses from the IMF using inverse transform sampling."""
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
