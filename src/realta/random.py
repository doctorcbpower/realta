import numpy as np


class RandomGenerator:
    """Random number generator utilizing NumPy's BitGenerator."""

    def __init__(self, seed: int):
        self.rng = np.random.default_rng(abs(seed))

    def random(self) -> float:
        """Generate a random number between 0 and 1."""
        return self.rng.random()

    def uniform(self, low: float = 0.0, high: float = 1.0) -> float:
        """Generate a uniform random number in [low, high)."""
        return self.rng.uniform(low, high)

    def exponential(self, scale: float = 1.0) -> float:
        """Generate an exponential random number."""
        return self.rng.exponential(scale)
