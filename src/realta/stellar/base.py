from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np


class LifetimeModel(ABC):
    @abstractmethod
    def lifetime(self, mass: float | np.ndarray) -> float | np.ndarray:
        """Stellar lifetime in Myr."""


class RemnantModel(ABC):
    @abstractmethod
    def remnant_mass(self, mass: float | np.ndarray) -> float | np.ndarray:
        """Remnant mass in solar masses."""


class IonisationModel(ABC):
    @abstractmethod
    def photon_rate(self, mass: float | np.ndarray) -> float | np.ndarray:
        """Log10 ionizing photons per unit solar mass."""
