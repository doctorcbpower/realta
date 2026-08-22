from __future__ import annotations

import importlib.resources
import logging
from pathlib import Path
import numpy as np
from typing import ClassVar

logger = logging.getLogger("realta")


class DataTable:
    """Base class for loading and interpolating tabulated data."""

    def __init__(self, data_dir: str = None):
        if data_dir is None:
            self.data_dir = importlib.resources.files("realta") / "data"
        else:
            self.data_dir = Path(data_dir)
        self.loaded = False

    def load(self):
        raise NotImplementedError


class LifetimeTable(DataTable):
    """Stellar lifetime data table."""

    METAL_FILES: ClassVar[dict[int, str]] = {
        1: "lifetimes_z0.dat",
        2: "lifetimes_z8e-3.dat",
        3: "lifetimes_z2e-2.dat",
    }

    def __init__(self, imetal: int = 2, data_dir: str = None):
        super().__init__(data_dir)
        self.imetal = imetal
        self.mass = np.array([])
        self.lifetime = np.array([])
        self.load()

    def load(self):
        filename = self.METAL_FILES.get(self.imetal, self.METAL_FILES[2])
        filepath = self.data_dir / filename

        if not filepath.exists():
            logger.warning(
                f"Lifetime data file {filepath} not found. Using placeholder data."
            )
            self._create_placeholder_data()
            return

        try:
            with open(filepath, "r") as f:
                lines = f.readlines()

            data_lines = lines[6:]
            masses = []
            lifetimes = []

            for line in data_lines:
                if line.strip() and not line.startswith("#"):
                    parts = line.split()
                    if len(parts) >= 2:
                        masses.append(float(parts[0]))
                        lifetimes.append(float(parts[1]))

            self.mass = np.log10(np.array(masses))
            self.lifetime = np.log10(np.array(lifetimes))
            self.loaded = True
        except (FileNotFoundError, ValueError, IndexError) as e:
            logger.error(f"Error loading lifetime data: {e}")
            self._create_placeholder_data()

    def _create_placeholder_data(self):
        masses = np.logspace(-1, 2, 100)
        lifetimes = 10.0 * masses ** (-2.5)
        self.mass = np.log10(masses)
        self.lifetime = np.log10(lifetimes)
        self.loaded = True

    def get_lifetime(self, star_mass: float) -> float:
        if not self.loaded or len(self.mass) == 0:
            self._create_placeholder_data()

        lmass = np.log10(max(star_mass, 1e-10))
        idx = np.searchsorted(self.mass, lmass, side="right") - 1
        idx = max(0, min(idx, len(self.mass) - 2))

        a = (self.lifetime[idx + 1] - self.lifetime[idx]) / (
            self.mass[idx + 1] - self.mass[idx]
        )
        b = self.lifetime[idx] - a * self.mass[idx]
        log_lifetime = a * lmass + b
        return 10.0**log_lifetime


class RemnantTable(DataTable):
    """Remnant mass data table."""

    def __init__(self, data_dir: str = None):
        super().__init__(data_dir)
        self.minit = np.array([])
        self.mfin = np.array([])
        self.load()

    def load(self):
        filepath = self.data_dir / "remnant_masses.dat"

        if not filepath.exists():
            logger.warning(
                f"Remnant data file {filepath} not found. Using placeholder data."
            )
            self._create_placeholder_data()
            return

        try:
            with open(filepath, "r") as f:
                lines = f.readlines()

            data_lines = lines[5:]
            minit = []
            mfin = []

            for line in data_lines:
                if line.strip() and not line.startswith("#"):
                    parts = line.split()
                    if len(parts) >= 2:
                        minit.append(float(parts[0]))
                        mfin.append(float(parts[1]))

            self.minit = np.log10(np.array(minit))
            self.mfin = np.log10(np.clip(np.array(mfin), 1e-10, None))
            self.loaded = True
        except (FileNotFoundError, ValueError, IndexError) as e:
            logger.error(f"Error loading remnant data: {e}")
            self._create_placeholder_data()

    def _create_placeholder_data(self):
        minit = np.logspace(0, 2, 100)
        mfin = np.zeros_like(minit)
        for i, m in enumerate(minit):
            if m < 8:
                mfin[i] = 1e-10
            elif m < 20:
                mfin[i] = 1.4
            elif m < 40:
                mfin[i] = 5.0
            else:
                mfin[i] = 10.0

        self.minit = np.log10(minit)
        self.mfin = np.log10(mfin)
        self.loaded = True

    def get_remnant_mass(self, star_mass: float) -> float:
        if not self.loaded or len(self.minit) == 0:
            self._create_placeholder_data()

        lmass = np.log10(max(star_mass, 1e-10))
        idx = np.searchsorted(self.minit, lmass, side="right") - 1
        idx = max(0, min(idx, len(self.minit) - 2))

        a = (self.mfin[idx + 1] - self.mfin[idx]) / (
            self.minit[idx + 1] - self.minit[idx]
        )
        b = self.mfin[idx] - a * self.minit[idx]
        log_mfin = a * lmass + b
        return 10.0**log_mfin


class IonizingPhotonTable(DataTable):
    """Ionizing photon data table."""

    MUNIT = 1.99e30
    MATOM = 1.67e-27

    def __init__(self, data_dir: str = None):
        super().__init__(data_dir)
        self.mstar = np.array([])
        self.ngamma = np.array([])
        self.load()

    def load(self):
        filepath = self.data_dir / "ionise.dat"

        if not filepath.exists():
            logger.warning(
                f"Ionizing photon data file {filepath} not found. Using placeholder data."
            )
            self._create_placeholder_data()
            return

        try:
            with open(filepath, "r") as f:
                lines = f.readlines()

            data_lines = lines[1:]
            mstar = []
            ngamma_raw = []

            for line in data_lines:
                if line.strip() and not line.startswith("#"):
                    parts = line.split()
                    if len(parts) >= 2:
                        mstar.append(float(parts[0]))
                        ngamma_raw.append(float(parts[1]))

            mstar_arr = np.array(mstar)
            ngamma_arr = np.clip(np.array(ngamma_raw), 1e-30, None)

            self.mstar = np.log10(mstar_arr)
            self.ngamma = (
                np.log10(ngamma_arr)
                + np.log10(mstar_arr)
                + np.log10(self.MUNIT)
                - np.log10(self.MATOM)
            )
            self.loaded = True
        except (FileNotFoundError, ValueError, IndexError) as e:
            logger.error(f"Error loading ionizing photon data: {e}")
            self._create_placeholder_data()

    def _create_placeholder_data(self):
        mstar = np.logspace(0, 2, 100)
        ngamma_raw = np.zeros_like(mstar)
        for i, m in enumerate(mstar):
            if m >= 8:
                ngamma_raw[i] = m**2 * 1e48
            else:
                ngamma_raw[i] = 1e-30

        self.mstar = np.log10(mstar)
        self.ngamma = (
            np.log10(ngamma_raw)
            + np.log10(mstar)
            + np.log10(self.MUNIT)
            - np.log10(self.MATOM)
        )
        self.loaded = True

    def get_ngamma(self, star_mass: float) -> float:
        if not self.loaded or len(self.mstar) == 0:
            self._create_placeholder_data()

        if star_mass < 8.0:
            return -10.0

        lmass = np.log10(max(star_mass, 1e-10))
        idx = np.searchsorted(self.mstar, lmass, side="right") - 1
        idx = max(0, min(idx, len(self.mstar) - 2))

        a = (self.ngamma[idx + 1] - self.ngamma[idx]) / (
            self.mstar[idx + 1] - self.mstar[idx]
        )
        b = self.ngamma[idx] - a * self.mstar[idx]
        return a * lmass + b
