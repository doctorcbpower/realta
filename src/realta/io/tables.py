from __future__ import annotations

import importlib.resources
import logging
from pathlib import Path
from typing import ClassVar

import numpy as np

logger = logging.getLogger("realta")


class DataTable:
    """Base class for loading and interpolating tabulated data."""

    def __init__(self, data_dir: str | None = None):
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

    def __init__(self, imetal: int = 2, data_dir: str | None = None):
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
        raw_idx = int(np.searchsorted(self.mass, lmass, side="right") - 1)
        idx = max(0, min(raw_idx, len(self.mass) - 2))

        a = (self.lifetime[idx + 1] - self.lifetime[idx]) / (
            self.mass[idx + 1] - self.mass[idx]
        )
        b = self.lifetime[idx] - a * self.mass[idx]
        log_lifetime = a * lmass + b
        return 10.0**log_lifetime


class RemnantTable(DataTable):
    """Remnant mass data table."""

    def __init__(self, data_dir: str | None = None):
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
        raw_idx = int(np.searchsorted(self.minit, lmass, side="right") - 1)
        idx = max(0, min(raw_idx, len(self.minit) - 2))

        a = (self.mfin[idx + 1] - self.mfin[idx]) / (
            self.minit[idx + 1] - self.minit[idx]
        )
        b = self.mfin[idx] - a * self.minit[idx]
        log_mfin = a * lmass + b
        return 10.0**log_mfin


class MSLuminosityTable(DataTable):
    """Population-total main-sequence bolometric luminosity vs time.

    Source: FSPS (Conroy, Gunn & White 2009; python-fsps wrapper by
    Foreman-Mackey et al., MIT licensed), Kroupa IMF, instantaneous-burst
    single stellar population. Generation script:
    notebooks_helper/generate_ms_luminosity_table.py.

    The tabulated values are baked to a fixed fiducial cluster mass of
    FIDUCIAL_CLUSTER_MASS_MSUN (1e6 Msun, matching Power et al. 2009's
    nominal model -- see each ms_lbol_*.dat file's header). get_lbol()
    takes the *actual* total mass formed in the population being
    simulated and rescales linearly: SSP bolometric luminosity scales
    linearly with total mass formed at fixed IMF, metallicity and age.
    Without this rescaling, comparing this table's output directly
    against a simulation run with a different `ntot`/`mmin`/`mmax` (and
    therefore a different total formed mass) silently mixes two
    inconsistent normalizations -- e.g. Realta's own ntot=100_000,
    mmin=0.1-100 Msun default population forms only ~8.9e4 Msun, not
    1e6, so the un-rescaled MS curve would be ~11x too bright relative
    to that population's actual HMXB/accretion luminosity.

    This is a genuinely different quantity from `IonizingPhotonTable`
    (a separate, currently-unused ionising-photon-budget estimate for
    M >= 8 Msun stars only) -- this table is the FULL population's
    total bolometric luminosity (all masses, all evolutionary phases in
    FSPS' isochrones). It is not part of the HMXB evolution loop and is
    not consumed by BinaryPopulation.evolve() -- it exists for
    reproducing Fig. 1 of Power et al. (2009) (the "Luminosity - MS
    lifetime" curve) and similar population-luminosity comparisons.

    Domain of validity: tabulated for 0.1-100 Myr only; get_lbol()
    returns 0.0 outside that range rather than extrapolating (a
    population's bolometric luminosity is not well described by a
    log-log linear extrapolation of a 64-point SSP track). imetal=1
    (Z=0) has no true zero-metallicity match in FSPS's isochrone
    libraries -- that table uses the lowest available FSPS metallicity
    as a documented proxy; see the file's own header and the generation
    script for the actual Z used.
    """

    METAL_FILES: ClassVar[dict[int, str]] = {
        1: "ms_lbol_z0.dat",
        2: "ms_lbol_z8e-3.dat",
        3: "ms_lbol_z2e-2.dat",
    }

    # Cluster mass baked into the tabulated FSPS values (see class
    # docstring and notebooks_helper/generate_ms_luminosity_table.py's
    # CLUSTER_MASS_MSUN). Must match the generation script exactly.
    FIDUCIAL_CLUSTER_MASS_MSUN = 1.0e6

    def __init__(self, imetal: int = 2, data_dir: str | None = None):
        super().__init__(data_dir)
        self.imetal = imetal
        self.log_age = np.array([])
        self.log_lbol = np.array([])
        self.load()

    def load(self):
        filename = self.METAL_FILES.get(self.imetal, self.METAL_FILES[2])
        filepath = self.data_dir / filename

        if not filepath.exists():
            logger.warning(
                f"MS luminosity data file {filepath} not found. "
                "get_lbol() will return 0.0 for all ages."
            )
            return

        try:
            with open(filepath, "r") as f:
                lines = f.readlines()

            data_lines = lines[3:]
            ages = []
            lbols = []

            for line in data_lines:
                if line.strip() and not line.startswith("#"):
                    parts = line.split()
                    if len(parts) >= 2:
                        ages.append(float(parts[0]))
                        lbols.append(float(parts[1]))

            self.log_age = np.log10(np.array(ages))
            self.log_lbol = np.array(lbols)
            self.loaded = True
        except (FileNotFoundError, ValueError, IndexError) as e:
            logger.error(f"Error loading MS luminosity data: {e}")

    def get_lbol(self, age_myr: float, total_mass_msun: float) -> float:
        """Total population bolometric luminosity at a given age, erg/s.

        total_mass_msun: total stellar mass actually formed in the
        population being simulated -- the sum of ALL sampled IMF masses
        (e.g. BinaryPopulation.total_mass_msun), not just the M >= mcut
        subset. Required to rescale the table's fiducial-1e6-Msun values
        to the population actually being compared against -- see the
        class docstring for why this matters.

        Returns 0.0 for ages outside the tabulated range (including
        age_myr <= 0) or if the table failed to load, rather than
        extrapolating -- see the class docstring.
        """
        if not self.loaded or len(self.log_age) < 2 or age_myr <= 0.0:
            return 0.0

        lage = np.log10(age_myr)
        if lage < self.log_age[0] or lage > self.log_age[-1]:
            return 0.0

        raw_idx = int(np.searchsorted(self.log_age, lage, side="right") - 1)
        idx = max(0, min(raw_idx, len(self.log_age) - 2))

        a = (self.log_lbol[idx + 1] - self.log_lbol[idx]) / (
            self.log_age[idx + 1] - self.log_age[idx]
        )
        b = self.log_lbol[idx] - a * self.log_age[idx]
        log_lbol = a * lage + b
        fiducial_lbol = 10.0**log_lbol
        return fiducial_lbol * (total_mass_msun / self.FIDUCIAL_CLUSTER_MASS_MSUN)


class UVLuminosityTable(DataTable):
    """Population-total far-UV (GALEX FUV, ~1528 A) luminosity vs time.

    Paper 1's L_UV(t) observable (docs/science/research-programme.md,
    Figs. 1-2) -- like `MSLuminosityTable`, this is this session's own
    addition, not part of either Power et al. paper. Source: FSPS,
    Kroupa IMF, instantaneous-burst single stellar population, GALEX
    FUV band. Generation script: scripts/generate_fuv_luminosities.py
    (band-choice rationale and the m_AB -> nu*L_nu conversion are
    documented there; band decision reviewed in
    docs/science/paper1-binary-interaction-proposal.md).

    Same fiducial-mass-then-rescale convention as `MSLuminosityTable`:
    tabulated values are baked to FIDUCIAL_CLUSTER_MASS_MSUN (1e6 Msun),
    and `get_luv()` rescales linearly to the actual population's
    `total_mass_msun`. Same domain-of-validity note applies (0.1-100
    Myr; no extrapolation; imetal=1 uses FSPS's lowest available
    metallicity as a documented Z=0 proxy).

    The `fuv_lbol_z*.dat` data files (generated via FSPS + SPS_HOME)
    now exist in `src/realta/data/` -- see docs/provenance.md Section
    7. This class still degrades gracefully exactly like
    `MSLuminosityTable` does for a missing file (e.g. a custom
    `data_dir` that doesn't have them): `get_luv()` returns 0.0 rather
    than raising.
    """

    METAL_FILES: ClassVar[dict[int, str]] = {
        1: "fuv_lbol_z0.dat",
        2: "fuv_lbol_z8e-3.dat",
        3: "fuv_lbol_z2e-2.dat",
    }

    # Must match scripts/generate_fuv_luminosities.py's CLUSTER_MASS_MSUN.
    FIDUCIAL_CLUSTER_MASS_MSUN = 1.0e6

    def __init__(self, imetal: int = 2, data_dir: str | None = None):
        super().__init__(data_dir)
        self.imetal = imetal
        self.log_age = np.array([])
        self.log_luv = np.array([])
        self.load()

    def load(self):
        filename = self.METAL_FILES.get(self.imetal, self.METAL_FILES[2])
        filepath = self.data_dir / filename

        if not filepath.exists():
            logger.warning(
                f"UV luminosity data file {filepath} not found. "
                "get_luv() will return 0.0 for all ages -- run "
                "scripts/generate_fuv_luminosities.py to generate it "
                "(requires FSPS + SPS_HOME)."
            )
            return

        try:
            with open(filepath, "r") as f:
                lines = f.readlines()

            data_lines = lines[3:]
            ages = []
            luvs = []

            for line in data_lines:
                if line.strip() and not line.startswith("#"):
                    parts = line.split()
                    if len(parts) >= 2:
                        ages.append(float(parts[0]))
                        luvs.append(float(parts[1]))

            self.log_age = np.log10(np.array(ages))
            self.log_luv = np.array(luvs)
            self.loaded = True
        except (FileNotFoundError, ValueError, IndexError) as e:
            logger.error(f"Error loading UV luminosity data: {e}")

    def get_luv(self, age_myr: float, total_mass_msun: float) -> float:
        """Total population far-UV luminosity at a given age, erg/s.

        Same rescaling/domain-of-validity behaviour as
        `MSLuminosityTable.get_lbol()` -- see that method's docstring.
        """
        if not self.loaded or len(self.log_age) < 2 or age_myr <= 0.0:
            return 0.0

        lage = np.log10(age_myr)
        if lage < self.log_age[0] or lage > self.log_age[-1]:
            return 0.0

        raw_idx = int(np.searchsorted(self.log_age, lage, side="right") - 1)
        idx = max(0, min(raw_idx, len(self.log_age) - 2))

        a = (self.log_luv[idx + 1] - self.log_luv[idx]) / (
            self.log_age[idx + 1] - self.log_age[idx]
        )
        b = self.log_luv[idx] - a * self.log_age[idx]
        log_luv = a * lage + b
        fiducial_luv = 10.0**log_luv
        return fiducial_luv * (total_mass_msun / self.FIDUCIAL_CLUSTER_MASS_MSUN)


class IonizingPhotonTable(DataTable):
    """Ionizing photon data table."""

    MUNIT = 1.99e30
    MATOM = 1.67e-27

    def __init__(self, data_dir: str | None = None):
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
        raw_idx = int(np.searchsorted(self.mstar, lmass, side="right") - 1)
        idx = max(0, min(raw_idx, len(self.mstar) - 2))

        a = (self.ngamma[idx + 1] - self.ngamma[idx]) / (
            self.mstar[idx + 1] - self.mstar[idx]
        )
        b = self.ngamma[idx] - a * self.mstar[idx]
        return a * lmass + b
