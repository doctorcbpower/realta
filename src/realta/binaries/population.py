import logging
import numpy as np

from realta.binaries.binary import Binary
from realta.config import SimulationConfig
from realta.imf.factory import get_imf
from realta.io.tables import IonizingPhotonTable, LifetimeTable, RemnantTable
from realta.random import RandomGenerator
from realta.xray.luminosity import XRayLuminosity

logger = logging.getLogger("realta")


class BinaryPopulation:
    """Manages a population of binary stars."""

    PFAC = 365.229126
    AFAC = 0.0193852859

    def __init__(self, config: SimulationConfig):
        self.config = config
        self.binaries: list[Binary] = []
        self.rng = RandomGenerator(config.iseed)
        self.imf = get_imf(config.imf_type)

        self.lifetime_table = LifetimeTable(config.imetal, config.data_dir)
        self.remnant_table = RemnantTable(config.data_dir)
        self.ionizing_table = IonizingPhotonTable(config.data_dir)

        self.xray_calc = XRayLuminosity(
            lxmin=10.0 ** (config.lxmin - np.log10(config.lunit)),
            lxmax=10.0 ** (config.lxmax - np.log10(config.lunit)),
            lunit=config.lunit,
        )

        self.generate_population()

    def generate_population(self):
        ntot = self.config.ntot
        mmin = self.config.mmin
        mmax = self.config.mmax
        mcut = self.config.mcut
        pmin = self.config.pmin
        pmax = self.config.pmax
        mcomp = self.config.mcomp
        fbin = self.config.fbin

        logger.info(f"Generating {ntot} stars with IMF type {self.config.imf_type}")

        masses = self.imf.sample(ntot, mmin, mmax, self.rng)

        primary_masses = []
        secondary_masses = []
        periods = []
        semi_major_axes = []

        nmass = 0
        nhmxb = 0

        for mass in masses:
            primary_masses.append(mass)

            if self.rng.random() <= fbin and mass >= mcut:
                nmass += 1
                log_period = np.log(pmax / pmin) * self.rng.random()
                period = np.exp(log_period)
                periods.append(period)

                if mcomp < 0:
                    companion_mass = mmin + (mass - mmin) * self.rng.random()
                else:
                    companion_mass = mcomp + (mass - mcomp) * self.rng.random()

                companion_mass = min(companion_mass, mass)
                secondary_masses.append(companion_mass)

                if companion_mass >= abs(mcomp):
                    nhmxb += 1

                a_val = (
                    self.AFAC
                    * mass ** (1.0 / 3.0)
                    * (1.0 + companion_mass / mass) ** (1.0 / 3.0)
                )
                a_val *= period ** (2.0 / 3.0)
                semi_major_axes.append(a_val)
            else:
                secondary_masses.append(0.0)
                periods.append(0.0)
                semi_major_axes.append(0.0)

        indices = np.argsort(-np.array(primary_masses))

        self.binaries = []
        for i in indices:
            binary = Binary(
                primary_mass=primary_masses[i],
                secondary_mass=secondary_masses[i],
                period=periods[i],
                a=semi_major_axes[i],
                index=i,
            )

            if primary_masses[i] >= self.config.mcut:
                binary.turnoff_time = self.lifetime_table.get_lifetime(
                    primary_masses[i]
                )

            self.binaries.append(binary)

        self.binaries = [b for b in self.binaries if b.primary_mass >= mcut]

        logger.info(
            f"Generated {len(self.binaries)} massive binaries, {nhmxb} HMXB progenitors"
        )

    def sort_by_turnoff_time(self):
        self.binaries.sort(key=lambda b: b.turnoff_time)
        for i, b in enumerate(self.binaries):
            b.index = i

    def evolve(self, tnow: float, dt: float) -> tuple[float, float, int, int]:
        lumx_tot = 0.0
        nphot_tot = 0.0
        nactive = 0
        ndead = 0

        self.sort_by_turnoff_time()

        for binary in self.binaries:
            binary.lum_xray = 0.0

            if tnow < binary.turnoff_time:
                continue

            if binary.turnoff_time == 0.0:
                ndead += 1
                continue

            if binary.nturn == 0:
                primary_remnant = self.remnant_table.get_remnant_mass(
                    binary.primary_mass
                )
                deltam = binary.primary_mass - primary_remnant
                floss = deltam / (binary.primary_mass + binary.secondary_mass)

                binary.primary_mass = primary_remnant
                binary.turnoff_time = self.lifetime_table.get_lifetime(
                    binary.secondary_mass
                )
                binary.nturn = 1
                nactive += 1

            elif binary.nturn == 1:
                secondary_remnant = self.remnant_table.get_remnant_mass(
                    binary.secondary_mass
                )
                deltam = binary.secondary_mass - secondary_remnant
                floss = deltam / (binary.primary_mass + binary.secondary_mass)

                binary.secondary_mass = secondary_remnant
                binary.turnoff_time = 0.0
                binary.nturn = 2
                ndead += 1
                continue

            if floss <= 0.5:
                binary.a *= deltam / (binary.primary_mass + binary.secondary_mass)
                binary.period = self.PFAC * np.sqrt(
                    binary.a**3 / (binary.primary_mass + binary.secondary_mass)
                )

            if (
                binary.turnoff_time > 0.0
                and floss <= 0.5
                and binary.secondary_mass >= abs(self.config.mcomp)
                and self.rng.random() <= self.config.fbin
            ):
                binary.lum_xray = self.xray_calc.get_lumx(
                    binary.primary_mass,
                    binary.secondary_mass,
                    binary.period,
                    binary.a,
                    iseed=None,
                    use_weibull=True,
                )

            lumx_tot += binary.lum_xray

            if binary.nturn == 0:
                ng1 = self.ionizing_table.get_ngamma(binary.primary_mass)
                ng2 = self.ionizing_table.get_ngamma(binary.secondary_mass)
                t1 = self.lifetime_table.get_lifetime(binary.primary_mass)
                t2 = self.lifetime_table.get_lifetime(binary.secondary_mass)

                nphot_tot += 10.0 ** (
                    ng1 + np.log10(dt / t1) + ng2 + np.log10(dt / t2) - 60
                )
            elif binary.nturn == 1:
                ng2 = self.ionizing_table.get_ngamma(binary.secondary_mass)
                t2 = self.lifetime_table.get_lifetime(binary.secondary_mass)

                nphot_tot += 10.0 ** (ng2 + np.log10(dt / t2) - 60)

        return lumx_tot, nphot_tot, nactive, ndead
