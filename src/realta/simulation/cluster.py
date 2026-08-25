from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import ClassVar

import numpy as np

from realta.binaries.population import BinaryPopulation
from realta.config import SimulationConfig, load_config
from realta.io.tables import IonizingPhotonTable, MSLuminosityTable, UVLuminosityTable

logger = logging.getLogger("realta")


class ClusterSimulation:
    """Main simulation orchestration engine."""

    IMF_MAP: ClassVar[dict[int, str]] = {
        1: "Salpeter",
        2: "Kroupa",
        3: "Chabrier",
    }

    # Julian year * 1e6 -- converts LifetimeTable's Myr to seconds for
    # the Q_H(m) rate calculation below (A3).
    MYR_TO_SECONDS = 3.1557e13

    def __init__(self, config: SimulationConfig | None = None):
        self.config = config if config is not None else load_config()
        self.population: BinaryPopulation | None = None
        # A2/A3 (docs/science/paper1-detailed-work-breakdown.md):
        # population-level L_bol/L_UV (FSPS SSP tables) and the
        # massive-star population's own ionizing-photon output,
        # wired into run()'s per-timestep output alongside
        # lumx_tot/nphot_tot -- see run()'s own comments for why these
        # live here rather than inside BinaryPopulation.evolve().
        self.ms_table = MSLuminosityTable(self.config.imetal, self.config.data_dir)
        self.uv_table = UVLuminosityTable(self.config.imetal, self.config.data_dir)
        self.ionizing_table = IonizingPhotonTable(self.config.data_dir)

    def initialize(self):
        logger.info("Initializing simulation...")
        self.population = BinaryPopulation(self.config)
        logger.info("Simulation initialized.")

    def run(self, output_dir: str = "output") -> list[dict]:
        if self.population is None:
            self.initialize()

        assert self.population is not None

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        self._write_initial_conditions(output_path)

        tmax = self.config.tmax
        dt = self.config.dt

        # Calculate discrete steps to avoid floating-point drift
        num_steps = math.floor(tmax / dt) + 1
        results = []

        logger.info(
            f"Starting time evolution to {tmax} Myr with dt={dt} Myr ({num_steps} steps)"
        )

        for step in range(num_steps):
            tnow = step * dt
            lumx_tot, nphot_tot, nactive, ndead = self.population.evolve(tnow, dt)

            # A2/A3: population-level L_bol/L_UV (FSPS SSP tables,
            # rescaled to this run's actual total_mass_msun -- see
            # MSLuminosityTable/UVLuminosityTable's own docstrings) and
            # the massive-star population's own ionizing-photon output
            # (independent of L_X/HMXB activity, unlike nphot_tot's
            # existing NPHOT_PER_LUMX-based accretion proxy -- see
            # _qh_ms_tot's docstring). lbol_tot/qh_tot follow the same
            # "population contribution + HMXB contribution" convention
            # already used for L_bol in scripts/run_paper1_experiment.py
            # (ms_lbol + lumx_tot); luv_tot is MS-only, no HMXB/accretion
            # UV model exists (same scope note as that script's own).
            # lumx_tot/nphot_tot themselves are unchanged -- this is
            # purely additive to the results dict.
            ms_lbol = self.ms_table.get_lbol(tnow, self.population.total_mass_msun)
            ms_luv = self.uv_table.get_luv(tnow, self.population.total_mass_msun)
            qh_ms_tot = self._qh_ms_tot(tnow)

            results.append(
                {
                    "time": tnow,
                    "lumx_tot": lumx_tot,
                    "nphot_tot": nphot_tot,
                    "nactive": nactive,
                    "ndead": ndead,
                    "lbol_tot": ms_lbol + lumx_tot,
                    "luv_tot": ms_luv,
                    "qh_tot": qh_ms_tot + nphot_tot,
                }
            )

        self._write_results(results, output_path)
        logger.info("Simulation complete.")
        return results

    def _qh_ms_tot(self, tnow: float) -> float:
        """Total ionizing-photon rate (photons/s) from currently-alive
        M >= 8 Msun massive stars in the population (A3, docs/science/
        paper1-detailed-work-breakdown.md) -- independent of L_X/HMXB
        activity, unlike the pre-existing `nphot_tot`
        (`BinaryPopulation.NPHOT_PER_LUMX * lumx_tot`, a fixed
        constant multiple of the X-ray luminosity, carrying zero
        independent information -- see docs/provenance.md Section 3).

        Source: `io/tables.py::IonizingPhotonTable.get_ngamma(mass)`,
        previously unused (see that class's docstring for why it's a
        "genuinely different quantity" from `MSLuminosityTable`).
        `get_ngamma()` returns `log10(N_gamma)`, the *total* number of
        ionizing photons emitted over the star's *whole* main-sequence
        lifetime (confirmed by its own MUNIT/MATOM conversion -- it
        multiplies the tabulated per-atom yield by the star's total
        baryon count -- and by a direct sanity check before adopting
        this interpretation: dividing by `LifetimeTable.get_lifetime()`
        converted to seconds gives Q_H(M) rates of ~2e47/s at 10 Msun
        up to ~8e49/s at 80 Msun, matching the well-known literature
        range for O/early-B main-sequence ionizing rates -- e.g. Vacca,
        Garmany & Shull 1996, ApJ 460, 914 -- to within the expected
        order of magnitude, not an exact-match requirement).

        A star is "currently alive" here if: it is `m1` and
        `nturn == 0` (hasn't had its own supernova yet), or it is `m2`
        (a real companion, `m2 >= 8`) and `tnow < t2_lifetime` (works
        correctly across RLOF-driven lifetime-clock resets, since
        `t2_lifetime` is always the secondary's current predicted
        death time, not a fixed formation-time value). Below 8 Msun a
        star contributes no ionizing photons at all (`get_ngamma`'s
        own physical cutoff, independent of `config.mcut`, which is a
        different, binary-formation-only threshold).
        """
        pop = self.population
        assert pop is not None
        total = 0.0

        alive_m1 = pop.m1[(pop.nturn == 0) & (pop.m1 >= 8.0)]
        alive_m2 = pop.m2[(pop.m2 >= 8.0) & (tnow < pop.t2_lifetime)]

        for m in np.concatenate([alive_m1, alive_m2]):
            lifetime_s = pop.lifetime_table.get_lifetime(m) * self.MYR_TO_SECONDS
            if lifetime_s > 0.0:
                total += 10.0 ** self.ionizing_table.get_ngamma(m) / lifetime_s

        return total

    def _get_imf_name(self) -> str:
        return self.IMF_MAP.get(
            self.config.imf_type, f"Custom_IMF_{self.config.imf_type}"
        )

    def _write_initial_conditions(self, output_dir: Path):
        if self.population is None:
            raise RuntimeError(
                "Population is not initialized. Call initialize() first."
            )

        pop = self.population
        m1 = pop.m1
        m2 = pop.m2
        period = pop.period
        a = pop.a
        turnoff_time = pop.turnoff_time
        t2_lifetime = pop.t2_lifetime
        remnant_table = pop.remnant_table

        imf_name = self._get_imf_name()
        filename = output_dir / f"{imf_name}.init.dat"

        with open(filename, "w") as f:
            f.write(f"# {imf_name} IMF\n")
            f.write(
                "# ntot (mmin,mmax,mcut)/Msol (pmin,pmax)/days (lxmin,lxmax)/ergs/s\n"
            )
            f.write(
                f"{self.config.ntot} {self.config.mmin} {self.config.mmax} "
                f"{self.config.mcut} {self.config.pmin} {self.config.pmax} "
                f"{self.config.lxmin} {self.config.lxmax}\n"
            )
            f.write("# n (m1,m2)/M* P/days a/AU (t1,t2)/Myrs (mr1,mr2)/M*\n")

            for i in range(len(m1)):
                t1 = turnoff_time[i]
                t2 = t2_lifetime[i]
                mr1 = remnant_table.get_remnant_mass(m1[i])
                mr2 = remnant_table.get_remnant_mass(m2[i]) if m2[i] > 0 else 0.0

                f.write(
                    f"{i + 1:9d} {m1[i]:12.4f} {m2[i]:12.4f} "
                    f"{period[i]:12.4f} {a[i]:12.4f} {t1:12.4f} {t2:12.4f} "
                    f"{mr1:12.4f} {mr2:12.4f}\n"
                )

        logger.info(f"Initial conditions written to {filename}")

    def _write_results(self, results: list[dict], output_dir: Path):
        imf_name = self._get_imf_name()
        filename = output_dir / f"{imf_name}.tevol.dat"

        with open(filename, "w") as f:
            f.write(f"# {imf_name} IMF\n")
            f.write("# ntot (mmin,mmax,mcut)/Msol (pmin,pmax)/days\n")
            f.write(
                f"{self.config.ntot} {self.config.mmin} {self.config.mmax} "
                f"{self.config.mcut} {self.config.pmin} {self.config.pmax}\n"
            )
            # lbol_tot/luv_tot/qh_tot appended, not inserted, so any
            # existing simple column-index parsing of the first five
            # columns is unaffected (A2/A3 additions) -- nothing else
            # in this repo currently reads this file back (confirmed
            # via grep before adding these).
            f.write(
                "# t/Myrs lx_tot/ergs nphot npop ndead "
                "lbol_tot/ergs luv_tot/ergs qh_tot/s^-1\n"
            )

            f.writelines(
                f"{r['time']:18.8e} {r['lumx_tot']:18.8e} "
                f"{r['nphot_tot']:18.8e} {r['nactive']:9d} {r['ndead']:9d} "
                f"{r['lbol_tot']:18.8e} {r['luv_tot']:18.8e} {r['qh_tot']:18.8e}\n"
                for r in results
            )

        logger.info(f"Results written to {filename}")
