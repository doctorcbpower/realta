import logging
from pathlib import Path
from typing import Dict, List, Optional

from realta.binaries.population import BinaryPopulation
from realta.config import SimulationConfig, load_config

logger = logging.getLogger("realta")


class ClusterSimulation:
    """Main simulation orchestration engine."""

    def __init__(self, config: Optional[SimulationConfig] = None):
        if config is None:
            self.config = load_config()
        else:
            self.config = config

        self.population: Optional[BinaryPopulation] = None

    def initialize(self):
        logger.info("Initializing simulation...")
        self.population = BinaryPopulation(self.config)
        logger.info("Simulation initialized.")

    def run(self, output_dir: str = "output") -> list[dict]:
        if self.population is None:
            self.initialize()

        Path(output_dir).mkdir(parents=True, exist_ok=True)
        self._write_initial_conditions(output_dir)

        tmax = self.config.tmax
        dt = self.config.dt
        tnow = 0.0
        results = []

        logger.info(f"Starting time evolution to {tmax} Myr with dt={dt} Myr")

        while tnow <= tmax:
            lumx_tot, nphot_tot, nactive, ndead = self.population.evolve(tnow, dt)

            results.append(
                {
                    "time": tnow,
                    "lumx_tot": lumx_tot,
                    "nphot_tot": nphot_tot,
                    "nactive": nactive,
                    "ndead": ndead,
                }
            )

            tnow += dt

        self._write_results(results, output_dir)
        logger.info("Simulation complete.")
        return results

    def _write_initial_conditions(self, output_dir: str):
        imf_name = {1: "Salpeter", 2: "Kroupa", 3: "Chabrier"}
        filename = Path(output_dir) / f"{imf_name[self.config.imf_type]}.init.dat"

        with open(filename, "w") as f:
            f.write(f"# {imf_name[self.config.imf_type]} IMF\n")
            f.write(
                "# ntot (mmin,mmax,mcut)/Msol (pmin,pmax)/days (lxmin,lxmax)/ergs/s\n"
            )
            f.write(
                f"{self.config.ntot} {self.config.mmin} {self.config.mmax} "
                f"{self.config.mcut} {self.config.pmin} {self.config.pmax} "
                f"{self.config.lxmin} {self.config.lxmax}\n"
            )
            f.write("# n (m1,m2)/M* P/days a/AU (t1,t2)/Myrs (mr1,mr2)/M*\n")

            for i, binary in enumerate(self.population.binaries):
                t1 = self.population.lifetime_table.get_lifetime(binary.primary_mass)
                t2 = self.population.lifetime_table.get_lifetime(binary.secondary_mass)
                mr1 = self.population.remnant_table.get_remnant_mass(
                    binary.primary_mass
                )
                mr2 = self.population.remnant_table.get_remnant_mass(
                    binary.secondary_mass
                )

                f.write(
                    f"{i + 1} {binary.primary_mass:12.4f} {binary.secondary_mass:12.4f} "
                    f"{binary.period:12.4f} {binary.a:12.4f} {t1:12.4f} {t2:12.4f} "
                    f"{mr1:12.4f} {mr2:12.4f}\n"
                )

        logger.info(f"Initial conditions written to {filename}")

    def _write_results(self, results: List[Dict], output_dir: str):
        imf_name = {1: "Salpeter", 2: "Kroupa", 3: "Chabrier"}
        filename = Path(output_dir) / f"{imf_name[self.config.imf_type]}.tevol.dat"

        with open(filename, "w") as f:
            f.write(f"# {imf_name[self.config.imf_type]} IMF\n")
            f.write("# ntot (mmin,mmax,mcut)/Msol (pmin,pmax)/days\n")
            f.write(
                f"{self.config.ntot} {self.config.mmin} {self.config.mmax} "
                f"{self.config.mcut} {self.config.pmin} {self.config.pmax}\n"
            )
            f.write("# t/Myrs lx_tot/ergs nphot npop ndead\n")

            for r in results:
                f.write(
                    f"{r['time']:18.8e} {r['lumx_tot']:18.8e} "
                    f"{r['nphot_tot']:18.8e} {r['nactive']:9d} {r['ndead']:9d}\n"
                )

        logger.info(f"Results written to {filename}")
