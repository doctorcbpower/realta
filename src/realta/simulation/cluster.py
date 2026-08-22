from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import ClassVar

from realta.binaries.population import BinaryPopulation
from realta.config import SimulationConfig, load_config

logger = logging.getLogger("realta")


class ClusterSimulation:
    """Main simulation orchestration engine."""

    IMF_MAP: ClassVar[dict[int, str]] = {
        1: "Salpeter",
        2: "Kroupa",
        3: "Chabrier",
    }

    def __init__(self, config: SimulationConfig | None = None):
        self.config = config if config is not None else load_config()
        self.population: BinaryPopulation | None = None

    def initialize(self):
        logger.info("Initializing simulation...")
        self.population = BinaryPopulation(self.config)
        logger.info("Simulation initialized.")

    def run(self, output_dir: str = "output") -> list[dict]:
        if self.population is None:
            self.initialize()

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

            results.append(
                {
                    "time": tnow,
                    "lumx_tot": lumx_tot,
                    "nphot_tot": nphot_tot,
                    "nactive": nactive,
                    "ndead": ndead,
                }
            )

        self._write_results(results, output_path)
        logger.info("Simulation complete.")
        return results

    def _get_imf_name(self) -> str:
        return self.IMF_MAP.get(
            self.config.imf_type, f"Custom_IMF_{self.config.imf_type}"
        )

    def _write_initial_conditions(self, output_dir: Path):
        imf_name = self._get_imf_name()
        filename = output_dir / f"{imf_name}.init.dat"
        pop = self.population

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

            # Iterate over vectorized array indices
            for i in range(len(pop.m1)):
                t1 = pop.turnoff_time[i]
                t2 = pop.t2_lifetime[i]
                mr1 = pop.remnant_table.get_remnant_mass(pop.m1[i])
                mr2 = (
                    pop.remnant_table.get_remnant_mass(pop.m2[i])
                    if pop.m2[i] > 0
                    else 0.0
                )

                f.write(
                    f"{i + 1:9d} {pop.m1[i]:12.4f} {pop.m2[i]:12.4f} "
                    f"{pop.period[i]:12.4f} {pop.a[i]:12.4f} {t1:12.4f} {t2:12.4f} "
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
            f.write("# t/Myrs lx_tot/ergs nphot npop ndead\n")

            f.writelines(
                f"{r['time']:18.8e} {r['lumx_tot']:18.8e} "
                f"{r['nphot_tot']:18.8e} {r['nactive']:9d} {r['ndead']:9d}\n"
                for r in results
            )

        logger.info(f"Results written to {filename}")
