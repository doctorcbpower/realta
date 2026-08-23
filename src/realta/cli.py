import argparse
import logging

from realta.config import load_config
from realta.simulation.cluster import ClusterSimulation

logger = logging.getLogger("realta")


def main():
    """Command-line interface for running realta."""
    parser = argparse.ArgumentParser(
        description="realta - High-Mass X-ray Binary Population Framework"
    )
    parser.add_argument(
        "--config",
        "-c",
        type=str,
        default=None,
        help="Path to configuration YAML file",
    )
    parser.add_argument(
        "--output", "-o", type=str, default="output", help="Output directory"
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Enable verbose logging"
    )

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    config = load_config(args.config)
    sim = ClusterSimulation(config)
    results = sim.run(args.output)

    print("\nSimulation complete!")
    print(f"Output directory: {args.output}")
    print(f"Total binaries: {len(sim.population.m1)}")
    print(f"Final time: {results[-1]['time']} Myr")
    print(f"Final X-ray luminosity: {results[-1]['lumx_tot']:.2e} erg/s")


if __name__ == "__main__":
    main()
