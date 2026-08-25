#!/usr/bin/env python3
"""Build Figure 3 -- the effect of mergers (docs/science/
research-programme.md, "Figure 3 -- Effect of mergers"; C3, docs/
science/paper1-detailed-work-breakdown.md).

Usage:
    python scripts/figure3_merger_effects.py \
        --config configs/figure3_merger_effects.yml \
        --output-dir output/figure3

Compares no-mergers / standard-mergers / enhanced-mergers
(configs/figure3_merger_effects.yml's own prescription mapping --
"non_interacting" / "standard_interaction" / "enhanced_mergers"),
showing both:
    - luminosity evolution: L_X(t), read directly from
      ClusterSimulation.run()'s own `lumx_tot` (A2/A3's per-timestep
      wiring), matching Figure 1/2's own convention.
    - compact-object formation: cumulative merger count over time,
      built from BinaryPopulation.merge_time (set for every merger,
      whether from the formation-time p_merge channel or an
      RLOF-classifier-driven IMMEDIATE_MERGER/CE-merge event -- see
      docs/provenance.md Section 6), and cumulative `ndead` (both
      stars now compact remnants) from `results` directly, since
      Realta has no explicit WD/NS/BH type classifier to draw a finer
      compact-object census from (a named scope limitation, not
      invented physics -- see this module's own docstring further
      down).
"""

from __future__ import annotations

import argparse
import dataclasses
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import yaml

from realta.config import SimulationConfig
from realta.simulation.cluster import ClusterSimulation

logger = logging.getLogger("realta")

PLOT_STYLE = {
    "font.family": "serif",
    "font.size": 12,
    "axes.labelsize": 14,
    "legend.fontsize": 9,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
}

PRESCRIPTION_LABELS = {
    "non_interacting": "No mergers",
    "standard_interaction": "Standard mergers",
    "enhanced_mergers": "Enhanced mergers",
}


def load_experiment_config(path: str | Path) -> tuple[dict, list[str]]:
    with open(path, "r") as f:
        experiment = yaml.safe_load(f) or {}

    if "base" not in experiment or "binary_prescriptions" not in experiment:
        raise ValueError(
            f"{path} must have top-level 'base' (SimulationConfig fields) "
            "and 'binary_prescriptions' (list of prescription names) keys."
        )

    base = experiment["base"]
    prescriptions = experiment["binary_prescriptions"]

    valid_fields = {field.name for field in dataclasses.fields(SimulationConfig)}
    unknown = sorted(set(base) - valid_fields - {"binary_prescription"})
    if unknown:
        raise ValueError(
            f"Unknown key(s) in {path}'s 'base' block: {unknown}. "
            f"Valid keys are: {sorted(valid_fields)}."
        )
    if "binary_prescription" in base:
        raise ValueError(
            "'base' must not set binary_prescription directly -- it is "
            "supplied per-entry from 'binary_prescriptions'."
        )

    return base, prescriptions


def run_variant(base: dict, prescription: str, output_dir: Path) -> dict:
    config = SimulationConfig(binary_prescription=prescription, **base)
    sim = ClusterSimulation(config)
    results = sim.run(output_dir=str(output_dir / prescription))

    time = np.array([r["time"] for r in results])
    merge_time = sim.population.merge_time  # NaN if never merged
    # Cumulative merger count at each output timestep -- merge_time is
    # an absolute time (Myr), same convention as `time` above.
    cumulative_mergers = np.array(
        [np.sum(np.nan_to_num(merge_time, nan=np.inf) <= t) for t in time]
    )

    return {
        "prescription": prescription,
        "label": PRESCRIPTION_LABELS.get(prescription, prescription),
        "time": time,
        "l_x": np.array([r["lumx_tot"] for r in results]),
        "ndead": np.array([r["ndead"] for r in results]),
        "cumulative_mergers": cumulative_mergers,
        "n_binaries": len(sim.population.m1),
    }


def make_figure3(all_results: list[dict], out_path: Path):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    ax_lx, ax_co = axes
    for r in all_results:
        ax_lx.plot(r["time"], r["l_x"], label=r["label"])
        ax_co.plot(r["time"], r["ndead"], label=f"{r['label']} (compact both)")
        ax_co.plot(
            r["time"],
            r["cumulative_mergers"],
            linestyle="--",
            label=f"{r['label']} (mergers)",
        )

    ax_lx.set_xlabel("Time (Myr)")
    ax_lx.set_ylabel(r"$L_X\ (\mathrm{erg/s})$")
    ax_lx.set_yscale("log")
    ax_lx.set_title("Luminosity evolution")
    ax_lx.legend()

    ax_co.set_xlabel("Time (Myr)")
    ax_co.set_ylabel("Cumulative count")
    ax_co.set_title("Compact-object formation")
    ax_co.legend(fontsize=7)

    fig.suptitle("Figure 3 -- Effect of mergers")
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Figure 3 written to {out_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config", default="configs/figure3_merger_effects.yml", type=Path
    )
    parser.add_argument("--output-dir", default="output/figure3", type=Path)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    plt.rcParams.update(PLOT_STYLE)

    base, prescriptions = load_experiment_config(args.config)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    all_results = [
        run_variant(base, prescription, args.output_dir)
        for prescription in prescriptions
    ]

    figures_dir = args.output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    make_figure3(all_results, figures_dir / "figure3_merger_effects.png")


if __name__ == "__main__":
    main()
