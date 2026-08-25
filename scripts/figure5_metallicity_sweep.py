#!/usr/bin/env python3
"""Build Figure 5 -- the metallicity sweep (docs/science/
research-programme.md, "Figure 5 -- Metallicity": repeat the principal
experiment for several metallicities, showing whether the binary
signature survives changes in Z; Figure 5, docs/science/
paper1-detailed-work-breakdown.md).

Usage:
    python scripts/figure5_metallicity_sweep.py \
        --config configs/figure5_metallicity_sweep.yml \
        --output-dir output/figure5

One config in, one figure out -- same loader pattern as
scripts/figure3_merger_effects.py and scripts/figure4_imf_binary_grid.py.
See configs/figure5_metallicity_sweep.yml for the config format and the
prescription/metallicity choices' own rationale.

Observable plotted: L_X/L_UV(t) (Figure 2's own central quantity, A2/A3's
per-timestep wiring), one panel per `imetal_values` entry, with both
`binary_prescriptions` overlaid per panel -- so the separation between
the two curves within a panel is "the binary signature", and comparing
that separation across panels is "whether it survives changes in Z".
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
    "single": "Single-star",
    "non_interacting": "Non-interacting binaries",
    "standard_interaction": "Standard interaction",
    "enhanced_interaction": "Enhanced interaction",
    "enhanced_mergers": "Enhanced mergers",
}

IMETAL_LABELS = {1: "Z = 0", 2: "Z = 0.008", 3: "Z = 0.02"}


def load_sweep_config(
    path: str | Path,
) -> tuple[dict, list[str], list[int]]:
    with open(path, "r") as f:
        experiment = yaml.safe_load(f) or {}

    required = {"base", "binary_prescriptions", "imetal_values"}
    missing = required - set(experiment)
    if missing:
        raise ValueError(f"{path} is missing required key(s): {sorted(missing)}")

    base = experiment["base"]
    valid_fields = {field.name for field in dataclasses.fields(SimulationConfig)}
    unknown = sorted(set(base) - valid_fields - {"binary_prescription", "imetal"})
    if unknown:
        raise ValueError(
            f"Unknown key(s) in {path}'s 'base' block: {unknown}. "
            f"Valid keys are: {sorted(valid_fields)}."
        )
    for forced in ("binary_prescription", "imetal"):
        if forced in base:
            raise ValueError(
                f"'base' must not set {forced!r} directly -- it is supplied "
                "per run from 'binary_prescriptions'/'imetal_values'."
            )

    return base, experiment["binary_prescriptions"], experiment["imetal_values"]


def run_variant(base: dict, prescription: str, imetal: int, output_dir: Path) -> dict:
    config = SimulationConfig(binary_prescription=prescription, imetal=imetal, **base)
    sim = ClusterSimulation(config)
    subdir = output_dir / f"imetal{imetal}" / prescription
    results = sim.run(output_dir=str(subdir))

    if not sim.uv_table.loaded:
        logger.warning(
            f"UVLuminosityTable did not load (imetal={imetal}) -- L_UV(t) "
            "will be 0.0, and this panel's L_X/L_UV will be degenerate."
        )

    return {
        "prescription": prescription,
        "label": PRESCRIPTION_LABELS.get(prescription, prescription),
        "imetal": imetal,
        "time": np.array([r["time"] for r in results]),
        "l_x": np.array([r["lumx_tot"] for r in results]),
        "l_uv": np.array([r["luv_tot"] for r in results]),
    }


def build_sweep(
    base: dict,
    binary_prescriptions: list[str],
    imetal_values: list[int],
    output_dir: Path,
) -> dict[int, list[dict]]:
    """Returns {imetal: [run_variant(...) for each prescription]}."""
    sweep = {}
    for imetal in imetal_values:
        sweep[imetal] = [
            run_variant(base, prescription, imetal, output_dir)
            for prescription in binary_prescriptions
        ]
    return sweep


def make_figure5(sweep: dict[int, list[dict]], out_path: Path):
    imetal_values = list(sweep.keys())
    n_panels = len(imetal_values)
    fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, 4.5), sharey=True)
    if n_panels == 1:
        axes = [axes]

    for ax, imetal in zip(axes, imetal_values):
        for result in sweep[imetal]:
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = np.where(
                    result["l_uv"] > 0.0, result["l_x"] / result["l_uv"], np.nan
                )
            ax.plot(result["time"], ratio, label=result["label"], linewidth=1.4)
        ax.set_yscale("log")
        ax.set_xlabel("Time (Myr)")
        ax.set_title(IMETAL_LABELS.get(imetal, f"imetal={imetal}"))
        ax.grid(True, which="both", linestyle=":", alpha=0.4)

    axes[0].set_ylabel(r"$L_X / L_{\mathrm{UV}}$")
    axes[0].legend(loc="best", frameon=False, fontsize=8)
    fig.suptitle("Figure 5 -- Metallicity sweep")
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Figure 5 written to {out_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config", default="configs/figure5_metallicity_sweep.yml", type=Path
    )
    parser.add_argument("--output-dir", default="output/figure5", type=Path)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    plt.rcParams.update(PLOT_STYLE)

    base, binary_prescriptions, imetal_values = load_sweep_config(args.config)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    sweep = build_sweep(base, binary_prescriptions, imetal_values, args.output_dir)

    figures_dir = args.output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    make_figure5(sweep, figures_dir / "figure5_metallicity_sweep.png")


if __name__ == "__main__":
    main()
