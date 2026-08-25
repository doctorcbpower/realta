#!/usr/bin/env python3
"""Build Figure 6 -- a first look at stochastic realisations (docs/
science/research-programme.md, "Figure 6 -- Stochastic realisations":
"For finite cluster masses, show distributions rather than means...
P(L_X/L_UV | M_cl, Z, t)... provides the bridge to Paper 2.").

Usage:
    python scripts/figure6_stochastic_realisations.py \
        --config configs/figure6_stochastic_realisations.yml \
        --output-dir output/figure6

Deliberately lightweight, per an explicit user request (chat,
2026-08-25): this is NOT Paper 2's parameter-sweep/experiment-runner
machinery (`docs/science/development-roadmap.md` item 4's `Event`/
`PopulationHistory` abstraction), which does not exist in Realta yet.
It just re-runs the SAME config `n_realizations` times per cluster-mass
proxy (`ntot_values`), varying only `iseed`, and reads L_X/L_UV at each
of `selected_ages` off each run -- the minimum needed to actually see
the shape of the distribution the research programme names, not a
full implementation of Paper 2's own scope.

`ntot` is a proxy for M_cl, not M_cl itself -- each realization's own
BinaryPopulation.total_mass_msun (already computed,
binaries/population.py) is read directly and reported as the group's
mean REALIZED cluster mass on the figure's axis labels, rather than
just labelling groups by their `ntot` input.

A realization with zero active HMXBs at a given age (L_X=0, a real,
common outcome for a small enough cluster -- exactly the point of this
figure) cannot contribute a finite log10(L_X/L_UV) point to the violin
-- its existence is reported instead as a "quiescent fraction"
annotation on each panel, not silently dropped.
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
    "axes.labelsize": 13,
    "legend.fontsize": 9,
    "xtick.labelsize": 10,
    "ytick.labelsize": 12,
}


def load_sweep_config(
    path: str | Path,
) -> tuple[dict, list[int], int, list[float], int]:
    with open(path, "r") as f:
        experiment = yaml.safe_load(f) or {}

    required = {"base", "ntot_values", "n_realizations", "selected_ages", "base_iseed"}
    missing = required - set(experiment)
    if missing:
        raise ValueError(f"{path} is missing required key(s): {sorted(missing)}")

    base = experiment["base"]
    valid_fields = {field.name for field in dataclasses.fields(SimulationConfig)}
    unknown = sorted(set(base) - valid_fields)
    if unknown:
        raise ValueError(
            f"Unknown key(s) in {path}'s 'base' block: {unknown}. "
            f"Valid keys are: {sorted(valid_fields)}."
        )
    for forced in ("ntot", "iseed"):
        if forced in base:
            raise ValueError(
                f"'base' must not set {forced!r} directly -- it is supplied "
                "per realization from 'ntot_values'/'base_iseed'."
            )

    return (
        base,
        experiment["ntot_values"],
        experiment["n_realizations"],
        experiment["selected_ages"],
        experiment["base_iseed"],
    )


def run_realization(base: dict, ntot: int, iseed: int, output_dir: Path) -> dict:
    config = SimulationConfig(ntot=ntot, iseed=iseed, **base)
    sim = ClusterSimulation(config)
    subdir = output_dir / f"ntot{ntot}" / f"seed{iseed}"
    results = sim.run(output_dir=str(subdir))
    return {
        "total_mass_msun": sim.population.total_mass_msun,
        "time": np.array([r["time"] for r in results]),
        "l_x": np.array([r["lumx_tot"] for r in results]),
        "l_uv": np.array([r["luv_tot"] for r in results]),
    }


def _nearest_index(time: np.ndarray, age: float) -> int:
    return int(np.argmin(np.abs(time - age)))


def build_distributions(
    base: dict,
    ntot_values: list[int],
    n_realizations: int,
    selected_ages: list[float],
    base_iseed: int,
    output_dir: Path,
) -> dict[int, dict]:
    """Returns {ntot: {"mean_mass": float, "ages": {age: {"log10_ratios":
    np.ndarray (finite-L_X realizations only), "quiescent_fraction":
    float, "n": int}}}}.
    """
    sweep: dict[int, dict] = {}
    for group_index, ntot in enumerate(ntot_values):
        masses = []
        ratios_by_age: dict[float, list[float]] = {age: [] for age in selected_ages}
        zero_count_by_age: dict[float, int] = {age: 0 for age in selected_ages}

        for k in range(n_realizations):
            iseed = base_iseed + group_index * 1000 + k
            realization = run_realization(base, ntot, iseed, output_dir)
            masses.append(realization["total_mass_msun"])
            for age in selected_ages:
                idx = _nearest_index(realization["time"], age)
                l_x = realization["l_x"][idx]
                l_uv = realization["l_uv"][idx]
                if l_x > 0.0 and l_uv > 0.0:
                    ratios_by_age[age].append(np.log10(l_x / l_uv))
                else:
                    zero_count_by_age[age] += 1

        sweep[ntot] = {
            "mean_mass": float(np.mean(masses)),
            "ages": {
                age: {
                    "log10_ratios": np.array(ratios_by_age[age]),
                    "quiescent_fraction": zero_count_by_age[age] / n_realizations,
                    "n": n_realizations,
                }
                for age in selected_ages
            },
        }
        logger.info(
            f"ntot={ntot}: mean realized cluster mass = "
            f"{sweep[ntot]['mean_mass']:.0f} Msun"
        )
    return sweep


def make_figure6(sweep: dict[int, dict], selected_ages: list[float], out_path: Path):
    ntot_values = list(sweep.keys())
    n_panels = len(selected_ages)
    fig, axes = plt.subplots(1, n_panels, figsize=(5.5 * n_panels, 5), sharey=True)
    if n_panels == 1:
        axes = [axes]

    positions = np.arange(1, len(ntot_values) + 1)
    labels = [
        f"$M_{{cl}}\\approx${sweep[ntot]['mean_mass']:.0f} $M_\\odot$"
        for ntot in ntot_values
    ]

    for ax, age in zip(axes, selected_ages):
        data = [sweep[ntot]["ages"][age]["log10_ratios"] for ntot in ntot_values]
        plot_positions = [
            pos for pos, arr in zip(positions, data, strict=True) if arr.size > 0
        ]
        plot_data = [arr for arr in data if arr.size > 0]
        if plot_data:
            parts = ax.violinplot(
                plot_data, positions=plot_positions, showmeans=True, showextrema=True
            )
            for body in parts["bodies"]:
                body.set_alpha(0.6)

        for pos, ntot in zip(positions, ntot_values, strict=True):
            q_frac = sweep[ntot]["ages"][age]["quiescent_fraction"]
            ax.annotate(
                f"quiescent: {q_frac:.0%}",
                xy=(pos, 1.0),
                xycoords=("data", "axes fraction"),
                ha="center",
                va="bottom",
                fontsize=8,
                rotation=0,
                annotation_clip=False,
            )

        ax.set_xticks(positions)
        ax.set_xticklabels(labels)
        ax.set_title(f"t = {age:.0f} Myr")
        ax.grid(True, axis="y", linestyle=":", alpha=0.4)

    axes[0].set_ylabel(r"$\log_{10}(L_X / L_{\mathrm{UV}})$")
    fig.suptitle(
        r"Figure 6 -- Stochastic realisations: $P(L_X/L_{\mathrm{UV}} \mid M_{cl}, Z, t)$"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Figure 6 written to {out_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config", default="configs/figure6_stochastic_realisations.yml", type=Path
    )
    parser.add_argument("--output-dir", default="output/figure6", type=Path)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    plt.rcParams.update(PLOT_STYLE)

    base, ntot_values, n_realizations, selected_ages, base_iseed = load_sweep_config(
        args.config
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)

    sweep = build_distributions(
        base, ntot_values, n_realizations, selected_ages, base_iseed, args.output_dir
    )

    figures_dir = args.output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    make_figure6(
        sweep, selected_ages, figures_dir / "figure6_stochastic_realisations.png"
    )


if __name__ == "__main__":
    main()
