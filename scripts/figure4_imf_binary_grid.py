#!/usr/bin/env python3
"""Build Figure 4 -- the (alpha_IMF, f_bin) degeneracy grid
(docs/science/research-programme.md, "Figure 4 -- IMF versus binary
degeneracy"; C2, docs/science/paper1-detailed-work-breakdown.md).

Usage:
    python scripts/figure4_imf_binary_grid.py \
        --config configs/figure4_imf_binary_grid.yml \
        --output-dir output/figure4

One config in, one figure out. Unlike scripts/run_paper1_experiment.py
(which sweeps binary_prescription), this sweeps two continuous
parameters directly: `imf_slope` (A4, imf_type forced to 1/Salpeter --
see imf/factory.py::get_imf's docstring for why the continuous slope
is Salpeter-only) and `binary_fraction` (A1). No interaction
prescription is varied here -- Figure 4 is specifically about the IMF/
binary-fraction degeneracy, not interaction physics -- so every grid
point uses plain fsur-based activation
(use_rlof_classifier/use_post_sn_rlof both left at their False
defaults).

Observable plotted: L_X/L_UV (matching Figure 2's own central
quantity), read directly from ClusterSimulation.run()'s own
`lumx_tot`/`luv_tot` output (A2/A3), at each of `selected_ages` --
the nearest available timestep to each requested age, not
interpolated (an explicit "nearest, not interpolated" convention like
the one used for Figure 2's slicing conventions elsewhere).
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


def load_grid_config(
    path: str | Path,
) -> tuple[dict, list[float], list[float], list[float]]:
    with open(path, "r") as f:
        experiment = yaml.safe_load(f) or {}

    required = {"base", "alpha_imf_values", "f_bin_values", "selected_ages"}
    missing = required - set(experiment)
    if missing:
        raise ValueError(f"{path} is missing required key(s): {sorted(missing)}")

    base = experiment["base"]
    valid_fields = {field.name for field in dataclasses.fields(SimulationConfig)}
    unknown = sorted(
        set(base) - valid_fields - {"imf_type", "imf_slope", "binary_fraction"}
    )
    if unknown:
        raise ValueError(
            f"Unknown key(s) in {path}'s 'base' block: {unknown}. "
            f"Valid keys are: {sorted(valid_fields)}."
        )
    for forced in ("imf_type", "imf_slope", "binary_fraction"):
        if forced in base:
            raise ValueError(
                f"'base' must not set {forced!r} directly -- it is supplied "
                "per grid point from alpha_imf_values/f_bin_values."
            )

    return (
        base,
        experiment["alpha_imf_values"],
        experiment["f_bin_values"],
        experiment["selected_ages"],
    )


def run_grid_point(
    base: dict, alpha_imf: float, f_bin: float, output_dir: Path
) -> list[dict]:
    config = SimulationConfig(
        imf_type=1, imf_slope=alpha_imf, binary_fraction=f_bin, **base
    )
    sim = ClusterSimulation(config)
    subdir = output_dir / f"alpha{alpha_imf:.3f}_fbin{f_bin:.3f}"
    return sim.run(output_dir=str(subdir))


def _nearest_step(results: list[dict], age: float) -> dict:
    return min(results, key=lambda r: abs(r["time"] - age))


def build_grid(
    base: dict,
    alpha_imf_values: list[float],
    f_bin_values: list[float],
    selected_ages: list[float],
    output_dir: Path,
) -> np.ndarray:
    """Returns an array of shape (len(selected_ages), len(f_bin_values),
    len(alpha_imf_values)) of L_X/L_UV, row-major in (age, f_bin, alpha_imf).
    """
    grid = np.full(
        (len(selected_ages), len(f_bin_values), len(alpha_imf_values)), np.nan
    )
    for i, alpha_imf in enumerate(alpha_imf_values):
        for j, f_bin in enumerate(f_bin_values):
            logger.info(f"Running alpha_imf={alpha_imf}, f_bin={f_bin}")
            results = run_grid_point(base, alpha_imf, f_bin, output_dir)
            for k, age in enumerate(selected_ages):
                r = _nearest_step(results, age)
                l_uv = r["luv_tot"]
                grid[k, j, i] = r["lumx_tot"] / l_uv if l_uv > 0.0 else np.nan
    return grid


def make_figure4(
    grid: np.ndarray,
    alpha_imf_values: list[float],
    f_bin_values: list[float],
    selected_ages: list[float],
    out_path: Path,
):
    n_panels = len(selected_ages)
    fig, axes = plt.subplots(
        1,
        n_panels + 1,
        figsize=(5 * n_panels + 1.0, 4.5),
        gridspec_kw={"width_ratios": [10] * n_panels + [0.5]},
    )
    panel_axes = axes[:-1] if n_panels > 1 else [axes[0]]
    cax = axes[-1]

    finite = grid[np.isfinite(grid)]
    vmin, vmax = (finite.min(), finite.max()) if finite.size else (0.0, 1.0)

    im = None
    for k, (ax, age) in enumerate(zip(panel_axes, selected_ages)):
        im = ax.imshow(
            grid[k],
            origin="lower",
            aspect="auto",
            extent=[
                min(alpha_imf_values),
                max(alpha_imf_values),
                min(f_bin_values),
                max(f_bin_values),
            ],
            vmin=vmin,
            vmax=vmax,
            cmap="viridis",
        )
        ax.set_xlabel(r"$\alpha_{\mathrm{IMF}}$")
        if k == 0:
            ax.set_ylabel(r"$f_{\mathrm{bin}}$")
        ax.set_title(f"t = {age:.0f} Myr")

    if im is not None:
        fig.colorbar(im, cax=cax, label=r"$L_X/L_{\mathrm{UV}}$")
    fig.suptitle(
        r"Figure 4 -- ($\alpha_{\mathrm{IMF}}$, $f_{\mathrm{bin}}$) degeneracy grid"
    )
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Figure 4 written to {out_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config", default="configs/figure4_imf_binary_grid.yml", type=Path
    )
    parser.add_argument("--output-dir", default="output/figure4", type=Path)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    plt.rcParams.update(PLOT_STYLE)

    base, alpha_imf_values, f_bin_values, selected_ages = load_grid_config(args.config)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    grid = build_grid(
        base, alpha_imf_values, f_bin_values, selected_ages, args.output_dir
    )

    figures_dir = args.output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    make_figure4(
        grid,
        alpha_imf_values,
        f_bin_values,
        selected_ages,
        figures_dir / "figure4_imf_binary_grid.png",
    )


if __name__ == "__main__":
    main()
