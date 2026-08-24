#!/usr/bin/env python3
"""Run Paper 1's basic experiment from a single YAML config and produce
Figures 1 and 2 (docs/science/research-programme.md, "Paper 1 --
X-ray fingerprints of massive-star multiplicity").

Usage:
    python scripts/run_paper1_experiment.py \
        --config configs/paper1_basic_experiment.yml \
        --output-dir output/paper1

One config in, both figures out -- no manual notebook wrangling. See
configs/paper1_basic_experiment.yml for the config format (a `base`
SimulationConfig block plus a `binary_prescriptions` sweep list; this
script's loader is deliberately specific to that shape, not a general
experiment-runner -- see the config file's own header comment).

Observables, exactly as currently implemented (see
docs/provenance.md Sections 4/6/7 for what each one is and isn't):
    L_bol(t) = MS bolometric luminosity (MSLuminosityTable) + L_X(t)
               -- matches notebooks/Power_etal_2009_Plots.ipynb's own
               "Total bolometric luminosity" convention.
    L_UV(t)  = MS far-UV luminosity only (UVLuminosityTable, GALEX FUV).
               No HMXB/accretion UV contribution is modelled -- Realta
               has no accretion spectral model that extends into the
               UV, only the X-ray band (xray/luminosity.py) and the
               X-ray-to-ionising-photon conversion
               (BinaryPopulation.NPHOT_PER_LUMX). Flagging this rather
               than inventing an unsourced accretion-UV correction.
    Q_H(t)   = nphot_tot, the HMXB ionising-photon rate only (Power et
               al. 2013 conversion, see provenance.md Section 3). No MS
               ionising-photon contribution is included: the existing
               notebook's MS ionising curve is an ad hoc Starburst99-
               shaped placeholder, not backed by any table
               (IonizingPhotonTable exists but estimates a different,
               currently-unused quantity -- see its class docstring in
               io/tables.py) -- not reused here for the same reason.
    L_X(t)   = lumx_tot, summed active-HMXB X-ray luminosity.

As of this writing UVLuminosityTable's data files
(src/realta/data/fuv_lbol_z*.dat) do not exist yet (see
docs/provenance.md Section 7) -- L_UV(t) and therefore Figure 2
(L_X/L_UV) will be degenerate (L_UV = 0 everywhere) until
scripts/generate_fuv_luminosities.py has been run and its output placed
in src/realta/data/. This script still produces both figures using
whatever UVLuminosityTable currently returns, and prints an explicit
warning if the table didn't load, rather than silently producing a
misleading Figure 2.
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
from realta.io.tables import MSLuminosityTable, UVLuminosityTable
from realta.simulation.cluster import ClusterSimulation

logger = logging.getLogger("realta")

PLOT_STYLE = {
    "font.family": "serif",
    "font.size": 12,
    "axes.labelsize": 14,
    "legend.fontsize": 9,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "figure.autolayout": True,
}

# Display order and labels for Paper 1's five basic-experiment variants
# (research-programme.md); see
# docs/science/paper1-binary-interaction-proposal.md for what each one
# means quantitatively.
PRESCRIPTION_LABELS = {
    "single": "Single-star",
    "non_interacting": "Non-interacting binaries",
    "standard_interaction": "Standard interaction",
    "enhanced_interaction": "Enhanced interaction",
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
    lumx_tot = np.array([r["lumx_tot"] for r in results])
    nphot_tot = np.array([r["nphot_tot"] for r in results])
    total_mass_msun = sim.population.total_mass_msun

    ms_table = MSLuminosityTable(imetal=config.imetal)
    uv_table = UVLuminosityTable(imetal=config.imetal)
    if not uv_table.loaded:
        logger.warning(
            f"UVLuminosityTable did not load (imetal={config.imetal}) -- "
            "L_UV(t) will be 0.0 for all t, and Figure 2 (L_X/L_UV) will "
            "be degenerate. Run scripts/generate_fuv_luminosities.py and "
            "place its output in src/realta/data/ to fix this."
        )

    ms_lbol = np.array([ms_table.get_lbol(t, total_mass_msun) for t in time])
    ms_luv = np.array([uv_table.get_luv(t, total_mass_msun) for t in time])

    return {
        "prescription": prescription,
        "label": PRESCRIPTION_LABELS.get(prescription, prescription),
        "time": time,
        "l_bol": ms_lbol + lumx_tot,
        "l_uv": ms_luv,
        "q_h": nphot_tot,
        "l_x": lumx_tot,
    }


def make_figure1(all_results: list[dict], out_path: Path):
    fig, axes = plt.subplots(2, 2, figsize=(11, 9), sharex=True)
    panels = [
        ("l_bol", r"$L_{\mathrm{bol}}\ (\mathrm{erg/s})$", axes[0, 0]),
        ("l_uv", r"$L_{\mathrm{UV}}\ (\mathrm{erg/s})$", axes[0, 1]),
        ("q_h", r"$Q_H\ (\mathrm{s}^{-1})$", axes[1, 0]),
        ("l_x", r"$L_X\ (\mathrm{erg/s})$", axes[1, 1]),
    ]

    for key, ylabel, ax in panels:
        for result in all_results:
            ax.plot(result["time"], result[key], label=result["label"], linewidth=1.4)
        ax.set_yscale("log")
        ax.set_xlim(0, all_results[0]["time"].max())
        ax.set_ylabel(ylabel)
        ax.grid(True, which="both", linestyle=":", alpha=0.4)

    axes[1, 0].set_xlabel("Time (Myr)")
    axes[1, 1].set_xlabel("Time (Myr)")
    axes[0, 0].legend(loc="upper right", frameon=False, fontsize=8)
    fig.suptitle("Figure 1 -- Population evolution (Paper 1 basic experiment)")
    fig.savefig(out_path)
    plt.close(fig)
    logger.info(f"Figure 1 written to {out_path}")


def make_figure2(all_results: list[dict], out_path: Path):
    fig, ax = plt.subplots(figsize=(7, 6))

    for result in all_results:
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(
                result["l_uv"] > 0.0, result["l_x"] / result["l_uv"], np.nan
            )
        ax.plot(result["time"], ratio, label=result["label"], linewidth=1.4)

    ax.set_yscale("log")
    ax.set_xlim(0, all_results[0]["time"].max())
    ax.set_xlabel("Time (Myr)")
    ax.set_ylabel(r"$L_X / L_{\mathrm{UV}}$")
    ax.legend(loc="best", frameon=False, fontsize=9)
    ax.grid(True, which="both", linestyle=":", alpha=0.4)
    ax.set_title("Figure 2 -- X-ray/UV evolution (central figure)")
    fig.savefig(out_path)
    plt.close(fig)
    logger.info(f"Figure 2 written to {out_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", default="configs/paper1_basic_experiment.yml", type=Path
    )
    parser.add_argument("--output-dir", default="output/paper1", type=Path)
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
    make_figure1(all_results, figures_dir / "figure1_population_evolution.png")
    make_figure2(all_results, figures_dir / "figure2_xray_uv_evolution.png")


if __name__ == "__main__":
    main()
