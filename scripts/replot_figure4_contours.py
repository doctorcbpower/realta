#!/usr/bin/env python3
"""Re-plot Figure 4 from the already-computed (alpha_IMF, f_bin) grid
output (output/figure4/alpha*_fbin*/Salpeter.tevol.dat) -- no rerun.

Adds to the original figure4_imf_binary_grid.py rendering:
  - contours of constant R_X (1, 3, 10) on each panel;
  - a marker at the fiducial model (alpha_IMF=2.35, f_bin=1.0);
  - the same shared colour scale across panels as the original script
    already used (one global vmin/vmax, one colorbar).
"""

from pathlib import Path
import re
import yaml
import numpy as np
import matplotlib.pyplot as plt

GRID_DIR = Path("output/figure4")
CONFIG = yaml.safe_load(open("configs/figure4_imf_binary_grid.yml"))
ALPHA_VALUES = CONFIG["alpha_imf_values"]
FBIN_VALUES = CONFIG["f_bin_values"]
SELECTED_AGES = CONFIG["selected_ages"]

FIDUCIAL_ALPHA = 2.35  # Salpeter canonical slope
FIDUCIAL_FBIN = 1.0  # config.py SimulationConfig default

pat = re.compile(r"alpha([0-9.]+)_fbin([0-9.]+)")


def load_point(d):
    data = np.loadtxt(d / "Salpeter.tevol.dat", comments="#", skiprows=3)
    return data[:, 0], data[:, 1], data[:, 6]  # t, lx, luv


def nearest(t, t0):
    return int(np.argmin(np.abs(t - t0)))


def main():
    points = {}
    for d in sorted(GRID_DIR.glob("alpha*_fbin*")):
        m = pat.match(d.name)
        if not m:
            continue
        alpha, fbin = round(float(m.group(1)), 3), round(float(m.group(2)), 3)
        t, lx, luv = load_point(d)
        points[(alpha, fbin)] = (t, lx, luv)

    print(f"Loaded {len(points)} grid points")

    n_age = len(SELECTED_AGES)
    n_f, n_a = len(FBIN_VALUES), len(ALPHA_VALUES)
    grid = np.full((n_age, n_f, n_a), np.nan)

    for i, alpha in enumerate(ALPHA_VALUES):
        for j, fbin in enumerate(FBIN_VALUES):
            key = (round(alpha, 3), round(fbin, 3))
            if key not in points:
                print("MISSING", key)
                continue
            t, lx, luv = points[key]
            for k, age in enumerate(SELECTED_AGES):
                idx = nearest(t, age)
                grid[k, j, i] = lx[idx] / luv[idx] if luv[idx] > 0 else np.nan

    finite = grid[np.isfinite(grid)]
    vmin, vmax = finite.min(), finite.max()
    print(f"Shared colour scale: vmin={vmin:.3f} vmax={vmax:.3f}")

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 12,
            "axes.labelsize": 14,
            "legend.fontsize": 9,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
        }
    )
    fig, axes = plt.subplots(
        1,
        n_age + 1,
        figsize=(5 * n_age + 1.0, 4.6),
        gridspec_kw={"width_ratios": [10] * n_age + [0.5]},
    )
    panel_axes = axes[:-1]
    cax = axes[-1]

    alpha_mesh, fbin_mesh = np.meshgrid(ALPHA_VALUES, FBIN_VALUES)
    contour_levels = [1, 3, 10]
    im = None
    for k, (ax, age) in enumerate(zip(panel_axes, SELECTED_AGES)):
        im = ax.imshow(
            grid[k],
            origin="lower",
            aspect="auto",
            extent=[
                min(ALPHA_VALUES),
                max(ALPHA_VALUES),
                min(FBIN_VALUES),
                max(FBIN_VALUES),
            ],
            vmin=vmin,
            vmax=vmax,
            cmap="viridis",
        )
        panel_finite = grid[k][np.isfinite(grid[k])]
        levels_here = [
            lv
            for lv in contour_levels
            if panel_finite.size and panel_finite.min() < lv < panel_finite.max()
        ]
        if len(levels_here) >= 1:
            zz = np.where(np.isfinite(grid[k]), grid[k], vmin)
            cs = ax.contour(
                alpha_mesh,
                fbin_mesh,
                zz,
                levels=levels_here,
                colors="white",
                linewidths=1.1,
                linestyles=["-", "--", ":"][: len(levels_here)],
            )
            ax.clabel(cs, inline=True, fontsize=8, fmt=lambda v: f"{v:g}")
        ax.plot(
            FIDUCIAL_ALPHA,
            FIDUCIAL_FBIN,
            marker="*",
            markersize=16,
            markerfacecolor="white",
            markeredgecolor="black",
            markeredgewidth=0.8,
            linestyle="none",
            zorder=5,
        )
        ax.set_xlabel(r"$\alpha_{\mathrm{IMF}}$")
        if k == 0:
            ax.set_ylabel(r"$f_{\mathrm{bin}}$")
        ax.set_title(f"t = {age:.0f} Myr")

    fig.colorbar(im, cax=cax, label=r"$L_X/L_{\mathrm{UV}}$")
    fig.suptitle(
        r"Figure 4 -- ($\alpha_{\mathrm{IMF}}$, $f_{\mathrm{bin}}$) degeneracy grid "
        r"(white $\star$: fiducial model; contours: $R_X=1,3,10$)",
        fontsize=13,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])

    out_dir = Path("output/figure4/figures")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "figure4_imf_binary_grid_contours.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    print("wrote", out_path)


if __name__ == "__main__":
    main()
