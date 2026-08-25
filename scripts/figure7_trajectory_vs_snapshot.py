#!/usr/bin/env python3
"""Figure 7 -- same single-age R_X snapshot, different R_X(t) trajectory.

Reuses the existing Figure 4 (alpha_IMF, f_bin) grid run output
(output/figure4/alpha*_fbin*/Salpeter.tevol.dat) -- no new simulations
are run. For each grid point we already have the full L_X(t)/L_UV(t)
time series, not just the four age-slices used by Figure 4 itself.

We pick a reference age t_ref, compute R_X(t_ref) for all 81 grid
points, and search for pairs whose R_X(t_ref) match closely (small
single-age separation) but whose full trajectories diverge
substantially elsewhere, quantified by a log-space trajectory
separation metric

    D_ij = sqrt( mean_k [ log10 R_X,i(t_k) - log10 R_X,j(t_k) ]^2 )

evaluated over t in [t_min, t_max] Myr (excluding non-positive R_X).
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

GRID_DIR = Path("output/figure4")
T_REF = 30.0  # Myr, reference "single-age snapshot"
SNAPSHOT_TOL = 0.05  # relative tolerance in R_X(t_ref) to call two points "matched"
T_MIN, T_MAX = 5.0, 100.0

pat = re.compile(r"alpha([0-9.]+)_fbin([0-9.]+)")


def load_point(d: Path):
    f = d / "Salpeter.tevol.dat"
    data = np.loadtxt(f, comments="#", skiprows=3)
    t = data[:, 0]
    lx = data[:, 1]
    luv = data[:, 6]
    return t, lx, luv


def nearest_idx(t, t0):
    return int(np.argmin(np.abs(t - t0)))


def main():
    points = []
    for d in sorted(GRID_DIR.glob("alpha*_fbin*")):
        m = pat.match(d.name)
        if not m:
            continue
        alpha, fbin = float(m.group(1)), float(m.group(2))
        t, lx, luv = load_point(d)
        with np.errstate(divide="ignore", invalid="ignore"):
            rx = np.where(luv > 0, lx / luv, np.nan)
        points.append({alpha: alpha, fbin: fbin, t: t, rx: rx, dir: d.name})

    print(f"Loaded {len(points)} grid points")

    for p in points:
        i = nearest_idx(p["t"], T_REF)
        p["rx_ref"] = p["rx"][i]
        p["t_ref_actual"] = p["t"][i]

    t_common = points[0]["t"]
    mask = (t_common >= T_MIN) & (t_common <= T_MAX)

    def log_traj(p):
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.log10(p["rx"][mask])

    for p in points:
        p["logrx_masked"] = log_traj(p)

    candidates = [p for p in points if np.isfinite(p["rx_ref"]) and p["rx_ref"] > 0.05]
    print(
        f"{len(candidates)} candidates with finite, non-negligible R_X({T_REF:.0f} Myr)"
    )

    pairs = []
    for i in range(len(candidates)):
        for j in range(i + 1, len(candidates)):
            a, b = candidates[i], candidates[j]
            rel_diff = abs(a["rx_ref"] - b["rx_ref"]) / max(a["rx_ref"], b["rx_ref"])
            if rel_diff > SNAPSHOT_TOL:
                continue
            param_dist = np.hypot(
                (a["alpha"] - b["alpha"]) / 1.4,
                (a["fbin"] - b["fbin"]) / 1.0,
            )
            if param_dist < 0.3:
                continue
            valid = np.isfinite(a["logrx_masked"]) & np.isfinite(b["logrx_masked"])
            if valid.sum() < 10:
                continue
            D = np.sqrt(
                np.mean((a["logrx_masked"][valid] - b["logrx_masked"][valid]) ** 2)
            )
            pairs.append((D, rel_diff, a, b))

    pairs.sort(key=lambda x: -x[0])
    print(f"{len(pairs)} matched-snapshot pairs found")
    for D, rel_diff, a, b in pairs[:10]:
        print(
            f"  D={D:.3f} dex  snapshot_rel_diff={rel_diff:.3f}  "
            f"A=(alpha={a['alpha']:.3f}, fbin={a['fbin']:.3f}, Rx_ref={a['rx_ref']:.3f})  "
            f"B=(alpha={b['alpha']:.3f}, fbin={b['fbin']:.3f}, Rx_ref={b['rx_ref']:.3f})"
        )

    used = set()
    chosen = []
    for D, rel_diff, a, b in pairs:
        key_a, key_b = a["dir"], b["dir"]
        if key_a in used or key_b in used:
            continue
        chosen.append((D, rel_diff, a, b))
        used.add(key_a)
        used.add(key_b)
        if len(chosen) == 3:
            break

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 12,
            "axes.labelsize": 14,
            "legend.fontsize": 9,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
        }
    )
    fig, axes = plt.subplots(
        1, len(chosen), figsize=(4.2 * len(chosen), 4.0), sharey=True
    )
    if len(chosen) == 1:
        axes = [axes]
    colors = ["#1f77b4", "#d62728"]
    for ax, (D, rel_diff, a, b) in zip(axes, chosen):
        for p, c, tag in ((a, colors[0], "A"), (b, colors[1], "B")):
            with np.errstate(divide="ignore", invalid="ignore"):
                rx = p["rx"]
            good = np.isfinite(rx) & (rx > 0)
            ax.plot(
                p["t"][good],
                rx[good],
                color=c,
                lw=1.8,
                label=rf"{tag}: $\alpha_{{\rm IMF}}={p['alpha']:.2f}$, $f_{{\rm bin}}={p['fbin']:.2f}$",
            )
        ax.axvline(T_REF, color="grey", ls=":", lw=1)
        ax.set_yscale("log")
        ax.set_xlabel("Time (Myr)")
        ax.set_title(
            f"$D={D:.2f}$ dex, snapshot diff $={rel_diff*100:.0f}\\%$", fontsize=11
        )
        ax.legend(loc="lower right", frameon=False)
        ax.grid(alpha=0.3, which="both", linestyle=":")
    axes[0].set_ylabel(r"$L_X/L_{\rm UV}$")
    fig.suptitle(
        rf"Figure 7 -- same $R_X(t={T_REF:.0f}\,{{\rm Myr}})$ snapshot, different trajectory",
        y=1.03,
    )
    fig.tight_layout()
    out_dir = Path("output/figure7")
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_path = out_dir / "figure7_trajectory_vs_snapshot.png"
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    print(f"Wrote {fig_path}")

    with open(out_dir / "figure7_pairs.txt", "w") as f:
        f.write(f"Reference age: {T_REF} Myr; snapshot tolerance: {SNAPSHOT_TOL}\n")
        f.write(
            f"Trajectory-distance window: [{T_MIN}, {T_MAX}] Myr, log10 R_X RMS\n\n"
        )
        for D, rel_diff, a, b in pairs[:15]:
            f.write(
                f"D={D:.3f} dex  snap_diff={rel_diff*100:.1f}%  "
                f"A=(alpha={a['alpha']:.3f}, fbin={a['fbin']:.3f}, Rx({T_REF:.0f})={a['rx_ref']:.3f})  "
                f"B=(alpha={b['alpha']:.3f}, fbin={b['fbin']:.3f}, Rx({T_REF:.0f})={b['rx_ref']:.3f})\n"
            )
    print(f"Wrote {out_dir / 'figure7_pairs.txt'}")


if __name__ == "__main__":
    main()
