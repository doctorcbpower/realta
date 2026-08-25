#!/usr/bin/env python3
"""Re-plot Figure 2 from already-computed Paper-1 output (no rerun),
grouping variants into a 'baseline' family (fsur-based, effectively
degenerate with non-interacting) and a 'delayed-interaction' family
(post-SN RLOF, wind-capture), and marking t_first,X.

Single-star is NOT plotted: L_X(t) = 0 identically for that variant
(confirmed directly from output/paper1/single/Kroupa.tevol.dat), so
log(L_X/L_UV) is undefined at every timestep and there is nothing to
draw -- it is noted in the caption instead of a phantom legend entry.
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

BASE = Path("output/paper1")
BASELINE = [
    "non_interacting",
    "standard_interaction",
    "enhanced_interaction",
    "enhanced_mergers",
]
DELAYED = ["post_sn_rlof", "wind_capture"]
LABELS = {
    "non_interacting": "Non-interacting",
    "standard_interaction": "Standard interaction",
    "enhanced_interaction": "Enhanced interaction",
    "enhanced_mergers": "Enhanced mergers",
    "post_sn_rlof": "Post-SN secondary RLOF",
    "wind_capture": "Wind-capture accretion",
}


def load(variant):
    f = BASE / variant / "Kroupa.tevol.dat"
    data = np.loadtxt(f, comments="#", skiprows=3)
    t, lx, luv = data[:, 0], data[:, 1], data[:, 6]
    with np.errstate(divide="ignore", invalid="ignore"):
        rx = np.where(luv > 0, lx / luv, np.nan)
    return t, rx, lx


# t_first,X from the non_interacting run (first nonzero L_X)
t0, rx0, lx0 = load("non_interacting")
nz = np.nonzero(lx0 > 0)[0]
t_first_x = t0[nz[0]] if len(nz) else None
print("t_first,X =", t_first_x, "Myr")

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
fig, ax = plt.subplots(figsize=(7.2, 6))

blues = plt.cm.Blues(np.linspace(0.45, 0.85, len(BASELINE)))
warm = ["#8c3b12", "#e8850c"]

for c, variant in zip(blues, BASELINE):
    t, rx, _ = load(variant)
    good = np.isfinite(rx) & (rx > 0)
    ax.plot(t[good], rx[good], color=c, lw=1.3, ls="-", label=LABELS[variant])

for c, variant in zip(warm, DELAYED):
    t, rx, _ = load(variant)
    good = np.isfinite(rx) & (rx > 0)
    ax.plot(t[good], rx[good], color=c, lw=2.2, ls="--", label=LABELS[variant])

if t_first_x is not None:
    ax.axvline(t_first_x, color="grey", lw=1, ls=":")
    ax.text(t_first_x + 1, 1.2e-2, r"$t_{\rm first,X}$", fontsize=10, color="grey")

ax.axvspan(20, 50, color="grey", alpha=0.08, lw=0)
ax.text(21, 6.0, "discriminating\nage range", fontsize=9, color="dimgrey")

ax.set_yscale("log")
ax.set_xlabel("Time (Myr)")
ax.set_ylabel(r"$L_X/L_{\rm UV}$")
ax.set_title("Figure 2 -- X-ray/UV evolution (central figure)")
ax.legend(loc="lower right", frameon=False, fontsize=8.5)
ax.grid(alpha=0.25, which="both", ls=":")
fig.tight_layout()

out = BASE / "figures" / "figure2_xray_uv_evolution_grouped.png"
fig.savefig(out, dpi=200)
print("wrote", out)
