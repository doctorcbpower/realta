#!/usr/bin/env python3
"""Generate a tabulated MS bolometric-luminosity-vs-time table via FSPS,
for use as Realta's `MSLuminosityTable` data (analogous to the existing
lifetimes_z*.dat / remnant_masses.dat tables in src/realta/data/).

Run this once, locally, wherever `fsps` + SPS_HOME are already set up --
it is NOT part of the Realta package and is not meant to run at import
time or in CI. Requires: pip install fsps (MIT), plus the FSPS Fortran/
data tree from https://github.com/cconroy20/fsps (also MIT) pointed to
by the SPS_HOME environment variable.

Produces one file per metallicity: ms_lbol_z<label>.dat, with columns
    age_myr    log10_lbol_total_erg_s
for an instantaneous-burst (SSP), Kroupa-IMF population of the given
total stellar mass, matching Power et al. (2009)'s cluster model
(M_cluster = 1e6 Msun).

IMPORTANT domain-of-validity note: FSPS's isochrone/track libraries do
not extend to true primordial (Z=0) metallicity -- there is no exact
match for Realta's imetal=1 (Z=0) case. This script uses FSPS's lowest
available metallicity as a documented proxy for Z=0 and prints exactly
what Z that turned out to be; review before trusting the imetal=1 output.
Z=0.008 (Realta's imetal=2) also isn't on FSPS's native grid but
zcontinuous=1 interpolates continuously in log(Z), so it is not a proxy
-- it's a real (interpolated) value.
"""

from __future__ import annotations

import numpy as np

import fsps

MSUN_LBOL_ERG_S = 3.828e33  # IAU nominal solar luminosity, erg/s
CLUSTER_MASS_MSUN = 1.0e6  # Power et al. (2009) fiducial cluster mass
AGE_MAX_MYR = 100.0

# Realta's three imetal options and their target Z. imetal=1 (Z=0) has no
# true match in FSPS -- see docstring above.
TARGETS = {
    1: ("z0", 0.0),
    2: ("z8e-3", 0.008),
    3: ("z2e-2", 0.02),
}


def make_table(imetal: int, label: str, z_target: float) -> None:
    sp = fsps.StellarPopulation(
        zcontinuous=1,  # continuously interpolate in metallicity
        sfh=0,  # SSP: instantaneous burst
        imf_type=2,  # Kroupa (2001) -- matches Realta's imf_type=2
        add_stellar_remnants=1,
    )

    if z_target <= 0.0:
        # No true Z=0 track library in FSPS. Use the lowest metallicity
        # on the underlying isochrone grid as a documented proxy rather
        # than silently substituting. `zlegend` reports that grid
        # regardless of zcontinuous mode, so no toggling is needed.
        z_actual = float(np.min(sp.zlegend))
        print(
            f"WARNING: imetal={imetal} (Z=0) has no FSPS equivalent; "
            f"using lowest available Z={z_actual:.2e} as a documented "
            f"proxy. Do not treat this as a true zero-metallicity result."
        )
    else:
        z_actual = z_target

    logzsol = np.log10(z_actual / sp.solar_metallicity)
    sp.params["logzsol"] = logzsol

    log_age = sp.log_age  # log10(age / yr), per-Msun-formed SSP grid
    log_lbol = sp.log_lbol  # log10(L_bol / Lsun), per-Msun-formed SSP grid

    age_myr = (10.0**log_age) / 1.0e6
    mask = (age_myr > 0.0) & (age_myr <= AGE_MAX_MYR)
    age_myr = age_myr[mask]
    log_lbol = log_lbol[mask]

    # Scale from "per 1 Msun formed" to the cluster's total stellar mass,
    # then to erg/s.
    log_lbol_total_erg_s = (
        log_lbol + np.log10(CLUSTER_MASS_MSUN) + np.log10(MSUN_LBOL_ERG_S)
    )

    order = np.argsort(age_myr)
    age_myr = age_myr[order]
    log_lbol_total_erg_s = log_lbol_total_erg_s[order]

    filename = f"ms_lbol_{label}.dat"
    with open(filename, "w") as f:
        f.write(
            f"# FSPS SSP bolometric luminosity vs time (Kroupa IMF, Z={z_actual:.4g})\n"
        )
        f.write(f"# Cluster mass: {CLUSTER_MASS_MSUN:.4g} Msun, instantaneous burst\n")
        f.write("# age_myr  log10_lbol_total_erg_s\n")
        for age, log_lbol in zip(age_myr, log_lbol_total_erg_s):
            f.write(f"{age:14.6f} {log_lbol:14.6f}\n")

    print(f"wrote {filename} ({len(age_myr)} rows, Z={z_actual:.4g})")


if __name__ == "__main__":
    for imetal, (label, z_target) in TARGETS.items():
        make_table(imetal, label, z_target)
