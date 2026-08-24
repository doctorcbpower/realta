#!/usr/bin/env python3
"""Generate a tabulated far-UV-luminosity-vs-time table via FSPS, for
Realta's Paper 1 `L_UV(t)` observable (Figs. 1-2, see
docs/science/research-programme.md). Analogous to
generate_ms_luminosities.py, and produces data in the same layout so it
can be loaded the same way (see io/tables.py::MSLuminosityTable).

Run this once, locally, wherever `fsps` + SPS_HOME are already set up
-- it is NOT part of the Realta package and is not meant to run at
import time or in CI. Requires: pip install fsps (MIT), plus the FSPS
Fortran/data tree from https://github.com/cconroy20/fsps (also MIT)
pointed to by the SPS_HOME environment variable.

Band choice: GALEX FUV (~1528 A rest-frame), not NUV -- reviewed and
decided in docs/science/paper1-binary-interaction-proposal.md's "UV
band decision" section. FUV is the standard massive-star/SFR UV tracer
(Kennicutt & Evans 2012) and is more sensitive to the O/B population
than NUV, which would also pick up a longer-lived A-star contribution
and dilute the timescale match to L_X/Q_H that Paper 1's Q2 depends on.

Produces one file per metallicity: fuv_lbol_z<label>.dat, with columns
    age_myr    log10_lfuv_total_erg_s
for an instantaneous-burst (SSP), Kroupa-IMF population of the given
total stellar mass, matching Power et al. (2009)'s cluster model
(M_cluster = 1e6 Msun) -- same fiducial-mass convention as
MSLuminosityTable, so a UVLuminosityTable consuming this data should
rescale linearly by actual total_mass_msun the same way
MSLuminosityTable.get_lbol() does.

IMPORTANT domain-of-validity note: FSPS's isochrone/track libraries do
not extend to true primordial (Z=0) metallicity -- there is no exact
match for Realta's imetal=1 (Z=0) case. This script uses FSPS's lowest
available metallicity as a documented proxy for Z=0 and prints exactly
what Z that turned out to be; review before trusting the imetal=1
output. Z=0.008 (Realta's imetal=2) also isn't on FSPS's native grid
but zcontinuous=1 interpolates continuously in log(Z), so it is not a
proxy -- it's a real (interpolated) value.

`L_FUV` here is `nu * L_nu` at the GALEX FUV pivot wavelength,
converted from FSPS's AB absolute magnitude via the standard
m_AB -> f_nu -> L_nu -> nu*L_nu chain -- see `fuv_luminosity_erg_s`
below for the exact conversion, so the value is directly comparable in
units (erg/s) to MSLuminosityTable's L_bol.
"""

from __future__ import annotations

import fsps
import numpy as np

CLUSTER_MASS_MSUN = 1.0e6  # Power et al. (2009) fiducial cluster mass
AGE_MAX_MYR = 100.0

# GALEX FUV effective/pivot wavelength (Morrissey et al. 2007), Angstrom.
FUV_PIVOT_WAVELENGTH_AA = 1528.0
SPEED_OF_LIGHT_CM_S = 2.99792458e10
AB_ZEROPOINT_ERG_S_CM2_HZ = 3.631e-20  # f_nu for m_AB = 0
PARSEC_CM = 3.0857e18
TEN_PC_CM = 10.0 * PARSEC_CM

# Realta's three imetal options and their target Z. imetal=1 (Z=0) has no
# true match in FSPS -- see docstring above.
TARGETS = {
    1: ("z0", 0.0),
    2: ("z8e-3", 0.008),
    3: ("z2e-2", 0.02),
}


def fuv_luminosity_erg_s(abs_mag_ab: np.ndarray) -> np.ndarray:
    """Convert FSPS AB absolute magnitude (per Msun formed) to nu*L_nu
    (erg/s, per Msun formed) at the GALEX FUV pivot frequency.

    m_AB -> f_nu at 10 pc -> L_nu = 4*pi*d^2*f_nu -> nu*L_nu.
    """
    f_nu_at_10pc = AB_ZEROPOINT_ERG_S_CM2_HZ * 10.0 ** (-0.4 * abs_mag_ab)
    l_nu = 4.0 * np.pi * TEN_PC_CM**2 * f_nu_at_10pc
    nu = SPEED_OF_LIGHT_CM_S / (FUV_PIVOT_WAVELENGTH_AA * 1.0e-8)
    return nu * l_nu


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
    abs_mag_fuv = sp.get_mags(bands=["galex_fuv"])[:, 0]  # AB, per Msun formed

    age_myr = (10.0**log_age) / 1.0e6
    mask = (age_myr > 0.0) & (age_myr <= AGE_MAX_MYR)
    age_myr = age_myr[mask]
    abs_mag_fuv = abs_mag_fuv[mask]

    lfuv_per_msun = fuv_luminosity_erg_s(abs_mag_fuv)
    log_lfuv_total_erg_s = np.log10(lfuv_per_msun) + np.log10(CLUSTER_MASS_MSUN)

    order = np.argsort(age_myr)
    age_myr = age_myr[order]
    log_lfuv_total_erg_s = log_lfuv_total_erg_s[order]

    filename = f"fuv_lbol_{label}.dat"
    with open(filename, "w") as f:
        f.write(
            f"# FSPS SSP GALEX-FUV (nu*L_nu, 1528 A) vs time "
            f"(Kroupa IMF, Z={z_actual:.4g})\n"
        )
        f.write(f"# Cluster mass: {CLUSTER_MASS_MSUN:.4g} Msun, instantaneous burst\n")
        f.write("# age_myr  log10_lfuv_total_erg_s\n")
        f.writelines(
            f"{age:14.6f} {log_lfuv:14.6f}\n"
            for age, log_lfuv in zip(age_myr, log_lfuv_total_erg_s)
        )

    print(f"wrote {filename} ({len(age_myr)} rows, Z={z_actual:.4g})")


if __name__ == "__main__":
    for imetal, (label, z_target) in TARGETS.items():
        make_table(imetal, label, z_target)
