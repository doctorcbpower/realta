"""Radiatively-driven (CAK) stellar wind physics for OB/supergiant
donors, feeding the wind-capture accretion module
(binaries/wind_capture.py).

Source: El Mellah & Casse (2017, MNRAS, arXiv:1609.01532), "A
numerical investigation of wind accretion in persistent Supergiant
X-ray Binaries I", building on Castor, Abbott & Klein (1975, CAK) and
Friend & Castor (1982, ApJ 261, 293). Both source papers pasted and
verified directly before implementing (not from memory), per this
project's established discipline for dense/precise physics.

Scope: this module gives the ISOLATED-star CAK wind (mass-loss rate
and velocity law) -- it does NOT include Friend & Castor's Roche-
potential modification of the wind's own acceleration (their eqs.
4/6-8), which would need a 1D ODE solve per binary. That was a
deliberate scope decision (chat, 2026-08-25): the isolated-star law,
evaluated at the orbital separation, is explicitly flagged as an
approximation at the call site (binaries/wind_capture.py) rather than
silently treated as the modified result -- see that module's own
docstring for how it's used and the accuracy consequence.

Unit-conversion discipline: every physical constant below is a
directly-cited cgs value, and every derived conversion factor
(L_Edd/c^2 in Msun/yr, etc.) is computed explicitly from them in code
-- not a memorized "standard" shortcut number. A first draft of this
module used a remembered "L_Edd(1 Msun)=3.2e4 Lsun" figure that turned
out to implicitly assume kappa_e=0.4 (pure hydrogen), inconsistent
with the kappa_e=0.34 (X=0.7) default actually used here, and a
remembered "L_Edd/c^2=1.5e-8 Msun/yr per Lsun" that was off by ~4
orders of magnitude (that number is closer to a *typical-efficiency
Eddington accretion rate*, a different quantity, not L_Edd/c^2
itself) -- both caught by direct numerical recomputation before this
version was written, not left in.
"""

from __future__ import annotations

# CODATA/IAU physical constants, cgs, cited directly rather than via
# a derived "standard" shortcut number.
G_CGS = 6.674e-8  # cm^3 g^-1 s^-2
C_CGS = 2.998e10  # cm/s
MSUN_G = 1.989e33  # g
RSUN_CM = 6.957e10  # cm, IAU 2015 nominal (same value used for
# RSUN_PER_AU in binaries/population.py)
LSUN_ERG_S = 3.828e33  # erg/s, IAU nominal
YEAR_S = 3.1557e7  # s (Julian year), matching
# simulation/cluster.py::ClusterSimulation.MYR_TO_SECONDS's own
# convention (that constant is this one * 1e6)

# Electron-scattering (Thomson) mass opacity for a fully-ionized
# solar-composition plasma, kappa_es = 0.2*(1+X) cm^2/g with X~0.7 the
# hydrogen mass fraction -- a standard stellar-structure result (e.g.
# Kippenhahn & Weigert, "Stellar Structure and Evolution"), not from
# either source paper specifically. Used only for the Eddington
# luminosity/factor helpers below.
ELECTRON_SCATTERING_OPACITY = 0.34  # cm^2/g

# G*Msun/Rsun in km^2/s^2 -- GM_sun/R_sun = 6.674e-8*1.989e33/6.957e10
# = 1.9081e15 cm^2/s^2 = 1.9081e5 km^2/s^2. Gives
# v_esc(1 Msun, 1 Rsun) = sqrt(2*this) = 617.8 km/s, matching the
# well-known solar escape velocity (~617.5-618 km/s) as a sanity
# check.
G_MSUN_RSUN_KM2_S2 = G_CGS * MSUN_G / RSUN_CM / 1.0e10  # cm^2/s^2 -> km^2/s^2


def eddington_luminosity(
    mass: float, kappa_e: float = ELECTRON_SCATTERING_OPACITY
) -> float:
    """Eddington luminosity (Lsun) for a star of `mass` (Msun),
    L_Edd = 4*pi*G*M*c/kappa_e -- standard formula, not paper-specific.
    Computed directly from cgs constants (see module docstring for why
    this replaced an earlier, incorrect memorized-shortcut version).
    """
    l_edd_erg_s = 4.0 * 3.141592653589793 * G_CGS * (mass * MSUN_G) * C_CGS / kappa_e
    return l_edd_erg_s / LSUN_ERG_S


def eddington_factor(
    mass: float, luminosity: float, kappa_e: float = ELECTRON_SCATTERING_OPACITY
) -> float:
    """Eddington factor Gamma = L/L_Edd (dimensionless), one of El
    Mellah & Casse's four shape parameters (their Sec. 2.2, `Gamma`).
    `luminosity` in Lsun (e.g. from
    stellar/main_sequence.py::ms_luminosity/hg_luminosity for the
    donor's own current luminosity).
    """
    return luminosity / eddington_luminosity(mass, kappa_e)


def wind_terminal_velocity(mass: float, radius: float, alpha: float) -> float:
    """CAK terminal wind velocity (km/s), El Mellah & Casse eq. (7):
    v_inf = 2.5 * v_esc * alpha/(1-alpha), with v_esc the star's own
    (unmodified, isolated-star) escape velocity. `alpha` is the CAK
    force multiplier (their Sec. 3.3(i): 0.45-0.65 for OB supergiants,
    citing Shimada et al. 1994).

    Uses eq. (7) -- their calibrated fit against "state-of-the-art
    hydrodynamical simulations" (Müller & Vink 2008; Noebauer & Sim
    2015), valid for effective temperatures beyond 21,000K (Vink et
    al. 1999) -- rather than the more basic point-source CAK form
    (their eq. 5, `v_inf = v_esc*sqrt(alpha/(1-alpha))`). Every
    realistic OB-supergiant HMXB donor this module targets sits above
    that temperature threshold, so eq. (7) is the appropriate choice
    for this project's actual use case -- confirmed directly: eq. (5)
    gives v_inf~=510 km/s for Vela X-1's own donor parameters (Friend
    & Castor Table 1: M=24 Msun, R=35 Rsun, alpha~0.5), well below
    the observed 700-1700 km/s range; eq. (7) gives ~1280 km/s,
    comfortably inside it. Not swapped in from memory -- caught by
    this direct numeric cross-check before finalizing the module, not
    assumed.
    """
    if not 0.0 < alpha < 1.0:
        raise ValueError(f"alpha must be in (0, 1), got {alpha}")
    v_esc = (2.0 * G_MSUN_RSUN_KM2_S2 * mass / radius) ** 0.5
    return 2.5 * v_esc * alpha / (1.0 - alpha)


def wind_velocity(r: float, radius: float, v_inf: float) -> float:
    """CAK wind velocity (km/s) at distance `r` (Rsun) from the
    star's centre, El Mellah & Casse eq. (7)'s calibrated profile:
    v(r) = v_inf * (1 - R*/r)^0.7 -- paired with
    `wind_terminal_velocity`'s own eq. (7) normalization (a different
    exponent, 0.7 not 0.5, from the basic point-source law, eq. 5 --
    the two are not meant to be mixed). `radius` is the star's own
    radius (Rsun). Isolated-star law -- see this module's docstring
    for the Roche-potential-modification scope cut.
    """
    if r <= radius:
        return 0.0
    return v_inf * (1.0 - radius / r) ** 0.7


def wind_mass_loss_rate(
    mass: float, luminosity: float, alpha: float, q_force: float, gamma: float
) -> float:
    """CAK wind mass-loss rate (Msun/yr), El Mellah & Casse eq. (20)
    (the Gayley 1995 Q-parametrization):

        Mdot = (1/(1+alpha))^(1/alpha) * [alpha/(1-alpha)]
               * Gamma * (Gamma*Q/(1-Gamma))^((1-alpha)/alpha)
               * L_Edd/c^2

    `q_force` is the Q force multiplier (their Sec. 3.3(ii): Q~900 for
    OB supergiants, expected to lie in 800-2000). `gamma` is the
    Eddington factor (`eddington_factor` above) -- must be < 1.
    `luminosity` (Lsun) is only used to derive L_Edd's normalization
    via `mass`, consistent with `gamma`'s own definition.
    """
    if not 0.0 < alpha < 1.0:
        raise ValueError(f"alpha must be in (0, 1), got {alpha}")
    if not 0.0 < gamma < 1.0:
        raise ValueError(f"gamma (Eddington factor) must be in (0, 1), got {gamma}")
    if q_force <= 0.0:
        raise ValueError(f"q_force must be > 0, got {q_force}")

    l_edd_erg_s = eddington_luminosity(mass) * LSUN_ERG_S
    l_edd_over_c2_g_s = l_edd_erg_s / C_CGS**2
    l_edd_over_c2_msun_per_yr = l_edd_over_c2_g_s * YEAR_S / MSUN_G

    prefactor = (1.0 / (1.0 + alpha)) ** (1.0 / alpha) * (alpha / (1.0 - alpha))
    q_term = (gamma * q_force / (1.0 - gamma)) ** ((1.0 - alpha) / alpha)
    return prefactor * gamma * q_term * l_edd_over_c2_msun_per_yr
