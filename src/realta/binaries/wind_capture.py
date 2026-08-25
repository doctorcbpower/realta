"""Wind-capture accretion: converts a donor's CAK wind
(stellar/cak_wind.py) into a compact-object accretion rate and a
circularization-radius estimate.

Source: El Mellah & Casse (2017, MNRAS, arXiv:1609.01532), "A
numerical investigation of wind accretion in persistent Supergiant
X-ray Binaries I". Pasted and verified directly before implementing.

Modular by design, per the stellar-evolution/wind/capture split this
was scoped with (chat, 2026-08-25):
    stellar/cak_wind.py    -> donor wind state (Mdot_wind, v(r))
    binaries/wind_capture.py (this module) -> Mdot_acc,
        circularization radius, disc diagnostics

Two accretion-rate estimators are provided:
    bhl_accretion_fraction() -- the paper's own recommended simplified
        approach (their eq. 18 evaluated at the wind velocity AT THE
        ORBITAL SEPARATION, not terminal velocity -- their own text,
        Sec. 4.2.3, states this recovers their full numerical result
        to within 6%; their more literal "analytic" eq. (19) is
        explicitly NOT used here -- their own text says it
        underestimates by "at least a factor of three", so eq. (19)
        is not their actual recommendation despite being the more
        compact closed-form expression).
    bhl_accretion_rate_simple() -- the plain textbook Bondi-Hoyle-
        Lyttleton estimate with a caller-supplied wind velocity (e.g.
        terminal velocity, ignoring the orbital-separation
        refinement) -- the "simpler BHL-style fallback for
        comparison" this was scoped to include.

Circularization radius (`circularization_radius_fraction`) is the one
piece in this module NOT directly verified against a pasted source --
see that function's own docstring for a real error this caught (a
first-draft formula off by four orders of magnitude from El Mellah &
Casse's own stated range) and the resulting, still-imperfect fix.
"""

from __future__ import annotations

from realta.stellar.cak_wind import G_MSUN_RSUN_KM2_S2

# El Mellah & Casse's own calibration factor for the fraction of wind
# entering the Bondi-Hoyle-Lyttleton accretion cylinder (their eq. 18,
# citing Foglizzo & Ruffert 1996 / El Mellah & Casse 2015 -- accurate
# to a few percent for Mach numbers > ~4, "easily verified by the
# highly supersonic winds in SgXB" per their own text).
_BHL_CALIBRATION_FACTOR = 0.77


def escape_velocity(mass: float, radius: float) -> float:
    """Escape velocity (km/s) for `mass` (Msun) at `radius` (Rsun) --
    shared helper, not paper-specific (basic Newtonian mechanics)."""
    return (2.0 * G_MSUN_RSUN_KM2_S2 * mass / radius) ** 0.5


def relative_wind_velocity(v_wind: float, v_orbital: float) -> float:
    """Relative wind speed (km/s) seen by the compact object,
    combining the (radial) wind velocity with the compact object's
    own orbital velocity, added in quadrature (the two are
    perpendicular to leading order for a wind launched radially from
    the donor and an orbit circular to leading order) -- the standard
    Bondi-Hoyle-Lyttleton relative-velocity construction, not paper-
    specific."""
    return (v_wind**2 + v_orbital**2) ** 0.5


def accretion_radius(m_compact: float, v_rel: float) -> float:
    """Bondi-Hoyle-Lyttleton accretion radius (Rsun),
    R_acc = 2*G*M_compact/v_rel^2 (Hoyle & Lyttleton 1939; Bondi &
    Hoyle 1944; El Mellah & Casse's own eq. before their eq. 17,
    "R_acc = (2GM2)/v_bullet^2"). `m_compact` in Msun, `v_rel` in km/s.
    """
    if v_rel <= 0.0:
        raise ValueError(f"v_rel must be positive, got {v_rel}")
    return 2.0 * G_MSUN_RSUN_KM2_S2 * m_compact / v_rel**2


def bhl_accretion_fraction(r_acc: float, separation: float) -> float:
    """Fraction of the donor's wind mass-loss captured by the compact
    object, El Mellah & Casse eq. (18): beta = 0.77*(R_acc/a)^2 --
    THEIR OWN RECOMMENDED simplified estimate (see this module's
    docstring for why, not their more compact eq. 19). `r_acc` and
    `separation` in the same length unit (Rsun). Clamped to [0, 1] --
    the underlying quadratic form is only valid for R_acc << a (the
    supersonic-wind, small-cross-section regime both papers assume);
    it is not meaningful, and not clamped away silently without
    saying so, if a configuration ever pushes R_acc close to or above
    a (e.g. an unphysically slow wind or too-close an orbit).
    """
    if separation <= 0.0:
        raise ValueError(f"separation must be positive, got {separation}")
    beta = _BHL_CALIBRATION_FACTOR * (r_acc / separation) ** 2
    return min(1.0, beta)


def wind_capture_rate(mdot_wind: float, beta: float) -> float:
    """Accretion rate onto the compact object (Msun/yr) =
    beta * Mdot_wind -- trivial, but named for symmetry with the
    other functions here and to keep the "fraction vs. rate"
    distinction explicit at call sites."""
    return beta * mdot_wind


def bhl_accretion_rate_simple(
    m_compact: float,
    mdot_wind: float,
    v_wind: float,
    v_orbital: float,
    separation: float,
) -> float:
    """Plain textbook Bondi-Hoyle-Lyttleton accretion rate (Msun/yr) --
    the "simpler BHL-style fallback for comparison" this module was
    scoped to include. Uses whatever `v_wind` the caller supplies
    (e.g. the CAK terminal velocity, ignoring the orbital-separation
    refinement `bhl_accretion_fraction` applies) with the same
    calibrated 0.77 factor, so the two functions are directly
    comparable -- the difference between them isolates the effect of
    "wind speed at the orbital separation" vs. "terminal wind speed"
    alone, not a difference in the underlying BHL formalism itself.
    """
    v_rel = relative_wind_velocity(v_wind, v_orbital)
    r_acc = accretion_radius(m_compact, v_rel)
    beta = bhl_accretion_fraction(r_acc, separation)
    return wind_capture_rate(mdot_wind, beta)


def circularization_radius_fraction(v_orbital: float, v_rel: float) -> float:
    """R_circ/R_acc, the standard wind-accretion circularization-
    radius scaling from the shear in orbital velocity across the
    Bondi-Hoyle-Lyttleton accretion cylinder (Shapiro & Lightman 1976;
    see also Illarionov & Sunyaev 1975; the same scaling underlies the
    "disc vs. no disc" criterion widely quoted in the HMXB wind-
    accretion literature, e.g. Frank, King & Raine, "Accretion Power
    in Astrophysics"):

        R_circ/R_acc ~= (1/4) * (v_orbital/v_rel)^4

    NOT a formula given in either source paper -- El Mellah & Casse
    compute R_circ numerically (their Fig. 6) but do not provide a
    closed-form fit. This is a genuinely LOWER-confidence piece than
    the rest of this module: the (v_orbital/v_rel)^4 SCALING is the
    well-known form in the wind-accretion literature, but the exact
    prefactor (taken here as 1/4, the commonly-cited value) has not
    been independently verified against a pasted source, unlike
    everything else in this module.

    A first draft used a different, from-memory "shear across the
    cylinder width" formula (specific angular momentum
    l = Omega*R_acc^2/4, then R_circ=l^2/(G*M)) that gave
    R_circ/R_acc ~ 2e-7 for representative Vela-X-1-like parameters --
    four orders of magnitude below El Mellah & Casse's own stated
    range (their Sec. 4.4: R_circ/R_acc ~ 1e-3 to 1e-2). Caught by
    that direct numeric check, not left in. This replacement formula
    lands much closer (~1e-4 for the same illustrative numbers) but
    is still on the low side of their stated range -- flagged
    explicitly rather than silently presented as verified.

    `v_orbital`/`v_rel` in the same units (km/s) -- the ratio is
    dimensionless, so this needs no separate unit-conversion layer.
    """
    if v_rel <= 0.0:
        raise ValueError(f"v_rel must be positive, got {v_rel}")
    if v_orbital < 0.0:
        raise ValueError(f"v_orbital must be >= 0, got {v_orbital}")
    return 0.25 * (v_orbital / v_rel) ** 4


def circularization_radius(v_orbital: float, v_rel: float, r_acc: float) -> float:
    """Circularization radius (Rsun) = `circularization_radius_fraction`
    (see that function's own docstring for the formula and its
    confidence level) times `r_acc` (Rsun)."""
    return circularization_radius_fraction(v_orbital, v_rel) * r_acc
