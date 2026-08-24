"""Giant-branch quantities needed as Hertzsprung-Gap boundary values.

Hurley, Pols & Tout (2000, MNRAS 315, 543), Section 5.1 (L_BGB, eq. 10)
and Section 5.2 (R_GB, eqs. 46-48). This module does NOT implement the
time-evolution of the giant branch itself (the M_c-L relation, eqs.
31-45, and the GB lifetime machinery) -- only the two static, purely
algebraic quantities HG's radius/luminosity endpoint needs:

    L_BGB(M, Z)       -- luminosity at the base of the GB (eq. 10)
    R_GB(M, L, Z)      -- GB radius as a function of L at fixed M (eq. 46)
    M_c,BGB(M)         -- core mass at the base of the GB (asymptotic
                           large-M approximation only -- see
                           `core_mass_bgb`'s own docstring)

Coefficients (a27-a32, b1-b7) were transcribed from Hurley et al.
(2000)'s Appendix A and verified by direct comparison against
user-pasted excerpts of the source PDF text -- see
docs/science/rlof-ce-classifier-proposal.md's "Update, 2026-08-24"
note for the two additional transcription errors this caught (a28's
eta exponent, and b4/b5's row-duplication) beyond the two already
caught in main_sequence.py's a40/a41.
"""

from __future__ import annotations

import numpy as np

# Appendix A rows, same (alpha, beta, gamma, eta, mu) convention as
# main_sequence.py's _A table: a_n(Z) = alpha + beta*zeta + gamma*zeta^2
# + eta*zeta^3 + mu*zeta^4, zeta = log10(Z/0.02).
_A = {
    27: (9.511033e1, 6.819618e1, -1.045625e1, -1.474939e1, 0.0),
    28: (3.113458e1, 1.012033e1, -4.650511, -2.463185, 0.0),
    "29p": (1.413057, 4.578814e-1, -6.850581e-2, -5.588658e-2, 0.0),
    30: (3.910862e1, 5.196646e1, 2.264970e1, 2.873680, 0.0),
    31: (4.597479, -2.855179e-1, 2.709724e-1, 0.0, 0.0),
    32: (6.682518, 2.827718e-1, -7.294429e-2, 0.0, 0.0),
}

# Fixed literal constants in eq. 10 (not Z-dependent).
_C2 = 9.301992
_C3 = 4.637345

# b1, b4-b7 rows, same convention (only 3 columns given for b1 in the
# source -- eta and mu are 0, the paper's own "blank means zero"
# convention).
_B = {
    1: (3.9700e-1, 2.8826e-1, 5.2930e-1, 0.0, 0.0),
    4: (9.960283e-1, 8.164393e-1, 2.38383, 2.223436, 8.638115e-1),
    5: (2.561062e-1, 7.072646e-2, -5.444596e-2, -5.798167e-2, -1.349129e-2),
    6: (1.157338, 1.467883, 4.299661, 3.130500, 6.992080e-1),
    7: (4.022765e-1, 3.05001e-1, 9.962137e-1, 7.914079e-1, 1.728098e-1),
}


def _zeta(z: float) -> float:
    return np.log10(z / 0.02)


def _a(n: int | str, zeta: float) -> float:
    alpha, beta, gamma, eta, mu = _A[n]
    return alpha + beta * zeta + gamma * zeta**2 + eta * zeta**3 + mu * zeta**4


def _b(n: int, zeta: float) -> float:
    alpha, beta, gamma, eta, mu = _B[n]
    return alpha + beta * zeta + gamma * zeta**2 + eta * zeta**3 + mu * zeta**4


def l_bgb(mass: float, z: float) -> float:
    """Luminosity at the base of the giant branch (Lsun). Eq. (10).

    L_BGB = (a27*M^a31 + a28*M^c2) / (a29 + a30*M^c3 + M^a32),
    c2 = 9.301992, c3 = 4.637345 (fixed constants, not Z-dependent).
    """
    zeta = _zeta(z)
    a27, a28, a30, a31, a32 = (_a(n, zeta) for n in (27, 28, 30, 31, 32))
    a29 = _a("29p", zeta) ** a32
    m = mass
    numerator = a27 * m**a31 + a28 * m**_C2
    denominator = a29 + a30 * m**_C3 + m**a32
    return numerator / denominator


def mass_radius_exponent(z: float) -> float:
    """x, the exponent to which R depends on M at constant L on the GB
    (R_GB proportional to M^-x). Eq. (47).
    """
    zeta = _zeta(z)
    return (
        0.30406 + 0.0805 * zeta + 0.0897 * zeta**2 + 0.0878 * zeta**3 + 0.0222 * zeta**4
    )


def r_gb(mass: float, luminosity: float, z: float) -> float:
    """GB radius (Rsun) as a function of L at fixed mass. Eqs. (46), (48).

    R_GB = A*(L^b1 + b2*L^b3), A = min(b4*M^-b5, b6*M^-b7).
    """
    zeta = _zeta(z)
    sigma = np.log10(z)

    b1 = min(0.54, _b(1, zeta))

    b2_raw = 10.0 ** (-4.6739 - 0.9394 * sigma)
    b2 = min(max(b2_raw, -0.04167 + 55.67 * z), 0.4771 - 9329.21 * z**2.94)

    b3_prime = max(-0.1451, -2.2794 - 1.5175 * sigma - 0.254 * sigma**2)
    b3 = 10.0**b3_prime
    if z > 0.004:
        b3 = max(b3, 0.7307 + 14265.1 * z**3.395)

    b4 = _b(4, zeta) + 0.1231572 * zeta**5
    b5 = _b(5, zeta)
    b6 = _b(6, zeta) + 0.01640687 * zeta**5
    b7 = _b(7, zeta)

    m = mass
    a_coeff = min(b4 * m**-b5, b6 * m**-b7)
    return a_coeff * (luminosity**b1 + b2 * luminosity**b3)


# eq. (44)'s own constants (Z-independent), read directly from the
# source text: "M_c,BGB = min[0.95*M_c,BAGB, (C + c1*M^c2)^(1/4)] ...
# The constants c1 = 9.20925e-5 and c2 = 5.402216 are independent of
# Z, so that for large enough M we have M_c,BGB ~= 0.098*M^1.35". Used
# here ONLY to derive the exact mass above which that asymptotic
# approximation exceeds the Chandrasekhar mass -- NOT to implement the
# full eq. (44) itself, which additionally needs `C` (the inverse of
# the full GB core-mass-luminosity relation, eqs. 31-43 -- a
# substantially larger addition than this module's scope; see
# docs/science/rlof-ce-classifier-proposal.md). Verified self-
# consistent against the paper's own stated asymptotic coefficients:
# c1**0.25 = 0.0979... ~= 0.098, c2/4 = 1.350554 ~= 1.35.
_EQ44_C1 = 9.20925e-5
_EQ44_C2 = 5.402216

# Chandrasekhar mass (Msun) -- duplicated from remnant.py's own
# M_CHANDRASEKHAR rather than imported, to avoid a circular import
# (remnant.py already imports main_sequence.py, which imports this
# module).
_M_CHANDRASEKHAR = 1.44

# Mass (Msun) above which 0.098*M^1.35 exceeds the Chandrasekhar mass,
# derived directly from eq. 44's own c1/c2 (not a separately-chosen
# cutoff): M_max = (M_ch^4 / c1)^(1/c2).
CORE_MASS_BGB_MAX_MASS = (_M_CHANDRASEKHAR**4 / _EQ44_C1) ** (1.0 / _EQ44_C2)


def core_mass_bgb(mass: float) -> float:
    """Core mass at the base of the giant branch (Msun), for
    intermediate-mass stars (M_HeF <= M < min(M_FGB,
    CORE_MASS_BGB_MAX_MASS)) -- the branch of Hurley et al. (2000)
    eq. (28) relevant to Realta's HG-donor mass range.

    NOT the full eq. (44) formula -- that needs `C`, the inverse of
    the full GB core-mass-luminosity relation (eqs. 31-43, needing a
    mass-dependent hydrogen-burning rate "constant" and piecewise
    D/p functions) and M_c,BAGB (eq. 66, coefficients b36-b38,
    separately verified but not usable without `C`). Implementing eq.
    44 properly is a substantially larger addition than this module's
    scope -- comparable in size to the MS or HG modules, not a
    coefficient tweak (see
    docs/science/rlof-ce-classifier-proposal.md's "Core mass"
    section). This function instead uses the Z-independent large-mass
    asymptotic limit Hurley et al. state directly in their Sec. 5.2
    text: "for large enough M we have M_c,BGB ~= 0.098*M^1.35,
    independent of Z."

    Raises ValueError for mass >= CORE_MASS_BGB_MAX_MASS (~7.3 Msun):
    above that mass the asymptotic approximation itself gives a
    super-Chandrasekhar "core mass" at the *end of HG* -- physically
    implausible (real HG/GB-base core masses for stars in this range
    are a few tenths of a solar mass; substantial core growth happens
    much later) and not something to silently clamp or guess past.
    This was found and corrected during implementation, not assumed --
    see docs/science/rlof-ce-classifier-proposal.md's "Core mass:
    implementation note" for the full account, including the earlier,
    now-superseded framing of this as valid "for M >= M_HeF".
    """
    if mass >= CORE_MASS_BGB_MAX_MASS:
        raise ValueError(
            f"mass={mass} Msun >= {CORE_MASS_BGB_MAX_MASS:.3f} Msun: "
            "the 0.098*M^1.35 asymptotic approximation for M_c,BGB is "
            "super-Chandrasekhar (physically implausible) above this "
            "mass. The full eq. (44) formula, which stays valid here, "
            "is not implemented -- see this function's docstring."
        )
    return 0.098 * mass**1.35
