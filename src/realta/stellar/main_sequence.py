"""Main-sequence radius, luminosity, lifetime and phase.

Hurley, Pols & Tout (2000, MNRAS 315, 543), Section 5.1 (main sequence,
eqs. 1-24) and Appendix A (coefficients a1-a81). Scope note: this
module implements ONLY the zero-age-main-sequence-to-terminal-main-
sequence phase (Hurley's stellar type k=0,1) plus the critical masses
and coefficient machinery it depends on.

Hertzsprung-Gap (k=2) and later phases are deliberately NOT
implemented here -- see docs/science/rlof-ce-classifier-proposal.md
and docs/provenance.md for why: HG's radius formula (eq. 27) requires
the giant-branch radius-luminosity relation (eq. 44-48), whose
coefficients (b1, b4-b7) could not be transcribed with confidence from
the source PDF (a repeated-digit anomaly was found across supposedly-
independent rows on re-reading -- a strong signal of a transcription
defect, not a real coincidence). Rather than hard-code a possibly-wrong
coefficient block into a codebase whose whole discipline is protecting
a trusted baseline, that phase is deferred until those coefficients are
re-verified. Callers must treat any t >= t_MS(mass, z) as "phase not
modelled by this module" -- see `phase()` below, which raises rather
than silently returning wrong physics.

Unit conventions (matching the paper): mass in Msun, radius in Rsun,
luminosity in Lsun, time/ages in Myr (the paper's own "unless otherwise
specified" convention, Section 5).
"""

from __future__ import annotations

import numpy as np

from realta.stellar import giant_branch
from realta.stellar.zams import zams_luminosity, zams_radius

# Appendix A (Hurley et al. 2000, pp. 566-568): each entry is the
# (alpha, beta, gamma, eta, mu) row for a_n(Z) = alpha + beta*zeta +
# gamma*zeta**2 + eta*zeta**3 + mu*zeta**4, zeta = log10(Z/0.02).
# A 0.0 in a slot means that table column was blank in the paper (the
# paper's own convention: "A blank entry in a table implies a zero
# value"). Only coefficients actually used by the MS-only scope of
# this module are included; GB/CHeB/AGB-only coefficients (a27-a32,
# b-series) are intentionally omitted -- see module docstring.
_A = {
    1: (1.593890e3, 2.053038e3, 1.231226e3, 2.327785e2, 0.0),
    2: (2.706708e3, 1.483131e3, 5.772723e1, 7.411230e1, 0.0),
    3: (1.466143e2, -1.048442e2, -6.795374e1, -1.391127e1, 0.0),
    4: (4.141990e-2, 4.564888e-2, 2.958542e-2, 5.571483e-3, 0.0),
    5: (3.426349e-1, 0.0, 0.0, 0.0, 0.0),
    6: (1.949814e1, 1.758178, -6.008212, -4.470533, 0.0),
    7: (4.903830, 0.0, 0.0, 0.0, 0.0),
    8: (5.212154e-2, 3.166411e-2, -2.750074e-3, -2.271549e-3, 0.0),
    9: (1.312179, -3.294936e-1, 9.231860e-2, 2.610989e-2, 0.0),
    10: (8.073972e-1, 0.0, 0.0, 0.0, 0.0),
    # a11 = a'11 * a14, a12 = a'12 * a14 -- primed rows below.
    "11p": (1.031538, -2.434480e-1, 7.732821, 6.460705, 1.374484),
    "12p": (1.043715, -1.577474e-1, -5.168234, -5.596506, -1.299394),
    13: (7.859875e2, -8.542048, -2.642511e1, -9.585707, 0.0),
    14: (3.858911e3, 2.459681e3, -7.630093e1, -3.486057e2, -4.861703e1),
    15: (2.888720e1, 2.952979e2, 1.850341e2, 3.797254e1, 0.0),
    16: (7.196580, 5.613746e-1, 3.805871e-1, 8.398728e-2, 0.0),
    "18p": (2.187715e-1, -2.154437, -3.768678, -1.975518, -3.021475e-1),
    "19p": (1.466440, 1.839725, 6.442199, 4.023635, 6.957529e-1),
    20: (2.652091e1, 8.178458e1, 1.156058e2, 7.633811e1, 1.950698e1),
    21: (1.472103, -2.947609, -3.312828, -9.945065e-1, 0.0),
    22: (3.071048, -5.679941e-1, -9.745523e-1, -3.594543e-1, 0.0),
    23: (2.617890, 1.019135, -3.292551e-1, -7.445123e-2, 0.0),
    24: (1.075567e-2, 1.773287e-2, 9.610479e-3, 1.732469e-3, 0.0),
    25: (1.476246, 1.899331, 1.195010, 3.035051e-1, 0.0),
    26: (5.502535, -6.601663e-2, 9.968207e-2, 3.599801e-2, 0.0),
    34: (1.910302e-1, 1.158624e-1, 3.348990e-2, 2.599706e-3, 0.0),
    35: (3.931056e-1, 7.277637e-2, -1.366593e-1, -4.508946e-2, 0.0),
    36: (3.267776e-1, 1.204424e-1, 9.988332e-2, 2.455361e-2, 0.0),
    37: (5.990212e-1, 5.570264e-2, 6.207626e-2, 1.777283e-2, 0.0),
    38: (7.330122e-1, 5.192827e-1, 2.316416e-1, 8.346941e-3, 0.0),
    39: (1.172768, -1.209262e-1, -1.193023e-1, -2.859837e-2, 0.0),
    40: (3.982622e-1, -2.296279e-1, -2.262539e-1, -5.219837e-2, 0.0),
    41: (3.571038, -2.223625e-2, -2.611794e-2, -6.359648e-3, 0.0),
    42: (1.9848, 1.1386, 3.5640e-1, 0.0, 0.0),
    43: (6.300e-2, 4.810e-2, 9.840e-3, 0.0, 0.0),
    44: (1.200, 2.450, 0.0, 0.0, 0.0),
    45: (2.321400e-1, 1.828075e-3, -2.232007e-2, -3.378734e-3, 0.0),
    46: (1.163659e-2, 3.427682e-3, 1.421393e-3, -3.710666e-3, 0.0),
    47: (1.048020e-2, -1.231921e-2, -1.686860e-2, -4.234354e-3, 0.0),
    48: (1.555590, -3.223927e-1, -5.197429e-1, -1.066441e-1, 0.0),
    49: (9.7700e-2, -2.3100e-1, -7.5300e-2, 0.0, 0.0),
    50: (2.4000e-1, 1.8000e-1, 5.9500e-1, 0.0, 0.0),
    51: (3.3000e-1, 1.3200e-1, 2.1800e-1, 0.0, 0.0),
    52: (1.1064, 4.1500e-1, 1.8000e-1, 0.0, 0.0),
    53: (1.1900, 3.7700e-1, 1.7600e-1, 0.0, 0.0),
    54: (3.855707e-1, -6.104166e-1, 5.676742, 1.060894e1, 5.284014),
    55: (3.579064e-1, -6.442936e-1, 5.494644, 1.054952e1, 5.280991),
    56: (9.587587e-1, 8.777464e-1, 2.017321e-1, 0.0, 0.0),
    58: (4.907546e-1, -1.683928e-1, -3.108742e-1, -7.202918e-2, 0.0),
    59: (4.537070e-1, -4.465455e-1, -1.612690e-1, -1.623246e-2, 0.0),
    60: (1.796220, 2.814020e-1, 1.423325, 3.421036e-1, 0.0),
    61: (2.256216, 3.773400e-1, 1.537867, 4.396373e-1, 0.0),
    62: (8.4300e-2, -4.7500e-2, -3.5200e-2, 0.0, 0.0),
    63: (7.3600e-2, 7.4900e-2, 4.4260e-2, 0.0, 0.0),
    64: (1.3600e-1, 3.5200e-2, 0.0, 0.0, 0.0),
    65: (1.564231e-3, 1.653042e-3, -4.439786e-3, -4.951011e-3, -1.216530e-3),
    66: (1.4770, 2.9660e-1, 0.0, 0.0, 0.0),
    67: (5.210157, -4.143695, -2.120870, 0.0, 0.0),
    68: (1.1160, 1.6600e-1, 0.0, 0.0, 0.0),
    69: (1.071489, -1.164852e-1, -8.623831e-2, -1.582349e-2, 0.0),
    70: (7.108492e-1, 7.935927e-1, 3.926983e-1, 3.622146e-2, 0.0),
    71: (3.478514, -2.585474e-2, -1.512955e-2, -2.833691e-3, 0.0),
    72: (9.132108e-1, -1.653695e-1, 0.0, 3.636784e-2, 0.0),
    73: (3.969331e-3, 4.539076e-3, 1.720906e-3, 1.897857e-4, 0.0),
    74: (1.600, 7.640e-1, 3.322e-1, 0.0, 0.0),
    75: (8.109e-1, -6.282e-1, 0.0, 0.0, 0.0),
    76: (1.192334e-2, 1.083057e-2, 1.230969, 1.551656, 0.0),
    77: (-1.668868e-1, 5.818123e-1, -1.105027e1, -1.668070e1, 0.0),
    78: (7.615495e-1, 1.068243e-1, -2.011333e-1, -9.371415e-2, 0.0),
    79: (9.409838, 1.522928, 0.0, 0.0, 0.0),
    80: (-2.7110e-1, -5.7560e-1, -8.3800e-2, 0.0, 0.0),
    81: (2.4930, 1.1475, 0.0, 0.0, 0.0),
}


def _zeta(z: float) -> float:
    return np.log10(z / 0.02)


def _a(n: int | str, zeta: float) -> float:
    alpha, beta, gamma, eta, mu = _A[n]
    return alpha + beta * zeta + gamma * zeta**2 + eta * zeta**3 + mu * zeta**4


def m_hook(z: float) -> float:
    """Initial mass above which a star develops a MS 'hook'. Eq. (1)."""
    zeta = _zeta(z)
    return 1.0185 + 0.16015 * zeta + 0.0892 * zeta**2


def t_bgb(mass: float, z: float) -> float:
    """Time to reach the base of the giant branch, Myr. Eq. (4)."""
    zeta = _zeta(z)
    a1, a2, a3, a4, a5 = (_a(n, zeta) for n in (1, 2, 3, 4, 5))
    m = mass
    return (a1 + a2 * m**4 + a3 * m**5.5 + m**7) / (a4 * m**2 + a5 * m**7)


def _mu_coefficient(mass: float, z: float) -> float:
    """Eq. (7): mu, used to define t_hook = mu * t_BGB."""
    zeta = _zeta(z)
    a6, a7, a8, a9, a10 = (_a(n, zeta) for n in (6, 7, 8, 9, 10))
    return max(0.5, 1.0 - 0.01 * max(a6 / mass**a7, a8 + a9 / mass**a10))


def _x_coefficient(z: float) -> float:
    """Eq. (6)."""
    zeta = _zeta(z)
    return max(0.95, min(0.95 - 0.03 * (zeta + 0.30103), 0.99))


def t_hook(mass: float, z: float) -> float:
    return _mu_coefficient(mass, z) * t_bgb(mass, z)


def t_ms(mass: float, z: float) -> float:
    """Main-sequence lifetime, Myr. Eq. (5)."""
    return max(t_hook(mass, z), _x_coefficient(z) * t_bgb(mass, z))


def l_tms(mass: float, z: float) -> float:
    """Luminosity at the end of the MS (Lsun). Eq. (8)."""
    zeta = _zeta(z)
    a14 = _a(14, zeta)
    a11 = _a("11p", zeta) * a14
    a12 = _a("12p", zeta) * a14
    a13, a15, a16 = (_a(n, zeta) for n in (13, 15, 16))
    m = mass
    return (a11 * m**3 + a12 * m**4 + a13 * m ** (a16 + 1.8)) / (
        a14 + a15 * m**5 + m**a16
    )


def _a17(z: float) -> float:
    """The M* boundary mass for eq. (9a)/(9b), via eq. numbered inline
    on p. 548 of Hurley et al. (2000): note this uses sigma = log10(Z),
    not zeta = log10(Z/0.02).
    """
    sigma = np.log10(z)
    log_a17 = max(
        0.097 - 0.1072 * (sigma + 3),
        max(0.097, min(0.1461, 0.1461 + 0.1237 * (sigma + 2))),
    )
    return min(1.6, max(1.4, 10.0**log_a17))


def r_tms(mass: float, z: float) -> float:
    """Radius at the end of the MS (Rsun). Eqs. (9a)/(9b)."""
    zeta = _zeta(z)
    a17 = _a17(z)
    m_star = a17 + 0.1
    c1 = -8.672073e-2

    def r_low(m: float) -> float:
        a18 = _a("18p", zeta) * _a(20, zeta)
        a19 = _a("19p", zeta) * _a(20, zeta)
        a20, a21, a22 = (_a(n, zeta) for n in (20, 21, 22))
        return (a18 + a19 * m**a21) / (a20 + m**a22)

    def r_high(m: float) -> float:
        a23, a24, a25, a26 = (_a(n, zeta) for n in (23, 24, 25, 26))
        return (c1 * m**3 + a23 * m**a26 + a24 * m ** (a26 + 1.5)) / (a25 + m**5)

    if mass <= a17:
        return r_low(mass)
    if mass >= m_star:
        return r_high(mass)
    # Straight-line interpolation between the two end-points (Hurley
    # et al. 2000, p. 548, "with straight-line interpolation to
    # connect equations (9a) and (9b) between the end-points").
    r_a17 = r_low(a17)
    r_mstar = r_high(m_star)
    frac = (mass - a17) / (m_star - a17)
    return r_a17 + frac * (r_mstar - r_a17)


def _delta_l(mass: float, z: float) -> float:
    """Luminosity perturbation, eq. (16)."""
    zeta = _zeta(z)
    hook = m_hook(z)
    a33 = min(1.4, 1.5135 + 0.3769 * zeta)
    a33 = max(0.6355 - 0.4192 * zeta, max(1.25, a33))
    a34, a35, a36, a37 = (_a(n, zeta) for n in (34, 35, 36, 37))

    def branch3(m: float) -> float:
        return min(a34 / m**a35, a36 / m**a37)

    if mass <= hook:
        return 0.0
    b = branch3(a33)
    if mass < a33:
        return b * ((mass - hook) / (a33 - hook)) ** 0.4
    return branch3(mass)


def _delta_r(mass: float, z: float) -> float:
    """Radius perturbation, eq. (17)."""
    zeta = _zeta(z)
    hook = m_hook(z)
    a38, a39, a40, a41 = (_a(n, zeta) for n in (38, 39, 40, 41))
    a42 = min(1.25, max(1.1, _a(42, zeta)))
    a43 = _a(43, zeta)
    a44 = min(1.3, max(0.45, _a(44, zeta)))

    def branch4(m: float) -> float:
        return (a38 + a39 * m**3.5) / (a40 * m**3 + m**a41) - 1.0

    if mass <= hook:
        return 0.0
    b = branch4(2.0)
    if mass <= a42:
        return a43 * ((mass - hook) / (a42 - hook)) ** 0.5
    if mass < 2.0:
        return a43 + (b - a43) * ((mass - a42) / (2.0 - a42)) ** a44
    return branch4(mass)


def _alpha_l(mass: float, z: float) -> float:
    """Eq. (19a)/(19b)."""
    zeta = _zeta(z)
    a48 = _a(48, zeta)

    def high(m: float) -> float:
        a45, a46, a47 = (_a(n, zeta) for n in (45, 46, 47))
        return (a45 + a46 * m**a48) / (m**0.4 + a47 * m**1.9)

    if mass >= 2.0:
        return high(mass)

    a49 = max(_a(49, zeta), 0.145)
    a50 = min(_a(50, zeta), 0.306 + 0.053 * zeta)
    a51 = min(_a(51, zeta), 0.3625 + 0.062 * zeta)
    a52 = max(_a(52, zeta), 0.9)
    if z > 0.01:
        a52 = min(a52, 1.0)
    a53 = max(_a(53, zeta), 1.0)
    if z > 0.01:
        a53 = min(a53, 1.1)

    if mass < 0.5:
        return a49
    if mass < 0.7:
        return a49 + 5.0 * (0.3 - a49) * (mass - 0.5)
    if mass < a52:
        return 0.3 + (a50 - 0.3) * (mass - 0.7) / (a52 - 0.7)
    if mass < a53:
        return a50 + (a51 - a50) * (mass - a52) / (a53 - a52)
    b = high(2.0)
    return a51 + (b - a51) * (mass - a53) / (2.0 - a53)


def _beta_l(mass: float, z: float) -> float:
    """Eq. (20)."""
    zeta = _zeta(z)
    a54, a55, a56 = (_a(n, zeta) for n in (54, 55, 56))
    a57 = min(1.4, 1.5135 + 0.3769 * zeta)
    a57 = max(0.6355 - 0.4192 * zeta, max(1.25, a57))

    beta = max(0.0, a54 - a55 * mass**a56)
    if mass > a57 and beta > 0.0:
        b = max(0.0, a54 - a55 * a57**a56)
        beta = max(0.0, b - 10.0 * (mass - a57) * b)
    return beta


def _alpha_r(mass: float, z: float) -> float:
    """Eq. (21a)/(21b)."""
    zeta = _zeta(z)
    a58, a59, a60, a61 = (_a(n, zeta) for n in (58, 59, 60, 61))

    def mid(m: float) -> float:
        return a58 * m**a60 / (a59 + m**a61)

    a62 = max(0.065, _a(62, zeta))
    a63 = _a(63, zeta)
    if z < 0.004:
        a63 = min(0.055, a63)
    a64 = max(0.091, min(0.121, _a(64, zeta)))
    a66 = max(_a(66, zeta), min(1.6, -0.308 - 1.046 * zeta))
    a66 = max(0.8, min(0.8 - 2.0 * zeta, a66))
    a67 = _a(67, zeta)
    a68 = max(0.9, min(_a(68, zeta), 1.0))
    if a68 > a66:
        a64 = mid(a66)
    a68 = min(a68, a66)

    if a66 <= mass <= a67:
        return mid(mass)
    if mass < 0.5:
        return a62
    if mass < 0.65:
        return a62 + (a63 - a62) * (mass - 0.5) / 0.15
    if mass < a68:
        return a63 + (a64 - a63) * (mass - 0.65) / (a68 - 0.65)
    if mass < a66:
        b = mid(a66)
        return a64 + (b - a64) * (mass - a68) / (a66 - a68)
    c = mid(a67)
    a65 = _a(65, zeta)
    return c + a65 * (mass - a67)


def _beta_r(mass: float, z: float) -> float:
    """Eq. (22a)/(22b). Returns beta_R = beta'_R - 1."""
    zeta = _zeta(z)
    a69, a70, a71 = (_a(n, zeta) for n in (69, 70, 71))

    def mid(m: float) -> float:
        return a69 * m**3.5 / (a70 + m**a71)

    a72 = _a(72, zeta)
    if z > 0.01:
        a72 = max(a72, 0.95)
    a73 = _a(73, zeta)
    a74 = max(1.4, min(_a(74, zeta), 1.6))

    if 2.0 <= mass <= 16.0:
        beta_prime = mid(mass)
    elif mass <= 1.0:
        beta_prime = 1.06
    elif mass < a74:
        beta_prime = 1.06 + (a72 - 1.06) * (mass - 1.0) / (a74 - 1.06)
    elif mass <= 2.0:
        b = mid(2.0)
        beta_prime = a72 + (b - a72) * (mass - a74) / (2.0 - a74)
    else:
        c = mid(16.0)
        beta_prime = c + a73 * (mass - 16.0)
    return beta_prime - 1.0


def _gamma_coefficient(mass: float, z: float) -> float:
    """Eq. (23)."""
    zeta = _zeta(z)
    a75 = max(1.0, min(_a(75, zeta), 1.27))
    a75 = max(a75, 0.6355 - 0.4192 * zeta)
    if mass > a75 + 0.1:
        return 0.0

    a76 = max(_a(76, zeta), -0.1015564 - 0.2161264 * zeta - 0.05182516 * zeta**2)
    a77 = max(
        -0.3868776 - 0.5457078 * zeta - 0.1463472 * zeta**2, min(0.0, _a(77, zeta))
    )
    a78 = max(0.0, min(_a(78, zeta), 7.454 + 9.046 * zeta))
    a79 = min(_a(79, zeta), max(2.0, -13.3 - 18.6 * zeta))
    a80 = max(0.0585542, _a(80, zeta))
    a81 = min(1.5, max(0.4, _a(81, zeta)))

    def low(m: float) -> float:
        return a76 + a77 * abs(m - a78) ** a79

    if mass <= 1.0:
        return low(mass)
    b = low(1.0)
    if mass <= a75:
        return b + (a80 - b) * ((mass - 1.0) / (a75 - 1.0)) ** a81
    c = a80 if a75 != 1.0 else b
    return c - 10.0 * (mass - a75) * c


def _tau1(t: float, hook_time: float) -> float:
    return min(1.0, t / hook_time)


def _tau2(t: float, hook_time: float, epsilon: float = 0.01) -> float:
    return max(0.0, min(1.0, (t - (1.0 - epsilon) * hook_time) / (epsilon * hook_time)))


def _eta_exponent(mass: float, z: float) -> float:
    """Eq. (18)."""
    if z > 0.0009:
        return 10.0
    if mass <= 1.0:
        return 10.0
    if mass >= 1.1:
        return 20.0
    return 10.0 + (20.0 - 10.0) * (mass - 1.0) / 0.1


def ms_luminosity(mass: float, z: float, t: float) -> float:
    """L_MS(mass, z, t) in Lsun. Eq. (12)."""
    tau = t / t_ms(mass, z)
    hook_time = t_hook(mass, z)
    tau1 = _tau1(t, hook_time)
    tau2 = _tau2(t, hook_time)

    l_zams = zams_luminosity(mass, z)
    alpha_l = _alpha_l(mass, z)
    beta_l = _beta_l(mass, z)
    eta = _eta_exponent(mass, z)
    delta_l = _delta_l(mass, z)

    log_ratio = np.log10(l_tms(mass, z) / l_zams)
    log_l = (
        alpha_l * tau
        + beta_l * tau**eta
        + (log_ratio - alpha_l - beta_l) * tau**2
        - delta_l * (tau1**2 - tau2**2)
    )
    return l_zams * 10.0**log_l


def ms_radius(mass: float, z: float, t: float) -> float:
    """R_MS(mass, z, t) in Rsun. Eq. (13), floored per eq. (24)."""
    tau = t / t_ms(mass, z)
    hook_time = t_hook(mass, z)
    tau1 = _tau1(t, hook_time)
    tau2 = _tau2(t, hook_time)

    r_zams = zams_radius(mass, z)
    alpha_r = _alpha_r(mass, z)
    beta_r = _beta_r(mass, z)
    gamma = _gamma_coefficient(mass, z)
    delta_r = _delta_r(mass, z)

    log_ratio = np.log10(r_tms(mass, z) / r_zams)
    log_r = (
        alpha_r * tau
        + beta_r * tau**10
        + gamma * tau**40
        + (log_ratio - alpha_r - beta_r - gamma) * tau**3
        - delta_r * (tau1**3 - tau2**3)
    )
    r = r_zams * 10.0**log_r

    # Eq. (24): low-mass degenerate floor (Pols et al. 1998, via Tout
    # et al. 1997). X = initial hydrogen abundance.
    x_hydrogen = 0.76 - 3.0 * z
    r_floor = 0.0258 * (1.0 + x_hydrogen) ** (5.0 / 3.0) * mass ** (-1.0 / 3.0)
    return max(r, r_floor)


def m_fgb(z: float) -> float:
    """Maximum initial mass for which helium ignites on the first
    giant branch (Msun). Eq. (3). Above this mass, a star has no GB
    phase at all -- it goes HG -> CHeB directly, and its HG endpoint
    values (L_EHG/R_EHG) come from helium-ignition quantities (L_HeI/
    R_HeI) that this module does not implement -- see `l_ehg`/`r_ehg`.
    """
    return 13.048 * (z / 0.02) ** 0.06 / (1.0 + 0.0012 * (0.02 / z) ** 1.27)


def m_hef(z: float) -> float:
    """Maximum initial mass for which helium ignites degenerately in a
    flash at the tip of the GB (Msun). Eq. (2). Below this mass, a
    star develops a degenerate core on the GB and its core mass at the
    end of HG comes from a different relation (M_c,GB(L_BGB)) than the
    intermediate-mass branch this module implements -- see
    `core_mass_ehg`.
    """
    zeta = _zeta(z)
    return 1.995 + 0.25 * zeta + 0.087 * zeta**2


def l_ehg(mass: float, z: float) -> float:
    """Luminosity at the end of the HG (Lsun) -- the start of the GB.

    L_EHG = L_BGB for M < M_FGB (eq. 10, via giant_branch.l_bgb).
    Raises ValueError for M >= M_FGB: those stars skip the GB entirely
    and need L_HeI (helium-ignition luminosity, part of the CHeB
    machinery this module does not implement) -- see module docstring
    and docs/science/rlof-ce-classifier-proposal.md.
    """
    if mass >= m_fgb(z):
        raise ValueError(
            f"mass={mass} Msun >= M_FGB({z})={m_fgb(z):.3f} Msun: this star "
            "has no GB phase (HG -> CHeB directly) and needs L_HeI, which "
            "is not implemented by this module."
        )
    return giant_branch.l_bgb(mass, z)


def r_ehg(mass: float, z: float) -> float:
    """Radius at the end of the HG (Rsun) -- the start of the GB.

    R_EHG = R_GB(M, L_BGB) for M < M_FGB (eq. 46, via
    giant_branch.r_gb). Raises ValueError for M >= M_FGB -- see
    `l_ehg`'s docstring; the same R_HeI gap applies to the radius.
    """
    if mass >= m_fgb(z):
        raise ValueError(
            f"mass={mass} Msun >= M_FGB({z})={m_fgb(z):.3f} Msun: this star "
            "has no GB phase (HG -> CHeB directly) and needs R_HeI, which "
            "is not implemented by this module."
        )
    return giant_branch.r_gb(mass, giant_branch.l_bgb(mass, z), z)


def core_mass_ehg(mass: float, z: float) -> float:
    """Core mass at the end of the HG (Msun) -- the start of the GB.
    Eq. (28), intermediate-mass branch only (M_HeF <= M < M_FGB) --
    the branch relevant to Realta's HG-donor mass range.

    M_c,EHG = M_c,BGB (giant_branch.core_mass_bgb -- itself an
    asymptotic approximation, see that function's docstring) for
    M_HeF <= M < M_FGB.

    Raises ValueError outside that bracket: M < M_HeF needs
    M_c,GB(L_BGB) (a degenerate-core low-mass relation this module
    does not implement -- rare for Realta's mcut=8 Msun default
    anyway), M >= M_FGB needs M_c,HeI (same CHeB gap as `l_ehg`/
    `r_ehg`).
    """
    if mass < m_hef(z):
        raise ValueError(
            f"mass={mass} Msun < M_HeF({z})={m_hef(z):.3f} Msun: this "
            "low-mass star develops a degenerate GB core and needs "
            "M_c,GB(L_BGB), which is not implemented by this module."
        )
    if mass >= m_fgb(z):
        raise ValueError(
            f"mass={mass} Msun >= M_FGB({z})={m_fgb(z):.3f} Msun: this star "
            "has no GB phase (HG -> CHeB directly) and needs M_c,HeI, which "
            "is not implemented by this module."
        )
    return giant_branch.core_mass_bgb(mass)


def _rho_coefficient(mass: float) -> float:
    """rho, the fraction of M_c,EHG the core mass starts at when HG
    begins. Eq. (29).
    """
    return (1.586 + mass**5.25) / (2.434 + 1.02 * mass**5.25)


def core_mass_hg(mass: float, z: float, t: float) -> float:
    """Core mass during the HG (Msun), for t_MS <= t < t_BGB. Eq. (30).

    Grows linearly in tau from rho*M_c,EHG at t_MS to M_c,EHG at
    t_BGB. Does not implement the mass-loss caveat in Hurley et al.'s
    text (taking the max of this value and the previous time-step's
    core mass) since this module does not model stellar-wind mass
    loss during the HG at all.
    """
    tau = (t - t_ms(mass, z)) / (t_bgb(mass, z) - t_ms(mass, z))
    m_c_ehg = core_mass_ehg(mass, z)
    rho = _rho_coefficient(mass)
    return ((1.0 - tau) * rho + tau) * m_c_ehg


def hg_luminosity(mass: float, z: float, t: float) -> float:
    """L_HG(mass, z, t) in Lsun, for t_MS <= t < t_BGB. Eq. (26)."""
    tau = (t - t_ms(mass, z)) / (t_bgb(mass, z) - t_ms(mass, z))
    l_tms_val = l_tms(mass, z)
    return l_tms_val * (l_ehg(mass, z) / l_tms_val) ** tau


def hg_radius(mass: float, z: float, t: float) -> float:
    """R_HG(mass, z, t) in Rsun, for t_MS <= t < t_BGB. Eq. (27)."""
    tau = (t - t_ms(mass, z)) / (t_bgb(mass, z) - t_ms(mass, z))
    r_tms_val = r_tms(mass, z)
    return r_tms_val * (r_ehg(mass, z) / r_tms_val) ** tau


def phase(mass: float, z: float, t: float) -> int:
    """Hurley stellar type k: 0 (deeply/fully convective MS, M<=0.7),
    1 (radiative-core MS, M>0.7), or 2 (Hertzsprung Gap).

    Identifies the evolutionary phase only -- it does NOT check
    M_FGB, so a mass >= M_FGB HG star still reports k=2 here even
    though `hg_luminosity`/`hg_radius` will raise for it (they need
    L_HeI/R_HeI, not implemented -- see those functions' docstrings).
    Callers that need an actual radius/luminosity, not just the phase
    label, must be prepared for that second raise.

    Raises ValueError for t >= t_BGB(mass, z) -- GB and later phases
    are out of scope for this module (see module docstring). Callers
    doing RLOF classification must catch this and treat it as "not
    modelled by this module" rather than guessing a phase.
    """
    lifetime = t_ms(mass, z)
    if t < lifetime:
        return 0 if mass <= 0.7 else 1
    if t < t_bgb(mass, z):
        return 2
    raise ValueError(
        f"t={t} Myr >= t_BGB({mass} Msun, Z={z})={t_bgb(mass, z)} Myr: "
        "star has reached the giant branch. GB and later phases are "
        "not implemented by this module -- see its docstring."
    )
