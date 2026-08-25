"""Zero-age main-sequence luminosity and radius.

Tout, Pols, Eggleton & Han (1996, MNRAS 281, 257), the paper Hurley,
Pols & Tout (2000) itself sources L_ZAMS/R_ZAMS from without reprinting
the coefficients (see that paper's Section 5, opening paragraph).
Coefficients transcribed directly from Tables 1 and 2 of the 1996
paper (page 258), supplied as a clean-text PDF -- high transcription
confidence, unlike some of the 2000 paper's later, denser appendix
tables (see docs/physics/stellar-tracks.md).

Equations (page numbers/eq numbers refer to Tout et al. 1996):

    L(M, Z) = (alpha*M^5.5 + beta*M^11)
              / (gamma + M^3 + delta*M^5 + epsilon*M^7
                 + zeta*M^8 + eta*M^9.5)                        (eq. 1)

    R(M, Z) = (theta*M^2.5 + iota*M^6.5 + kappa*M^11
               + lambda_*M^19 + mu*M^19.5)
              / (nu + xi*M^2 + omicron*M^8.5 + M^18.5
                 + pi_*M^19.5)                                  (eq. 2)

Each coefficient c is itself a quartic in log10(Z/0.02) (eq. 3 for the
luminosity coefficients, eq. 4 for radius -- same functional form):

    c(Z) = a + b*zeta + c2*zeta**2 + d*zeta**3 + e*zeta**4,
    zeta = log10(Z / 0.02)

except nu (radius), which the paper fixes at its Z=0.02 value for all
metallicities (Section 3: "the coefficient nu was fixed at its optimum
value at Z = 0.02 ... To avoid it[large coefficient swings] ... the
coefficient nu was fixed").

Valid for M in [0.1, 100] Msun, Z in [0.0001, 0.03] (paper's own
stated domain); this module does not enforce those bounds -- callers
(io/tables or the SSE main-sequence module) are expected to keep
callers within Realta's own IMF mass range, which is already inside
this domain for realistic configs.
"""

from __future__ import annotations

import numpy as np

# Table 1 (Tout et al. 1996, p. 258): luminosity coefficients alpha..eta,
# each row (a, b, c, d, e) for eq. 3.
_L_COEFFICIENTS = {
    "alpha": (0.39704170, -0.32913574, 0.34776688, 0.37470851, 0.09011915),
    "beta": (8.52762600, -24.41225973, 56.43597107, 37.06152575, 5.45624060),
    "gamma": (0.00025546, -0.00123461, -0.00023246, 0.00045519, 0.00016176),
    "delta": (5.43288900, -8.62157806, 13.44202049, 14.51584135, 3.39793084),
    "epsilon": (5.56357900, -10.32345224, 19.44322980, 18.97361347, 4.16903097),
    "zeta": (0.78866060, -2.90870942, 6.54713531, 4.05606657, 0.53287322),
    "eta": (0.00586685, -0.01704237, 0.03872348, 0.02570041, 0.00383376),
}

# Table 2 (Tout et al. 1996, p. 258): radius coefficients theta..pi,
# each row (a, b, c, d, e) for eq. 4. nu's row is (a, 0, 0, 0, 0) --
# fixed at its Z=0.02 value, per the paper's own Section 3.
_R_COEFFICIENTS = {
    "theta": (1.71535900, 0.62246212, -0.92557761, -1.16996966, -0.30631491),
    "iota": (6.59778800, -0.42450044, -12.13339427, -10.73509484, -2.51487077),
    "kappa": (10.08855000, -7.11727086, -31.67119479, -24.24848322, -5.33608972),
    "lambda_": (1.01249500, 0.32699690, -0.00923418, -0.03876858, -0.00412750),
    "mu": (0.07490166, 0.02410413, 0.07233664, 0.03040467, 0.00197741),
    "nu": (0.01077422, 0.0, 0.0, 0.0, 0.0),
    "xi": (3.08223400, 0.94472050, -2.15200882, -2.49219496, -0.63848738),
    "omicron": (17.84778000, -7.45345690, -48.96066856, -40.05386135, -9.09331816),
    "pi_": (0.00022582, -0.00186899, 0.00388783, 0.00142402, -0.00007671),
}


def _zeta(z: float) -> float:
    return np.log10(z / 0.02)


def _coefficient(row: tuple[float, float, float, float, float], zeta: float) -> float:
    a, b, c, d, e = row
    return a + b * zeta + c * zeta**2 + d * zeta**3 + e * zeta**4


def zams_luminosity(mass: float, z: float) -> float:
    """L_ZAMS(mass, z) in Lsun. Tout et al. (1996), eqs. 1 and 3."""
    zeta = _zeta(z)
    c = {name: _coefficient(row, zeta) for name, row in _L_COEFFICIENTS.items()}
    m = mass
    numerator = c["alpha"] * m**5.5 + c["beta"] * m**11
    denominator = (
        c["gamma"]
        + m**3
        + c["delta"] * m**5
        + c["epsilon"] * m**7
        + c["zeta"] * m**8
        + c["eta"] * m**9.5
    )
    return numerator / denominator


def zams_radius(mass: float, z: float) -> float:
    """R_ZAMS(mass, z) in Rsun. Tout et al. (1996), eqs. 2 and 4."""
    zeta = _zeta(z)
    c = {name: _coefficient(row, zeta) for name, row in _R_COEFFICIENTS.items()}
    m = mass
    numerator = (
        c["theta"] * m**2.5
        + c["iota"] * m**6.5
        + c["kappa"] * m**11
        + c["lambda_"] * m**19
        + c["mu"] * m**19.5
    )
    denominator = (
        c["nu"] + c["xi"] * m**2 + c["omicron"] * m**8.5 + m**18.5 + c["pi_"] * m**19.5
    )
    return numerator / denominator
