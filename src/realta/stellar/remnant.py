"""Compact-remnant radius -- used as a stellar core's radius (R_c1) in
the common-envelope binding-energy formula (Hurley, Tout & Pols 2002,
eq. 69).

Hurley, Pols & Tout (2000, MNRAS 315, 543), Section 6.2.1 (white
dwarfs, eq. 91) and Section 6.3 ("Small-envelope behaviour and hot
subdwarfs"): the core radius `R_c` of a star still on the HG or GB is
defined as the radius that star's core would have if it lost its
envelope immediately -- `R_c = R_ZHe(M_c)` (naked-helium-star radius,
eq. 78) for `M < M_HeF`, `R_c = R_WD(M_c)` (white-dwarf mass-radius
relation, eq. 91) otherwise.

Sanity-checked before implementation, not just transcribed: this is
essentially the standard non-relativistic Chandrasekhar mass-radius
relation. `white_dwarf_radius(0.6)` gives ~0.0125 Rsun (~8700 km),
matching the well-known real value for a 0.6 Msun white dwarf (Sirius
B territory) -- this, combined with the formula being clean text (not
a dense appendix table), is why this one wasn't put through another
paste-verification round the way the a/b coefficient tables were.
"""

from __future__ import annotations

from realta.stellar import main_sequence

# Chandrasekhar mass, fixed for normal (non-cataclysmic-variable)
# composition -- mu_e ~= 2, so M_ch = 5.8/mu_e^2 ~= 1.44 Msun (HTP02
# Sec. 6.2.1). "so it is composition-dependent... but for low-mass
# stars in cataclysmic variables mu_e ~= 2, except... so we use
# M_ch = 1.44 at all times."
M_CHANDRASEKHAR = 1.44

# Neutron-star radius (Rsun), used only as a floor in eq. 91 -- HTP02's
# own stated value, "simply set to 10 km."
R_NEUTRON_STAR = 1.4e-5


def white_dwarf_radius(mass: float) -> float:
    """R_WD(mass) in Rsun. Eq. (91).

    R_WD = max[R_NS, 0.0115*sqrt((M_ch/M)^(2/3) - (M/M_ch)^(2/3))].
    """
    scaled = (M_CHANDRASEKHAR / mass) ** (2.0 / 3.0) - (mass / M_CHANDRASEKHAR) ** (
        2.0 / 3.0
    )
    return max(R_NEUTRON_STAR, 0.0115 * scaled**0.5)


def core_radius(core_mass: float, donor_mass: float, z: float) -> float:
    """Core radius R_c1 (Rsun) for a donor still on the HG or GB, for
    use in the CE binding-energy formula.

    R_c = R_WD(core_mass) for donor_mass >= M_HeF -- the only regime
    reachable by this codebase's current core-mass tracking scope (see
    main_sequence.py::core_mass_ehg, which is itself restricted to
    M_HeF <= M < M_FGB). Raises ValueError for donor_mass < M_HeF:
    that regime uses R_ZHe(core_mass) (naked-helium-star ZAMS radius,
    eq. 78), which is not implemented here.
    """
    m_hef = main_sequence.m_hef(z)
    if donor_mass < m_hef:
        raise ValueError(
            f"donor_mass={donor_mass} Msun < M_HeF({z})={m_hef:.3f} Msun: "
            "core radius for this regime uses R_ZHe(core_mass) "
            "(naked-helium-star radius, eq. 78), which is not "
            "implemented here."
        )
    return white_dwarf_radius(core_mass)
