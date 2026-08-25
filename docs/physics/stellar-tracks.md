# Stellar radius/luminosity tracks

Analytic stellar tracks, used internally by the RLOF/CE classifier
([`rlof-classifier.md`](rlof-classifier.md)) to compute donor
radius/luminosity vs. age. **Used internally only** — Realta's
reported `L_bol`/`L_UV` come from `MSLuminosityTable`/
`UVLuminosityTable` ([`observables.md`](observables.md)), not from
this module.

## Scope

Stellar types k=0 (fully convective MS, M ≤ 0.7 Msun), k=1
(radiative-core MS, M > 0.7 Msun), and k=2 (Hertzsprung Gap, M <
`M_FGB` for radius/luminosity; additionally `M_HeF <= M <
CORE_MASS_BGB_MAX_MASS` (~7.3 Msun) for core mass/radius). The true
giant branch (k=3) and later phases are out of scope:
`phase()` raises `ValueError` for `t >= t_BGB`; `l_ehg`/`r_ehg`/
`hg_radius`/`hg_luminosity`/`core_mass_ehg` raise for HG stars with
`M >= M_FGB` (no GB phase modelled at those masses); `core_mass_ehg`
also raises for `M < M_HeF` (needs a different, degenerate-core
relation); `core_mass_bgb`/`core_mass_ehg`/`core_mass_hg`/
`core_radius` raise for `M >= CORE_MASS_BGB_MAX_MASS`. Callers must
treat a raise as "not modelled", not a bug.

## Zero-age main sequence

`L_ZAMS(M,Z)`, `R_ZAMS(M,Z)` — Tout, Pols, Eggleton & Han (1996, MNRAS
281, 257), eqs. 1-4, Tables 1-2.

Implementation: `stellar/zams.py::zams_luminosity`, `zams_radius`.
Tests: `tests/test_hurley_main_sequence.py`.

## Main sequence (k=0,1)

`t_MS`, `L_MS(t)`, `R_MS(t)` — Hurley, Pols & Tout (2000, MNRAS 315,
543), Sec. 5.1 (eqs. 1-24), Appendix A coefficients a1-a81.

Implementation: `stellar/main_sequence.py`.
Tests: `tests/test_hurley_main_sequence.py` — solar-calibration bounds,
MS-lifetime bounds, monotonic radius/luminosity growth, exact
convergence to `R_TMS`/`L_TMS` at `t -> t_MS`, the `phase()` scope
guard, the low-mass degenerate-radius floor (eq. 24).

## Giant-branch base (endpoint values only)

`L_BGB(M,Z)`, `R_GB(M,L,Z)`, mass-radius exponent `x` — Hurley et al.
(2000), eq. 10 (a27-a32) and eqs. 46-48 (b1-b7). Used as the HG's
endpoint boundary values, not for time-evolution along the GB itself.

Implementation: `stellar/giant_branch.py::l_bgb`, `r_gb`,
`mass_radius_exponent`.
Tests: `tests/test_giant_branch.py` — solar-calibration bound for
`L_BGB`, monotonic growth with mass, cross-validation of `r_gb()`
against Hurley et al.'s own illustrative Z=0.02 formula (`R_GB ≈
1.1*M^-0.3*(L^0.4 + 0.383*L^0.76)`, Sec. 5.2), agreeing to ~10-20%.

## Hertzsprung Gap (k=2, M < M_FGB)

`L_HG(t)`, `R_HG(t)` — Hurley et al. (2000), Sec. 5.1.2, eqs. 25-30.

Implementation: `stellar/main_sequence.py::hg_luminosity`, `hg_radius`,
`l_ehg`, `r_ehg`, `m_fgb`, `phase()`.
Tests: `tests/test_hertzsprung_gap.py` — monotonic radius growth,
exact convergence to `R_TMS`/`R_EHG` (and `L_TMS`/`L_EHG`) at both
endpoints, `phase()` reporting k=2 during HG, `M_FGB` sanity bound
(10-16 Msun at solar Z).

## HG core mass (M_HeF ≤ M < CORE_MASS_BGB_MAX_MASS)

`M_c,HG(t)` — Hurley et al. (2000), eqs. 2 (`M_HeF`), 28-30
(core-mass growth). `M_c,BGB` uses the paper's stated large-mass
asymptotic limit `0.098*M^1.35` (`= c1^0.25 * M^(c2/4)`, `c1 =
9.20925e-5`, `c2 = 5.402216` from eq. 44), not the full eq. 44
core-mass-luminosity relation. This exceeds the Chandrasekhar mass
above `M = (M_ch^4/c1)^(1/c2) ≈ 7.317` Msun — `core_mass_bgb` raises
above that mass (`CORE_MASS_BGB_MAX_MASS`) rather than returning an
unphysical value.

Implementation: `stellar/main_sequence.py::m_hef`, `core_mass_ehg`,
`core_mass_hg`, `_rho_coefficient`; `stellar/giant_branch.py::
core_mass_bgb`, `CORE_MASS_BGB_MAX_MASS`.
Tests: `tests/test_hg_core_mass.py`.

## Core radius

`R_c1` — Hurley et al. (2000), Sec. 6.2.1 (white dwarfs, eq. 91) and
Sec. 6.3 (`R_c = R_WD(M_c)` for `M >= M_HeF`; `R_c = R_ZHe(M_c)`
otherwise, not implemented — out of scope).

Implementation: `stellar/remnant.py::white_dwarf_radius`,
`core_radius`, `M_CHANDRASEKHAR`, `R_NEUTRON_STAR`.
Tests: `tests/test_remnant.py` — 0.6 Msun white dwarf real-value
cross-check (~0.0115-0.013 Rsun), mass-radius direction, shrinkage
toward the Chandrasekhar mass, neutron-star-radius floor, end-to-end
`core_mass_hg -> core_radius` chain.
