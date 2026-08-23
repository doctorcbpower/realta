# Provenance: Power et al. → Realta

This is the traceability chain the development brief asks for (§19):
**paper assumption/equation → Python implementation → test coverage.**
It exists so that anyone touching the physics can check whether a
change preserves the published baseline without re-deriving it from
scratch, and so gaps in test coverage are visible in one place rather
than scattered across module docstrings.

## References

- **Power, C., Wynn, G. A., Combet, C., Wilkinson, M. I. (2009)**,
  "Primordial globular clusters, X-ray binaries and cosmological
  reionization", *MNRAS*, 395, 2, 1146-1152.
  [arXiv:0902.1897](https://arxiv.org/abs/0902.1897) ·
  [DOI: 10.1111/j.1365-2966.2009.14628.x](https://doi.org/10.1111/j.1365-2966.2009.14628.x) ·
  [ADS](https://ui.adsabs.harvard.edu/abs/2009MNRAS.395.1146P/abstract)

  The primary reference for Realta's population model: coeval star
  formation, the IMF choices, the sudden-mass-loss survival criterion,
  `f_sur`, and the peaked/Eddington-limited X-ray luminosity
  distribution all originate here (mostly Sec. 2.1-2.2). Realta is a
  Python port of the ~20-year-old Fortran Monte Carlo code this paper's
  results were produced with.

- **Power, C., James, G., Combet, C., Wynn, G. (2013)**, "Feedback from
  High-mass X-ray Binaries on the High-redshift Intergalactic Medium:
  Model Spectra", *ApJ*, 764, 1, 76.
  [arXiv:1211.5854](https://arxiv.org/abs/1211.5854) ·
  [DOI: 10.1088/0004-637X/764/1/76](https://doi.org/10.1088/0004-637X/764/1/76) ·
  [ADS](https://ui.adsabs.harvard.edu/abs/2013ApJ...764...76P/abstract)

  Direct follow-up to the 2009 paper, using the same Monte Carlo
  population model, that derives model HMXB spectra (Cygnus X-1-like,
  black-body + power-law accretion states) and ionizing-photon output.
  Relevant to Realta's X-ray-to-ionizing-photon conversion (Section 3
  below) -- that conversion's physical basis lives here, not in the
  2009 paper.

Status column: `done` = implemented and matches the paper as far as
verified; `flagged` = a known discrepancy or open question, not yet
resolved (see the linked docstring for detail); `untested` = correct as
far as manual verification goes, but has no automated test pinning its
numeric output; `pinned (integration)` = exercised and value-checked by
`tests/test_regression.py::test_reference_cluster_run_pinned_trajectory`,
which pins the *combined* end-to-end trajectory (`lumx_tot`/`nphot_tot`/
`nactive`/`ndead` vs. time, plus `total_mass_msun`/binary count) for two
fixed config+seed cases (`fsur=1.0` and `fsur=0.5`) -- a change to any
row marked this way would move the pinned numbers and fail that test,
but the row is not tested in isolation the way a unit test would;
`unit-tested (phase)` = exercised in isolation by
`tests/test_evolve.py`, which constructs a `BinaryPopulation` and
directly overwrites its internal state to set up one exact scenario per
`evolve()` phase (survival, death, activation gate), independent of
what a real IMF-sampled population would or wouldn't reach.

---

## 1. Population generation

| Power et al. (2009) reference | Realta implementation | Status |
|---|---|---|
| IMF sampling, all three forms (Salpeter, Kroupa, Chabrier -- Sec. 2.2) | `imf/salpeter.py`, `imf/kroupa.py`, `imf/chabrier.py`, `imf/factory.py::get_imf` | done, **unit-tested** — `tests/test_imf.py` covers all three (CDF boundary + monotonicity + sampling bounds). Kroupa additionally **pinned (integration)**: `n_massive`/`total_mass_msun` are direct output of Kroupa sampling in the regression configs. `config.mmin` default is 0.1 Msun (the practical stellar lower-mass cutoff) -- see `imf/kroupa.py`'s class docstring. |
| Sec. 2.1: all massive stars (M ≥ `mcut`) are in binaries at formation | `binaries/population.py::generate_population` — every `m1 >= mcut` star unconditionally gets `m2`/`period`/`a` | pinned (integration) — `n_massive` (251) is a direct pinned value |
| Companion mass distribution, flat between `mcmpct` and `m1` (Sec. 2.1) | `binaries/population.py::generate_population` (`m2 = cfg.mcomp + (m1 - cfg.mcomp) * rng.random()`, clipped to `m1`) | pinned (integration) — `m2` feeds `floss`/`mtot` in `evolve()` Phase 1, so a change here would move the pinned trajectory, but `m2` itself is not directly asserted |
| Orbital period, log-flat between `pmin`/`pmax` (Sec. 2.1) | `binaries/population.py::generate_population` | **unit-tested** — `tests/test_population_generation.py::test_orbital_period_bounds_and_log_uniform_distribution` checks exact bounds plus the log-uniform shape statistically (fixed seed). Outside the integration pin's reach (`period` never feeds `floss`), so this test is the *only* coverage. |
| Semi-major axis from Kepler's third law (Sec. 2.1) | `binaries/population.py::generate_population` (`AFAC` constant) | **unit-tested** — `test_semi_major_axis_matches_reference_formula` recomputes `a` independently from `m1`/`m2`/`period` and asserts an exact match (`rtol=1e-12`), and pins the exact value of `AFAC`. Note: `AFAC` is not independently re-derivable to full precision from Kepler's third law alone (a from-scratch physical derivation gives a ~1% different prefactor) — a units/precision-convention detail, not a bug. |
| MS lifetime vs. mass, tabulated | `io/tables.py::LifetimeTable`, `data/lifetimes_z*.dat` (data table headers cite Schaerer et al. 1993) | pinned (integration) — `turnoff_time` gates the SN1 trigger, so the `ndead`/`nactive` trajectory is sensitive to it |

## 2. Supernova transitions and HMXB activation

| Power et al. (2009) reference | Realta implementation | Status |
|---|---|---|
| Sec. 2.1: sudden-mass-loss survival criterion (binary stays bound iff `floss <= 0.5`) | `binaries/population.py::evolve`, Phase 1 | **unit-tested (phase)** — `test_evolve_phase1_survival_criterion` hand-constructs two systems straddling the 0.5 threshold and asserts the exact survive/disrupt split. Also pinned (integration) in both `fsur` cases; I verified sensitivity directly (perturbed the `0.5` threshold, confirmed the pinned test fails, reverted). |
| Sec. 2.1: `f_sur` — probability a surviving binary becomes an active HMXB ("if the binary remains bound, then it has a probability of f_sur that it will evolve into a HMXB") | `config.py::SimulationConfig.fsur`, applied in `binaries/population.py::evolve` Phase 1. See `config.py`'s field docstring for the naming history (this port's own port-fidelity fix, not a paper discrepancy). | **unit-tested (phase)** — `test_evolve_phase1_fsur_partial_activation` (statistical, fixed seed, ~10σ tolerance) directly exercises the `fsur<1` rejection branch and confirms it both activates and rejects a real fraction; I verified sensitivity by hardcoding the gate to always-true and confirming the test fails, then reverted. Also pinned (integration) at both `fsur=1.0` and `fsur=0.5` in `tests/test_regression.py`. |
| Remnant mass vs. progenitor mass, tabulated | `io/tables.py::RemnantTable`, `data/remnant_masses.dat` | pinned (integration) — feeds `floss` directly, for the mass range sampled in the regression configs. Its *value* for a given progenitor mass is not independently tested (`test_evolve_phase1_survival_criterion` stubs it out entirely to isolate the survival-criterion logic from the table lookup). |
| Secondary SN → system marked dead | `binaries/population.py::evolve`, Phase 2 | **unit-tested (phase)** — `test_evolve_phase2_death` hand-constructs a mid-lifetime active HMXB and asserts the exact death transition (remnant mass applied, `turnoff_time`/`is_survived`/`lum_xray` cleared, `nturn=2`). Also pinned (integration) via `ndead`. |

## 3. X-ray luminosity and ionising photons

| Power et al. reference | Realta implementation | Status |
|---|---|---|
| Per-HMXB X-ray luminosity draw, Weibull-like/peaked, Eddington-limited, `L_X^peak ≈ 1e38 erg/s` (2009, Sec. 2.2) | `xray/luminosity.py::XRayLuminosity` (`distribution="weibull"`, Realta's default) | pinned (integration) — the pinned configs use the weibull distribution, and I confirmed sensitivity directly (doubling the `lumx_tot` accumulation broke the pin, then reverted). This was the site of this session's RNG-seeding bug fix (see `xray/luminosity.py` module docstring). |
| Flat log-uniform X-ray luminosity draw, offered as an alternative to the peaked distribution above | `xray/luminosity.py::XRayLuminosity` (`distribution="uniform"`), `config.py::xray_distribution` | **unit-tested** — `tests/test_xray_luminosity.py` covers bounds, log-uniform shape, and a real asymmetry this pass surfaced: the "uniform" branch is **not** Eddington-limited (no rejection against `L_Edd`), unlike "weibull" — `test_uniform_distribution_is_not_eddington_limited_unlike_weibull` documents and locks in this behavior. Not Realta's default. |
| X-ray luminosity persists once assigned, until SN2 death | `binaries/population.py::evolve` — L_X drawn once at SN1 in Phase 1, zeroed at SN2 death, **not** redrawn per timestep | **unit-tested (phase)** — `test_evolve_phase1_survival_criterion` asserts `lum_xray` is set (once, at activation) and stays nonzero only for the activated survivor; `test_evolve_phase2_death` asserts it is zeroed exactly at SN2 death, not before. This was one of this session's four originally-diagnosed porting bugs. |
| Ionising photon rate from X-ray luminosity, via a model spectral shape (2013, Cygnus X-1-like black-body + power-law accretion states — this is the paper the conversion's physical basis actually comes from, not the 2009 paper) | `binaries/population.py::NPHOT_PER_LUMX` class constant | pinned (integration) — `nphot_tot` is pinned directly, and is a fixed multiple of `lumx_tot`, so this constant's value is confirmed unchanged at every pinned timestep |

## 4. Fig. 1 reproduction (this session's own addition, not in either paper's original code)

| Quantity | Source | Realta implementation | Status |
|---|---|---|---|
| Main-sequence bolometric luminosity vs. time | Not part of the ported Monte Carlo model — the 2009 paper's Fig. 1 MS curve was almost certainly generated with a population-synthesis code of the day (Starburst99 being the standard circa 2009) | `io/tables.py::MSLuminosityTable`, sourced from FSPS (`notebooks_helper/generate_ms_luminosity_table.py`) | done, untested, **and flagged**: FSPS vs. Starburst99 track differences (Padova vs. Geneva stellar evolution) are a documented literature-level systematic of up to a factor of a few in this exact age window (3–20 Myr) — see this session's SB99-vs-FSPS research summary. Not a bug, but worth stating explicitly wherever this reproduction is presented. |
| MS/HMXB luminosity mass-normalization | n/a (this session's own bug, introduced when the FSPS table was added, not present in the original model) | `binaries/population.py::BinaryPopulation.total_mass_msun`, `io/tables.py::MSLuminosityTable.get_lbol(age_myr, total_mass_msun)` | `total_mass_msun` itself is **pinned (integration)** (17891.377637950565 Msun in both regression cases). `get_lbol()`'s rescaling formula is now **unit-tested** — `tests/test_ms_luminosity_table.py` checks exact linear scaling with mass (2x mass -> 2x luminosity, at several ages), zero-mass/zero-luminosity, the no-extrapolation domain boundary, and a value cross-check against `total_mass_msun`'s own pinned figure. All three metallicity tables (`imetal=1,2,3`) are loaded and checked, not just the default. |

## 5. Random number generation

| Power et al. reference | Realta implementation | Status |
|---|---|---|
| Single seeded random stream drives the whole simulation | `binaries/population.py::BinaryPopulation.np_rng` (`np.random.default_rng(config.iseed)`), threaded through population generation, SN-survival, HMXB activation, **and** `XRayLuminosity.get_lumx(rng=...)` | **tested** — `tests/test_regression.py::test_run_is_deterministic_for_fixed_seed` directly asserts two runs of the same config produce bit-identical output. This was this session's other RNG bug fix (X-ray draws previously used the unseeded global `np.random`, breaking reproducibility). |

---

## Known gaps

Closed this session:

- ~~No numeric regression test~~ — `tests/test_regression.py` now has
  `test_reference_cluster_run_pinned_trajectory` (pins `lumx_tot`/
  `nphot_tot`/`nactive`/`ndead` at 7 timesteps, plus `total_mass_msun`/
  binary count, for a fixed config+seed) and
  `test_run_is_deterministic_for_fixed_seed` (asserts run-to-run
  reproducibility directly, rather than relying on it being true because
  the RNG-seeding row above says so). Sensitivity was hand-verified by
  deliberately perturbing the survival threshold and the `lumx_tot`
  accumulation and confirming each broke the pinned test, then
  reverting — see the two rows in Section 2/3 above with that note.
- ~~Salpeter and Chabrier IMFs are untested~~ — `tests/test_imf.py` now
  covers all three IMFs (CDF boundary, monotonicity, sampling bounds).
- ~~No unit-test layer for individual `evolve()` phases~~ —
  `tests/test_evolve.py` now hand-constructs an exact scenario for each
  phase (survival/disruption split, secondary-SN death, X-ray
  luminosity draw-once-and-persist) by overwriting `BinaryPopulation`'s
  internal state directly and stubbing `remnant_table.get_remnant_mass`,
  rather than relying on what a real IMF-sampled population happens to
  reach. Each of these was individually confirmed to actually fail when
  the corresponding logic is broken (see the "I verified..." notes on
  the relevant rows above), not just written and trusted.
- ~~`fsur < 1` not exercised by any automated test~~ — closed two ways:
  `test_evolve_phase1_fsur_partial_activation` (unit-level, statistical,
  fixed seed) isolates the activation gate itself; a second pinned
  integration case (`fsur=0.5` in `REGRESSION_CASES`) confirms the
  gate's effect survives through the full end-to-end trajectory
  alongside the existing `fsur=1.0` case.
- ~~Orbital period/semi-major axis, the `uniform` X-ray distribution,
  and `MSLuminosityTable.get_lbol()`'s rescaling formula untested~~ —
  the three gaps the pinned integration test structurally cannot reach
  (nothing about them feeds `lumx_tot`/`nphot_tot`/`nactive`/`ndead`)
  now each have dedicated unit tests: `tests/test_population_generation.py`,
  `tests/test_xray_luminosity.py`, `tests/test_ms_luminosity_table.py`.
  Sensitivity was hand-verified for the semi-major-axis and MS-rescaling
  tests specifically (broke each, confirmed failure, reverted) the same
  way as the other pinned/unit tests in this document. Surfaced one new,
  real finding along the way: the "uniform" X-ray distribution has no
  Eddington-limit rejection, unlike "weibull" — see Section 3 above.
- ~~The Kroupa `mmin`-handling question~~ — resolved by decision:
  `config.mmin` default changed from 0.01 to 0.1 Msun (the practical
  stellar lower-mass cutoff -- objects below ~0.08 Msun are substellar).
  This excludes the Kroupa IMF's shallowest (beta0=0.3) segment from
  being sampled by default. `config.yml` and `SimulationConfig`'s
  default both updated; `tests/test_regression.py`'s pinned values
  were regenerated against the new default (`total_mass_msun` moved
  from 15836.45 to 17891.38 Msun, `n_massive` from 225 to 251, for the
  same `ntot`/`iseed` -- higher `mmin` means fewer low-mass stars
  compete for the same star-count budget, so more of it lands above
  `mcut`) and reverified deterministic/sensitive the same way as every
  other pinned value in this document.

Still open:

- **No `examples/` directory** (brief §26) walking through a full
  run end-to-end outside of the paper-reproduction notebook.
