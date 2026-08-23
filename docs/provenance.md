# Provenance: Power et al. (2009) → reference Fortran → Realta

This is the traceability chain the development brief asks for (§19):
**paper assumption/equation → reference Fortran implementation → Python
implementation → test coverage.** It exists so that anyone touching the
physics can check whether a change preserves the Power et al. (2009)
baseline (MNRAS 395, 1146) without re-deriving it from scratch, and so
gaps in test coverage are visible in one place rather than scattered
across module docstrings.

Reference Fortran source: `gc_hmxbs/*.f` (not part of this repository;
the ~20-year-old Monte Carlo model Realta is ported from). File names
below refer to that source tree.

Status column: `done` = implemented and matches the reference as far as
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

| Paper / assumption | Reference Fortran | Realta implementation | Status |
|---|---|---|---|
| IMF sampling, all three forms (Salpeter/Kroupa/Chabrier) | `make_stars.f` (driver), `salpeter.f`, `kroupa.f`, `log_normal_IMF.f` | `imf/salpeter.py`, `imf/kroupa.py`, `imf/chabrier.py`, `imf/factory.py::get_imf` | done, **unit-tested** — `tests/test_imf.py` now covers all three (CDF boundary + monotonicity + sampling bounds), not just Kroupa. Kroupa additionally **pinned (integration)**: `n_massive`/`total_mass_msun` are direct output of Kroupa sampling for `REGRESSION_CONFIG`. Kroupa remains **flagged** — see `imf/kroupa.py` class docstring for the `mmin`-handling discrepancy (only manifests if `config.mmin != 0.01`, not exercised by the current pin). |
| Sec. 2.1: all massive stars (M ≥ `mcut`) are in binaries at formation (`fpbin=1.0`) | `make_stars.f` | `binaries/population.py::generate_population` — every `m1 >= mcut` star unconditionally gets `m2`/`period`/`a` | pinned (integration) — `n_massive` (225) is a direct pinned value |
| Companion mass distribution, flat between `mcmpct` and `m1` | `make_stars.f` | `binaries/population.py::generate_population` (`m2 = cfg.mcomp + (m1 - cfg.mcomp) * rng.random()`, clipped to `m1`) | pinned (integration) — `m2` feeds `floss`/`mtot` in `evolve()` Phase 1, so a change here would move the pinned trajectory, but `m2` itself is not directly asserted |
| Orbital period, log-flat between `pmin`/`pmax` | `make_stars.f` | `binaries/population.py::generate_population` | **unit-tested** — `tests/test_population_generation.py::test_orbital_period_bounds_and_log_uniform_distribution` checks exact bounds plus the log-uniform shape statistically (fixed seed). Still outside the integration pin's reach (`period` never feeds `floss`), so this test is the *only* coverage. |
| Semi-major axis from Kepler's third law | `make_stars.f` | `binaries/population.py::generate_population` (`AFAC` constant) | **unit-tested** — `test_semi_major_axis_matches_reference_formula` recomputes `a` independently from `m1`/`m2`/`period` via the reference's exact formula and asserts an exact match (`rtol=1e-12`), and pins `AFAC == 0.0193852859` against the literal value in `main.f`. Note: `AFAC` is the reference's own hardcoded constant, not independently re-derived from Kepler's third law -- a from-scratch physical derivation gives a ~1% different prefactor, which is the reference's own precision/convention, not a bug (see the test's docstring). |
| MS lifetime vs. mass, tabulated (Schaerer et al. 1993-derived) | `lifetime.f`, reads `lifetimes_z*.dat` | `io/tables.py::LifetimeTable`, same bundled `.dat` files | pinned (integration) — `turnoff_time` gates the SN1 trigger, so the `ndead`/`nactive` trajectory is sensitive to it |

## 2. Supernova transitions and HMXB activation

| Paper / assumption | Reference Fortran | Realta implementation | Status |
|---|---|---|---|
| Sec. 2.1: sudden-mass-loss survival criterion (binary stays bound iff `floss <= 0.5`) | `main.f` | `binaries/population.py::evolve`, Phase 1 | **unit-tested (phase)** — `test_evolve_phase1_survival_criterion` hand-constructs two systems straddling the 0.5 threshold and asserts the exact survive/disrupt split. Also pinned (integration) in both `fsur` cases; I verified sensitivity directly (perturbed the `0.5` threshold, confirmed the pinned test fails, reverted). |
| Sec. 2.1: `f_sur` — probability a surviving binary becomes an active HMXB | `main.f` (`ran3(iseed).le.fbin`) | `config.py::SimulationConfig.fsur`, applied in `binaries/population.py::evolve` Phase 1. **Renamed** from the reference's confusingly-named `fbin` — see `config.py` field docstring for the full naming history/rationale (this session's own fix, not a paper discrepancy). | **unit-tested (phase)** — `test_evolve_phase1_fsur_partial_activation` (statistical, fixed seed, ~10σ tolerance) directly exercises the `fsur<1` rejection branch and confirms it both activates and rejects a real fraction; I verified sensitivity by hardcoding the gate to always-true and confirming the test fails, then reverted. Also pinned (integration) at both `fsur=1.0` and `fsur=0.5` in `tests/test_regression.py`. |
| Remnant mass vs. progenitor mass, tabulated | `get_mremnant.f`, reads remnant mass table | `io/tables.py::RemnantTable`, `data/remnant_masses.dat` | pinned (integration) — feeds `floss` directly, for the mass range sampled in the regression configs. Its *value* for a given progenitor mass is not independently tested (`test_evolve_phase1_survival_criterion` stubs it out entirely to isolate the survival-criterion logic from the table lookup). |
| Secondary SN → system marked dead | `main.f` | `binaries/population.py::evolve`, Phase 2 | **unit-tested (phase)** — `test_evolve_phase2_death` hand-constructs a mid-lifetime active HMXB and asserts the exact death transition (remnant mass applied, `turnoff_time`/`is_survived`/`lum_xray` cleared, `nturn=2`). Also pinned (integration) via `ndead`. |

## 3. X-ray luminosity and ionising photons

| Paper / assumption | Reference Fortran | Realta implementation | Status |
|---|---|---|---|
| Per-HMXB X-ray luminosity draw, Weibull-like/peaked, Eddington-limited, `L_X^peak ≈ 1e38 erg/s` | `get_lumx.f` | `xray/luminosity.py::XRayLuminosity` (`distribution="weibull"`, Realta's default) | pinned (integration) — `REGRESSION_CONFIG` uses the weibull distribution, and I confirmed sensitivity directly (doubling the `lumx_tot` accumulation broke the pin, then reverted). This was the site of this session's RNG-seeding bug fix (see `xray/luminosity.py` module docstring). |
| `get_lumx.f`'s dead "uniform" branch (only reachable if `iseed == -1`, a debug sentinel never hit by `main.f`) | `get_lumx.f` | `xray/luminosity.py::XRayLuminosity` (`distribution="uniform"`), `config.py::xray_distribution` | **unit-tested** — `tests/test_xray_luminosity.py` covers bounds, log-uniform shape, and a real asymmetry this pass surfaced: the "uniform" branch is **not** Eddington-limited (no `lumx <= ledd` rejection), unlike "weibull" — `test_uniform_distribution_is_not_eddington_limited_unlike_weibull` documents and locks in this behavior. Still **not the reference default**. |
| X-ray luminosity persists once assigned, until SN2 death | (implicit in `main.f`'s single-draw-per-lifetime structure) | `binaries/population.py::evolve` — L_X drawn once at SN1 in Phase 1, zeroed at SN2 death, **not** redrawn per timestep | **unit-tested (phase)** — `test_evolve_phase1_survival_criterion` asserts `lum_xray` is set (once, at activation) and stays nonzero only for the activated survivor; `test_evolve_phase2_death` asserts it is zeroed exactly at SN2 death, not before. This was one of this session's four originally-diagnosed porting bugs. |
| Ionising photon rate from X-ray luminosity (spectral-shape conversion, 13.6 eV–1500 eV vs. 1e6 eV reference) | `main.f` (`lum_xray(i)*(6.2415e11/13.6)*alog(1500/13.6)/alog(1e6/13.6)`) | `binaries/population.py::NPHOT_PER_LUMX` class constant | pinned (integration) — `nphot_tot` is pinned directly, and is a fixed multiple of `lumx_tot`, so this constant's value is confirmed unchanged at every pinned timestep |
| `get_ngamma.f`'s MS-based ionising photon estimate | `get_ngamma.f` | **not ported** — dead in the reference too as far as this port is concerned; do not confuse with `IonizingPhotonTable`, which is a different (M≥8 Msun ionising-budget) quantity, also not currently wired into `evolve()` | not applicable |

## 4. Fig. 1 reproduction (this session's own addition, not in the original Fortran)

| Quantity | Source | Realta implementation | Status |
|---|---|---|---|
| Main-sequence bolometric luminosity vs. time | Not in `gc_hmxbs` — the original paper's Fig. 1 MS curve was almost certainly generated with a population-synthesis code of the day (Starburst99 being the standard circa 2009), external to the ported Fortran | `io/tables.py::MSLuminosityTable`, sourced from FSPS (`notebooks_helper/generate_ms_luminosity_table.py`) | done, untested, **and flagged**: FSPS vs. Starburst99 track differences (Padova vs. Geneva stellar evolution) are a documented literature-level systematic of up to a factor of a few in this exact age window (3–20 Myr) — see this session's SB99-vs-FSPS research summary. Not a bug, but worth stating explicitly wherever this reproduction is presented. |
| MS/HMXB luminosity mass-normalization | n/a (this session's own bug, introduced when the FSPS table was added, not present in the reference) | `binaries/population.py::BinaryPopulation.total_mass_msun`, `io/tables.py::MSLuminosityTable.get_lbol(age_myr, total_mass_msun)` | `total_mass_msun` itself is **pinned (integration)** (15836.446390529192 Msun in both regression cases). `get_lbol()`'s rescaling formula is now **unit-tested** — `tests/test_ms_luminosity_table.py` checks exact linear scaling with mass (2x mass -> 2x luminosity, at several ages), zero-mass/zero-luminosity, the no-extrapolation domain boundary, and a value cross-check against `total_mass_msun`'s own pinned figure. All three metallicity tables (`imetal=1,2,3`) are loaded and checked, not just the default. |

## 5. Random number generation

| Paper / assumption | Reference Fortran | Realta implementation | Status |
|---|---|---|---|
| Single seeded stream drives the whole simulation (`ran3`, Numerical Recipes) | `ran3.f`, seeded from `iseed` in `main.f` | `binaries/population.py::BinaryPopulation.np_rng` (`np.random.default_rng(config.iseed)`), threaded through population generation, SN-survival, HMXB activation, **and** `XRayLuminosity.get_lumx(rng=...)` | **tested** — `tests/test_regression.py::test_run_is_deterministic_for_fixed_seed` directly asserts two runs of the same config produce bit-identical output. This was this session's other RNG bug fix (X-ray draws previously used the unseeded global `np.random`, breaking reproducibility). |

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

Still open:

- **The Kroupa `mmin`-handling discrepancy** (§1 above) is flagged, not
  resolved — needs a decision on which behaviour (Fortran's hardcoded
  0.01, or this port's caller-supplied `mmin`) is actually intended
  before `mmin` is treated as a free parameter. Not exercised by the
  current pin (`mmin=0.01` throughout).
- **No `examples/` directory** (brief §26) walking through a full
  run end-to-end outside of the paper-reproduction notebook.
