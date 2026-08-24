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
| Continuous single power-law IMF slope, `config.imf_slope` (A4, `docs/science/paper1-detailed-work-breakdown.md` -- not from Power et al., needed for Figure 4's `(alpha_IMF, f_bin)` degeneracy grid) | `imf/factory.py::get_imf`'s `slope` parameter, overriding `SalpeterIMF`'s own default `alpha=2.35`; `config.imf_slope: float \| None`, `None` reproducing the baseline exactly | **unit-tested, sensitivity-verified, one pre-existing edge case newly guarded**: `tests/test_imf.py` confirms the override reaches `SalpeterIMF.alpha`, that `None` matches the unmodified default, that it is silently ignored for Kroupa/Chabrier (no single power-law slope to override there), the expected direction (a shallower slope samples a higher mean mass), and a population-level numeric regression pin (`total_mass_msun`/`n_massive` for `imf_slope=1.8`). Deliberately Salpeter-only, not a new IMF family -- see `get_imf`'s own docstring. Found, while adding this, that `SalpeterIMF.cdf()`'s denominator (`mmax**beta - mmin**beta`, `beta = 1 - alpha`) is identically `0.0` at `alpha=1.0` -- a pre-existing singularity, unreachable before `imf_slope` made `alpha` user-configurable (the `2.35` default never hits it). Guarded in `config.py::SimulationConfig.__post_init__` (`imf_slope` must be `> 0` and `!= 1.0`) rather than fixed at the `SalpeterIMF` level (a log-form CDF for that case is not implemented) -- `tests/test_binary_prescriptions.py`/`test_imf.py` cover both the config-level rejection and the raw, unguarded `ZeroDivisionError` `get_imf()` itself would still hit if called directly with `slope=1.0`. |
| Sec. 2.1: all massive stars (M ≥ `mcut`) are in binaries at formation | `binaries/population.py::generate_population` — every `m1 >= mcut` star unconditionally gets `m2`/`period`/`a` (default `binary_fraction=1.0` -- see next row) | pinned (integration) — `n_massive` (251) is a direct pinned value |
| Independently configurable `binary_fraction`/`mass_ratio_distribution`/`period_distribution` (A1, `docs/science/paper1-detailed-work-breakdown.md` -- roadmap item 7, not from Power et al.) | `config.py::SimulationConfig.binary_fraction`/`mass_ratio_distribution`/`period_distribution`; `binaries/population.py::generate_population` | **unit-tested, sensitivity-verified, baseline confirmed untouched**: `tests/test_binary_sampling_distributions.py` covers `binary_fraction=1.0` (default) still pairing every star, `binary_fraction=0.0` giving all-placeholder (`m2=0`/`period=0`/`a=0`, still tracked through `m1` -- generalizing, not replacing, the `"single"` prescription's own array-emptying shortcut) companions, an intermediate value matching the configured fraction statistically, `mass_ratio_distribution="flat_q"` (uniform in `q=m2/m1`) vs. the default's different shape, `period_distribution="log_normal"` staying within `[pmin, pmax]` and centring on the pmin/pmax-derived `mu` (a *named, non-literature-sourced* simplification -- `mu`/`sigma` derived directly from `pmin`/`pmax`, truncated via `scipy.stats.truncnorm` so it can never sample outside the configured bounds), and a combined non-default pinned regression case. I verified sensitivity directly twice: (1) forced the new per-star companion-assignment Bernoulli draw to always fire, even at `binary_fraction=1.0` -- confirmed `tests/test_regression.py`'s pinned baseline trajectory breaks exactly as expected (the RNG stream shifts), reverted; (2) removed the RLOF-classifier precompute loop's new `self.m2[i] <= 0.0` skip guard -- confirmed a `binary_fraction<1` run crashes with `ValueError: mass_ratio must be positive, got 0.0` (a real divide-by-zero/inf-`q1` risk `find_rlof_onset`/`roche_lobe_radius` have no guard against on their own, since every existing unit test for them only ever passes a real, positive companion mass), reverted. The `has_companion` mask is also applied to the pre-SN merger channel's `did_merge` eligibility (a no-companion star's placeholder `period=0` would otherwise trivially satisfy `periods < p_merge_max_period`). |
| Companion mass distribution, flat between `mcmpct` and `m1` (Sec. 2.1) | `binaries/population.py::generate_population` (`m2 = cfg.mcomp + (m1 - cfg.mcomp) * rng.random()`, clipped to `m1`) -- `mass_ratio_distribution="uniform"` (default, see above) | pinned (integration) — `m2` feeds `floss`/`mtot` in `evolve()` Phase 1, so a change here would move the pinned trajectory, but `m2` itself is not directly asserted |
| Orbital period, log-flat between `pmin`/`pmax` (Sec. 2.1) | `binaries/population.py::generate_population` -- `period_distribution="log_uniform"` (default, see above) | **unit-tested** — `tests/test_population_generation.py::test_orbital_period_bounds_and_log_uniform_distribution` checks exact bounds plus the log-uniform shape statistically (fixed seed). Outside the integration pin's reach (`period` never feeds `floss`), so this test is the *only* coverage. |
| Semi-major axis from Kepler's third law (Sec. 2.1) | `binaries/population.py::generate_population` (`AFAC` constant) | **unit-tested** — `test_semi_major_axis_matches_reference_formula` recomputes `a` independently from `m1`/`m2`/`period` and asserts an exact match (`rtol=1e-12`), and pins the exact value of `AFAC`. Note: `AFAC` is not independently re-derivable to full precision from Kepler's third law alone (a from-scratch physical derivation gives a ~1% different prefactor) — a units/precision-convention detail, not a bug. `self.a` (this constant's output, threaded everywhere in `BinaryPopulation`) is in **AU**, not Rsun — see the "real bug found and fixed" row in Section 10 below for a case where that unit mattered and was originally missed. |
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
| L_bol/L_UV wired into `ClusterSimulation.run()`'s own per-timestep output (A2, `docs/science/paper1-detailed-work-breakdown.md`) | n/a -- Realta-specific integration | `simulation/cluster.py::ClusterSimulation.run()` (`lbol_tot`/`luv_tot` keys, alongside the pre-existing `lumx_tot`/`nphot_tot`) | **unit-tested**: `tests/test_cluster_simulation_observables.py` confirms both keys are present and finite for every timestep, and that they exactly match `ms_table.get_lbol()`/`uv_table.get_luv()` evaluated independently -- these used to be recomputed *after* a run completed, in `scripts/run_paper1_experiment.py`, duplicating logic; that script now reads them straight from `results` instead (`run_variant()` simplified accordingly, `MSLuminosityTable`/`UVLuminosityTable` imports removed there as now-unused). `lbol_tot = ms_lbol + lumx_tot`; `luv_tot` is MS-only (no accretion-UV model exists, unchanged scope note). |
| Independent Q_H(t) from the massive-star population (A3, `docs/science/paper1-detailed-work-breakdown.md` -- previously `nphot_tot = NPHOT_PER_LUMX * lumx_tot`, a fixed constant multiple of L_X carrying zero independent information, see Section 3) | Massive-star ionizing-photon output, per star, summed over currently-alive `M >= 8` Msun stars | `simulation/cluster.py::ClusterSimulation._qh_ms_tot`, using the previously-unused `io/tables.py::IonizingPhotonTable.get_ngamma`; `qh_tot = qh_ms_tot + nphot_tot` in `run()`'s results | **unit-tested, sensitivity-verified, one real gap found and fixed along the way**: `tests/test_cluster_simulation_observables.py` confirms (1) a calibration check -- `get_ngamma(m)` interpreted as the *total* photon count over the star's whole MS lifetime (confirmed via its own MUNIT/MATOM baryon-count conversion), divided by `LifetimeTable.get_lifetime(m)` in seconds, lands within the well-known literature range for O/early-B ionizing rates at 10/40/80 Msun (Vacca, Garmany & Shull 1996, ApJ 460, 914) -- an order-of-magnitude sanity check, not an exact-match requirement; (2) `qh_tot` is not a fixed multiple of `lumx_tot` (confirms real independent variation, the actual thing A3 fixes); (3) a star stops contributing once its own `nturn`/lifetime clock says it has died, and a companion contributes correctly while alive across the whole `t2_lifetime` window (both hand-constructed scenarios, following `tests/test_evolve.py`'s own technique). I verified sensitivity directly: removed the `nturn==0` guard on the primary's contribution, confirmed both of those last two tests fail exactly as expected, reverted. Measured negligible runtime overhead on the full Paper 1 pipeline (~49s either way, `ntot=100_000` x 5 prescriptions x 201 steps) before finalizing the per-star-loop approach over a vectorized alternative. **Real gap found and fixed**: the `"single"` `binary_prescription` used to empty `BinaryPopulation.m1` entirely (`n_massive=0`), a shortcut that was harmless while only `L_X`/HMXB-related quantities read `m1` (correctly zero for single stars either way) -- but `_qh_ms_tot` also reads `m1`/`nturn` to know which massive stars exist, so with an empty array it silently reported `Q_H=0` for single-star populations, which is physically wrong (massive single stars still ionize). Migrated `"single"` onto the same `has_companion=False` mechanism A1 already built for `binary_fraction<1` (`m1` stays populated, `m2=0`) -- `L_X`/HMXB activation stay exactly zero either way (`m2=0` blocks that unconditionally), only `Q_H`/`L_bol` tracking is fixed. `tests/test_binary_prescriptions.py::test_single_prescription_suppresses_binary_formation` updated to assert the new (correct) `m1`-populated/`m2=0` state instead of the old empty-array one; `tests/test_cluster_simulation_observables.py::test_single_prescription_gives_nonzero_qh` pins the fix directly. |

## 6. Paper 1 binary-interaction and merger prescriptions

**Not paper-derived.** Power et al. (2009/2013) contain no mass-transfer,
common-envelope, or merger physics at all -- these rows document a new,
Realta-specific parameterization added to support Paper 1's basic
experiment (`docs/science/research-programme.md`), reviewed and
accepted in `docs/science/paper1-binary-interaction-proposal.md`
(chat, 2026-08-24). Every parameter below defaults to a value that
makes it a no-op, so the pre-existing `fsur`-only baseline (Sections
1-3 above) is reproduced exactly unless `binary_prescription` is
explicitly set away from `"non_interacting"` -- confirmed directly: the
full pre-existing test suite (`tests/test_regression.py`'s pinned
trajectories included) passes unchanged after this addition.

| Prescription/mechanism | Realta implementation | Status |
|---|---|---|
| `binary_prescription="single"` -- no companion assigned to any M >= mcut star, no HMXB channel | `binaries/population.py::generate_population` (`n_massive` forced to 0) | **unit-tested** -- `tests/test_binary_prescriptions.py::test_single_prescription_suppresses_binary_formation` |
| Pre-SN merger channel (`p_merge`/`p_merge_max_period`/`f_merge`) -- eligible short-period binaries merge at formation, `m1 -> m1 + f_merge*m2`, `m2 -> 0`, lifetime recomputed for the merged (rejuvenated) mass. **No longer auto-enabled by `enhanced_mergers`** (see the reconciliation note below) -- remains available as an explicit, independent override. | `binaries/population.py::generate_population` | **unit-tested** -- `test_merger_channel_folds_companion_into_primary_and_disables_hmxb`, `test_merger_channel_respects_period_threshold`. I verified sensitivity directly (stopped the merger code from zeroing `m2`, confirmed the mechanism test fails, reverted). Minimal `did_merge`/`merge_time` bookkeeping is recorded for a future Figure 3 but no `Event`/`PopulationHistory` abstraction is built yet (see `docs/science/development-roadmap.md` item 4). |
| Config validation (`binary_prescription` enum, non-negative/`[0,1]` bounds on the four new parameters) | `config.py::SimulationConfig.__post_init__` | **unit-tested** -- `test_invalid_binary_prescription_rejected`, `test_out_of_range_interaction_params_rejected` |
| No-op guarantee for the default prescription (`interaction_boost=1.0`, `p_merge=0.0`, `use_rlof_classifier=False` resolved automatically) | `config.py::_PRESCRIPTION_DEFAULTS`, `SimulationConfig.__post_init__` | **unit-tested** -- `test_default_prescription_matches_baseline_exactly`; also implicitly reconfirmed by every pre-existing pinned/unit test in this document continuing to pass unmodified. |

### Reconciliation with the physics-based RLOF classifier (2026-08-24)

`standard_interaction`/`enhanced_interaction`/`enhanced_mergers` now
drive their behaviour through the real Roche-lobe-overflow classifier
(Section 10) rather than purely through the placeholder parameters
above -- see `docs/science/rlof-ce-classifier-proposal.md` "Decision
3". This is a genuine, deliberate scientific-behaviour change for
these three prescriptions specifically -- NOT for `"single"` or
`"non_interacting"`, which remain exactly as pinned (confirmed: the
full pre-existing test suite, including `tests/test_regression.py`'s
pinned trajectories, passes unchanged).

| What changed | Old (pre-2026-08-24) | New | Realta implementation | Status |
|---|---|---|---|---|
| `standard_interaction`/`enhanced_interaction`: when `interaction_boost` applies | Unconditionally, to every surviving binary regardless of interaction history | Only to binaries the RLOF classifier found underwent stable mass transfer on the MS (`use_rlof_classifier` now implied `True` by these prescriptions) -- a binary that never interacted uses plain `fsur` | `binaries/population.py::evolve`, Phase 1 (`had_stable_mt` gate) | **unit-tested, sensitivity-verified** -- `tests/test_binary_prescriptions.py::test_interaction_boost_raises_activation_fraction_for_stable_mt_systems` (boost applies when flagged as interacted) and `test_interaction_boost_not_applied_to_never_interacted_systems` (plain `fsur` when not) are a matched pair; I verified sensitivity directly for both (forced `had_stable_mt=True` unconditionally, confirmed the "not applied" test fails as expected; reverted the boost condition to always-True path already covered by the first test's own sensitivity check from the wiring session). |
| `enhanced_mergers`: merger driver | Independent random formation-time draw (`p_merge=0.2`/`p_merge_max_period=10.0`/`f_merge=0.5`), unrelated to any real physics | `q_crit_ms=0.4` (lower than HTP02's own fiducial 0.695 -- see `binaries/interaction.py::Q_CRIT_MS`), `use_rlof_classifier=True` -- more RLOF-ing systems are classified as dynamically-unstable immediate mergers by the real classifier. `p_merge` defaults to 0 for this prescription now (available as an independent override, not auto-enabled) | `config.py::_PRESCRIPTION_DEFAULTS["enhanced_mergers"]` | **unit-tested** -- covered transitively by `tests/test_rlof_classifier.py`/`test_rlof_wiring.py`'s `q_crit_ms`-override tests (the mechanism is generic, not prescription-specific); no dedicated `enhanced_mergers`-specific numeric pin exists yet since this milestone didn't require reproducing a specific published `enhanced_mergers` figure. |
| Genuine numpy bug found and fixed during this reconciliation | n/a | `np.full(n, some_str_enum_member, dtype=object)` silently truncates/corrupts the fill value (confirmed: results in a plain, truncated `str` that fails equality against the real enum member, despite the array's dtype correctly showing `object`); per-element assignment and per-element scalar comparison are unaffected, but a **vectorized** `array == RLOFOutcome.X` comparison is also silently broken (returns all-`False` even when every element genuinely equals it) -- a real footgun for any future diagnostic code counting outcomes. Fixed by switching to list-based array construction (`np.array([...], dtype=object)`) and documented directly in `RLOFOutcome`'s own docstring. | `binaries/population.py::generate_population`, `binaries/interaction.py::RLOFOutcome` | **unit-tested** -- `tests/test_rlof_wiring.py::test_rlof_outcome_array_construction_avoids_np_full_corruption`. I verified sensitivity directly: reverted to `np.full(...)`, confirmed the test fails with the exact truncated-string symptom, reverted back. |

### Paper 1 config: `pmin` raised for this experiment only (2026-08-24)

Found while running `scripts/run_paper1_experiment.py` end-to-end for
the first time: `configs/paper1_basic_experiment.yml` inherited the
global `pmin=0.1` days default (`config.py`/`config.yml`, itself a
pre-existing Power et al. 2009 value, unchanged by this session). At
that `pmin`, Kepler's third law puts most 8-100 Msun binaries' orbital
separations below the stars' own ZAMS radii -- born already in
contact. This was invisible before the RLOF classifier existed (the
old `fsur`-only model never compared separation against stellar
radius), and is a separate finding from the AU/Rsun units bug above
(confirmed independently: raising `pmin` alone, with the units bug
still present, did not fix the all-merger degeneracy -- see the
git history around this fix). With the units bug fixed, `pmin=1.0` day
(raised from `configs/paper1_basic_experiment.yml` only, not the
global default -- confirmed the full pre-existing test suite, which
uses the global default, is unaffected) gives a real, non-degenerate
mix of RLOF outcomes. See that config file's own comment for the
verification detail.

**Residual, understood degeneracy, root cause found (2026-08-24)**:
`standard_interaction`/`enhanced_interaction`/`enhanced_mergers`
produce bit-identical Figure 1/2 output in this config -- confirmed by
diffing their `.tevol.dat` files directly (zero differences across all
201 timesteps). Two separate causes, both now understood as expected
physics, not bugs:

1. `enhanced_mergers`'s lower `q_crit_ms=0.4` produces an *identical*
   RLOF-outcome distribution to the `q_crit_ms=0.695` default --
   checked directly by re-running `find_rlof_onset` at both thresholds
   for every generated binary; zero outcomes differ. This is the
   already-documented `find_rlof_onset` emergent property (Section 10
   above): the automatically-selected donor is almost always the more
   massive star (`q1 > 1`), already far above either threshold, so
   lowering `q_crit_ms` has no binaries left to reclassify here.

2. `interaction_boost` (`standard_interaction=1.5` vs
   `enhanced_interaction=3.0`) never has anything to act on: it only
   applies to binaries the classifier found underwent
   `STABLE_MASS_TRANSFER`, and **`STABLE_MASS_TRANSFER` cannot fire at
   all in the current design, for any config, not just this one**.
   Traced this precisely (all 51 `STABLE_MASS_TRANSFER`-classified
   binaries in this run, e.g. `m1=54.0, m2=2.6 Msun`): `generate_population`
   always enforces `m2 <= m1`, and `STABLE_MASS_TRANSFER` requires
   `q1 = donor/companion < q_crit < 1`, so the stable-MT donor is
   *always* `m2`, the lighter star, by construction. But Realta only
   tracks one explosion clock (`turnoff_time`), computed from `m1`'s
   own (Schaerer) lifetime -- and `m1`, being far more massive, always
   explodes first (e.g. `turnoff_time=4.4 Myr` for that 54 Msun
   primary vs. `donor_tbgb=480 Myr` for its 2.6 Msun companion's own
   Hurley/Tout MS+HG duration). `nturn` flips to 1 long before the
   donor's own predicted `rlof_time`, and Phase 0 requires
   `nturn==0`, so the event is permanently gated out.

   This is *not* the two lifetime prescriptions disagreeing with each
   other for a given star -- directly compared Schaerer vs. Hurley/
   Tout `t_bgb` across 8-100 Msun at Z=0.008 and they agree to within
   ~5% throughout (ratio 0.95-1.04, and Hurley/Tout is even slightly
   *longer* above ~80 Msun). It is a structural/physical fact about
   the mass hierarchy: pre-SN "stable mass transfer with the lighter
   star as donor" requires the *lighter* star to evolve off the MS on
   its own (hundred-Myr-scale) timescale, which is always far longer
   than the *heavier* primary's own (few-Myr-scale) pre-SN lifetime --
   so the primary has essentially always already exploded before this
   channel could ever become physically reachable. `IMMEDIATE_MERGER`/
   `COMMON_ENVELOPE` don't have this problem: their donor is (per the
   same emergent property) typically `m1` itself, so both the RLOF
   clock and the SN clock are the same star's, which agree to ~5% as
   shown above.

   The astrophysically important, missing piece this exposes is not a
   fix to the pre-SN classifier, but a **structurally different,
   not-yet-modelled channel**: the secondary star's *own* later
   Roche-lobe overflow onto the now-compact primary (Case B/C mass
   transfer onto a compact object, wind-accretion vs. RLOF) -- the
   actual dominant real-world HMXB-formation pathway once the primary
   has already collapsed. Phase 0 was scoped this session to pre-SN
   interaction between two still-evolving stars only (`nturn==0`,
   matching HTP02's own CE-eligible donor list); a post-SN
   secondary-RLOF channel (`nturn==1`, donor = the still-evolving
   secondary, companion = a compact remnant) is new scope, deliberately
   not attempted here, and belongs on the roadmap for whenever the
   mass/timescale regime it addresses is actually the one under study
   -- see the "Known gaps" note below.

## 7. Fig. 1/2 UV observable (Paper 1, not in either paper's original code)

| Quantity | Source | Realta implementation | Status |
|---|---|---|---|
| Population far-UV luminosity vs. time, `L_UV(t)` | GALEX FUV (~1528 A), via FSPS `get_mags(bands=["galex_fuv"])` -- band choice documented and reviewed in `docs/science/paper1-binary-interaction-proposal.md`'s "UV band decision" section | `scripts/generate_fuv_luminosities.py` (data-generation script, mirrors `scripts/generate_ms_luminosities.py`), `io/tables.py::UVLuminosityTable` (mirrors `MSLuminosityTable`'s fiducial-mass-then-rescale convention exactly), wired into `scripts/run_paper1_experiment.py` | **done** -- `fuv_lbol_z*.dat` (all three metallicities) generated and placed in `src/realta/data/`; `UVLuminosityTable` loads them and is unit-tested against the real data (`tests/test_uv_luminosity_table.py::test_all_three_metallicity_tables_load_and_scale`, mirroring `MSLuminosityTable`'s equivalent test) plus a synthetic-table test isolating the rescaling arithmetic. `scripts/run_paper1_experiment.py` runs end-to-end and produces both figures with real `L_UV(t)`. |

## 9. Stellar radius/luminosity module (RLOF/CE classifier prerequisite)

Not part of the Power et al. baseline -- new physics added to support
Paper 1's "enhanced interaction"/"enhanced mergers" prescriptions with
real Roche-lobe-overflow physics in place of the illustrative
`interaction_boost`/`p_merge` placeholders in Section 6 above. See
`docs/science/rlof-ce-classifier-proposal.md` for the full design
record, scope decisions and the transcription-risk finding below.
**Used internally only** -- Realta's reported `L_bol`/`L_UV` continue
to come from `MSLuminosityTable`/`UVLuminosityTable` (Section 4/7
above), not from this module (see the proposal doc's "Decision 1").

| Quantity | Source | Realta implementation | Status |
|---|---|---|---|
| Zero-age main-sequence luminosity/radius, `L_ZAMS(M,Z)`/`R_ZAMS(M,Z)` | Tout, Pols, Eggleton & Han (1996, MNRAS 281, 257), eqs. 1-4, Tables 1-2 | `stellar/zams.py::zams_luminosity`, `zams_radius` | **unit-tested** -- `tests/test_hurley_main_sequence.py::test_zams_luminosity_radius_solar_calibration` (1 Msun/Z=0.02 close to but below 1 Lsun/Rsun, as physically expected for the ZAMS vs. the Sun's present, evolved state). Coefficients transcribed from a clean-text (non-image) PDF read -- high transcription confidence, unlike the 2000 paper's denser appendix tables. |
| Main-sequence lifetime, luminosity, radius vs. age, `t_MS`/`L_MS(t)`/`R_MS(t)` (stellar types k=0,1) | Hurley, Pols & Tout (2000, MNRAS 315, 543), Section 5.1 (eqs. 1-24), Appendix A coefficients a1-a81 | `stellar/main_sequence.py` | **unit-tested, self-consistency-pinned, one real bug caught and fixed**: `tests/test_hurley_main_sequence.py` covers solar-calibration sanity bounds, MS-lifetime sanity bounds, monotonic radius/luminosity growth across the full mass range up to and beyond `mcut`, exact convergence to `R_TMS`/`L_TMS` at t->t_MS, the phase() scope guard, and the low-mass degenerate radius floor (eq. 24). **A genuine coefficient transcription error was caught during this session**: `a40`'s gamma exponent was originally read as -2 instead of -1, and `a41`'s alpha exponent as -1 instead of 0 (10x too small) -- both in the Delta_R "hook" perturbation block (eq. 17). This caused `ms_radius()` to collapse to near-planet-size (0.02-0.05 Rsun) for 5/20 Msun stars through most of the MS, recovering the correct value only in the last ~1% of the lifetime -- caught by the monotonicity test, root-caused, and fixed by cross-checking a user-supplied copy-pasted excerpt of the source table against the original image-based read. I verified sensitivity directly: reverted the fix, confirmed `test_ms_radius_and_luminosity_monotonically_increase` fails for M>=3 Msun exactly as expected, then re-applied the fix. Pinned self-consistency values (`test_pinned_values_z_solar`) are this implementation's own output, not an independently-computed reference -- see "Still open" below. |
| Giant-branch base luminosity/radius-vs-L relation, `L_BGB(M,Z)`/`R_GB(M,L,Z)` -- used as HG's endpoint boundary values, not for time-evolution along the GB itself | Hurley et al. (2000), eq. 10 (a27-a32) and eqs. 46-48 (b1-b7, mass-radius exponent x) | `stellar/giant_branch.py::l_bgb`, `r_gb`, `mass_radius_exponent` | **unit-tested, three more real transcription errors caught and fixed**: `tests/test_giant_branch.py` covers a solar-calibration sanity bound for `L_BGB` (a few Lsun, matching the textbook expectation for the Sun's eventual BGB point), monotonic growth with mass, and -- the key check -- cross-validation of `r_gb()` against Hurley et al.'s own illustrative Z=0.02 example formula (`R_GB ~= 1.1*M^-0.3*(L^0.4+0.383*L^0.76)`, stated directly in their Section 5.2 text). This caught: (1)/(2) `b4`'s gamma/eta/mu had been duplicated from `b5`'s row, and `b5`'s own alpha exponent was wrong (both caught via user-pasted excerpts, before this row made it into code); (3) an entire `b2`-clamping step (`b2 = min[max(b2, -0.04167+55.67Z), 0.4771-9329.21Z^2.94]`) had been dropped from the implementation entirely, with the retained constant also mis-read as `-0.14167`. The illustrative-formula cross-check disagreed by up to 14x before this fix and agrees to ~10-20% after (consistent with "simplified approximation vs. full fit", not a remaining bug) -- see `docs/science/rlof-ce-classifier-proposal.md`'s "Update, implementation session" note for the full account. I verified sensitivity directly: reverted the `b2` clamp, confirmed the illustrative-formula test fails with a ~1.4x-too-small radius, reverted back. `a27-a32` (`L_BGB`) were separately verified against user-pasted excerpts and caught one more error (`a28`'s eta exponent, `e-2` instead of `e0`, 100x off) before any code was written against them. |
| Hertzsprung-Gap luminosity/radius vs. age, `L_HG(t)`/`R_HG(t)` (stellar type k=2, M < M_FGB only) | Hurley et al. (2000), Section 5.1.2, eqs. 25-30 | `stellar/main_sequence.py::hg_luminosity`, `hg_radius`, `l_ehg`, `r_ehg`, `m_fgb`, updated `phase()` | **unit-tested**: `tests/test_hertzsprung_gap.py` covers monotonic radius growth across the HG (the defining physical feature -- rapid expansion crossing the gap), exact convergence to `R_TMS`/`R_EHG` (and `L_TMS`/`L_EHG`) at both endpoints, `phase()` correctly reporting k=2 during HG and raising past `t_BGB` (GB itself remains out of scope), a sanity bound on `M_FGB` (10-16 Msun at solar Z), and that `l_ehg`/`r_ehg`/`hg_radius` raise rather than silently guess for M >= M_FGB (those stars skip the GB entirely and need L_HeI/R_HeI, which are CHeB-phase quantities this module does not implement). `phase()` itself still reports k=2 for M >= M_FGB stars (it identifies the evolutionary type, not whether this module can compute a radius for it) -- callers needing an actual radius must be prepared for the second raise. |
| HG core mass vs. age, `M_c,HG(t)` (intermediate-mass branch only, `M_HeF <= M < min(M_FGB, CORE_MASS_BGB_MAX_MASS)`, i.e. roughly 2-7.3 Msun at solar Z) -- the CE-prerequisite this module was extended for | Hurley et al. (2000), eqs. 2 (`M_HeF`), 28-30 (core-mass growth); `M_c,BGB` itself via the paper's own stated large-M asymptotic limit, NOT the full eq. 44 (see status column) | `stellar/main_sequence.py::m_hef`, `core_mass_ehg`, `core_mass_hg`, `_rho_coefficient`; `stellar/giant_branch.py::core_mass_bgb`, `CORE_MASS_BGB_MAX_MASS` | **unit-tested; one approximation caught being wrong and corrected mid-implementation, not just flagged in advance**: `tests/test_hg_core_mass.py` covers monotonic core-mass growth with initial mass, `M_HeF`'s solar-calibration sanity value (~2 Msun), exact convergence to `rho*M_c,EHG`/`M_c,EHG` at the start/end of HG, that `core_mass_ehg` raises outside the `M_HeF <= M < M_FGB` bracket, and (the key addition) that `core_mass_bgb` stays sub-Chandrasekhar within its supported range and raises above `CORE_MASS_BGB_MAX_MASS`. **What happened**: `M_c,BGB` uses `0.098*M^1.35` (Z-independent), the large-mass asymptotic limit Hurley et al. state directly in their Sec. 5.2 text, adopted to avoid the full eq. 44 (which needs the entire GB core-mass-luminosity relation, eqs. 31-43 -- a materially larger addition than a coefficient block, discovered only when eq. 44 was actually re-read for `core_radius`, not when the approximation was first adopted). Initially accepted as "valid for M >= M_HeF" with no upper bound; implementing `core_radius` immediately surfaced a real problem -- at M=10 Msun the approximation gives a super-Chandrasekhar (2.17 Msun) "core mass", which is physically implausible for a star still at the end of HG (real HG/GB-base core masses in this range are a few tenths of a solar mass) and broke `white_dwarf_radius()` outright (NaN under the sqrt). Root cause identified precisely: `0.098*M^1.35` is exactly `c1^0.25*M^(c2/4)` using eq. 44's own `c1=9.20925e-5`, `c2=5.402216` (confirmed self-consistent: `c1**0.25~=0.098`, `c2/4=1.350554~=1.35`), and exceeds `M_ch=1.44` above `M=(M_ch^4/c1)^(1/c2)~=7.317` Msun -- derived directly from those constants, not a separately chosen cutoff. `core_mass_bgb` now raises above that mass rather than returning the wrong value, the same "raise, don't guess" idiom already used for `M_FGB`/`M_HeF`/`t_BGB` elsewhere in this module. I verified sensitivity directly: replaced the linear eq. 30 growth with a constant (endpoint-convergence tests fail as expected) and separately disabled the new cap (the raise-above-cutoff test fails as expected, reproducing the original NaN bug); both reverted. |
| Core radius, `R_c1` -- the CE binding-energy formula's other core-structure input alongside core mass | Hurley et al. (2000), Section 6.2.1 (white dwarfs, eq. 91) and Section 6.3 (core-radius definition for HG/GB stars: `R_c = R_WD(M_c)` for `M >= M_HeF`, `R_c = R_ZHe(M_c)` otherwise -- not implemented, out of scope) | `stellar/remnant.py::white_dwarf_radius`, `core_radius`, `M_CHANDRASEKHAR`, `R_NEUTRON_STAR` | **unit-tested; not put through a paste-verification round, and why**: `tests/test_remnant.py` covers a real-value cross-check (0.6 Msun white dwarf, ~0.0115-0.013 Rsun, matching the well-known Sirius-B-territory value), the defining mass-radius *direction* (heavier white dwarfs are smaller, not larger -- an easy sign error to miss), shrinkage toward zero approaching the Chandrasekhar mass, the neutron-star-radius floor, and an end-to-end `core_mass_hg -> core_radius` chain for a real HG donor. Unlike the a/b coefficient tables, this formula is clean text (not a dense appendix table) and was cross-checked against a well-known real value *before* writing any code -- both factors together were judged sufficient confidence without a separate user-paste round; the 0.6 Msun check above is that same cross-check, now pinned as a regression test. `core_radius` raises for `M < M_HeF` (the `R_ZHe` branch, not implemented) exactly like `core_mass_ehg` does, and inherits `core_mass_hg`'s upper bound (`CORE_MASS_BGB_MAX_MASS`) transitively through whatever core mass is passed in. |

**Scope note**: only stellar types k=0 (fully/deeply convective MS,
M<=0.7), k=1 (radiative-core MS, M>0.7) and k=2 (Hertzsprung Gap,
M < M_FGB for radius/luminosity, and additionally
M_HeF <= M < CORE_MASS_BGB_MAX_MASS (~7.3 Msun) for core mass/radius
specifically) are implemented. The true giant branch (k=3) and later
phases remain out of scope -- `phase()` raises `ValueError` for any
`t >= t_BGB`, and `l_ehg`/`r_ehg`/`hg_radius`/`hg_luminosity`/
`core_mass_ehg` raise for HG stars with M >= M_FGB specifically (no GB
phase at all for those masses; they'd need L_HeI/R_HeI/M_c,HeI
instead), `core_mass_ehg` additionally raises for M < M_HeF (needs a
different, degenerate-core low-mass relation, M_c,GB(L_BGB), not
implemented -- rare for Realta's mcut=8 Msun default), and
`core_mass_bgb`/`core_mass_ehg`/`core_mass_hg`/`core_radius` raise for
M >= CORE_MASS_BGB_MAX_MASS even when M < M_FGB (the asymptotic
core-mass approximation itself breaks down there -- see the core-mass
row above). Callers must treat these as "not modelled" rather than
guessing. See the proposal doc for the full transcription-verification
history.

## 10. RLOF outcome classifier (Stage 2 of the interaction module)

Not part of the Power et al. baseline. Given a binary's masses,
separation and metallicity, and the donor's age, classifies whether
Roche-lobe overflow is occurring and what the outcome is, for MS
donors only (see Section 9's scope note; HG+ donors return
`PHASE_NOT_MODELLED`, not a guess).

**Wired into `BinaryPopulation` as a new, opt-in event**
(`config.use_rlof_classifier`, default `False`) -- see
`config.py::use_rlof_classifier`'s field docstring. It ADDS a new
"Phase 0" MS-RLOF event in `evolve()`, alongside the existing
`fsur`/`interaction_boost` gate. As of 2026-08-24 it is also
reconciled with `interaction_boost`/`enhanced_mergers` -- see Section
6's "Reconciliation" subsection for the old-vs-new behaviour and the
proposal doc's "Decision 3" for the rationale. The onset time is
precomputed once per binary in `generate_population()`
via root-finding (`interaction.py::find_rlof_onset`), mirroring how
`turnoff_time` is precomputed from `LifetimeTable` -- this keeps
`evolve()`'s per-timestep cost identical to before (a cheap numpy
comparison), rather than calling the classifier itself every step for
every binary (`ntot=100,000` default makes that infeasible).

| Quantity | Source | Realta implementation | Status |
|---|---|---|---|
| Roche-lobe radius, `R_L1(a, q1)` | Eggleton (1983) fit, as used by Hurley, Tout & Pols (2002, MNRAS 329, 897), eq. 53 | `binaries/interaction.py::roche_lobe_radius` | **unit-tested** -- `tests/test_rlof_classifier.py::test_roche_lobe_radius_equal_mass_sanity_value` checks the well-known q=1 literature value (R_L1/a~0.379), plus a separation-scaling check and a non-positive-mass-ratio rejection test. |
| MS-donor RLOF classification (detached / stable mass transfer / immediate merger) | HTP02 Sec. 2.6.4 (`q_crit_ms=0.695` for dynamically-unstable MS donors, extended from k1=0 to k1=1 as a named simplification -- see `docs/science/rlof-ce-classifier-proposal.md` "Decision 2") and Sec. 2.7.1 (MS donors are not CE-eligible, so dynamical instability merges directly rather than forming a CE) | `binaries/interaction.py::classify_rlof` | **unit-tested** -- `tests/test_rlof_classifier.py` covers all three reachable outcomes with concrete mass/separation scenarios, the `PHASE_NOT_MODELLED` guard past `t_MS`, and that overriding `q_crit_ms` actually changes the classification (not a dead parameter). I verified sensitivity directly: removed the `q1 > q_crit_ms` branch, confirmed two tests fail exactly as expected, reverted. |
| HG-donor RLOF classification (detached / stable mass transfer / common envelope) | HTP02 eqs. 56-57 (GB `q_crit`, reused for HG donors) and Sec. 2.7.1 (HG donors ARE CE-eligible, unlike MS donors) -- the reuse itself is Zuo & Li (2014, MNRAS 442, 1980)'s eq. 1, citing Shao & Li (in prep.); HTP02's own HG treatment is a crude fixed `q_crit=4`, which it calls "rather approximate" | `binaries/interaction.py::hg_q_crit`, extended `classify_rlof` | **unit-tested, sensitivity-verified, now reachable through the full pipeline**: `tests/test_hg_ce_classifier.py` covers `hg_q_crit` matching a direct formula evaluation and being materially below HTP02's crude HG default (confirms the refinement does something, not a no-op), the key behavioural difference from MS donors (unstable HG donor -> `COMMON_ENVELOPE`, not `IMMEDIATE_MERGER`), stable mass transfer and detached cases, that a donor still identified as HG (`phase()` returns k=2) but with mass outside the core-mass-tracking range correctly falls back to `PHASE_NOT_MODELLED` rather than crashing, and (since `find_rlof_onset` was extended to search HG -- see the row below) end-to-end confirmation that the root-finder actually reaches this classification, not just that it exists in isolation. I verified sensitivity directly: removed the `q1 > q_crit` branch, confirmed the CE-outcome test fails exactly as expected, reverted. The CE outcome itself (survive vs. merge, mass/orbit update) is not resolved here -- `COMMON_ENVELOPE` is a terminal classification with no consequence model yet (see `evolve()`'s Phase 0, Section 12), same status as the still-unimplemented energy-balance solve. |
| Immediate-merger mass combination | Not from HTP02 -- that paper gives an explicit envelope mass-loss prescription only for common-envelope mergers (eqs. 69-77), not for a direct/dynamical MS-MS collision, the only merger channel reachable by this module's current scope. Conservative merging (no mass loss) is used as the natural default absent an explicit prescription -- a named simplification, not a citation. | `binaries/interaction.py::merge_stellar_masses` | **unit-tested** -- trivial addition check. Flagged explicitly in the function's own docstring as not paper-derived. |
| Per-binary RLOF onset time (which star donates first, and when), across MS and HG | Root-find `R_donor(t) = R_L1` once per binary via bisection, using `ms_radius()`'s confirmed MS-wide monotonicity and (since 2026-08-24) `hg_radius()`'s confirmed HG-wide monotonicity, continuous at the MS/HG boundary | `binaries/interaction.py::find_rlof_onset` (renamed from `find_ms_rlof_onset` -- the old name became inaccurate once HG search was added) | **unit-tested, one emergent finding documented, one real bug caught downstream**: `tests/test_rlof_classifier.py` covers a wide-binary (through both MS and HG) never-interacts case, a correctly-identified MS crossing time/donor (cross-checked directly against `ms_radius()`/`roche_lobe_radius()` at the returned time), correct donor selection when the *lighter*-labelled star (m2) is actually the physical donor, the born-overflowing (t=0) edge case, and a documented emergent property: because Eggleton's `R_L1/a` increases monotonically with the donor's own mass ratio, `IMMEDIATE_MERGER` dominates over `STABLE_MASS_TRANSFER` for automatically-selected MS donors in practice (a parameter sweep during development found no natural stable-MT case) -- not a bug, see the module docstring. `tests/test_hg_ce_classifier.py` separately covers both `COMMON_ENVELOPE` and `STABLE_MASS_TRANSFER` being reached via HG-phase root-finding, and a separation that stays detached through the MS but is crossed once the donor expands during HG (confirms the HG search actually does something, not a no-op). **Real bug found and fixed downstream, in `evolve()`'s consequence code** (Section 12): the `STABLE_MASS_TRANSFER` branch unconditionally called `ms_radius()` for the donor's radius, which silently breaks (`brentq` bracket-sign error) once a stable-MT crossing can legitimately happen during HG -- fixed by checking `phase()` at `rlof_time` and using `hg_radius()` for HG donors; see `tests/test_rlof_wiring.py::test_evolve_applies_stable_mass_transfer_for_hg_donor`, sensitivity-verified (reverted the fix, confirmed the test fails with the same class of error, reverted back). |
| Wiring into `BinaryPopulation` (opt-in `config.use_rlof_classifier`, precomputed once per binary, processed in `evolve()`'s new "Phase 0") | n/a -- Realta-specific integration, not itself a citable physics choice | `config.py::use_rlof_classifier`/`q_crit_ms`, `binaries/population.py::generate_population`/`evolve` | **unit-tested, sensitivity-verified, baseline confirmed untouched**: `tests/test_rlof_wiring.py` covers the disabled-by-default inert state, bit-identical output between a config that doesn't mention `use_rlof_classifier` and one that sets it `False` explicitly, graceful Z=0 (imetal=1) handling (warns, skips, doesn't crash), that a realistic massive-star population actually produces some finite RLOF times, that `evolve()` correctly applies an immediate merger (mass combination, `m2` zeroed, `did_merge`/`merge_time` set, lifetime clock reset from the merge time) and does not reprocess it on a later timestep, and (Section 12 below) that stable mass transfer is correctly applied and not reprocessed either. I verified sensitivity directly: removed the `m2` zeroing on merge, confirmed the merger-application test fails exactly as expected, reverted. The full pre-existing test suite passes unchanged, confirming the opt-in flag does not perturb the pinned baseline. |
| **Real bug found and fixed, AU/Rsun unit mismatch** | n/a -- a units-consistency bug in the wiring, not a citable physics choice | `binaries/population.py::BinaryPopulation.RSUN_PER_AU`, applied at every call into `binaries/interaction.py` from `generate_population`/`evolve` | **found by running the full Paper 1 pipeline end-to-end** (`scripts/run_paper1_experiment.py`) for the first time this session, not by unit tests -- every unit test for `interaction.py` was written with hand-picked, self-consistent Rsun-scale separations, so none of them could have caught a units mismatch at the `population.py` boundary. `self.a` (`AFAC`/`PFAC`, confirmed via an Earth-Sun sanity check: `AFAC`'s formula gives `a~=0.99` for `M=1 Msun, P=365.25 days`, only sensible as ~1 AU, and `PFAC=365.229126` is essentially a sidereal year in days -- i.e. these constants encode Kepler's third law in the standard AU/Msun/year convention) is in **AU**, but every stellar-radius function this session added (`ms_radius`/`hg_radius`/`core_radius`, from Hurley et al. 2000/2002) returns **Rsun**, and `roche_lobe_radius`/`classify_rlof`/`find_rlof_onset`/`apply_stable_mass_transfer`/`apply_common_envelope` all compare `separation` directly against those Rsun-scale radii with no conversion. Effect before the fix: every donor looked ~215x closer to its Roche lobe than it really is, so running the actual Paper 1 config produced `IMMEDIATE_MERGER` for essentially 100% of massive binaries in every RLOF-classifier-enabled prescription, making Figure 2 (the central figure) completely degenerate (`L_X=0`/`Q_H=0` for `standard_interaction`/`enhanced_interaction`/`enhanced_mergers` across the whole 100 Myr run). Fixed by converting `self.a[i]` to Rsun (`* RSUN_PER_AU`) at each `interaction.py` call site and converting returned separations back to AU before storing -- `self.a` itself stays in AU everywhere else (the pre-existing SN1 mass-loss orbit-widening code, and every already-pinned regression value, is untouched). `RSUN_PER_AU = 215.032` (1 au = 1.495978707e11 m, IAU 2012 exact definition; R_sun = 6.957e8 m, IAU 2015 nominal). Confirmed the fix by re-running the full pipeline: `standard_interaction`/`enhanced_interaction`/`enhanced_mergers` now show a real mix of `DETACHED`/`STABLE_MASS_TRANSFER`/`IMMEDIATE_MERGER`/`PHASE_NOT_MODELLED` outcomes and non-zero `L_X`/`Q_H`, and Figure 2 shows genuine separation between `non_interacting` and the RLOF-classifier prescriptions. Three `tests/test_rlof_wiring.py` tests that hand-constructed `pop.a` directly at Rsun-scale values needed updating to set AU-scale values instead (dividing by `RSUN_PER_AU`), with assertions correspondingly converted back for comparison -- the physics/outcome each test checks is unchanged, only the raw `pop.a` input/assertion values. |

## 12. Stable mass-transfer consequence model

Not an HTP02 prescription -- see Section 10's framing and the proposal
doc's "Decision" on instantaneous vs. rate-integrated treatment: HTP02
rate-integrates mass transfer via Kelvin-Helmholtz/nuclear time-scales
(eqs. 58-61), which does not fit Realta's instantaneous-event
architecture (SN and merger events are already both instantaneous
state changes at a precomputed time). Applied instantaneously instead,
at the same precomputed `rlof_time` used for the merger channel.

| Quantity | Source | Realta implementation | Status |
|---|---|---|---|
| Instantaneous conservative mass transfer to the new detachment point | Not HTP02 -- a named simplification: solves for the mass Δm (donor->companion, conservative) such that the widened orbit's Roche-lobe radius exactly equals the donor's current (unchanged-during-the-instant) radius. Orbital widening itself is standard two-body angular-momentum conservation (`a_f = a_i*(M1i*M2i/(M1f*M2f))^2`), not an HTP02-specific result. Only reachable for `donor_mass < companion_mass` (the only case `classify_rlof()` labels `STABLE_MASS_TRANSFER`) -- an initial design ("transfer until masses equalize") was caught and corrected before implementation because it moved mass in the physically wrong direction for that regime; see the proposal doc's "Decision 3" history. | `binaries/interaction.py::apply_stable_mass_transfer`, `_widened_separation` | **unit-tested, sensitivity-verified**: `tests/test_rlof_classifier.py` confirms mass conservation, correct direction (donor gets lighter, not heavier), orbit widening, and that the new Roche-lobe radius exactly matches the donor's radius at the new separation/mass-ratio (the self-consistent detachment point this function is defined to reach) -- plus a rejection test for the wrong-direction call. I verified sensitivity directly: broke `_widened_separation` to not widen the orbit, confirmed the test fails with a `brentq` bracket-sign error (no valid detachment point exists without widening), reverted. |
| Wiring into `evolve()`'s Phase 0 | n/a | `binaries/population.py::evolve` | **unit-tested, sensitivity-verified, one bug found and fixed once HG search was added**: `tests/test_rlof_wiring.py` confirms `evolve()` actually applies the transfer (mass conservation, correct direction, orbit widening) and does not reprocess an already-applied event on a later timestep. Both stars' lifetime clocks (`turnoff_time`/`t2_lifetime`) are reset from `tnow` at their new masses -- the same full-reset simplification used for the merger channel, not partial (Tout et al. 1997/Brček et al. 2026) rejuvenation. I verified sensitivity directly: removed the orbital-separation update, confirmed the widening assertion fails, reverted. Originally computed the donor's radius via `ms_radius()` unconditionally; once `find_rlof_onset` could also find a stable-MT crossing during HG, this broke (see Section 10's `find_rlof_onset` row) -- fixed by checking `main_sequence.phase()` at `rlof_time` and dispatching to `ms_radius()`/`hg_radius()` accordingly. |
| `COMMON_ENVELOPE` outcome in `evolve()`'s Phase 0 | HTP02 eqs. 69-73 (energy-balance solve, see Section 12a below) | `binaries/population.py::evolve`, `binaries/interaction.py::apply_common_envelope` | **unit-tested, sensitivity-verified, wired into both consequence branches**: `tests/test_rlof_wiring.py` confirms `evolve()` calls `apply_common_envelope` and applies the result -- survival strips the donor to its core mass, leaves the companion untouched, tightens the orbit to `a_f`, and resets both lifetime clocks; merger routes the donor's *core* mass (not its pre-CE full mass) plus the companion through the same `did_merge`/`merge_time`/lifetime-reset pathway as `IMMEDIATE_MERGER`. A separate test confirms the `config.alpha_ce`/`lambda_ce` override actually reaches `apply_common_envelope` inside `evolve()` (a low `alpha_ce` forces a merger that would otherwise survive at the default), not just that the config field itself is set. |

## 12a. Common-envelope energy-balance solve (survive vs. merge)

HTP02 Sec. 2.7.1, eqs. 69-73 (the alpha-lambda formalism) -- see the
proposal doc's "CE alpha-lambda: implementation outline" section for
the full derivation and the coalescence-check reasoning (whichever of
`a_f` or the companion/core's own Roche-lobe-filling separation `a_L`
is reached at the *larger* separation happens first during inspiral).

| Quantity | Source | Realta implementation | Status |
|---|---|---|---|
| Envelope binding energy, initial/final orbital energy, `a_f` | HTP02 eqs. 69-72 (`G` cancels exactly -- every term is linear in it, confirmed algebraically before implementing, so no unit-consistency risk) | `binaries/interaction.py::apply_common_envelope` | **unit-tested, sensitivity-verified**: `tests/test_ce_energy_balance.py` confirms the typical mid-HG-donor merge outcome, that the same code path CAN return `survives=True` for a hand-picked small-envelope-fraction case (rules out "always merges" being a structural bug), that `new_donor_mass` equals the core mass in both branches, that a higher `alpha_ce` favours survival (larger `a_f`, the physically expected direction), and that flipping `lambda_ce`'s magnitude changes `a_f` (confirms it's actually wired into the energy balance, not silently ignored). |
| Coalescence check (`a_L` via the Eggleton fit on the bare core's and companion's actual radii) | HTP02, the paragraph introducing eq. 73 | `binaries/interaction.py::apply_common_envelope` | **unit-tested**: `tests/test_ce_energy_balance.py::test_apply_common_envelope_uses_core_radius_for_the_donor_core` confirms the *compact* core radius (`remnant.core_radius`, white-dwarf-like), not the donor's pre-CE giant radius, feeds this check. |
| Scope gap: eqs. 74-77 not implemented | HTP02 eqs. 74-77 (partial-envelope-retention Newton-Raphson solve for the merged star's mass) | n/a | **documented gap, not an oversight**: needs `R_i`, "the radius the system would have if it were to coalesce immediately," which HTP02 does not define operationally in a way this module could implement without further study. The merge branch instead assumes full envelope loss (bare-core + companion conservative merge via `merge_stellar_masses`), the same simplification already used for `IMMEDIATE_MERGER`. |
| `alpha_CE`/`lambda_CE` values and their exposure as overridable config fields | Zuo & Li (2014, MNRAS 442, 1980) for `alpha_CE=0.9` (their own basic-model value, middle of their HMXB-population-calibrated 0.8-1.0 preferred range); HTP02 eq. 69 for `lambda_CE=0.5` (its own fixed value, explicitly flagged there as not a true constant) | `binaries/interaction.py::ALPHA_CE`/`LAMBDA_CE`; `config.py::SimulationConfig.alpha_ce`/`lambda_ce` (default to the module constants when unset, matching the `q_crit_ms` pattern) | **unit-tested**: `tests/test_rlof_wiring.py` confirms the config defaults match the module constants, that overriding them is honoured and validated (`> 0`), and that an override actually propagates into `evolve()`'s CE consequence, not just the config field. |
| **Finding, verified not a bug**: realistic mid-HG donors generically merge | n/a -- a numerical/physical consequence of the equations above, not a separate citation | n/a | Every tested mid-HG donor/companion combination merges (`E_bind,i` large relative to available orbital energy, since the donor still carries most of its mass in an extended envelope, forcing `a_f` tighter than the companion's own Roche-lobe-filling separation). Confirmed genuine physics, not a sign/arithmetic error, by constructing a small-envelope-fraction case that does survive (same code path) -- see the proposal doc's write-up. Matches the literature's general expectation that HG-donor CE is merger-prone (HTP02 Sec. 2.7.1 itself; population-synthesis codes such as StarTrack/COMPAS often treat HG-donor CE as a forced merger by convention). |

## 13. Random number generation

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
- ~~No automated test exercises the RLOF-classifier pipeline or
  `scripts/run_paper1_experiment.py` end to end~~ -- this was the
  actual gap that let the AU/Rsun units bug (Section 10 above) ship
  undetected: every existing `interaction.py` unit test uses hand-
  picked, self-consistent Rsun-scale separations, so none of them
  could have caught a units mismatch at the `population.py` boundary;
  it was only found by manually running the script. Closed by
  `tests/test_paper1_pipeline_regression.py`: a numeric regression pin
  (population summary stats, RLOF-outcome-distribution counts, and a
  per-timestep trajectory) for a `standard_interaction` run at
  Paper-1-relevant scale (`mcut=8`, `pmin=1.0`, matching
  `configs/paper1_basic_experiment.yml`), plus a smoke test that loads
  `scripts/run_paper1_experiment.py`'s own functions directly and
  confirms an RLOF-classifier prescription produces non-zero
  `L_X`/`Q_H` and that both figure files are actually written -- the
  specific symptom the units bug produced (all-zero `L_X`/`Q_H` via
  near-total `IMMEDIATE_MERGER`) is exactly what would fail first if
  that conversion were silently removed again.

Still open:

- **No `examples/` directory** (brief §26) walking through a full
  run end-to-end outside of the paper-reproduction notebook.
- **Figure 3 (mergers vs. compact-object formation) not built** -- only
  the minimal `did_merge`/`merge_time` event bookkeeping it will need
  is in place (Section 6). Deliberately deferred, along with the full
  `Event`/`PopulationHistory` abstraction from
  `docs/science/development-roadmap.md` item 4, per the Paper 1
  implementation prompt's stated scope.
- **Figures 4-6** (IMF/binary degeneracy grid, metallicity sweep,
  stochastic realisations) are out of scope for this milestone --
  Figure 6 in particular is explicitly the bridge to Paper 2's
  parameter-sweep/experiment-runner machinery (Phase 3), not something
  to build here.
- ~~RLOF classifier not wired into the simulation loop~~ -- now wired
  as an opt-in event (Section 10), and now reconciled with
  `standard_interaction`/`enhanced_interaction`/`enhanced_mergers`
  (Section 6's "Reconciliation" subsection) rather than left
  independent.
- ~~Stable mass-transfer consequence model not implemented~~ -- now
  implemented as an instantaneous, conservative transfer to the new
  detachment point (Section 12 below), not HTP02's rate-integrated
  eqs. 58-61 -- see `docs/science/rlof-ce-classifier-proposal.md`'s
  "Decision" on instantaneous vs. rate-integrated treatment.
- ~~`classify_rlof` MS-donor-only~~ -- now also classifies HG donors
  (Section 10 above: `COMMON_ENVELOPE` outcome, `hg_q_crit`), correctly
  distinguishing them from MS donors' `IMMEDIATE_MERGER`.
- ~~`find_rlof_onset` MS-only, HG path unreachable through the
  pipeline~~ -- now also root-finds HG radius crossings and is wired
  all the way through `evolve()` (Section 10/12 above); one real bug
  found and fixed in the process (`evolve()`'s stable-MT branch was
  unconditionally using `ms_radius()`).
- ~~CE alpha-lambda formalism (eqs. 69-73) not implemented~~ -- now
  implemented (Section 12a above): the energy-balance solve, the
  coalescence check, and the survive/merge consequence model, all
  wired into `evolve()`. Eqs. 74-77 (partial-envelope-retention
  Newton-Raphson solve for the merged mass) remain a documented,
  deliberate gap -- see Section 12a's "scope gap" row. `C` (the inverse
  of the full GB core-mass-luminosity relation) would only be needed
  for a faithful eq. 44 above `CORE_MASS_BGB_MAX_MASS` -- still
  deferred, not blocking, since CE only needs core mass/radius within
  the already-supported range.
- **Brček, Hirai, Mandel & Lower (2026) not available** -- the task
  brief names this paper for post-interaction MS core-mass/radius
  response and rejuvenation; it has not been supplied. What's actually
  implemented (Section 12) is a simpler full-reset simplification
  (both stars' lifetime clocks restart from the event time at their
  new masses), not even HTP02's own cited Tout et al. (1997) partial
  fractional-age rejuvenation -- that remains an explicit extension
  point, as does the fuller Brček et al. treatment once available.
- **Hovis-Afflerbach et al. (2025) stripped-donor properties** --
  named in the task brief as a downstream consequence, not required
  for this task. No interface stub exists yet; left as a named future
  extension point per the task's own scope note.
- **Post-SN secondary Roche-lobe overflow not modelled (found
  2026-08-24, see Section 6's "Residual, understood degeneracy" note
  above)** -- Phase 0's RLOF classifier only checks pre-SN interaction
  between two still-evolving stars (`nturn==0`), matching HTP02's own
  CE-eligible donor list. Because `generate_population` always enforces
  `m2 <= m1`, and `STABLE_MASS_TRANSFER` requires the donor to be the
  *lighter* star, that channel's donor is always `m2` -- whose own
  pre-SN evolutionary timescale is always far longer than `m1`'s
  (the only star whose SN Realta's `turnoff_time` tracks), so `m1`
  has essentially always already exploded before this pre-SN channel
  becomes reachable. This makes `STABLE_MASS_TRANSFER`, and by
  extension `interaction_boost`, structurally unable to fire for any
  realistic massive-star binary in the current design -- confirmed,
  not just suspected: zero of 51 `STABLE_MASS_TRANSFER`-classified
  binaries in the Paper 1 config ever got processed, across the whole
  100 Myr run.

  **This needs to be implemented properly once work moves to mass
  scales and timescales where it is the relevant physics**: the
  astrophysically dominant real HMXB-formation channel is the
  secondary's *own*, later Roche-lobe overflow onto the by-then-compact
  primary (Case B/C mass transfer onto a neutron star/black hole,
  wind-accretion vs. RLOF) -- not pre-SN mass transfer between two
  still-live stars. That needs a genuinely new channel: donor = the
  still-evolving secondary, companion = a compact remnant (not another
  main-sequence/HG star), gated on `nturn==1` rather than `nturn==0`,
  with its own Roche-lobe/accretion treatment (a compact object has no
  stellar radius to overflow *from*, and accretion onto it is a
  different physical regime from CE onto a live companion). Not
  attempted this session -- deliberately out of scope for the pre-SN
  classifier work done here, and only worth building when a study's
  mass range/timescale actually needs it (e.g. lower-mass secondaries
  over longer baselines, where this channel would dominate over what's
  modelled today).
