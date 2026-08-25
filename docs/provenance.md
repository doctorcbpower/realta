# Provenance: Power et al. → Realta

Traceability for the ported baseline: paper equation/assumption →
Python implementation → test coverage. Physics added beyond the
Power et al. papers is documented separately under `docs/physics/`
(see `docs/physics/README.md` for an index) and in
`docs/known-gaps.md` for what remains unimplemented.

## References

- **Power, C., Wynn, G. A., Combet, C., Wilkinson, M. I. (2009)**,
  "Primordial globular clusters, X-ray binaries and cosmological
  reionization", *MNRAS*, 395, 2, 1146-1152.
  [arXiv:0902.1897](https://arxiv.org/abs/0902.1897) ·
  [DOI: 10.1111/j.1365-2966.2009.14628.x](https://doi.org/10.1111/j.1365-2966.2009.14628.x) ·
  [ADS](https://ui.adsabs.harvard.edu/abs/2009MNRAS.395.1146P/abstract)

  Primary reference for the population model: coeval star formation,
  IMF choices, the sudden-mass-loss survival criterion, `f_sur`, and
  the peaked/Eddington-limited X-ray luminosity distribution (Sec.
  2.1-2.2). Realta is a Python port of the Fortran Monte Carlo code
  this paper's results were produced with.

- **Power, C., James, G., Combet, C., Wynn, G. (2013)**, "Feedback from
  High-mass X-ray Binaries on the High-redshift Intergalactic Medium:
  Model Spectra", *ApJ*, 764, 1, 76.
  [arXiv:1211.5854](https://arxiv.org/abs/1211.5854) ·
  [DOI: 10.1088/0004-637X/764/1/76](https://doi.org/10.1088/0004-637X/764/1/76) ·
  [ADS](https://ui.adsabs.harvard.edu/abs/2013ApJ...764...76P/abstract)

  Source of the X-ray-to-ionising-photon conversion (Section 3):
  model HMXB spectra (black-body + power-law accretion states) and
  ionising-photon output, from the same Monte Carlo population model.

Status legend: `pinned` = value(s) pinned in
`tests/test_regression.py::test_reference_cluster_run_pinned_trajectory`
(`lumx_tot`/`nphot_tot`/`nactive`/`ndead` vs. time, `total_mass_msun`,
binary count, for two fixed config+seed cases); `unit` = exercised in
isolation by a dedicated unit test; `phase` = exercised by
`tests/test_evolve.py`, which constructs a `BinaryPopulation` and sets
internal state directly to isolate one `evolve()` phase.

---

## 1. Population generation

| Reference | Implementation | Tests |
|---|---|---|
| IMF sampling — Salpeter, Kroupa, Chabrier (Sec. 2.2) | `imf/salpeter.py`, `imf/kroupa.py`, `imf/chabrier.py`, `imf/factory.py::get_imf` | `tests/test_imf.py` (CDF boundary, monotonicity, sampling bounds). Kroupa `n_massive`/`total_mass_msun`: pinned. `config.mmin = 0.1` Msun (practical hydrogen-burning lower bound). |
| Sec. 2.1: every star with M ≥ `mcut` is placed in a binary at formation | `binaries/population.py::generate_population` | pinned (`n_massive`) |
| Companion mass, flat between `mcomp` and `m1` (Sec. 2.1) | `binaries/population.py::generate_population` | pinned (`m2` feeds `evolve()` Phase 1) |
| Orbital period, log-flat between `pmin`/`pmax` (Sec. 2.1) | `binaries/population.py::generate_population` | `tests/test_population_generation.py::test_orbital_period_bounds_and_log_uniform_distribution` |
| Semi-major axis, Kepler's third law (Sec. 2.1) | `binaries/population.py::generate_population` (`AFAC`) — internal units **AU** | `tests/test_population_generation.py::test_semi_major_axis_matches_reference_formula`. `AFAC = ` Kepler's third law in the AU/Msun/year convention (`PFAC` is a sidereal year in days). |
| MS lifetime vs. mass, tabulated (Schaerer et al. 1993) | `io/tables.py::LifetimeTable`, `data/lifetimes_z*.dat` | pinned (gates the SN1 trigger) |

## 2. Supernova transitions and HMXB activation

| Reference | Implementation | Tests |
|---|---|---|
| Sec. 2.1: sudden-mass-loss survival criterion, binary stays bound iff `floss <= 0.5` | `binaries/population.py::evolve`, Phase 1 | `phase`, pinned |
| Sec. 2.1: `f_sur`, probability a surviving binary becomes an active HMXB | `config.py::SimulationConfig.fsur`, `binaries/population.py::evolve` Phase 1 | `phase` (`test_evolve_phase1_fsur_partial_activation`), pinned (`fsur=1.0`, `fsur=0.5`) |
| Remnant mass vs. progenitor mass, tabulated | `io/tables.py::RemnantTable`, `data/remnant_masses.dat` | pinned (feeds `floss`) |
| Secondary SN → system marked dead | `binaries/population.py::evolve`, Phase 2 | `phase` (`test_evolve_phase2_death`), pinned (`ndead`) |

## 3. X-ray luminosity and ionising photons

| Reference | Implementation | Tests |
|---|---|---|
| Per-HMXB X-ray luminosity draw, peaked/Weibull-shaped, Eddington-limited, `L_X^peak ≈ 1e38 erg/s` (Sec. 2.2) | `xray/luminosity.py::XRayLuminosity` (`distribution="weibull"`, default) | pinned |
| Flat log-uniform alternative draw (not in Power et al.; offered as a comparison mode) | `xray/luminosity.py::XRayLuminosity` (`distribution="uniform"`), `config.py::xray_distribution` | `tests/test_xray_luminosity.py`. Not Eddington-limited, unlike `"weibull"`. |
| X-ray luminosity persists once drawn, until SN2 | `binaries/population.py::evolve` — drawn once at SN1 (Phase 1), zeroed at SN2 death, never redrawn | `phase` |
| Ionising photon rate from X-ray luminosity, via a model spectral shape (2013, Cygnus X-1-like black-body + power-law accretion states) | `binaries/population.py::NPHOT_PER_LUMX` | pinned (`nphot_tot` is a fixed multiple of `lumx_tot`) |

## 4. Random number generation

| Reference | Implementation | Tests |
|---|---|---|
| Single seeded stream drives the whole simulation | `binaries/population.py::BinaryPopulation.np_rng` (`np.random.default_rng(config.iseed)`), threaded through population generation, SN survival, HMXB activation, and `XRayLuminosity.get_lumx(rng=...)` | `tests/test_regression.py::test_run_is_deterministic_for_fixed_seed` |
