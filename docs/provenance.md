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
far as manual verification goes, but has no automated regression test
pinning its numeric output.

---

## 1. Population generation

| Paper / assumption | Reference Fortran | Realta implementation | Status |
|---|---|---|---|
| IMF sampling, all three forms (Salpeter/Kroupa/Chabrier) | `make_stars.f` (driver), `salpeter.f`, `kroupa.f`, `log_normal_IMF.f` | `imf/salpeter.py`, `imf/kroupa.py`, `imf/chabrier.py`, `imf/factory.py::get_imf` | Salpeter, Chabrier: done, **untested** (no test file). Kroupa: done but **flagged** — see `imf/kroupa.py` class docstring for an `mmin`-handling discrepancy that only manifests if `config.mmin != 0.01`. Tested: `tests/test_imf.py` (Kroupa CDF + sampling bounds only). |
| Sec. 2.1: all massive stars (M ≥ `mcut`) are in binaries at formation (`fpbin=1.0`) | `make_stars.f` | `binaries/population.py::generate_population` — every `m1 >= mcut` star unconditionally gets `m2`/`period`/`a` | done, untested |
| Companion mass distribution, flat between `mcmpct` and `m1` | `make_stars.f` | `binaries/population.py::generate_population` (`m2 = cfg.mcomp + (m1 - cfg.mcomp) * rng.random()`, clipped to `m1`) | done, untested |
| Orbital period, log-flat between `pmin`/`pmax` | `make_stars.f` | `binaries/population.py::generate_population` | done, untested |
| Semi-major axis from Kepler's third law | `make_stars.f` | `binaries/population.py::generate_population` (`AFAC` constant) | done, untested |
| MS lifetime vs. mass, tabulated (Schaerer et al. 1993-derived) | `lifetime.f`, reads `lifetimes_z*.dat` | `io/tables.py::LifetimeTable`, same bundled `.dat` files | done, untested |

## 2. Supernova transitions and HMXB activation

| Paper / assumption | Reference Fortran | Realta implementation | Status |
|---|---|---|---|
| Sec. 2.1: sudden-mass-loss survival criterion (binary stays bound iff `floss <= 0.5`) | `main.f` | `binaries/population.py::evolve`, Phase 1 | done, untested |
| Sec. 2.1: `f_sur` — probability a surviving binary becomes an active HMXB | `main.f` (`ran3(iseed).le.fbin`) | `config.py::SimulationConfig.fsur`, applied in `binaries/population.py::evolve` Phase 1. **Renamed** from the reference's confusingly-named `fbin` — see `config.py` field docstring for the full naming history/rationale (this session's own fix, not a paper discrepancy). | done, untested |
| Remnant mass vs. progenitor mass, tabulated | `get_mremnant.f`, reads remnant mass table | `io/tables.py::RemnantTable`, `data/remnant_masses.dat` | done, untested |
| Secondary SN → system marked dead | `main.f` | `binaries/population.py::evolve`, Phase 2 | done, untested |

## 3. X-ray luminosity and ionising photons

| Paper / assumption | Reference Fortran | Realta implementation | Status |
|---|---|---|---|
| Per-HMXB X-ray luminosity draw, Weibull-like/peaked, Eddington-limited, `L_X^peak ≈ 1e38 erg/s` | `get_lumx.f` | `xray/luminosity.py::XRayLuminosity` (`distribution="weibull"`, Realta's default) | done, untested — this was the site of this session's RNG-seeding bug fix (see `xray/luminosity.py` module docstring) |
| `get_lumx.f`'s dead "uniform" branch (only reachable if `iseed == -1`, a debug sentinel never hit by `main.f`) | `get_lumx.f` | `xray/luminosity.py::XRayLuminosity` (`distribution="uniform"`), `config.py::xray_distribution` | done, but **not the reference default** — Realta exposes it as an opt-in config choice rather than reproducing it as reachable-but-dead code |
| X-ray luminosity persists once assigned, until SN2 death | (implicit in `main.f`'s single-draw-per-lifetime structure) | `binaries/population.py::evolve` — L_X drawn once at SN1 in Phase 1, zeroed at SN2 death, **not** redrawn per timestep | done, untested — this was one of this session's four originally-diagnosed porting bugs |
| Ionising photon rate from X-ray luminosity (spectral-shape conversion, 13.6 eV–1500 eV vs. 1e6 eV reference) | `main.f` (`lum_xray(i)*(6.2415e11/13.6)*alog(1500/13.6)/alog(1e6/13.6)`) | `binaries/population.py::NPHOT_PER_LUMX` class constant | done, untested |
| `get_ngamma.f`'s MS-based ionising photon estimate | `get_ngamma.f` | **not ported** — dead in the reference too as far as this port is concerned; do not confuse with `IonizingPhotonTable`, which is a different (M≥8 Msun ionising-budget) quantity, also not currently wired into `evolve()` | not applicable |

## 4. Fig. 1 reproduction (this session's own addition, not in the original Fortran)

| Quantity | Source | Realta implementation | Status |
|---|---|---|---|
| Main-sequence bolometric luminosity vs. time | Not in `gc_hmxbs` — the original paper's Fig. 1 MS curve was almost certainly generated with a population-synthesis code of the day (Starburst99 being the standard circa 2009), external to the ported Fortran | `io/tables.py::MSLuminosityTable`, sourced from FSPS (`notebooks_helper/generate_ms_luminosity_table.py`) | done, untested, **and flagged**: FSPS vs. Starburst99 track differences (Padova vs. Geneva stellar evolution) are a documented literature-level systematic of up to a factor of a few in this exact age window (3–20 Myr) — see this session's SB99-vs-FSPS research summary. Not a bug, but worth stating explicitly wherever this reproduction is presented. |
| MS/HMXB luminosity mass-normalization | n/a (this session's own bug, introduced when the FSPS table was added, not present in the reference) | `binaries/population.py::BinaryPopulation.total_mass_msun`, `io/tables.py::MSLuminosityTable.get_lbol(age_myr, total_mass_msun)` | done, untested (no regression test), verified manually via notebook re-execution |

## 5. Random number generation

| Paper / assumption | Reference Fortran | Realta implementation | Status |
|---|---|---|---|
| Single seeded stream drives the whole simulation (`ran3`, Numerical Recipes) | `ran3.f`, seeded from `iseed` in `main.f` | `binaries/population.py::BinaryPopulation.np_rng` (`np.random.default_rng(config.iseed)`), threaded through population generation, SN-survival, HMXB activation, **and** `XRayLuminosity.get_lumx(rng=...)` | done — this was this session's other RNG bug fix (X-ray draws previously used the unseeded global `np.random`, breaking reproducibility) |

---

## Known gaps (not yet closed)

- **No numeric regression test.** Nothing above has automated coverage
  pinning actual `lumx_tot`/`nphot_tot`/`nactive`/`ndead` trajectories
  for a fixed seed — `tests/test_regression.py` only checks the
  simulation doesn't crash or produce NaNs. This is the single highest-
  value gap: without it, none of the "done" rows above are protected
  from silent drift by a future change. (Audit §32 item 5 / suggested
  next steps item 6.)
- **Salpeter and Chabrier IMFs are untested** — only Kroupa has a test
  (`tests/test_imf.py`), despite all three being equally reachable via
  `config.imf_type`.
- **The Kroupa `mmin`-handling discrepancy** (§1 above) is flagged, not
  resolved — needs a decision on which behaviour (Fortran's hardcoded
  0.01, or this port's caller-supplied `mmin`) is actually intended
  before `mmin` is treated as a free parameter.
- **No `examples/` directory** (brief §26) walking through a full
  run end-to-end outside of the paper-reproduction notebook.
