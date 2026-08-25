# RLOF outcome classifier

Classifies whether Roche-lobe overflow occurs and its outcome, given a
binary's masses, separation, metallicity, and the donor's age. MS
donors only ([`stellar-tracks.md`](stellar-tracks.md)'s scope note);
HG+ donors beyond that scope return `PHASE_NOT_MODELLED`.

Wired in as an opt-in event (`config.use_rlof_classifier`, default
`False`) — a new "Phase 0" MS/HG-RLOF event in `evolve()`, alongside
the existing `fsur`/`interaction_boost` gate. The onset time is
precomputed once per binary in `generate_population()` via
root-finding (`interaction.py::find_rlof_onset`), mirroring
`turnoff_time`'s precomputation from `LifetimeTable` — this keeps
`evolve()`'s per-timestep cost independent of the classifier.

Internal separations are stored in **AU** (`binaries/population.py`);
every function in this module expects **Rsun**. Call sites convert via
`BinaryPopulation.RSUN_PER_AU` (`= 215.032`, from the IAU 2012 AU and
IAU 2015 nominal solar radius).

## Roche-lobe radius

`R_L1(a, q1)` — Eggleton (1983) fit, as used by Hurley, Tout & Pols
(2002, MNRAS 329, 897), eq. 53.

Implementation: `binaries/interaction.py::roche_lobe_radius`.
Tests: `tests/test_rlof_classifier.py` — q=1 literature value
(`R_L1/a ≈ 0.379`), separation scaling, non-positive-mass-ratio
rejection.

## MS-donor classification

Detached / stable mass transfer / immediate merger. `q_crit_ms =
0.695` for dynamically-unstable MS donors (HTP02 Sec. 2.6.4, k1=0
extended to k1=1 as a named simplification). MS donors are not
CE-eligible (HTP02 Sec. 2.7.1), so dynamical instability merges
directly rather than forming a common envelope.

Implementation: `binaries/interaction.py::classify_rlof`.
Tests: `tests/test_rlof_classifier.py` — all three outcomes, the
`PHASE_NOT_MODELLED` guard past `t_MS`, `q_crit_ms` override.

## HG-donor classification

Detached / stable mass transfer / common envelope. GB `q_crit`
(HTP02 eqs. 56-57) reused for HG donors — the reuse follows Zuo & Li
(2014, MNRAS 442, 1980) eq. 1 (citing Shao & Li, in prep.); HTP02's
own HG treatment is a fixed `q_crit = 4`. HG donors ARE CE-eligible
(HTP02 Sec. 2.7.1), unlike MS donors.

Implementation: `binaries/interaction.py::hg_q_crit`, `classify_rlof`.
Tests: `tests/test_hg_ce_classifier.py`.

## Immediate-merger mass combination

HTP02 gives an explicit envelope mass-loss prescription only for
common-envelope mergers (eqs. 69-77), not for a direct MS-MS
collision. Conservative merging (no mass loss) is used as the default
absent an explicit prescription.

Implementation: `binaries/interaction.py::merge_stellar_masses`.

## RLOF onset time

Root-finds `R_donor(t) = R_L1` once per binary via bisection, using
`ms_radius()`'s MS-wide monotonicity and `hg_radius()`'s HG-wide
monotonicity (continuous at the MS/HG boundary).

Implementation: `binaries/interaction.py::find_rlof_onset`.
Tests: `tests/test_rlof_classifier.py`, `tests/test_hg_ce_classifier.py`.

### Donor-selection property

Eggleton's `R_L1/a` increases monotonically with the donor's own mass
ratio, so the automatically-selected donor is almost always the more
massive star (`q1 > 1`). Consequences:

- `IMMEDIATE_MERGER` dominates over `STABLE_MASS_TRANSFER` for
  MS donors in practice.
- `q_crit_ms` has limited leverage on the outcome distribution — see
  [`interaction-prescriptions.md`](interaction-prescriptions.md)'s
  `STABLE_MASS_TRANSFER` structural-limit note.
- Cross-checked against Xu et al. (2025, A&A 704, A218,
  arXiv:2503.23876)'s SMC statistics (`M1=5-100` Msun, `q=0.3-0.95`,
  `P=1-3162` d: 8% post-mass-transfer, 7% merger observed). Realta
  gives 1.6% post-mass-transfer / 29.5% merger for the same selection
  window (`scripts/xu2025_smc_crosscheck.py`) — merger-heavy relative
  to observations, traced to this donor-selection property (for that
  `q` window, `q1 = M_donor/M_companion` is always > 1, comfortably
  above `Q_CRIT_MS`). `~22%` of the sample falls to
  `PHASE_NOT_MODELLED` (GB/AGB donors beyond MS/HG scope). Not
  compensated for by adjusting `Q_CRIT_MS`/`ALPHA_CE`/`LAMBDA_CE`
  against this one external comparison.

## Wiring

`config.py::use_rlof_classifier`/`q_crit_ms`,
`binaries/population.py::generate_population`/`evolve`. Requires
`imetal = 2` or `3` (Hurley/Tout formulae are undefined at Z=0); a run
with `use_rlof_classifier=True` and `imetal=1` logs a warning and
skips RLOF classification.

Tests: `tests/test_rlof_wiring.py`.
