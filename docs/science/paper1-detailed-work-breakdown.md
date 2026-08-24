# Paper 1 detailed work breakdown

A detailed, ordered implementation prompt for Paper 1's remaining work,
organized into workstreams by actual dependency rather than as one flat
sequential list. Written so it can be pasted into a fresh session with
no other context. This extends `docs/science/paper1-implementation-prompt.md`
(the original scoping prompt) with a concrete breakdown once the
codebase's actual gaps had been checked directly rather than assumed --
see that file for the original framing, and
`docs/science/paper1-implementation-prompt.md`'s Item 2 for the
interaction-model design question this breakdown's Workstream B answers.

---

You're working on Realta, a modular stellar/binary population-synthesis
Python framework (repo already checked out) built on a ~20-year-old
Fortran Monte Carlo model of HMXBs in globular clusters (Power et al.
2009, MNRAS 395, 1146, arXiv:0902.1897). Before touching any code, read:

- docs/provenance.md -- paper-equation -> implementation -> test
  traceability, and the pinned-regression discipline this project
  follows for every physics change.
- docs/science/research-programme.md, the "Paper 1" section -- the
  scientific target all of this serves (working title: "The X-ray
  fingerprints of massive-star multiplicity").
- docs/science/development-roadmap.md -- target architecture; this task
  implements a deliberately narrow slice of it.
- docs/science/paper1-implementation-prompt.md -- the original Paper 1
  scoping prompt. This task supersedes/extends it with a concrete,
  ordered work breakdown.

Governing principles (a "Development and Scientific Software Brief"
applies throughout): preserve the Power et al. (2009) baseline exactly;
never silently change scientific behaviour; flag ambiguity rather than
resolve it yourself; explain any nontrivial physics or design decision
before implementing it; avoid over-engineering; keep the public API
small; treat external tools (MESA, FSPS, etc.) as optional; never
commit or push without being explicitly asked; prefer small,
incremental, reviewable changes.

## Current state (verified 2026-08-24, not assumed -- re-check before starting)

**This section was substantially out of date as of 2026-08-24** -- it
was written before an intervening session implemented most of
Workstream B and part of Workstream C below. Re-verified directly
against the codebase (not from memory) on that date; see
`docs/provenance.md` Sections 6, 9, 10, 12, 12a and
`docs/science/rlof-ce-classifier-proposal.md` for the full
implementation history this section summarizes.

**Done, beyond the original "25 passing tests" baseline**:

- Workstream B is essentially complete. B1 (Hurley, Pols & Tout 2000
  radius/phase module: `src/realta/stellar/{zams,main_sequence,
  giant_branch,remnant}.py`, MS+HG scope, a documented, deliberate
  reduction from full GB/AGB) is implemented and unit-tested; its own
  "flag before implementing" question (does Hurley's L(t) replace the
  FSPS-sourced `L_bol`, or stay internal-only?) was already resolved
  and documented as "internal only" (`rlof-ce-classifier-proposal.md`
  "Decision 1"). B2 (interaction outcome classifier:
  `binaries/interaction.py::classify_rlof`/`hg_q_crit`, HTP02 eqs.
  56-57/2.6.4/2.7.1, `Q_CRIT_MS`/`ALPHA_CE`/`LAMBDA_CE` as named,
  citable, overridable constants/config fields) is implemented and
  unit-tested. B3 (consequence model:
  `apply_stable_mass_transfer`/`apply_common_envelope`, HTP02 eqs.
  69-73) is implemented for the stable-MT and CE-survive/merge
  outcomes. B4 (binary-prescription variants:
  `standard_interaction`/`enhanced_interaction`/`enhanced_mergers` via
  `config.py`'s `_PRESCRIPTION_DEFAULTS`) is implemented and wired into
  `evolve()`. B5's minimal event bookkeeping (`did_merge`/`merge_time`)
  exists, as B5 itself specifies (not the full Event taxonomy).
- Most of C1: a single YAML config
  (`configs/paper1_basic_experiment.yml` +
  `scripts/run_paper1_experiment.py`) produces Figure 1 and Figure 2
  in one reproducible invocation, now numerically regression-pinned
  (`tests/test_paper1_pipeline_regression.py`) and confirmed
  non-degenerate after a real AU/Rsun units bug (found by actually
  running the pipeline end-to-end, not by any unit test -- see
  `docs/provenance.md` Section 6/10) was found and fixed.

**Update, 2026-08-24**: A1 and A4 are now done -- see their own
sections below for the implementation writeup, and
`docs/provenance.md` Section 1 for the full traceability. `"single"`/
`"non_interacting"` still come from the discrete `binary_prescription`
enum, unchanged -- A1's `binary_fraction` generalizes rather than
replaces that mechanism (a `binary_fraction<1` star is still tracked
through `m1`, just with placeholder `m2=0`/`period=0`/`a=0`, unlike
`"single"`'s own array-emptying shortcut).

**Update, 2026-08-24**: A2 and A3 are now also done -- Workstream A is
complete. See their own sections above for the implementation
writeup, and `docs/provenance.md` Section 4 for the full
traceability.

**Still not done**:

- B2's Xu et al. (2025) SMC-statistics cross-check: not done -- the
  classifier's population-level outcome fractions have never been
  compared against that paper's 8%-post-MT/7%-merger figures.
- B3's Brček, Hirai, Mandel & Lower (2026) rejuvenation and
  Hovis-Afflerbach et al. (2025) stripped-donor properties: not
  implemented (as this document itself allows -- "leave an explicit
  extension point ... even if not implemented now"). What exists
  instead is a simpler full lifetime-clock-reset simplification,
  explicitly documented as not the Brček treatment.
- C2 (Figure 4, IMF/binary grid): not started -- blocked on A1 and A4.
- C3 (Figure 3, mergers): not started -- blocked on B5's figure-side
  work (the bookkeeping exists, the figure does not).
- Figure 5 (metallicity sweep): not built as an actual figure yet,
  though the underlying 3-preset `imetal` system this item relies on
  is unchanged/available.

**New finding, not on this document's original list**: the same
session that closed out B1-B4/C1 also found that
`STABLE_MASS_TRANSFER`/`interaction_boost` cannot structurally fire
for the mass regime Paper 1's basic experiment covers -- the donor for
that channel is always the lighter, much-slower-evolving companion
(`generate_population` enforces `m2 <= m1`), so the heavier primary
has essentially always already exploded before the channel becomes
reachable (confirmed: zero of 51 `STABLE_MASS_TRANSFER`-classified
binaries in the pinned Paper 1 config were ever processed across a
full 100 Myr run). This is not covered by A1-C3 above -- see
`docs/science/paper1-followup-prompt.md` for the full writeup and a
proposed scope (a post-SN secondary-Roche-lobe-overflow channel onto
the by-then-compact primary). That document and this one currently
overlap/aren't reconciled with each other; worth resolving into a
single source of truth before further planning.

## Workstream A -- independent of the interaction-model decisions

Can proceed now, in any order, without waiting on the Hurley/Tout/
Pols mass-transfer-stability and CE alpha/lambda parameters.

### A1. Binary sampling distributions (roadmap item 7) -- DONE (2026-08-24)

Make binary fraction, mass-ratio distribution, and period
distribution independently configurable, rather than the current
hard-wired 100%-fraction / uniform-companion-mass / log-uniform-
period scheme. This alone unlocks two of Paper 1's five basic-
experiment variants for free: "single-star populations"
(binary_fraction=0) and "non-interacting binaries"
(binary_fraction>0, no interaction model applied -- i.e. today's
existing default behaviour, just made explicit and selectable
rather than being the only option). Extend config.py /
config.yml accordingly. Pin regression values for at least one
non-default binary_fraction and one non-default period/mass-ratio
choice.

Implemented as three new `SimulationConfig` fields:
`binary_fraction` (default `1.0`, a per-star Bernoulli draw skipped
entirely at the default so the RNG stream/baseline is untouched -- see
below for why this generalizes, rather than replaces, the `"single"`
prescription's own mechanism), `mass_ratio_distribution` (`"uniform"`
default / `"flat_q"`), `period_distribution` (`"log_uniform"` default
/ `"log_normal"`, generic pmin/pmax-derived parameters, per the
2026-08-24 chat decision -- no literature source). See
`docs/provenance.md` Section 1 for the full writeup, including two
sensitivity checks (the RNG-skip claim, and a real divide-by-zero risk
in the RLOF classifier for no-companion stars, both found and guarded
against while implementing this).

### A2. L_UV(t) observable -- DONE (2026-08-24)

Wire an FSPS UV band into the actual simulation run loop (not
just the notebook). FLAG BEFORE IMPLEMENTING: confirm the specific
band (FUV? NUV? which FSPS filter?) -- this directly determines
Fig 1 and Fig 2 (L_X/L_UV is the paper's central quantity), so
get the definition confirmed rather than picking one
unilaterally. Extend MSLuminosityTable or add a parallel table as
appropriate; wire into ClusterSimulation.run()'s per-timestep
output alongside lumx_tot/nphot_tot.

Band question was already resolved earlier this session (GALEX FUV --
see `docs/provenance.md` Section 7). What was still missing: `L_bol`/
`L_UV` were only ever computed *after* a run completed, in
`scripts/run_paper1_experiment.py`, not inside `ClusterSimulation.run()`
itself. Now computed in `run()`'s own per-timestep loop
(`lbol_tot`/`luv_tot` keys), alongside `lumx_tot`/`nphot_tot`; the
script simplified to read them from `results` directly instead of
recomputing. See `docs/provenance.md` Section 4 for the full writeup.

### A3. Independent Q_H(t) -- DONE (2026-08-24)

Replace the current fixed-constant-times-L_X placeholder with an
actual ionizing-photon-rate calculation from the massive-star
population (src/realta/io/tables.py already has an unused
IonizingPhotonTable -- check whether it's suitable, or whether it
needs a genuine per-star Q_H(m, age) source). This is required for
Q4 ("which observables best distinguish IMF and binary physics")
and R_Q(t) to be meaningful -- right now Q_H is degenerate with
L_X by construction. Not strictly required for Fig 2 (L_X/L_UV
only), but shouldn't be left long past the MVP.

`IonizingPhotonTable` turned out to be suitable: `get_ngamma(m)` is a
*total* integrated photon count over the star's whole MS lifetime
(confirmed by its own MUNIT/MATOM baryon-count conversion, and by a
literature sanity check against Vacca, Garmany & Shull 1996 before
adopting that interpretation), divided by `LifetimeTable.get_lifetime(m)`
to get a genuine per-star Q_H(m) rate. `ClusterSimulation._qh_ms_tot`
sums this over currently-alive M >= 8 Msun stars each timestep;
`qh_tot = qh_ms_tot + nphot_tot` in `run()`'s results (`nphot_tot`
itself, the pre-existing HMXB accretion proxy, is unchanged). Found
and fixed one real, newly-exposed gap along the way: the `"single"`
prescription used to empty `BinaryPopulation.m1` entirely, which
silently zeroed `Q_H` for single-star populations (harmless before,
since only `L_X`/HMXB quantities read `m1`) -- migrated onto A1's
`has_companion=False` mechanism instead, fixing it. See
`docs/provenance.md` Section 4 for the full writeup.

### A4. Continuous IMF slope parameter -- DONE (2026-08-24)

Extend the IMF interface (src/realta/imf/) so the power-law slope
is a free continuous parameter rather than a choice between three
named IMFs. Required for Figure 4's (alpha_IMF, f_bin) degeneracy
grid. Keep the existing Kroupa/Salpeter/Chabrier presets working
unchanged -- this is an additive interface extension, not a
replacement.

Implemented as `config.imf_slope: float | None`, overriding
`SalpeterIMF`'s own `alpha` (default `2.35`) via
`imf/factory.py::get_imf`'s new `slope` parameter -- Salpeter-only,
deliberately not extended to Kroupa (no single slope to sweep given
its multi-segment break structure) or Chabrier (no power-law slope at
all). See `docs/provenance.md` Section 1 for the full writeup,
including a pre-existing `SalpeterIMF.cdf()` singularity at
`alpha=1.0` found and guarded against while adding this.

Regression-test and provenance-document each of A1-A4 as it lands,
per the existing discipline -- do not batch them into one untested
change.

## Workstream B -- the interaction-physics stack

Sequential; each step depends on the previous one. Step B2's specific
parameter choices are still being scoped by the user (a conversation
with Jarrod Hurley is in progress) -- implement using the
literature-standard Hurley, Tout & Pols (2002) default values,
explicitly labeled and cited as provisional/reviewable, rather than
blocking on that conversation finishing. Do not treat those defaults
as final without confirmation.

### B1. Radius/phase module

Implement Hurley, Pols & Tout (2000, MNRAS 315, 543) -- the SSE
analytic fitting formulae for R(t), L(t), core mass, and
evolutionary phase as a function of (mass, metallicity, age).
Self-contained, dependency-free (closed-form polynomials).

FLAG BEFORE IMPLEMENTING: Realta's existing L_bol comes from an
FSPS-sourced table, a different stellar-model family from
Hurley's fits. Decide and document explicitly whether Hurley's
L(t) is used only internally (for radius/phase bookkeeping and
RLOF timing) while FSPS remains the source of the reported
L_bol/L_UV (from A2), or whether this migrates the primary
luminosity source too. Record the decision in docs/provenance.md
-- do not let the two silently diverge.

### B2. Interaction outcome classifier

Given (M1, M2, a, Z) at the moment B1's R(t) meets the donor's
Roche lobe radius (Eggleton 1983 fit), classify: no interaction /
stable mass transfer / common envelope survives / common envelope
merges / immediate contact merger. Use Hurley, Tout & Pols (2002,
MNRAS 329, 897)'s mass-transfer-stability criteria and CE
alpha-lambda energetics as the baseline.

Cite the specific q_crit values and alpha/lambda choices used, as
a named and clearly revisable constant set (e.g. a small
dataclass or config block), not hard-coded numbers buried in
logic. Cross-check the resulting population-level outcome
fractions against Xu et al. (2025, A&A 704, A218,
arXiv:2503.23876)'s SMC statistics (8% post-mass-transfer, 7%
merger, for M1=5-100 Msun, q=0.3-0.95, P=1-3162 d) as a sanity
check for a comparable setup -- not an exact-match requirement.

### B3. Consequence model

For stable-MT and CE-survive outcomes, update (M1, M2, a) via
mass/angular-momentum conservation and the alpha-lambda CE
formalism. For MS mass gainers specifically, apply Brček, Hirai,
Mandel & Lower (2026, ApJ 1002, 78, arXiv:2512.13838) -- as
implemented in COMPAS -- for core-mass/radius response and
rejuvenation, rather than restarting B1's single-star track from
scratch. Leave an explicit extension point for stripped-donor
properties (Hovis-Afflerbach et al. 2025, A&A 697, A239,
arXiv:2412.05356) even if not implemented now.

### B4. Binary-prescription variants

Parameterize B2-B3 into the paper's actual named comparison
variants: "standard binary interaction," "enhanced interaction,"
"enhanced massive-star mergers" (the two variants that don't need
any of this -- single-star, non-interacting -- come from A1
instead). This is a distinct integration/config step on top of
B1-B3, not automatic once they exist. Likely implemented as
tunable knobs on B2/B3 (e.g. scaled interaction/merger
probabilities) rather than entirely separate physics per variant.

### B5. Merger/event tracking

Log B3's merger outcomes as discrete events (minimal -- do not
build the full Event taxonomy from development-roadmap.md item 4,
just enough to drive Figure 3: luminosity evolution and
compact-object formation compared across no-mergers / standard /
enhanced-mergers).

## Workstream C -- integration

Depends on enough of A and B being done to be meaningful.

### C1. YAML experiment config + Figure 1 / Figure 2

One YAML config selecting IMF (incl. A4's slope), binary
prescription (A1 + B4's variants), and metallicity, producing
Figure 1 (L_bol/L_UV/Q_H/L_X vs t for several binary models) and
Figure 2 (L_X/L_UV vs age across prescriptions -- the paper's
central figure) in one reproducible invocation. This is
achievable once A1, A2, and at minimum one real interaction
variant from B4 exist -- it does not need A3/A4/B5.

### C2. Figure 4 (IMF vs binary degeneracy grid)

Needs A4 (continuous IMF slope) and A1 (binary fraction) -- does
not need workstream B at all, since the grid is over
(alpha_IMF, f_bin), not interaction prescription. Could ship
before B is finished.

### C3. Figure 3 (effect of mergers)

Needs the full B stack including B5's event tracking.

Figure 5 (metallicity) is close to achievable now with the existing
3-preset imetal system once C1 exists -- treat as a cheap add-on to
C1, not a separate block, unless a smoother Z sweep is wanted later.
Figure 6 (stochastic realisations) is explicitly out of scope for
this task -- it's the stated bridge to Paper 2 in
docs/science/research-programme.md.

## Testing and documentation discipline (applies to every block above)

Follow tests/test_regression.py and tests/test_evolve.py's existing
pattern: pin exact values for every new physics path, verify
sensitivity by deliberately breaking the code and confirming the test
fails, then revert. Extend docs/provenance.md with new rows citing
the specific paper/section/equation for whatever gets implemented --
author-year is not enough, cite the actual criteria/formula used. The
existing fsur-based Power et al. (2009) path (imf_type=2, mmin=0.1,
etc.) must remain bit-identical by default -- the 25 currently passing
tests should be unaffected unless new config options are explicitly
turned on.

Work incrementally, explain before implementing anything nontrivial,
flag rather than silently resolve any remaining ambiguity, and do not
commit or push without being asked.
