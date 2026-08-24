# Proposal: stellar-radius module and RLOF/CE interaction classifier

Status: **Stage 1 (radius/phase module), Stage 2's classifier
(`binaries/interaction.py::classify_rlof`/`find_ms_rlof_onset`), the
stable-mass-transfer consequence model
(`apply_stable_mass_transfer`), and all of it wired into
`BinaryPopulation` (opt-in `config.use_rlof_classifier`, processed in
`evolve()`'s Phase 0) are implemented and tested, for MS donors only --
and reconciled with the `enhanced_interaction`/`enhanced_mergers`
prescriptions (see "Decision 3" below).** The radius-module
prerequisite for CE is now also done: `stellar/giant_branch.py`
(`l_bgb`, `r_gb`, `mass_radius_exponent`) and `stellar/
main_sequence.py`'s HG extension (`hg_luminosity`, `hg_radius`,
`l_ehg`, `r_ehg`, `m_fgb`, updated `phase()`) are implemented and
tested, for k=0,1,2 (MS+HG) with M < M_FGB -- three more real
transcription errors were caught and fixed during this implementation
pass (see "Update, implementation session" below). HG core-mass
tracking (`M_HeF <= M < min(M_FGB, CORE_MASS_BGB_MAX_MASS)`, roughly
2-7.3 Msun at solar Z) is also done -- `stellar/
main_sequence.py::core_mass_ehg`/`core_mass_hg`/`m_hef`,
`stellar/giant_branch.py::core_mass_bgb` -- using the paper's own
stated large-mass asymptotic limit for `M_c,BGB` rather than the full
eq. 44. **Core radius is also done** -- `stellar/remnant.py::
white_dwarf_radius`/`core_radius`, using HTP02's white-dwarf mass-
radius relation (eq. 91), cross-checked against a well-known real
value (0.6 Msun WD, ~0.0125 Rsun) before writing any code. Both core-
mass and core-radius work needed a real correction during
implementation, not just an advance flag -- see "Core mass:
implementation note" below for the mass cap this surfaced.
`alpha_CE=0.9` (Zuo & Li 2014, HMXB-specific constraint) and
`lambda_CE=0.5` (HTP02's fixed value, entanglement caveat documented)
have been decided. `classify_rlof` now also handles HG donors --
`binaries/interaction.py::hg_q_crit` (HTP02's own GB q_crit formula,
reused for HG per Zuo & Li 2014's eq. 1) and a new `COMMON_ENVELOPE`
outcome, correctly distinguished from MS donors' `IMMEDIATE_MERGER`
(HG donors are CE-eligible, MS donors are not). **The root-finder is
now also wired all the way through**: `find_rlof_onset` (renamed from
`find_ms_rlof_onset`) searches both the MS and HG radius tracks, and
`BinaryPopulation.evolve()` reaches `COMMON_ENVELOPE` and HG-donor
`STABLE_MASS_TRANSFER` through a real population run, not just direct
`classify_rlof` calls -- this surfaced and fixed one real bug (see
"Trigger integration" below): `evolve()`'s stable-MT consequence code
was unconditionally using `ms_radius()`, which silently breaks for an
HG donor. **Not yet implemented**: the CE energy-balance solve itself
(survive vs. merge) and the consequence model -- `COMMON_ENVELOPE` is
currently a terminal, explicit no-op in `evolve()` (Section "CE
alpha-lambda: implementation outline" below still describes what's
needed there). This document is the written record the
task's two "stop and flag before implementing" points required, plus
the scope changes forced by transcription-risk findings during
implementation.

## Problem

Paper 1's "enhanced interaction"/"enhanced massive-star mergers" binary
prescriptions (see `docs/science/paper1-binary-interaction-proposal.md`)
were originally implemented as an illustrative, not-physically-derived
parameterization (`interaction_boost`, `p_merge`) layered on the
existing `fsur` gate. This task replaces that placeholder with real
stellar-radius and Roche-lobe-overflow physics: given a binary's
(M1, M2, a, Z) at some time, determine whether the donor overflows its
Roche lobe, and if so, whether the outcome is stable mass transfer, a
common envelope that survives, a common envelope that merges, or an
immediate/contact merger.

## Decision 1: Hurley R(t)/L(t) scope -- internal only

Hurley's SSE radius/luminosity fits are used **internally only**, for
per-star radius bookkeeping and RLOF timing. Realta's reported
`L_bol`/`L_UV` observables continue to come from `MSLuminosityTable`/
`UVLuminosityTable` (FSPS, population-level SSP) exactly as before --
see `docs/provenance.md` Section 4/7. The two model families are
expected to diverge somewhat (different physics, different purposes:
per-star mechanistic bookkeeping vs. population-synthesis reported
luminosity) and that divergence is not something this work tries to
reconcile. A future migration to a per-star-summed reported luminosity
would be a separate, larger architectural decision needing its own
proposal.

## Decision 2: q_crit / CE alpha-lambda values, deterministic classifier

Adopt Hurley, Tout & Pols (2002)'s own fiducial prescriptions directly,
exposed as named, overridable `SimulationConfig` parameters (same
pattern as `interaction_boost` etc.), not buried constants:

- Mass-transfer stability: `zeta_ad` vs `zeta_L` comparison (HTP02
  Section 2.6). For this implementation's donor scope (MS only, see
  below): `q_crit_ms = 0.695` (HTP02 Section 2.6.4, "Dynamical mass
  transfer from low-mass main-sequence stars" -- stated there
  specifically for deeply-convective k1=0 donors; **extended here, as
  a named simplification, to radiative k1=1 donors too**, since HTP02
  does not give a separate explicit dynamical-instability criterion for
  radiative MS donors and a full `zeta_ad(M, Mc)` treatment is out of
  this task's scope).
- Common envelope: `lambda_CE = 0.5` (HTP02 Section 2.7.1, eqs. 69-77 --
  the paper's own stated typical value, explicitly flagged there as
  uncertain: "lambda ... is probably not a constant"). `alpha_CE`
  **superseded below** (see "Literature findings" and "Decided") --
  HTP02's own generic "`alpha_CE ~= 1` is used" was the original
  placeholder here; now replaced with `alpha_CE=0.9` from Zuo & Li
  (2014)'s HMXB-specific constraint, a materially better-targeted
  source than HTP02's own offhand remark.
- Classifier form: **deterministic threshold**, not tabulated/
  probabilistic. Xu et al. (2025)'s SMC grid is used only as a
  population-level qualitative sanity check (compare Realta's aggregate
  post-MT/merger fractions against their fiducial ~8%/~7%), not as a
  lookup table -- Realta does not have that grid's data file, and
  building one would be its own undertaking (same pattern as the
  FSPS-table generation scripts).

## Scope change: ZAMS + MS only (not ZAMS + MS + HG as originally agreed)

The original plan (agreed earlier in this session) was to implement
Hurley (2000) through the Hertzsprung Gap (k=0,1,2), covering Case A
(MS) and Case B (HG) RLOF -- the two dominant close-massive-binary
channels. This had to be cut back to **MS only (k=0,1)** during
implementation:

HG's radius formula (eq. 27) requires the giant-branch radius-
luminosity relation (eq. 44-48) to compute its endpoint boundary
condition `R_EHG`, which needs coefficients b1, b4-b7. On a careful
re-read of that specific table, coefficients from supposedly-
independent least-squares fits shared identical values to 6-7
significant figures across different rows (b4/b5's gamma-eta-mu block,
b5/b7's beta) -- a near-certain sign of a transcription defect in
reading that particular dense table, not a real coincidence. Per this
project's "never fabricate data" discipline, that block was not
hard-coded. **HG (k=2) is deferred until those coefficients are
independently re-verified** (see `src/realta/stellar/main_sequence.py`'s
module docstring).

Practical consequence: Stage 2's classifier can currently only reach a
real outcome for MS-donor RLOF (k1 in {0,1}). Any binary whose donor
has evolved past the MS at the moment of RLOF gets an explicit
`"phase_not_modelled"` outcome (see `main_sequence.phase()`, which
raises rather than guessing) -- not silently wrong physics.

## Transcription-risk finding worth recording

During implementation, `ms_radius()` was found to collapse to
near-planet-size (0.02-0.05 Rsun) for 5 and 20 Msun stars for most of
the main sequence, only recovering the correct terminal-MS radius in
the last ~1 per cent of the lifetime -- physically nonsensical, and
specifically in Realta's HMXB-relevant mass range (mcut=8 Msun
default). Root cause: two genuine coefficient transcription errors in
the `Delta_R` perturbation block (`a40`'s gamma exponent was `-2`
instead of `-1`; `a41`'s alpha exponent was `-1` instead of `0`, i.e.
10x too small). Caught by comparing a fresh image-based read against a
user-supplied copy-pasted excerpt of the same table, and confirmed via
`tests/test_hurley_main_sequence.py`'s monotonicity tests (which fail
loudly if this regresses). This is the concrete justification for the
"never fabricate/trust unverified transcription" discipline applied
throughout this session -- an initial "high confidence" image read was
still wrong.

**Update, 2026-08-24**: the same paste-and-compare process was applied
to the GB radius-luminosity coefficients (b1, b4-b7) and the L_BGB
coefficients (a27-a32) that were originally left unverified when HG
was deferred. Two more real errors were caught: `b4`'s gamma/eta/mu had
been duplicated from `b5`'s row (explaining the original "repeated
digits across rows" red flag exactly), `b5`'s own alpha exponent was
also wrong (`e-2` instead of `e-1`), and `a28`'s eta exponent was
`e-2` instead of `e0` -- 100x off. All values below are now
user-verified against the source PDF text directly, not an image read:

```
b1: (alpha=3.97e-1, beta=2.8826e-1, gamma=5.293e-1, eta=0, mu=0)
b4: (9.960283e-1, 8.164393e-1, 2.38383, 2.223436, 8.638115e-1)
b5: (2.561062e-1, 7.072646e-2, -5.444596e-2, -5.798167e-2, -1.349129e-2)
b6: (1.157338, 1.467883, 4.299661, 3.1305, 6.99208e-1)
b7: (4.022765e-1, 3.05001e-1, 9.962137e-1, 7.914079e-1, 1.728098e-1)

a27: (9.511033e1, 6.819618e1, -1.045625e1, -1.474939e1, 0)
a28: (3.113458e1, 1.012033e1, -4.650511, -2.463185, 0)  # eta corrected
a'29: (1.413057, 4.578814e-1, -6.850581e-2, -5.588658e-2, 0)  # a29 = a'29^a32
a30: (3.910862e1, 5.196646e1, 2.264970e1, 2.873680, 0)
a31: (4.597479, -2.855179e-1, 2.709724e-1, 0, 0)  # 3 values only, no eta
a32: (6.682518, 2.827718e-1, -7.294429e-2, 0, 0)  # 3 values only, no eta
```

Each row is `(alpha, beta, gamma, eta, mu)` for
`a_n(Z) = alpha + beta*zeta + gamma*zeta^2 + eta*zeta^3 + mu*zeta^4`,
`zeta = log10(Z/0.02)`, matching `main_sequence.py`'s existing `_A`
table convention. Modifier formulas for b1-b7 -- **update, implementation
session**: re-reading this block directly (not from memory) caught a
third real error here, beyond the table-row ones above: an entire
`b2` clamping step had been dropped, and the constant in it was
mis-read (`-0.14167` instead of `-0.04167`). Confirmed fixed by
cross-checking `r_gb()`'s output against Hurley et al.'s own
illustrative Z=0.02 example formula (`R_GB ~= 1.1*M^-0.3*(L^0.4 +
0.383*L^0.76)`, stated in their Section 5.2 text) -- computing the
corrected `b2` at Z=0.02 gives ~0.383, matching that example's
coefficient exactly; `r_gb()`'s output vs. the illustrative formula
now agrees to ~10-20% (consistent with "simplified illustrative
approximation vs. full fit", not a remaining bug), where it disagreed
by up to 14x before the fix.

```
b1 = min(0.54, b1)
b2 = 10^(-4.6739 - 0.9394*sigma)          # sigma = log10(Z), NOT zeta
b2 = min[max(b2, -0.04167 + 55.67*Z), 0.4771 - 9329.21*Z^2.94]
b'3 = max(-0.1451, -2.2794 - 1.5175*sigma - 0.254*sigma^2)
b3 = 10^b'3
b3 = max(b3, 0.7307 + 14265.1*Z^3.395)    # for Z > 0.004
b4 = b4 + 0.1231572*zeta^5                 # applied on top of the table value
b6 = b6 + 0.01640687*zeta^5
```

`R_GB(M, L) = A*(L^b1 + b2*L^b3)`, `A = min(b4*M^-b5, b6*M^-b7)`
(eq. 46); mass-radius exponent `x = 0.30406 + 0.0805*zeta +
0.0897*zeta^2 + 0.0878*zeta^3 + 0.0222*zeta^4` (eq. 47) -- this `x` is
also the one used directly in Zuo & Li (2014)'s HG/GB `q_crit` formula
above.

## Decision 3: prescription wiring -- additive, not a replacement (yet)

**Status: reconciled, 2026-08-24 -- see docs/provenance.md Section 6's
"Reconciliation" subsection for the implemented old-vs-new mapping.**
`interaction_boost` now applies only to binaries the classifier found
actually underwent stable mass transfer (not unconditionally);
`enhanced_mergers` now drives its merger rate through a lowered
`q_crit_ms=0.4` on the real classifier rather than the independent
`p_merge` random draw, which is no longer auto-enabled for that
prescription (still available as an explicit override). `p_merge`/
`interaction_boost` themselves are not deleted -- they remain usable
knobs, just no longer the primary driver for these three prescriptions.

The task's own framing calls this classifier "the physics gap
`enhanced_interaction`/`enhanced_mergers` depend on," implying it
should eventually *replace* those prescriptions' illustrative
`interaction_boost`/`p_merge` placeholders (see
`paper1-binary-interaction-proposal.md`). But those placeholders are
already tested, pinned, and are the current definition of those
prescriptions -- used in the Paper 1 figures already generated (see
`scripts/run_paper1_experiment.py`). Swapping them wholesale now would
silently change scientific behaviour for existing configs.

Chosen instead: a new, independent, opt-in config field
(`use_rlof_classifier: bool = False`). When enabled, it adds a genuine
new MS-RLOF event (stable mass transfer or immediate merger) to
`BinaryPopulation.evolve()`, running alongside -- not instead of -- the
existing `fsur`/`interaction_boost` gate, and independent of
`binary_prescription`. `enhanced_interaction`/`enhanced_mergers`
continue to mean exactly what they meant before. Reconciling the two
(having those prescriptions actually consume this classifier's output
instead of the placeholder parameters) remains a deliberate, separate,
later decision -- not resolved here.

## Stage 2 design (implemented for MS donors; wiring also done)

Given (M1, M2, a, Z) and the current age of each star:

1. Compute each star's Roche-lobe radius via the Eggleton (1983) fit
   (HTP02 eq. 53): `R_L1/a = 0.49 q1^(2/3) / (0.6 q1^(2/3) +
   ln(1+q1^(1/3)))`, `q1 = M1/M2`.
2. Compute the donor's actual radius via `main_sequence.ms_radius()`
   (donor must be MS-phase; otherwise return `"phase_not_modelled"`).
3. If `R_donor < R_L1`: `"detached"` (no interaction).
4. Else (RLOF): if `q1 > q_crit_ms` (0.695): `"immediate_merger"`
   (HTP02 Section 2.6.4 -- MS donors are not in HTP02's CE-eligible
   donor list, Section 2.7.1, so dynamical MS mass transfer merges
   directly rather than forming a CE).
5. Else: `"stable_mass_transfer"`.
6. (HG donors, q_crit=4 leading to CE, deferred with phase 2 above.)

Consequence model (mass/orbit update) -- **implemented, instantaneous**
(see "Decision: instantaneous vs. rate-integrated" below for why):

- Stable MT: **decided, and corrected once before implementation.**
  The first version of this decision ("transfer until masses
  equalize", `delta_m = (M_donor - M_companion)/2`) was accepted, but
  turned out to be physically backwards: `classify_rlof()` only labels
  `STABLE_MASS_TRANSFER` when `q1 <= q_crit_ms < 1`, i.e. the donor is
  *already* the lighter star -- mass in RLOF always flows donor ->
  companion, so continued transfer moves the mass ratio further from
  equality, never toward it. Caught and corrected before writing any
  code. The implemented rule instead root-solves the mass Δm
  transferred (conservative: total mass and orbital angular momentum
  both conserved, standard `a_f = a_i*(M1i*M2i/(M1f*M2f))^2`, not a
  citation-requiring result) such that the widened orbit's Roche-lobe
  radius exactly equals the donor's *current* radius -- the physical
  point of new detachment, treating the donor's radius as fixed during
  the instantaneous jump (no post-mass-loss stellar-structure response
  is modelled -- that is exactly the Brček, Hirai, Mandel & Lower
  (2026) concern named in the task brief, not available this session).
  See `binaries/interaction.py::apply_stable_mass_transfer`.
  Both stars' lifetime clocks are reset from the transfer time at
  their new masses (the same full-reset simplification used for the
  merger channel) rather than partial (Tout et al. 1997/Brček et al.
  2026) rejuvenation -- an explicit extension point, not implemented.
- Immediate merger: conservative mass combination
  (`merge_stellar_masses`), restricted to the k1,k2 in {0,1} (MS+MS)
  subset relevant to this scope -- not HTP02's full Table 2 stellar-
  type combination logic (that table spans k=0-14, most of which is
  out of this module's scope).
- Stripped-donor properties (Hovis-Afflerbach et al. 2025): explicit
  extension point/interface stub only, per the task's own scope note --
  not implemented.

## Decision: instantaneous vs. rate-integrated mass transfer

HTP02 rate-integrates mass transfer via Kelvin-Helmholtz/nuclear
time-scales (eqs. 58-61) -- a genuinely different kind of physics from
anything else in `evolve()`, which is entirely event-based (SN and
merger both happen instantaneously at a precomputed time). Chosen:
**instantaneous**, for architectural consistency and because rate-
integration would be a materially larger addition than this milestone
calls for -- reviewed and accepted in chat before implementation.

## Validation plan

- Component tests for the classifier (fixed M1/M2/a/Z scenarios spanning
  detached / stable-MT / immediate-merger), mirroring
  `tests/test_evolve.py`'s phase-isolation style.
- Population-level check: run a Kroupa/SMC-like config through the
  classifier at every RLOF-eligible timestep and compare the aggregate
  post-MT and merger fractions against Xu et al. (2025)'s fiducial SMC
  numbers (~8 per cent post-MT OB stars, ~7 per cent merger products)
  as a qualitative sanity check, not an exact-match requirement.
- `docs/provenance.md` gets new rows citing HTP02 section/equation
  numbers directly, following the existing discipline.

## CE alpha-lambda: implementation outline (not yet implemented)

Written 2026-08-24, to scope the next piece of work and to identify
exactly what a literature search needs to cover before implementing.
CE requires more prerequisite work than the classifier itself did, since
HTP02's CE-eligible donors are giant/HG-only (Sec. 2.7.1) and Stage 1's
radius module currently stops at the MS (k=0,1) -- see "Scope change"
above.

### 1. Blocking prerequisite: post-MS radius + core-mass module

**HG radius: done.** `stellar/giant_branch.py`/`stellar/
main_sequence.py`'s HG extension are implemented and tested (M <
M_FGB) -- see the "Update, implementation session" note above for the
coefficient-verification history (a27-a32, b1-b7, and a dropped `b2`
clamp, all caught and fixed).

**HG core mass: done**, for `M_HeF <= M < min(M_FGB,
CORE_MASS_BGB_MAX_MASS)` -- roughly 2-7.3 Msun at solar Z, narrower
than the `M_HeF` to `M_FGB` range originally assumed (see "Core mass:
implementation note" below for why) -- `core_mass_ehg`/`core_mass_hg`/
`m_hef`/`core_mass_bgb`.

**HG core radius: also done** -- `remnant.py::white_dwarf_radius`/
`core_radius`, HTP02's white-dwarf mass-radius relation (eq. 91),
cross-checked against a well-known real value before writing code
(0.6 Msun WD, ~0.0125 Rsun) rather than put through a paste-
verification round -- see docs/provenance.md Section 9's core-radius
row for why that was judged sufficient here specifically.

CE can only trigger for a donor with a "dense core" (k1 in
{2,3,4,5,6,8,9} -- HG through AGB and naked-helium phases). Both
prerequisite quantities the energy-balance formula (eq. 69) needs,
`M_c1` and `R_c1`, are now available within the supported mass range.

Practically: the scope call remains HG-only donors,
`M_HeF <= M < CORE_MASS_BGB_MAX_MASS`, for a first CE pass -- the true
GB (eqs. 31-45, full core-mass-luminosity evolution, needed both for
`M >= CORE_MASS_BGB_MAX_MASS` and for donors that reach the actual
giant branch) is a separate, larger addition, not something to fold in
here without its own scope decision, the same way Stage 1 made this
call for the MS/HG boundary.

### Core mass: implementation note (2026-08-24)

`M_c,EHG` (eq. 28, the anchor `core_mass_hg` grows toward across HG)
branches by mass. For Realta's HG-relevant range (`M_HeF~2` to
`M_FGB~13` Msun at solar Z, covering essentially all HG stars above
the default `mcut=8`), the applicable branch is `M_c,BGB` -- core mass
at the base of the giant branch for intermediate-mass stars.

`M_c,BGB`'s own full formula (eq. 44) is not self-contained: it needs
`M_c,BAGB` (core mass at the base of the AGB, eq. 66), which needs new,
unverified coefficients (`b36-b38`), plus two further constants that
confusingly reuse eq. 10's `c1`/`c2` names for unrelated values.
Implementing eq. 44 faithfully would mean another coefficient-
verification round before any code.

**Decided (chat, 2026-08-24)**: use the paper's own stated large-mass
asymptotic limit instead -- "for large enough M we have `M_c,BGB ~=
0.098*M^1.35`, independent of Z" (Hurley et al. 2000, Sec. 5.2 text).
Z-independent, needs no new coefficients. Flagged explicitly in
`giant_branch.py::core_mass_bgb`'s docstring: this is NOT the full eq.
44 formula, and how good the approximation is at the *lower* end of
Realta's HG-relevant range (`M_HeF~2` Msun) is not independently
quantified -- the paper does not state a precise "large enough"
threshold.

**Update, implementing `core_radius` (2026-08-24)**: re-reading eq. 44
directly (needed to check what `core_radius` would require) showed
this framing had understated the gap -- `C` in eq. 44 is not "two more
constants," it's the inverse of the *entire* GB core-mass-luminosity
relation (eqs. 31-43: piecewise `D`/`p`, a mass-dependent hydrogen-
burning rate, time-integration) -- a materially larger addition than a
coefficient block, comparable in size to the MS or HG modules
themselves. That re-read also surfaced a real, concrete problem with
the accepted approximation, not just the abstract "unquantified at low
M" caveat above: `core_mass_bgb(10.0)` returns 2.17 Msun, exceeding
the Chandrasekhar mass, at a phase (end of HG) where real stars in
this range have core masses of a few tenths of a solar mass. This
broke `white_dwarf_radius()` outright (`NaN` under the square root)
rather than merely being imprecise.

Root cause pinned down exactly: `0.098*M^1.35` is `c1^0.25*M^(c2/4)`
using eq. 44's own `c1=9.20925e-5`, `c2=5.402216` (self-consistency
confirmed: `c1**0.25~=0.098`, `c2/4=1.350554~=1.35`), and this
exceeds `M_ch=1.44` above `M=(M_ch^4/c1)^(1/c2)~=7.317` Msun --
derived directly from those constants, not a separately chosen
cutoff. **Decided (chat, 2026-08-24)**: cap `core_mass_bgb` at that
mass (`giant_branch.CORE_MASS_BGB_MAX_MASS`), raising above it rather
than returning a wrong number -- the same idiom already used for
`M_FGB`/`M_HeF`/`t_BGB` elsewhere in this codebase, chosen over either
building the full eq. 44 machinery now or silently clamping to a
plausible-looking but still-wrong value. Realta's actually-supported
HG-donor core-mass/core-radius range is therefore
`M_HeF (~2) <= M < CORE_MASS_BGB_MAX_MASS (~7.3)` Msun at solar Z, not
the full `M_HeF` to `M_FGB (~13)` range originally assumed.

### 2. The energy balance itself (HTP02 Sec. 2.7.1, eqs. 69-77)

Given a donor with core mass `M_c1`, core radius `R_c1`, envelope mass
`M_env1 = M1 - M_c1`, at CE onset:

- Envelope binding energy: `E_bind,i = -G(M1*M_env1/(lambda*R1) +
  M2*M'_env2/R2)`, HTP02's fiducial `lambda=0.5`.
- Initial orbital energy: `E_orb,i = -G*M_c1*M'_c2/(2*a_i)`.
- Efficiency: `E_bind,f - E_bind,i = alpha_CE*(E_orb,f - E_orb,i)`,
  HTP02's fiducial `alpha_CE~=1`.
- Solved via Newton-Raphson for the final mass `M_f`, using the
  giant's mass-radius power-law index `x` (from `R ~ M^-x`, part of the
  GB M_c-L relation) to get `E_bind,f`.
- Survive vs. merge: compare the resulting `a_f` against the sum of
  the two (post-CE) stellar/core radii -- overlap means the cores
  merge instead, routed through HTP02's Table 2 stellar-type
  combination matrix (already transcribed during this session's
  reading of the 2002 paper, see the earlier "Common envelopes,
  coalescence and collisions" notes).

### 3. What to search the literature for

Both `alpha_CE` and (especially) `lambda_CE` are explicitly flagged as
uncertain in HTP02 itself ("it is probably not a constant... generally
alpha_CE ~= 1 is used"; lambda's fiducial 0.5 is a single fixed number
standing in for something the literature treats as
structure/mass/evolutionary-state-dependent). Concretely worth
checking:

- Whether to keep `lambda_CE` fixed at 0.5 (HTP02's own simplification)
  or adopt a mass/evolutionary-phase-dependent `lambda(M, R, k)`
  prescription from later work (this is exactly the kind of thing
  Jarrod Hurley's own more recent work, or the broader post-2002 CE
  literature, would address -- the task's own brief anticipated this
  coming up).
- Whether `alpha_CE~1` is still the standard default choice, or
  whether more recent population-synthesis codes (e.g. COMPAS, the
  same code referenced for the Brček et al. 2026 MS-response
  treatment) use a different fiducial value or a distribution over it.
- Any updated giant-branch core-mass/core-radius fits that would
  affect the b1/b4-b7 coefficient re-verification in item 1 above --
  worth checking whether a cleaner source than the original 2000 paper's
  dense appendix table exists (a review, an errata, or a maintained
  reference implementation).

Whatever is found should be adopted the same way `Q_CRIT_MS` and the
now-decided `alpha_CE=0.9`/`lambda_CE=0.5` values (below) already are:
named, overridable `SimulationConfig` fields, not buried constants,
with the specific source cited directly in `binaries/interaction.py` and
`docs/provenance.md`.

### 3a. Literature findings (2026-08-24): Zuo & Li (2014, MNRAS 442, 1980)

"On the common envelope efficiency" -- directly on-target: constrains
`alpha_CE` specifically from HMXB populations (matching simulated
`L_X` vs. displacement-from-star-cluster distributions against real
data from three starburst galaxies -- M82, NGC 1569, NGC 5253; source
data from Kaaret et al. 2004), the same kind of coeval-starburst HMXB
population Realta models. Findings relevant to this module:

- **`alpha_CE ~ 0.8-1.0` preferred; `alpha_CE < ~0.4` excluded**
  (their Table 3/4, Fig. 1-2 -- models with `alpha_CE>=0.8` reproduce
  the observed displacement distribution, `alpha_CE<0.4` clearly
  fails). A strong, HMXB-specific candidate default -- better targeted
  than a generic literature value, since it comes from fitting the
  same kind of population Realta is built to model.
- **Explicit caveat from the authors**: this range is entangled with
  their specific `lambda` treatment -- "due to an ambiguous definition
  for the core boundary in the literature, the used lambda here still
  carries almost two orders of magnitude uncertainty, which may
  translate directly to the expected value of alpha_CE." Not an
  independent alpha_CE measurement; a paired (alpha_CE, lambda)
  choice, conditional on the lambda prescription below.
- **A low-risk HG-donor `q_crit` refinement** (their eq. 1, citing
  Shao & Li, in prep.): `q_crit = [1.67 - x + 2(M_c1/M1)^5] / 2.13`
  for HG (and GB/AGB) donors -- this is *the same formula* HTP02
  already gives for GB donors (its own eq. 56-57, which this document
  already has verified coefficients for), simply extended to HG donors
  instead of HTP02's own admittedly crude fixed `q_crit=4` for that
  phase (HTP02 itself calls this "rather approximate," inviting future
  calibration). Since the formula and its coefficients are already
  trusted from HTP02, adopting this for HG donors needs no new
  coefficient-transcription risk -- just applying an existing,
  verified formula to a phase HTP02 didn't originally apply it to.
- **Where variable lambda actually comes from**: Loveridge, van der
  Sluys & Kalogera (2011, ApJ 743, 49) -- cited as the source of
  envelope-binding-energy fitting formulae with implicitly
  mass/phase-dependent lambda, built on Dewi & Tauris (2000)'s
  core-boundary definition (core = mass below the point where
  hydrogen fraction X=10 per cent). This is the concrete next paper to
  read if a variable lambda(M, phase) is wanted instead of HTP02's
  fixed 0.5 -- not yet read/verified in this session.
- Zuo & Li also use a more detailed envelope-binding-energy
  prescription than HTP02's eq. 69 -- their eq. 4 adds an internal-
  energy term (`alpha_in` times the integrated internal energy, not
  just gravitational potential energy), following van der Sluys,
  Verbunt & Pols (2006). This is a materially bigger addition (needs
  detailed stellar-structure internal-energy profiles, not just
  fitting-formula radius/mass) -- noted here as a further-out option,
  not part of the near-term plan.

**Decided (chat, 2026-08-24)**: `lambda_CE=0.5` fixed (HTP02's own
value), `alpha_CE=0.9` as the single-point default -- Zuo & Li's own
basic-model choice (their model A09/M1, explicitly identified in their
paper as "our basic model according to our calculations"), sitting in
the middle of their empirically-preferred `0.8-1.0` range. `alpha_CE`
should be exposed as a named, overridable `SimulationConfig` field
(same pattern as `q_crit_ms`) so the `0.8-1.0` range can be explored as
a sensitivity range, not just the single point default. The
lambda-alpha_CE entanglement caveat above still applies and should be
carried into `docs/provenance.md`'s citation for this choice rather
than presented as an independent, unconditional measurement. Loveridge
et al. (2011) is deliberately NOT being read for this pass -- fixed
lambda accepted, per this decision, not because the variable-lambda
question is resolved.

### 4. Consequence model

New orbit (`a_f` from the energy-balance solve, a genuinely different
mechanism from `apply_stable_mass_transfer`'s conservative-widening
formula, not a generalization of it), new stellar type per Table 2,
envelope stripped from both stars. Needs its own instantaneous-vs-
rate-integrated consideration too, though CE itself is already treated
as effectively instantaneous (dynamical-timescale) in HTP02, so that
question is likely moot here in a way it wasn't for stable MT.

### 5. Trigger integration -- done, for HG donors (2026-08-24)

`classify_rlof()` now branches on donor phase (`binaries/
interaction.py`): MS (k=0,1) uses the flat `q_crit_ms=0.695` and
merges directly when unstable; HG (k=2) uses `hg_q_crit()` (HTP02's
own GB q_crit formula, eqs. 56-57, reused per Zuo & Li 2014's eq. 1 --
not HTP02 Sec. 2.6.1's `zeta_ad`/`zeta_eq` comparison as originally
guessed here, since the reused GB formula already gives an equivalent,
better-targeted result without needing that separate machinery) and
returns `COMMON_ENVELOPE` when unstable, since HG donors are CE-
eligible. GB (k=3) and later donor phases remain unreachable
(`main_sequence.phase()` raises past `t_BGB`).

**Also done (2026-08-24)**: the root-finder itself, renamed
`find_rlof_onset` (from `find_ms_rlof_onset`, now inaccurate), extended
to also search the HG radius track (`hg_radius`, monotonic across HG
the same way `ms_radius` is monotonic across the MS -- confirmed in
`tests/test_hertzsprung_gap.py`, with continuity at the MS/HG boundary
by construction of eq. 27) so a binary whose donor crosses its Roche
lobe after leaving the MS but before reaching the GB is now found.
Outcome determination after finding a crossing now checks
`main_sequence.phase()` at the crossing time to pick `q_crit_ms` vs.
`hg_q_crit(...)`, correctly returning `COMMON_ENVELOPE`/
`STABLE_MASS_TRANSFER` for HG donors and `PHASE_NOT_MODELLED` when the
crossing is computable (M < M_FGB) but the outcome isn't (M >=
CORE_MASS_BGB_MAX_MASS, where `hg_q_crit` needs core mass).

Wiring this all the way through surfaced one real bug, not just a
missing feature: `BinaryPopulation.evolve()`'s Phase 0
`STABLE_MASS_TRANSFER` consequence code unconditionally called
`ms_radius()` for the donor's radius. That's correct as long as
stable-MT crossings could only happen during the MS -- once
`find_rlof_onset` could also find one during HG, it silently broke
(`apply_stable_mass_transfer`'s `brentq` bracket-sign error, since
`ms_radius()` gives the wrong value entirely for an HG-phase star).
Fixed by checking `phase()` at `rlof_time` and dispatching to
`ms_radius()`/`hg_radius()` accordingly -- see
`tests/test_rlof_wiring.py::test_evolve_applies_stable_mass_transfer_for_hg_donor`,
sensitivity-verified.

### 2/4, done (2026-08-24): the energy-balance solve and consequence model

`apply_common_envelope` (`binaries/interaction.py`) implements HTP02
eqs. 69-73 for an HG donor with an MS companion. Given the donor's
core mass/radius (`core_mass_hg`/`remnant.core_radius`) and envelope
mass (`M1 - M_c1`):

- `E_bind,i/G = -(1/lambda_CE)*(M1*M_env1/R1)` (eq. 69; the companion's
  own envelope term vanishes for an MS "effective core," `M'_env2=0`).
- `E_orb,i/G = -(1/2)*M_c1*M'_c2/a_i`, `E_orb,f/G = E_bind,i/(G*alpha_CE)
  + E_orb,i/G` (eqs. 70-71), solved directly for `a_f` (eq. 72) -- `G`
  cancels exactly, since every term above is linear in it, so no
  unit-consistency risk from mixing `G`'s value in.
- Coalescence check: as the CE inspiral proceeds, the orbit shrinks
  monotonically from `a_i`; either it reaches `a_f` (full envelope
  ejection) or the bare core/companion first fill their own Roche
  lobes (via the Eggleton fit on their *actual* radii) at some larger
  separation `a_L` -- whichever happens at the larger separation
  happens first. `a_f > a_L` -> survives (new orbit `a_f`, donor
  stripped to `M_c1`); `a_f <= a_L` -> merge.

**Scope gap, not an oversight**: eqs. 74-77 (the partial-envelope-
retention Newton-Raphson solve for the merged star's final mass) are
NOT implemented -- they need `R_i`, "the radius the system would have
if it were to coalesce immediately," which HTP02 does not define
operationally in a way this module could implement without further,
separate study. The merge branch instead assumes full envelope loss at
coalescence (a bare-core + companion conservative merge via the
existing `merge_stellar_masses`), the same "no partial retention"
simplification already used for the MS-MS `IMMEDIATE_MERGER` case.

**Finding from development, verified not to be a bug**: every tested
realistic mid-HG donor/companion combination merges (`survives=False`)
-- e.g. `donor=5, companion=3, a=8 Msun/Rsun` at mid-HG age. Traced
this down before treating it as done: a mid-HG donor still carries most
of its mass in an extended, diffuse envelope (e.g. `M_env1=4.1` of `M1
=5` Msun), so `E_bind,i` is large relative to the available orbital
energy, forcing `a_f` to an unrealistically tight separation --
smaller than the actual stellar radii allow before the companion (or
core) already fills its own Roche lobe. Confirmed this is genuine
physics, not an arithmetic/sign error, by hand-picking a small-
envelope-fraction case (`M_env1=0.1` of `M1=5` Msun) and showing the
same code path DOES return `survives=True` there -- see
`tests/test_ce_energy_balance.py::test_apply_common_envelope_survives_when_envelope_is_small`.
This matches the literature's general expectation that HG-donor CE is
merger-prone (HTP02 Sec. 2.7.1 itself flags HG donors as the case where
"a common envelope and spiral-in are the only outcome" is likely;
population-synthesis codes such as StarTrack/COMPAS often treat
HG-donor CE as a forced merger by convention) -- not a defect in this
implementation.

`alpha_ce`/`lambda_ce` are now exposed as overridable
`SimulationConfig` fields (`config.py`), matching the `q_crit_ms`
pattern, defaulting to the module constants (`ALPHA_CE=0.9`,
`LAMBDA_CE=0.5`) when unset; no `binary_prescription` currently varies
them (unlike `q_crit_ms`). Wired into `BinaryPopulation.evolve()`'s
Phase 0: `COMMON_ENVELOPE` now calls `apply_common_envelope` and
applies its consequence -- survive updates `m1`/`m2`/`a` and resets
both lifetime clocks (the same full-reset simplification used
elsewhere in this module); merge routes through the same
`did_merge`/`merge_time`/lifetime-reset pathway as `IMMEDIATE_MERGER`,
merging the donor's core mass (not its pre-CE full mass) with the
companion. See `tests/test_ce_energy_balance.py` and
`tests/test_rlof_wiring.py`'s new CE-wiring tests, both sensitivity-
checked (deliberate `alpha_ce`/`lambda_ce` perturbation confirmed to
change the result in the expected direction).

### Units bug, found and fixed (2026-08-24): self.a is AU, this module assumed Rsun

Found by running `scripts/run_paper1_experiment.py` end-to-end for the
first time -- not by any unit test, since every test for this module
was written with hand-picked, internally-consistent Rsun-scale
separations. `BinaryPopulation.self.a` (`AFAC`/`PFAC`) is in **AU**:
confirmed via an Earth-Sun sanity check (`AFAC`'s formula gives
`a~=0.99` for `M=1 Msun, P=365.25 days`, only sensible as ~1 AU, and
`PFAC=365.229126` is essentially a sidereal year in days -- these
constants encode Kepler's third law in the standard AU/Msun/year
convention, the same one Power et al.'s original Fortran uses). But
every stellar-radius function this session added returns **Rsun**, and
`roche_lobe_radius`/`classify_rlof`/`find_rlof_onset`/
`apply_stable_mass_transfer`/`apply_common_envelope` all compare
`separation` directly against those Rsun-scale radii, with no
conversion. Effect: every donor looked ~215x closer to its Roche lobe
than it really is, so the actual Paper 1 config produced
`IMMEDIATE_MERGER` for ~100% of massive binaries in every RLOF-
classifier-enabled prescription, making Figure 2 completely
degenerate.

Fixed in `binaries/population.py`: a new `RSUN_PER_AU = 215.032`
constant (1 au = 1.495978707e11 m, IAU 2012; R_sun = 6.957e8 m, IAU
2015 nominal), applied at every call from `population.py` into
`interaction.py` (`self.a[i] * RSUN_PER_AU` going in, `new_a_rsun /
RSUN_PER_AU` coming back out). `self.a` itself stays in AU everywhere
else -- the pre-existing SN1 mass-loss orbit-widening code and every
already-pinned regression value are untouched. See
`docs/provenance.md`'s Section 10 for the full writeup, including the
three `tests/test_rlof_wiring.py` tests that needed their hand-
constructed `pop.a` values converted to AU to match.

A second, independent finding surfaced by the same end-to-end run:
even with the units fixed, `configs/paper1_basic_experiment.yml`'s
inherited `pmin=0.1` days puts most massive binaries in contact at
birth (Kepler's third law again) -- raised to `pmin=1.0` day for that
config specifically (not the global default), per the "Paper 1 config"
note in `docs/provenance.md` Section 6.
