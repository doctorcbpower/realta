# Case-A mass transfer for q1 > 1 (donor heavier than companion) — proposal

**Status: not implemented.** This is a fix proposal and literature audit,
not a design in progress or a description of current behaviour. Read-only
audit performed 2026-09-27/28; nothing below has been coded. Tracked as a
known gap (`docs/known-gaps.md`); deferred to future work — does not block
the current paper submission (which drops `enhanced_interaction`,
`enhanced_mergers` and `wind_capture` from its main text instead).

## The gap

[`interaction-prescriptions.md`](interaction-prescriptions.md)'s
"`STABLE_MASS_TRANSFER` donor-mass structural limit" and
[`rlof-classifier.md`](rlof-classifier.md)'s "Donor-selection property"
already document the symptom: `find_rlof_onset` selects whichever star
overflows first, which for realistic massive-star binaries is almost
always the heavier one, so `q1 = M_donor/M_companion >= 1` for the
selected donor in practice (confirmed directly:
`tests/test_rlof_classifier.py::test_find_rlof_onset_favours_the_more_massive_star_as_donor`
found no counterexample in a parameter sweep; a population-level check
against the paper1 base population, seed 42, found 0/1114 systems with a
star-2 MS donor selected earlier than star 1).

This audit traces the exact mechanism, not just the symptom.
`classify_rlof`'s MS branch (`interaction.py:156-162`) is not itself
asymmetric — it returns `STABLE_MASS_TRANSFER` for a low-q1 donor
regardless of which star that is. The asymmetry is enforced one level
down: `apply_stable_mass_transfer` (`interaction.py:349-354`) raises
`ValueError` for `donor_mass >= companion_mass`, and every config in the
codebase keeps `q_crit_ms <= 1` (`config.py` validation permits this but
nothing sets it otherwise). So a primary-donor RLOF event — donor mass
>= companion mass by construction, since `generate_population` clamps
`m2 <= m1` — is forced to `IMMEDIATE_MERGER` regardless of the numeric
`q_crit_ms` value. Raising `q_crit_ms` cannot fix this; the solver
downstream cannot process a stable outcome for q1 >= 1 at all.

## Literature basis

**Ge et al. series — does not cover realta's OB-star range.**
- Ge, Webbink, Han & Chen (2010, ApJ 717, 724 — Paper I): methods only,
  explicitly defers results to future papers; no q_crit(M,R) fit.
- Ge, Webbink, Chen & Han (2015, ApJ 812, 40 — Paper II): full
  0.1-100 Msun range at Z=0.02 (realta's own mass range) — but gives
  tables/figures only (their Table 3, Figs. 8-9), no closed-form fit.
  Modern implementations (COMPAS's `GE`/`GE_IC` prescriptions)
  interpolate this grid directly rather than fitting it.
- Ge, Tout, Chen, Sarkar, Walton & Han (2023, ApJ, "Criteria for
  Dynamical Timescale Mass Transfer of Metal-poor Intermediate-mass
  Stars", arXiv:2302.00183): does give a closed-form fit — q~_ad = a +
  b*log10(M/Msun) + c*log10(R/Rsun), coefficients (Z=0.001:
  a=2.88940,b=-2.46266,c=2.80378; Z=0.02: a=3.09721,b=-3.05344,
  c=3.24722), verified directly from the paper text — but the paper's
  own stated scope is **1.6 <= M/Msun <= 10**, their own IMXB/HMXB
  donor-mass boundary (following Tauris & van den Heuvel 2006) drawn at
  10 Msun. realta's donors run up to ~99 Msun with `mcut=8`. **Do not
  use this fit** — it would mean extrapolating roughly an order of
  magnitude past its calibrated range and outside its own stated
  regime (IMXBs, not HMXBs), despite matching the functional form
  someone might reasonably guess at from the literature.

**Wellstein, Langer & Braun (2001, A&A 369, 939) + Pols (1994) — the
applicable source, verified from primary text.**
- 74 detailed Case A/B evolutionary tracks, primaries **12-25 Msun**,
  secondaries 6-24 Msun — inside realta's OB range.
- Direct quote: "Case A systems are prone to contact due to reverse
  mass transfer during or after the primary's main sequence phase, all
  systems obtain contact for **initial mass ratios below ~0.65**, with
  a merger as the likely outcome" (q = M2/M1 at formation).
- Cross-checked within the same paper: Pols (1994) found Case A
  primaries **8-16 Msun** (starting exactly at realta's own `mcut`)
  avoid contact for **q >~ 0.7**, P0 > 1.6 d — a second, independent
  study at essentially the same number.
- Assumptions match realta's own design constraints: conservative
  transfer (beta=1) assumed for contact-free systems; accretor is a
  normal (not yet compact) star throughout Case A — the MS-MS regime
  needed here.
- Independent corroboration of the event-ordering problem below: the
  paper states they "found it relevant for contact and supernova order
  in Case A systems, particularly for the highest considered primary
  masses" — the literature already recognises SN order can reverse in
  Case A, not a hypothetical edge case invented for this audit.
- **Caveat.** This is an *initial* (formation-time) mass-ratio
  threshold from full evolutionary tracks — "does this system, over
  its whole pre-contact evolution, avoid contact" — not an instantaneous
  q1(t_rlof) criterion the way `classify_rlof` currently evaluates
  things. Rough cross-check: q(initial)~0.65 corresponds to
  q1(donor/accretor)~1/0.65~1.5, broadly consistent with the ~3-4
  literature consensus for radiative-donor q_crit that the Ge et al.
  (2023) paper itself cites (Hjellming 1989; Kalogera & Webbink 1996) —
  two independent methodologies landing in a compatible range.

## The evolutionary-event-ordering problem

`self.turnoff_time` (star-1's clock) and `self.t2_lifetime` (star-2's
clock) are each correctly recomputed per-star after a mass-transfer
event (`population.py:427-465`) — this respects fixed identity/updated
state. But the code that *fires* a supernova never compares the two
clocks: `sn1_mask` (`population.py:552-553`) is unconditionally
`(nturn==0) & (tnow>=self.turnoff_time)`, always processing `self.m1[i]`
as the exploding star; `sn2_mask` (`population.py:778`) is gated on
`nturn==1`, i.e. structurally cannot fire until "SN1" (always m1) has
already happened. If a stable q1>1 transfer reverses the mass ratio and
the originally-lighter star ends up with the smaller absolute explosion
time, the current code would still wait for the originally-heavier
star's (now later) explosion — silently getting the supernova order
backwards. Not reachable today (q1>1 stable transfer cannot occur), but
must be fixed *before or alongside* any q1>1 extension, not left as a
stale downstream assumption.

## Proposed fix

1. **New criterion, additive, not a redefinition of `q_crit_ms`.**
   Evaluate at population-generation time (before RLOF timing is
   computed at all): if q(initial) = m2/m1 (secondary/primary, at
   formation) >= a threshold, permit the RLOF classifier to return
   `STABLE_MASS_TRANSFER` (via new q1>1 mass-ratio-reversal machinery)
   instead of unconditional `IMMEDIATE_MERGER` for whichever star
   overflows first; below threshold, keep existing `IMMEDIATE_MERGER`
   behaviour (matches the literature — those do go into contact).
   Existing `q_crit_ms` (compared against instantaneous q1, q1<1
   branch) is untouched.
2. **Single constant threshold, not a mass-dependent fit** — matches
   the "deliberately simplified, no paper-specific hard-coded mass
   range" design goal, and doesn't manufacture precision the
   literature doesn't support. Suggested default **0.65** (Wellstein,
   Langer & Braun 2001); **0.7** (Pols 1994) is the more conservative
   alternative. Open decision, not made here — see below.
3. **Fix evolutionary event ordering first or alongside** (see above)
   — `sn1_mask`/`sn2_mask` need to dispatch on whichever of
   `turnoff_time`/`t2_lifetime` is actually smaller, not on array
   position.
4. **Extend `apply_stable_mass_transfer`'s root search past the
   mass-ratio-reversal point.** Analytically (from the existing,
   already-generic `_widened_separation` formula): parametrizing by
   mass transferred, the orbit shrinks from RLOF onset to a minimum
   exactly at mass-ratio equality, then widens past it, diverging as
   the original donor's mass -> 0. A root should exist past the
   reversal point but this is derived, not numerically tested — the
   search bracket and root uniqueness need verifying, not assumed from
   the existing q1<1 case.
5. **Generalise the post-SN-RLOF/wind-capture/`fsur`-activation
   blocks** (`population.py:610-617`, `695-774`) to operate on
   "whichever star has already exploded" rather than hardcoded
   `self.m1`/`self.m2`.
6. `interaction_boost` needs no change to its own mechanism
   (`population.py:599-609` already does the right thing) — it is
   inert only because zero pre-SN1 `STABLE_MASS_TRANSFER` events
   currently exist. Once (1)-(3) produce real ones, it activates
   automatically. Whether its existing illustrative values (1.5
   standard, 3.0 enhanced) still give a non-degenerate separation once
   real cases exist is untested, not guaranteed by this fix alone.

## Tests needed before trusting this

Unit level, spanning q1 = 1.05, 1.2, 1.5, 2, 3: stability classification
against the adopted threshold; mass and angular-momentum conservation
(exact); orbital separation strictly decreasing before reversal,
strictly increasing after, minimum exactly at mass-ratio equality;
stellar identity preserved (no m1/m2 swap, only masses change); correct
evolutionary ordering post-reversal (construct a case where reversal
flips which star has the smaller absolute lifetime, confirm the correct
star explodes next — the most important new test, directly exercising
the event-ordering fix); correct rejuvenation of the gainer (existing
`rejuvenate_ms_gainer` should not need to change); correct downstream
HMXB donor/accretor role assignment after a reversed system's first SN.

Population-level diagnostic (not a success criterion): the paper1 base
population's classification should no longer be saturated at 100%
merger for every q1>=1 primary-donor MS event.

## Open decisions

Single constant (0.65 vs 0.7) vs a mass-interpolated threshold between
the two source studies; whether `interaction_boost`'s existing 1.5/3.0
values should be revisited once real stable-MT systems exist.
