# Proposal: binary-interaction and merger parameterization for Paper 1

Status: **reviewed and accepted** (chat, 2026-08-24) -- illustrative
defaults, clearly flagged; GALEX FUV for `L_UV`. This document is the
record required by the Paper 1 implementation prompt's item 2 ("stop
and flag before implementing... produce a short written proposal...
and get it reviewed") before any interaction/merger physics is added.

## Problem

Paper 1's basic experiment (`docs/science/research-programme.md`,
"Basic experiment") compares five binary-prescription variants at fixed
IMF:

1. single-star populations;
2. non-interacting binaries;
3. standard binary interaction;
4. enhanced interaction;
5. enhanced massive-star mergers.

Realta currently has **no** mass-transfer, common-envelope, or merger
physics -- only the Power et al. (2009) `floss`-based sudden-mass-loss
survival criterion and the `fsur` HMXB-activation gate (see
`docs/provenance.md` Sections 1-2). None of variants 3-5 map onto
anything in the existing code, and the paper this is reproducing
predates and doesn't itself specify this physics -- so unlike every
other row in `docs/provenance.md`, these have no citable source. They
are new, Realta-specific parameterizations, and must be flagged as such
everywhere they appear (provenance table, config docstrings, plots).

The explicit constraint from the brief and the implementation prompt:
this is **not** a full binary-evolution code. The parameterization
below is deliberately the smallest addition that gives Paper 1's five
variants distinguishable `L_X(t)`/`L_UV(t)` trajectories, layered on
top of the existing `floss`/`fsur` gate rather than replacing it.

## Proposed mapping

| Variant | Realta mechanism | New parameters |
|---|---|---|
| 1. Single-star | No binary assignment for `M >= mcut`; no HMXB channel at all. `L_X(t) = 0` identically. | none |
| 2. Non-interacting binaries | Exactly today's model, unchanged. This is the calibration anchor -- the existing `floss <= 0.5` survival criterion is a dynamical (sudden-mass-loss) effect, not a mass-transfer interaction, so "non-interacting" is what Realta already does. | none |
| 3. Standard interaction | Multiplicative boost on the activation gate at primary SN: `fsur_eff = min(1, fsur * interaction_boost)`. Represents mass-transfer episodes that raise the fraction of survivors observed as active HMXBs (spin-up, tighter post-interaction orbits, etc.) without modelling the transfer itself. | `interaction_boost` (illustrative default: 1.5) |
| 4. Enhanced interaction | Same mechanism, larger boost. | `interaction_boost` (illustrative default: 3.0) |
| 5. Enhanced mergers | New pre-SN merger channel (see below). Merged systems exit the binary/HMXB channel entirely and fold into the single-star luminosity budget (no HMXB possible from a merged system in this minimal model). | `p_merge`, `p_merge_max_period`, `f_merge` (illustrative defaults below) |

## Merger channel (variant 5)

At binary formation (equivalently, before any SN physics runs -- this
happens once, at `generate_population()` time, since nothing about it
depends on evolved state), binaries with `period < p_merge_max_period`
are each given one Bernoulli draw at probability `p_merge`. A "merge"
means:

- `m1_new = m1 + f_merge * m2` (mass lost during the merger, e.g. to a
  common-envelope ejection, is `(1 - f_merge) * m2`);
- the system's lifetime is recomputed for `m1_new` via the existing
  `LifetimeTable` (this is the standard "rejuvenation" effect of a
  merger -- a merged star behaves roughly like a fresh star of its new,
  higher mass, not like the aged primary it used to be);
- the system is flagged as merged and removed from the binary/HMXB
  channel -- `m2` is zeroed, `nturn` is set so Phase 1/2 of `evolve()`
  never re-triggers SN bookkeeping for the (now nonexistent) secondary;
  the merged star still explodes on its own recomputed lifetime, but as
  a single star, contributing to `L_bol`/`L_UV` only (via
  `MSLuminosityTable`, which is population-total and doesn't
  distinguish single vs. merged-single stars) and never to `L_X`.

Illustrative defaults: `p_merge = 0.2` for short-period systems,
`p_merge_max_period = 10` days, `f_merge = 0.5`. These are **not**
derived from any source -- they exist to give variant 5 a
distinguishable trajectory for a first Paper 1 pass and should be
revisited once real Fig. 1/2 output exists to compare against.

Minimal merger-event bookkeeping (`did_merge: bool[]`, `merge_time:
float[]` on `BinaryPopulation`) is added now because Figure 3 (out of
scope for this milestone) will need it later, and it costs nothing to
record at merger time. No `Event`/`PopulationHistory` abstraction is
built for it yet -- see `docs/science/development-roadmap.md` item 4;
that generalization is deferred until Figure 3 is actually attempted.

## Configuration surface

New `SimulationConfig` field: `binary_prescription: str`, one of
`"single"`, `"non_interacting"` (default -- preserves current
behaviour exactly), `"standard_interaction"`, `"enhanced_interaction"`,
`"enhanced_mergers"`. Each non-default value implies a fixed set of the
parameters above (not independently configurable per-run in this first
pass -- avoids a combinatorial YAML surface before there's a reason for
one). `interaction_boost`/`p_merge`/`p_merge_max_period`/`f_merge`
remain present as individually overridable `SimulationConfig` fields
(with the prescription's implied values as defaults) so they can be
recalibrated later without code changes.

## What this does NOT do

- No mass-transfer/RLOF/common-envelope physics is modelled explicitly
  -- `interaction_boost` is a pure activation-probability multiplier,
  not a rate calculation.
- No orbital-period or mass-ratio dependence in the interaction boost
  itself (only the merger channel is period-dependent, via
  `p_merge_max_period`).
- No compact-object type distinction (BH vs NS) is introduced by any of
  this -- out of scope per the implementation prompt.
- Figure 3 (mergers vs. compact-object formation) is not built here;
  only the event bookkeeping it will eventually need is added.

## UV band decision (implementation prompt item 3)

`L_UV(t)` uses GALEX FUV (~1528 A), via FSPS `get_mags(bands=
["galex_fuv"])`, generated the same way as the existing
`MSLuminosityTable` (`scripts/generate_ms_luminosities.py` ->
`scripts/generate_fuv_luminosities.py`, same FSPS/SPS_HOME requirement,
same per-metallicity fiducial-mass-then-rescale convention). FUV is the
standard massive-star/SFR UV tracer (Kennicutt & Evans 2012) and is
more sensitive to the O/B population than NUV, which would pick up a
longer-lived A-star contribution and dilute the timescale match to
`L_X`/`Q_H` that Paper 1's Q2 depends on. Like `MSLuminosityTable`,
this is flagged in `docs/provenance.md` as this session's own addition,
not part of either paper's original code, with FSPS-vs-Starburst99
systematics carrying the same caveat already recorded for the bolometric
table.
