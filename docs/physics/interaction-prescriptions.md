# Binary-interaction prescriptions

Power et al. (2009/2013) contain no mass-transfer, common-envelope, or
merger physics. `config.binary_prescription` selects one of five named
parameter sets built on top of the RLOF classifier
([`rlof-classifier.md`](rlof-classifier.md)) and mass-transfer
consequence models ([`mass-transfer.md`](mass-transfer.md)). The
default (`"non_interacting"`) is a no-op relative to the ported
baseline.

| Prescription | `use_rlof_classifier` | `interaction_boost` | `q_crit_ms` | Notes |
|---|---|---|---|---|
| `single` | — | — | — | No companion assigned to any `M >= mcut` star; no HMXB channel at all (`binaries/population.py::generate_population`). |
| `non_interacting` | `False` | `1.0` | — | Baseline `f_sur`-only activation. |
| `standard_interaction` | `True` | `1.5` | `0.695` (default) | |
| `enhanced_interaction` | `True` | `3.0` | `0.695` (default) | |
| `enhanced_mergers` | `True` | `1.0` | `0.4` | Lower `q_crit_ms` than HTP02's own fiducial 0.695 (`binaries/interaction.py::Q_CRIT_MS`) — more RLOF-ing systems classified as dynamically-unstable mergers. |

Defaults resolved in `config.py::_PRESCRIPTION_DEFAULTS`.

## `interaction_boost`

Multiplies `fsur` for binaries the RLOF classifier found underwent
stable mass transfer on the MS (`had_stable_mt` gate,
`binaries/population.py::evolve` Phase 1): `fsur_eff = min(1, fsur *
interaction_boost)`. Not applied to binaries that never interacted,
which use plain `fsur`. `interaction_boost = 1.0` for every
prescription except `standard_interaction`/`enhanced_interaction`.

## Pre-SN merger channel

`config.p_merge` / `p_merge_max_period` / `f_merge` — an independent,
formation-time draw: eligible short-period binaries merge at
formation (`m1 -> m1 + f_merge*m2`, `m2 -> 0`), lifetime recomputed
for the merged mass. Not tied to `binary_prescription`; available as
an explicit override on any prescription (`enhanced_mergers` drives
its mergers through the RLOF classifier's `q_crit_ms` instead, see
above).

Implementation: `binaries/population.py::generate_population`.

## Reconciliation with the RLOF classifier

`standard_interaction`/`enhanced_interaction`/`enhanced_mergers` drive
their behaviour through the physics-based RLOF classifier rather than
independent placeholder parameters. `interaction_boost` only applies
to binaries the classifier found underwent stable mass transfer; a
binary that never interacted uses plain `fsur`.

## `STABLE_MASS_TRANSFER` donor-mass structural limit

`generate_population` enforces `m2 <= m1`, and
`classify_rlof`'s `STABLE_MASS_TRANSFER` requires the donor's mass
ratio `q1 = M_donor/M_companion < q_crit < 1` — so the stable-MT donor
is always `m2`, the lighter star. Realta tracks one explosion clock
(`turnoff_time`, from `m1`'s lifetime only); since `m1` is more
massive it explodes first, so `nturn` flips to 1 before the donor's
own predicted `rlof_time`, and Phase 0 (gated on `nturn == 0`) is
permanently closed out. `IMMEDIATE_MERGER`/`COMMON_ENVELOPE` do not
have this problem — their donor is typically `m1` itself, so the RLOF
clock and SN clock track the same star.

This makes `interaction_boost` structurally inert for any realistic
massive-star binary under the current phase structure — see
[`post-sn-rlof.md`](post-sn-rlof.md) for the channel that covers the
astrophysically dominant regime this exposes (post-SN secondary RLOF
onto the compact primary).

Config validation: `config.py::SimulationConfig.__post_init__`
(`binary_prescription` enum, `[0,1]`/non-negative bounds on the four
interaction parameters).

Tests: `tests/test_binary_prescriptions.py`, `tests/test_rlof_wiring.py`.
