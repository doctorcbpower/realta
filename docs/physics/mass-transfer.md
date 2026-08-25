# Mass-transfer consequence models

Consequence models for the two non-merger RLOF outcomes from
[`rlof-classifier.md`](rlof-classifier.md): stable mass transfer and
common envelope. Applied instantaneously at the precomputed
`rlof_time`, not HTP02's rate-integrated Kelvin-Helmholtz/nuclear
timescale treatment (eqs. 58-61) — Realta's SN and merger events are
already instantaneous state changes, so mass transfer follows the
same architecture.

## Stable mass transfer

Instantaneous conservative mass transfer to the new detachment point:
solves for `Δm` (donor → companion, conservative) such that the
widened orbit's Roche-lobe radius exactly equals the donor's current
radius. Orbital widening is standard two-body angular-momentum
conservation, `a_f = a_i * (M1i*M2i/(M1f*M2f))^2`. Only reachable for
`donor_mass < companion_mass` (the only case `classify_rlof()` labels
`STABLE_MASS_TRANSFER`).

Implementation: `binaries/interaction.py::apply_stable_mass_transfer`,
`_widened_separation`.
Wiring: `binaries/population.py::evolve`, Phase 0.
Tests: `tests/test_rlof_classifier.py` (mass conservation, direction,
orbit widening, Roche-lobe self-consistency at the new separation),
`tests/test_rlof_wiring.py`.

The donor's lifetime clock is a full reset from `tnow` at its new
mass (no verified response prescription exists for a mass-losing
donor). The companion's clock uses rejuvenation (below), not a reset.

## MS-mass-gainer rejuvenation

Tout, Aarseth, Pols & Eggleton (1997, MNRAS 291, 732), Sec. 5.1, eq.
(41):

```
t' = (mu/mu') * (tau'_MS/tau_MS) * t
```

`mu' = mu` (mass ratio collapses to 1, i.e. simple fractional-age
preservation) for `0.3 < M/Msun < 1.3` (radiative core); `mu = M`
(old mass), `mu' = M'` (new mass) otherwise (convective core) — the
surviving `M/M'` factor gives extra rejuvenation beyond fractional-age
preservation for a convective-core gainer. Remaining fraction is
clamped to `[1e-6, 1.0]`.

Applies only when the companion is genuinely MS-phase at `rlof_time`
(`classify_rlof()` places no phase constraint on the companion for
`STABLE_MASS_TRANSFER`); falls back to the full-reset simplification
otherwise.

Implementation: `binaries/interaction.py::rejuvenate_ms_gainer`, wired
into `evolve()`'s `STABLE_MASS_TRANSFER` branch.
Tests: `tests/test_ms_gainer_rejuvenation.py`.

## Common envelope

`COMMON_ENVELOPE` outcomes resolve via the alpha-lambda energy-balance
solve, HTP02 Sec. 2.7.1, eqs. 69-73 (`G` cancels — every term is
linear in it):

- Envelope binding energy and initial/final orbital energy determine
  `a_f`.
- Coalescence check: whichever of `a_f` or the companion/core's own
  Roche-lobe-filling separation `a_L` (Eggleton fit, on the bare
  core's and companion's actual radii) is reached at the *larger*
  separation happens first during inspiral.
- Survival: donor stripped to its core mass, orbit tightened to `a_f`;
  only the donor's own lifetime clock is reset (the companion is mass-
  unaffected, so its clock is left untouched — not reset, not
  rejuvenated).
- Merger: the donor's core mass (not its pre-CE full mass) plus the
  companion merge via the same `merge_stellar_masses` pathway as
  `IMMEDIATE_MERGER`.

`alpha_CE = 0.9` (Zuo & Li 2014, MNRAS 442, 1980, their own basic-model
value, middle of their HMXB-population-calibrated 0.8-1.0 range).
`lambda_CE = 0.5` (HTP02 eq. 69's own fixed value, explicitly flagged
there as not a true constant). Both overridable
(`config.alpha_ce`/`lambda_ce`).

**Not implemented**: HTP02 eqs. 74-77 (partial-envelope-retention
Newton-Raphson solve for the merged star's mass) — needs `R_i`, "the
radius the system would have if it were to coalesce immediately",
which HTP02 does not define operationally. The merger branch assumes
full envelope loss instead.

For MS/HG-range donors, `E_bind,i` is generically large relative to
available orbital energy (the donor still carries most of its mass in
an extended envelope), so realistic mid-HG-donor CE generically
merges — matching the literature's general expectation for HG-donor
CE (HTP02 Sec. 2.7.1; StarTrack/COMPAS often treat it as a forced
merger by convention).

Implementation: `binaries/interaction.py::apply_common_envelope`.
Wiring: `binaries/population.py::evolve`, Phase 0.
Tests: `tests/test_ce_energy_balance.py`, `tests/test_rlof_wiring.py`.
