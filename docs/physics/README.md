# Physics beyond the Power et al. baseline

None of this is in Power et al. (2009/2013) — see `docs/provenance.md`
for the ported baseline itself. Each file below covers one physics
addition: literature source, equations used, implementation location,
and test coverage.

- [`binary-sampling.md`](binary-sampling.md) — binary fraction, mass-ratio
  and period distributions, continuous IMF slope.
- [`interaction-prescriptions.md`](interaction-prescriptions.md) —
  `binary_prescription` phenomenology (interaction boost, mergers) and
  its reconciliation with the RLOF classifier.
- [`stellar-tracks.md`](stellar-tracks.md) — Hurley/Tout main-sequence,
  Hertzsprung-gap, and core-mass/radius tracks (RLOF/CE prerequisite).
- [`rlof-classifier.md`](rlof-classifier.md) — Roche-lobe-overflow
  outcome classification for MS and HG donors.
- [`mass-transfer.md`](mass-transfer.md) — stable mass-transfer
  consequence model, MS-mass-gainer rejuvenation, common-envelope
  energy balance.
- [`post-sn-rlof.md`](post-sn-rlof.md) — secondary Roche-lobe overflow
  onto the post-SN compact primary.
- [`wind-capture.md`](wind-capture.md) — CAK-driven wind mass loss and
  Bondi-Hoyle-Lyttleton capture accretion.
- [`observables.md`](observables.md) — L_bol/L_UV/Q_H population-level
  observables (FSPS-tabulated).

Open items and scope limits are tracked in `docs/known-gaps.md`, not
repeated in each file. Figure-generation scripts are catalogued in
`docs/figures.md`.
