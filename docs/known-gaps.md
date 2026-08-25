# Known gaps

- **No `examples/` directory** walking through a full run end-to-end
  outside of the paper-reproduction notebook.
- **Paper 2 parameter-sweep/experiment-runner machinery** — the
  `Event`/`PopulationHistory` abstraction (`docs/science/
  development-roadmap.md` item 4) does not exist. Figure 6's
  stochastic-realisations script (`docs/figures.md`) is a standalone
  repeat-seed runner, not this machinery. Only the minimal
  `did_merge`/`merge_time` event bookkeeping Figure 3 needs exists.
- **Brček, Hirai, Mandel & Lower (2026) rejuvenation** — not
  available; MS-mass-gainer rejuvenation uses Tout et al. (1997)
  instead (`docs/physics/mass-transfer.md`).
- **Hovis-Afflerbach et al. (2025) stripped-donor properties** — no
  interface stub exists. Post-SN RLOF and stable mass transfer both
  leave the donor's post-transfer envelope/structure unmodelled beyond
  the mass change itself.
- **CE eqs. 74-77** (partial-envelope-retention Newton-Raphson solve
  for the merged mass) — not implemented; the merger branch assumes
  full envelope loss (`docs/physics/mass-transfer.md`). The full GB
  core-mass-luminosity relation (eq. 44) this would need above
  `CORE_MASS_BGB_MAX_MASS` is also not implemented.
- **No compact-object-type (WD/NS/BH) census** — `RemnantTable`
  returns mass only, no type label.
- **No accretion-UV spectral model** — `L_UV(t)` is MS-only; no HMXB
  contribution.
