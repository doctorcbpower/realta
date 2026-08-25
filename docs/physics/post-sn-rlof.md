# Post-SN secondary Roche-lobe overflow

Not from Power et al. The astrophysically dominant real HMXB-formation
channel: the secondary's *own*, later Roche-lobe overflow onto the
already-compact primary (Case B/C mass transfer onto a neutron star or
black hole) — as opposed to `f_sur`'s phenomenological wind-fed-HMXB
approximation, and structurally distinct from the pre-SN RLOF channel
in [`rlof-classifier.md`](rlof-classifier.md) (see
[`interaction-prescriptions.md`](interaction-prescriptions.md)'s
`STABLE_MASS_TRANSFER` donor-mass structural-limit note for why the
pre-SN channel cannot cover this regime).

Opt-in: `config.use_post_sn_rlof`, default `False`. Gated on
`nturn == 1` (primary already compact, secondary not yet exploded),
checked live every timestep against the secondary's own Roche lobe
(no root-finder needed — `nturn == 1` is already the correct gate,
unlike Phase 0's precompute-once pattern). Uses the same Eggleton
Roche-lobe machinery as Phase 0
(`roche_lobe_radius`, `main_sequence.ms_radius`/`hg_radius`/`phase`).
On trigger, activation is certain (an `L_X` draw via the existing
`xray_calc`), not a stochastic `f_sur`-style draw.

**Deliberately minimal scope**: a single RLOF-only trigger. No
secondary mass-loss/envelope-stripping consequence, no compact-primary
mass growth. Requires `imetal = 2` or `3` (same Z=0 restriction as the
RLOF classifier).

Implementation: `binaries/population.py::evolve`, "Phase 1.5" (between
the SN1 and SN2 phases).
Tests: `tests/test_post_sn_rlof.py`.
