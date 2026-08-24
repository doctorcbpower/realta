# Paper 1 implementation prompt

A self-contained prompt for kicking off the implementation work behind
Paper 1's first end-to-end experiment (see
`docs/science/research-programme.md`). Written so it can be pasted into
a fresh session with no other context. Item 2 below is deliberately left
open -- the quantitative treatment of binary interaction / mergers needs
a written proposal before any of that physics is implemented.

---

You're working on Realta, a modular stellar/binary population-synthesis
Python framework (repo already checked out) built on a ~20-year-old
Fortran Monte Carlo model of HMXBs in globular clusters (Power et al.
2009, MNRAS 395, 1146, arXiv:0902.1897). Before touching any code, read:

- docs/provenance.md -- paper-equation -> implementation -> test
  traceability, and the discipline this project follows for pinning
  numeric regression values.
- docs/science/research-programme.md, the "Paper 1" section -- the
  scientific target this work serves.
- docs/science/development-roadmap.md -- the target architecture; you
  are implementing a deliberately narrow slice of it, not all 28 items.

Governing principles (a "Development and Scientific Software Brief"
applies throughout): preserve the Power et al. (2009) baseline exactly;
never silently change scientific behaviour; flag ambiguity rather than
resolve it yourself; explain any nontrivial physics or design decision
before implementing it; avoid over-engineering; keep the public API
small; treat external tools (MESA, FSPS, etc.) as optional; never
commit or push without being explicitly asked; prefer small,
incremental, reviewable changes.

## Goal

Reproduce, from a single YAML experiment config, Paper 1's two central
figures:

- Figure 1: L_bol(t), L_UV(t), Q_H(t), L_X(t) for several binary
  models, fixed IMF.
- Figure 2 (the central figure): L_X/L_UV vs age, compared across
  binary-interaction prescriptions.

"Reproducible from a YAML config" means: one config file in, the exact
figure out, no manual notebook wrangling in between.

## Scope (narrowed Phase 1 + Phase 2 of development-roadmap.md -- not
the full roadmap)

1. Extend the existing BinaryPopulation / ClusterSimulation /
   SimulationConfig classes into whatever minimal Population /
   PopulationHistory shape is actually load-bearing for this milestone
   (the new observables, plus merger events since Figure 3 needs them).
   Do not build the full Population/Event/environment architecture
   described in the roadmap -- only what this milestone needs. Prefer
   renaming/extending the existing classes over building a parallel
   abstraction.

2. **Stop and flag before implementing.** Paper 1's basic experiment
   compares 5 binary-prescription variants -- single-star, non-
   interacting binaries, standard binary interaction, enhanced
   interaction, enhanced massive-star mergers. Realta currently has
   *no* mass-transfer, common-envelope, or merger physics at all --
   only the Power et al. (2009) fsur-based survival/HMXB-activation
   gate. Before writing any of this, produce a short written proposal
   for what "standard interaction," "enhanced interaction," and
   "enhanced mergers" mean quantitatively (e.g. as parameterized rates
   layered on the existing survival/activation gate, not a full
   binary-evolution code), and get it reviewed. This is exactly the
   kind of ambiguous scientific choice the brief says to surface, not
   resolve unilaterally.

3. Add L_UV(t) as an observable. Check what's already available via
   the FSPS-sourced tables in src/realta/io/tables.py
   (MSLuminosityTable, and the currently-unused IonizingPhotonTable) --
   determine whether a UV band is already tabulated or needs adding,
   and confirm the specific band definition (FUV? NUV? which FSPS
   filter?) before implementing, since both figures depend on it.

4. Wire one YAML experiment config (extending config.yml's existing
   pattern) that selects IMF + binary-prescription variant and runs
   everything Figure 1/2 need in a single reproducible invocation.

5. Build Figure 1 and Figure 2 using the existing plotting conventions.

6. Follow the exact regression-testing discipline already established
   in tests/test_regression.py and tests/test_evolve.py: pin exact
   values for every new physics path, verify sensitivity by
   deliberately breaking and reverting the code, and extend
   docs/provenance.md with new rows citing the relevant Paper 1
   sections/prescriptions.

7. Do not perturb the existing default behaviour. The current
   fsur-based Power et al. (2009) path (imf_type=2, mmin=0.1, etc.)
   must stay bit-identical to what's already pinned -- the 25 currently
   passing tests should still pass unchanged unless the new
   binary-prescription selector is explicitly turned on.

## Out of scope for this milestone

Defer to later phases per development-roadmap.md: the full Event
taxonomy beyond what Figure 3 needs, the general
experiment-runner/parameter-sweep engine (that's Phase 3, for Paper 2),
FSPS integration beyond the one UV band needed here, the
compact-object/accretion/environment abstractions (the Investigation),
and the response-function/SFH-convolution machinery (Paper 3). Figure 6
(stochastic realisations) is explicitly the bridge to Paper 2 -- note
it, don't build it here.

Work incrementally, explain before implementing anything nontrivial,
and do not commit or push without being asked.
