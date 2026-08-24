# Paper 1 follow-up prompt: from "figures render" to "the physics is right"

**Update, 2026-08-25**: the RLOF-only version of this prompt's core
ask is DONE -- `config.use_post_sn_rlof` (default `False`), a new
"Phase 1.5" in `binaries/population.py::evolve` gated on `nturn==1`,
using the same Eggleton Roche-lobe machinery already built for Phase 0.
See `docs/provenance.md` Section 2 for the full writeup (confirmed
genuinely additive on a realistic population: `lumx_tot` +29% at
t=100 Myr with the channel on vs. off). Deliberately minimal scope,
per this prompt's own item 1's proposal step -- see that writeup for
exactly what was cut (no wind-accretion trigger, no consequence model
for the secondary's mass/envelope or the compact primary's mass).

**Update, 2026-08-25: wind-accretion physics is now implemented as
standalone modules, but NOT yet wired into `evolve()`.** The user
provided both source papers directly (El Mellah & Casse 2017,
arXiv:1609.01532; Friend & Castor 1982, ApJ 261, 293) and a detailed
spec. Both pasted and read section-by-section before implementing --
see `stellar/cak_wind.py` (CAK wind mass-loss rate/velocity law) and
`binaries/wind_capture.py` (BHL-style accretion rate + circularization
radius), and `docs/provenance.md` Section 15 for the full writeup,
including three real transcription/derivation errors caught and fixed
during implementation (two memorized-constant errors in the Eddington-
luminosity calculation, one velocity-law mismatch, all caught by direct
numerical cross-checks against Friend & Castor's own Vela X-1 data
before being accepted). The circularization-radius formula is the one
piece NOT verified against either paper (neither gives a closed-form
fit) -- flagged as the lowest-confidence part of this addition.

Deliberately NOT wired into `BinaryPopulation.evolve()` yet -- built
as standalone, independently-tested physics modules only, per the
user's own "keep the prescription modular" framing. What's needed for
that next step, not yet decided:
- Does wind accretion replace, gate, or run alongside the existing
  `fsur`/`use_post_sn_rlof` channels for a given binary?
- How does `L_X` get computed from `Mdot_acc` -- a new
  `eta*Mdot_acc*c^2`-style conversion (per Friend & Castor eq. 12-13),
  or fed through the existing `xray_calc` luminosity-draw machinery?
- Gating would be similar to the RLOF channel (`nturn==1`) but without
  any Roche-lobe check -- just "secondary is MS-phase and above some
  wind-relevant mass threshold" (`stellar/cak_wind.py`'s functions
  need mass/luminosity/radius/alpha/Q/Gamma per donor -- luminosity
  and radius are already available via
  `main_sequence.ms_luminosity`/`ms_radius`; `alpha`/`Q`/`Gamma`
  are the four CAK/Eddington shape parameters and would need sensible
  defaults or new config fields, following the `Q_CRIT_MS`/`ALPHA_CE`
  precedent of named, overridable constants).

The original scoping text below (for what wind accretion would need,
written before the papers were provided) is now superseded by the
actual implementation above -- kept for the historical record of what
was anticipated vs. what the papers actually gave.

A future session should get the actual
source text for Vink et al. (2001) (and whichever Bondi-Hoyle-
Lyttleton reference is preferred) pasted and verified before
implementing, the same way HTP02/Tout et al. (1997) were for this
session's other additions.

---

A self-contained prompt for the next phase of Paper 1 work. Written so
it can be pasted into a fresh session with no other context. It picks
up after `docs/science/paper1-implementation-prompt.md`'s scope is
fully delivered (Figures 1/2 reproduce from a single YAML config, both
now numerically pinned and non-degenerate -- see
`docs/provenance.md` Sections 6, 10, and the "Known gaps" entry on the
RLOF-classifier pipeline regression tests). That earlier prompt's goal
was "the two figures render, from one config, with real machinery
behind them." This one's goal is "the dominant physics for the mass
regime the figures actually cover is in the model" -- a different,
harder bar.

---

You're working on Realta, a modular stellar/binary population-synthesis
Python framework (repo already checked out) built on a ~20-year-old
Fortran Monte Carlo model of HMXBs in globular clusters (Power et al.
2009, MNRAS 395, 1146, arXiv:0902.1897). Before touching any code, read:

- `docs/provenance.md` -- paper-equation -> implementation -> test
  traceability, and the discipline this project follows for pinning
  numeric regression values. Read Sections 6, 10, 12, 12a in full, and
  the "Known gaps" section at the bottom (both "Closed this session"
  and "Still open") -- they describe exactly where the previous phase
  left off and why.
- `docs/science/rlof-ce-classifier-proposal.md` -- the RLOF/CE
  classifier's full design history, including the "Units bug" and
  emergent-donor-selection findings.
- `docs/science/research-programme.md`, the "Paper 1" section -- the
  scientific target this work serves.
- `docs/science/development-roadmap.md` -- the target architecture;
  this prompt is Phase 2's remainder, not the full 28-item roadmap.

Governing principles (a "Development and Scientific Software Brief"
applies throughout, unchanged from the previous prompt): preserve the
Power et al. (2009) baseline exactly; never silently change scientific
behaviour; flag ambiguity rather than resolve it yourself; explain any
nontrivial physics or design decision before implementing it; avoid
over-engineering; keep the public API small; treat external tools
(MESA, FSPS, etc.) as optional; never commit or push without being
explicitly asked; prefer small, incremental, reviewable changes.

## Why this prompt exists

Running `scripts/run_paper1_experiment.py` end-to-end (not just unit
tests -- see `docs/provenance.md`'s "no automated test exercises the
RLOF pipeline end to end" entry, now closed by
`tests/test_paper1_pipeline_regression.py`) surfaced a real, structural
finding, not a bug: for the mass range this experiment actually probes
(`mcut=8` Msun primaries paired with companions from `mcomp=0.5` Msun
up), the RLOF/CE classifier's pre-SN `STABLE_MASS_TRANSFER` channel can
essentially never fire. `generate_population` always enforces
`m2 <= m1`, and that outcome requires the donor to be the *lighter*
star -- so its donor is always the low-mass companion, whose own
pre-SN main-sequence+HG lifetime (hundreds of Myr for a 2-7 Msun star)
is always far longer than the massive primary's (a few Myr). The
primary has essentially always already exploded, flipping `nturn` to 1
and permanently gating Phase 0 out (`nturn==0` required), before this
channel becomes reachable. Confirmed directly, not just suspected: zero
of 51 `STABLE_MASS_TRANSFER`-classified binaries in the pinned Paper 1
config were ever processed across a full 100 Myr run.

This means `interaction_boost` (`standard_interaction`/
`enhanced_interaction`) currently has **no measurable effect on Figure
1/2's output** -- not because it's wired wrong, but because the
physical channel it's meant to boost doesn't exist yet for this mass
regime. The astrophysically dominant real HMXB-formation channel --
the *secondary* star's own, later Roche-lobe overflow onto the
by-then-compact primary (Case B/C mass transfer onto a neutron
star/black hole) -- is not modelled at all. Phase 0 was deliberately
scoped to pre-SN interaction between two still-live stars only,
matching HTP02's own CE-eligible donor list; this prompt is about
building the missing post-SN piece.

## Goal

Make `standard_interaction`/`enhanced_interaction`/`enhanced_mergers`
actually diverge from `non_interacting` for the *right physical
reason* -- genuine post-SN mass transfer driving HMXB activation, not
merely a pre-SN merger/detached-outcome split (which is all that
currently differentiates them). Concretely: extend `BinaryPopulation`
with a second Roche-lobe-overflow channel, evaluated after the primary
has already become a compact remnant (`nturn==1`), where the *still-
evolving secondary* can fill its Roche lobe and transfer mass onto (or
wind-feed) the compact primary.

## Scope

1. **Stop and flag before implementing.** This is new physics, not an
   extension of an existing prescription -- the brief's "flag
   ambiguity" principle applies especially hard here. Before writing
   any code, produce a short written proposal (following the same
   pattern as `docs/science/paper1-binary-interaction-proposal.md` and
   the RLOF/CE proposal doc) covering at minimum:
   - What triggers this channel: does the secondary need to be
     classified via the same `classify_rlof`-style Roche-lobe check
     (now against a *compact* primary's mass/radius-less point source,
     not another star's stellar radius), or is a coarser criterion
     (e.g. period below some threshold, following Power et al. 2009's
     own `fsur`/`floss` spirit) more appropriate given what's already
     in the model?
   - Wind-accretion vs. Roche-lobe accretion: real HMXBs span both
     regimes (persistent wind-fed vs. transient RLOF-fed). Does this
     milestone need to distinguish them, or is a single simplified
     "secondary interacts" gate enough for Paper 1's purposes (check
     `docs/science/research-programme.md`'s stated scope for Paper 1
     specifically before assuming detailed accretion physics is
     needed)?
   - How this interacts with the *existing* `fsur`-based HMXB
     activation gate (Power et al. 2009, Sec. 2.1) -- is the new
     channel a *replacement* for `fsur`, an additional *gate* on top of
     it, or a *separate* activation pathway alongside it? This is
     exactly the kind of quantitative-prescription decision the
     original Paper 1 prompt's item 2 flagged for `standard_interaction`
     etc. -- treat it with the same care.
   - What a "compact primary's Roche lobe" even means for the Eggleton
     fit (`binaries/interaction.py::roche_lobe_radius`) currently used
     -- it's donor-radius-agnostic (just needs `separation`/`q1`), so
     the formula itself likely still applies; the donor-selection and
     phase-tracking logic (`find_rlof_onset`, `classify_rlof`) will
     need real changes to handle "companion has no stellar radius,
     donor is whatever `main_sequence`/`giant_branch` phase the
     secondary has reached by then."

2. **Reconcile, or explicitly decide not to reconcile, the two
   independent stellar-lifetime prescriptions.** `turnoff_time`
   (Schaerer et al. 1993, pre-existing) drives the primary's SN;
   `t_ms`/`t_bgb` (Hurley et al. 2000, this session's addition) drives
   RLOF timing. They agree to within ~5% for the *same* star (already
   quantified, see `docs/provenance.md` Section 6), so this is a
   secondary concern relative to item 1's structural mass-hierarchy
   issue -- but a post-SN channel needs the *secondary's* lifetime
   clock too, and right now `t2_lifetime` (Schaerer-based) and the
   Hurley/Tout functions are two independent sources for what should
   be the same star's evolutionary state. Decide whether the new
   channel uses `t2_lifetime` (consistent with the existing SN-timing
   convention) or the Hurley/Tout functions (consistent with the RLOF
   classifier), and document the choice.

3. **Extend `BinaryPopulation`'s Phase structure.** The current
   `evolve()` has Phase 0 (pre-SN RLOF), Phase 1 (primary SN), Phase 2
   (secondary SN), Phase 3 (aggregate observables). A new post-SN RLOF
   channel most naturally sits between Phase 1 and Phase 2 -- gated on
   `nturn==1` (primary already compact, secondary not yet exploded).
   Reuse `apply_stable_mass_transfer`/`apply_common_envelope`'s
   *pattern* (instantaneous, energy/angular-momentum-conservative
   simplifications, not rate-integrated) rather than inventing a third
   style, unless the proposal in item 1 finds a specific reason not to.

4. **Consequence model.** What happens to the compact primary's mass
   when it accretes (does it grow? is there an Eddington cap, reusing
   `xray/luminosity.py`'s existing Eddington-limited draw machinery?),
   and what happens to the secondary (does it lose its envelope the
   same way an HG donor does in the existing CE model, or is a wind-fed
   regime non-destructive to the donor entirely)? This determines
   whether `merge_stellar_masses`-style bookkeeping applies here too,
   or whether a compact-primary accretion event needs its own,
   separate consequence function.

5. **Follow the exact regression-testing discipline already
   established** (`tests/test_regression.py`,
   `tests/test_paper1_pipeline_regression.py`,
   `tests/test_ce_energy_balance.py` for the pattern): pin exact
   values for every new physics path, verify sensitivity by
   deliberately breaking and reverting the code, and extend
   `docs/provenance.md` with new rows citing the relevant sections.
   Specifically add a regression case that would have caught the
   finding this prompt is about -- i.e. one that fails loudly if
   `interaction_boost` ever again has zero measurable effect on
   `standard_interaction` vs. `enhanced_interaction` output.

6. **Do not perturb existing default behaviour.** The current
   `fsur`-based Power et al. (2009) path, and every already-pinned
   RLOF-classifier regression value, must stay bit-identical unless the
   new post-SN channel is explicitly opted into (a new config flag,
   following the `use_rlof_classifier` pattern -- default `False`/a
   no-op for every prescription unless a decision in item 1 says
   otherwise).

## Out of scope for this milestone

- Distinguishing wind-fed vs. Roche-lobe-fed HMXBs at the spectral/
  luminosity level (i.e. changing `xray/luminosity.py`'s own
  luminosity-draw physics) -- unless item 1's proposal finds Paper 1
  specifically needs this distinction, treat the new channel as
  feeding the *same* existing X-ray luminosity draw, just gated by a
  different (physically motivated) activation criterion.
- Full binary-evolution accretion-disk/spin physics.
- Figure 3 (mergers vs. compact-object formation), Figures 4-6, the
  full `Event`/`PopulationHistory` abstraction from
  `docs/science/development-roadmap.md` item 4, an `examples/`
  directory, Brček, Hirai, Mandel & Lower (2026) rejuvenation physics,
  and Hovis-Afflerbach et al. (2025) stripped-donor properties --  all
  already-named, deliberately deferred items from the original Paper 1
  prompt and `docs/provenance.md`'s "Known gaps" section; nothing here
  changes their status.

Work incrementally, explain before implementing anything nontrivial,
and do not commit or push without being asked.
