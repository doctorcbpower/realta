#!/usr/bin/env python3
"""B2 (docs/science/paper1-detailed-work-breakdown.md): cross-check
Realta's RLOF/CE classifier's population-level outcome fractions
against Xu et al. (2025, A&A 704, A218, arXiv:2503.23876)'s SMC
statistics -- 8% post-mass-transfer, 7% merger, for a sample selected
with M1=5-100 Msun, q=0.3-0.95, P=1-3162 d. Explicitly a sanity check,
not an exact-match requirement (per the task's own wording) -- Realta
and Xu et al. are different codes with different physics, evolved to
different endpoints (see the caveats printed below).

Methodology (a named, documented choice, not itself from the paper):
Realta has no generation-time q_min/q_max sampling knob, so this
script generates a broad population (mcut=5, mmax=100, pmin=1,
pmax=3162, mass_ratio_distribution="flat_q" -- i.e. q ~ Uniform(0,1)
before filtering) and then POST-FILTERS to Xu et al.'s exact selection
window (0.3 <= q <= 0.95), which is a more faithful reproduction of
"a comparable setup" than trying to force the generation-time
distribution to match. Outcome classification uses
find_rlof_onset()'s precomputed classification at generation time;
COMMON_ENVELOPE-classified systems are additionally resolved via
apply_common_envelope() to split them into "survives" (counted as
post-mass-transfer, alongside STABLE_MASS_TRANSFER) vs. "merges"
(counted as merger, alongside IMMEDIATE_MERGER) -- this is not
optional bookkeeping, it is what "merger" vs. "post-mass-transfer"
actually means physically. Fractions are reported relative to the
full q/P-selected sample (including DETACHED/never-interacts and
PHASE_NOT_MODELLED systems in the denominator), matching how such
population-synthesis fractions are conventionally reported.

Known, deliberate mismatches with Xu et al.'s actual setup (not
attempting to eliminate these -- this is a sanity check):
- Metallicity: Realta only has three metallicity presets (Z=0, 0.008,
  0.02). SMC's true metallicity is Z~0.002-0.004 (~0.1-0.2 Zsun) --
  imetal=2 (Z=0.008, ~0.4 Zsun) is the closer of the two usable
  presets, but is not a true SMC metallicity match.
- Realta's classifier scope: MS/HG donors only (see
  binaries/interaction.py's module docstring) -- GB/AGB donors, which
  a real SMC population over Xu et al.'s full mass/period range would
  include, fall through as PHASE_NOT_MODELLED here rather than being
  classified at all.
- Realta's CE consequence model omits HTP02 eqs. 74-77 (partial
  envelope retention) -- see docs/provenance.md Section 12a.
"""

from __future__ import annotations

import numpy as np

from realta.binaries.interaction import RLOFOutcome, apply_common_envelope
from realta.binaries.population import BinaryPopulation
from realta.config import SimulationConfig

Z_IMETAL2 = 0.008  # see module docstring's metallicity caveat


def main(ntot: int = 200_000):
    config = SimulationConfig(
        ntot=ntot,
        imf_type=2,
        mmin=0.1,
        mmax=100.0,
        mcut=5.0,
        pmin=1.0,
        pmax=3162.0,
        mcomp=0.1,
        mass_ratio_distribution="flat_q",
        period_distribution="log_uniform",
        use_rlof_classifier=True,
        imetal=2,
        iseed=42,
    )
    pop = BinaryPopulation(config)

    q = pop.m2 / np.where(pop.m1 > 0, pop.m1, np.nan)
    selected = (q >= 0.3) & (q <= 0.95)
    n_selected = int(np.sum(selected))
    print(
        f"Generated {len(pop.m1)} massive binaries; {n_selected} pass "
        f"Xu et al.'s q in [0.3, 0.95] selection window."
    )

    n_post_mt = 0
    n_merger = 0
    n_detached = 0
    n_not_modelled = 0

    for i in np.where(selected)[0]:
        outcome = pop.rlof_outcome[i]
        if outcome == RLOFOutcome.STABLE_MASS_TRANSFER:
            n_post_mt += 1
        elif outcome == RLOFOutcome.IMMEDIATE_MERGER:
            n_merger += 1
        elif outcome == RLOFOutcome.DETACHED:
            n_detached += 1
        elif outcome == RLOFOutcome.PHASE_NOT_MODELLED:
            n_not_modelled += 1
        elif outcome == RLOFOutcome.COMMON_ENVELOPE:
            donor_is_1 = pop.rlof_donor_is_star1[i]
            donor_mass = pop.m1[i] if donor_is_1 else pop.m2[i]
            companion_mass = pop.m2[i] if donor_is_1 else pop.m1[i]
            separation_rsun = pop.a[i] * BinaryPopulation.RSUN_PER_AU
            survives, _, _, _ = apply_common_envelope(
                donor_mass, companion_mass, separation_rsun, Z_IMETAL2, pop.rlof_time[i]
            )
            if survives:
                n_post_mt += 1
            else:
                n_merger += 1

    frac_post_mt = n_post_mt / n_selected
    frac_merger = n_merger / n_selected

    print()
    print(f"{'Outcome':<20} {'count':>8} {'fraction':>10}")
    print(f"{'post-mass-transfer':<20} {n_post_mt:>8} {frac_post_mt:>10.3%}")
    print(f"{'merger':<20} {n_merger:>8} {frac_merger:>10.3%}")
    print(f"{'detached':<20} {n_detached:>8} {n_detached / n_selected:>10.3%}")
    print(
        f"{'not_modelled':<20} {n_not_modelled:>8} {n_not_modelled / n_selected:>10.3%}"
    )
    print()
    print("Xu et al. (2025) SMC statistics: post-mass-transfer=8%, merger=7%")
    print(
        f"Realta (this setup):             post-mass-transfer={frac_post_mt:.1%}, "
        f"merger={frac_merger:.1%}"
    )


if __name__ == "__main__":
    main()
