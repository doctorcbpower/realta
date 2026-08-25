"""Roche-lobe-overflow (RLOF) outcome classifier.

Stage 2 of the RLOF/CE module described in
docs/science/rlof-ce-classifier-proposal.md. Given a binary's masses,
separation and metallicity, and the donor's age, this classifies
whether Roche-lobe overflow is occurring and, if so, what the outcome
is -- for MS (k=0,1) and HG (k=2, within `stellar/main_sequence.py`'s
supported mass range) donors. GB (k=3) and later phases remain out of
scope; see that module's docstring.

Source: Hurley, Tout & Pols (2002, MNRAS 329, 897), hereafter HTP02.

Scope-driven simplification: HTP02's common-envelope-eligible donor
list (Section 2.7.1: "at the onset of RLOF where mass is transferred
from a giant (k1 in {2,3,4,5,6,8,9}) ... to a main-sequence star") does
NOT include MS donors (k1 in {0,1}). A dynamically unstable MS donor
therefore merges directly (Section 2.6.4: "mass transfer to a
companion proceeds dynamically if q1>0.695 ... only a single star
remains") rather than forming a common envelope. HG donors ARE
CE-eligible, so `classify_rlof` returns COMMON_ENVELOPE (not
IMMEDIATE_MERGER) for a dynamically unstable HG donor, using
`hg_q_crit` (HTP02's own GB q_crit formula, reused for HG per Zuo &
Li 2014). The CE outcome itself (survive vs. merge, and the resulting
mass/orbit) is not resolved here -- that needs the alpha-lambda
energy-balance solve (HTP02 eqs. 69-77), not yet implemented; see the
proposal doc.

Emergent finding (`find_rlof_onset` below): because the Eggleton
Roche-lobe fraction R_L1/a increases monotonically with the donor's
own mass ratio q1, `IMMEDIATE_MERGER` dominates in practice over
`STABLE_MASS_TRANSFER` for automatically-detected MS-MS donors -- the
star that reaches its (proportionally larger) Roche lobe first is
almost always the more massive one, so q1 > 1 > Q_CRIT_MS for the
detected donor in most configurations. `classify_rlof()` itself still
returns `STABLE_MASS_TRANSFER` correctly when a low-q1 donor is
specified directly; this is about which donor `find_rlof_onset`
actually selects, not a bug in the classification itself -- see
`tests/test_rlof_classifier.py::test_find_rlof_onset_favours_the_more_massive_star_as_donor`.
"""

from __future__ import annotations

from enum import Enum

import numpy as np
from scipy.optimize import brentq

from realta.stellar import giant_branch, main_sequence, remnant

# HTP02 Section 2.6.4: critical mass ratio (q1 = M_donor/M_companion)
# above which mass transfer from a main-sequence donor proceeds on a
# dynamical time-scale, ending in a single merged star. Stated there
# specifically for deeply-convective k1=0 donors; extended here, as a
# named simplification, to radiative k1=1 donors too -- see the
# proposal doc's "Decision 2" for why. Exposed as a parameter (not a
# buried constant) so it can be overridden/recalibrated.
Q_CRIT_MS = 0.695

# Common-envelope efficiency and binding-energy parameter (HTP02 eqs.
# 69-73). Decided (chat, 2026-08-24) after a literature review -- see
# docs/science/rlof-ce-classifier-proposal.md's "Literature findings"
# note: ALPHA_CE=0.9 is Zuo & Li (2014, MNRAS 442, 1980)'s own
# basic-model value, sitting in the middle of their HMXB-population-
# calibrated 0.8-1.0 preferred range (a materially better-targeted
# source than HTP02's own generic "alpha_CE~=1 is used"). LAMBDA_CE=0.5
# is HTP02's own fixed value (eq. 69) -- explicitly flagged there as
# not a true constant, and its accepted value here is entangled with
# Zuo & Li's alpha_CE range (their own caveat: differing core-boundary
# definitions can shift lambda by up to two orders of magnitude, which
# would shift the implied alpha_CE too). Both exposed as named,
# overridable parameters, not buried constants.
ALPHA_CE = 0.9
LAMBDA_CE = 0.5


class RLOFOutcome(str, Enum):
    """Roche-lobe-overflow classification outcome.

    numpy gotcha, confirmed empirically during this module's
    development: scalar comparison (`some_outcome == RLOFOutcome.X`,
    including `array[i] == RLOFOutcome.X` for a single element pulled
    out of an object-dtype array) works correctly, but a *vectorized*
    comparison of a whole numpy object array against a bare member
    (`array_of_outcomes == RLOFOutcome.X`) silently returns all-False
    even when every element genuinely equals it -- numpy appears to
    broadcast the str-Enum member as if it were itself a
    sequence/array rather than a scalar. Use a per-element comparison
    (list/generator comprehension, or a plain Python loop) instead of
    a vectorized one when checking an array of RLOFOutcome values --
    see `binaries/population.py`'s own per-element loops, which are
    unaffected, and `tests/test_rlof_wiring.py` for the pattern to
    use instead of `array == RLOFOutcome.X`.
    """

    DETACHED = "detached"
    STABLE_MASS_TRANSFER = "stable_mass_transfer"
    IMMEDIATE_MERGER = "immediate_merger"
    COMMON_ENVELOPE = "common_envelope"
    PHASE_NOT_MODELLED = "phase_not_modelled"


def roche_lobe_radius(separation: float, mass_ratio: float) -> float:
    """Roche-lobe radius of the star with mass ratio q1 = M1/M2.

    Eggleton (1983) fit, as used by HTP02 eq. (53):
        R_L1/a = 0.49 q1^(2/3) / (0.6 q1^(2/3) + ln(1 + q1^(1/3)))
    Accurate to within 1 per cent for 0 < q1 < infinity (HTP02's own
    stated accuracy for eq. 53). `mass_ratio` = M_this_star/M_other_star
    -- call with q2 = M2/M1 to get the companion's Roche-lobe radius
    instead (HTP02's own note directly below eq. 53).
    """
    if mass_ratio <= 0.0:
        raise ValueError(f"mass_ratio must be positive, got {mass_ratio}")
    q_cbrt = mass_ratio ** (1.0 / 3.0)
    q_two_thirds = mass_ratio ** (2.0 / 3.0)
    return (
        separation * 0.49 * q_two_thirds / (0.6 * q_two_thirds + np.log(1.0 + q_cbrt))
    )


def hg_q_crit(donor_mass: float, z: float, donor_age: float) -> float:
    """Critical mass ratio (q1 = M_donor/M_companion) above which mass
    transfer from an HG donor is dynamically unstable, forming a
    common envelope rather than transferring mass stably.

        q_crit = [1.67 - x + 2*(M_c1/M1)^5] / 2.13

    where x is the GB mass-radius exponent (giant_branch.py, eq. 47)
    and M_c1 is the donor's core mass. This is HTP02's own GB q_crit
    formula (eqs. 56-57), reused here for HG donors following Zuo & Li
    (2014, MNRAS 442, 1980)'s eq. 1, which cites Shao & Li (in prep.)
    for exactly this extension -- HTP02 itself uses a crude fixed
    `q_crit=4` for HG donors, which it calls "rather approximate."
    Since the formula and its constants (1.67, 2.13, the exponent 5)
    are already trusted from HTP02, this needed no new coefficient-
    verification round -- see
    docs/science/rlof-ce-classifier-proposal.md's literature-findings
    note. Unlike Q_CRIT_MS, this is not exposed as a single overridable
    scalar (it is itself a formula, not one fixed number) -- override
    `alpha_CE`/`lambda_CE`/`Q_CRIT_MS`-style parameters do not apply
    here.

    Propagates whatever ValueError `main_sequence.core_mass_hg` raises
    (e.g. mass outside the M_HeF <= M < CORE_MASS_BGB_MAX_MASS range)
    -- callers must be prepared for that, the same as for
    `hg_radius`/`hg_luminosity`.
    """
    core_mass = main_sequence.core_mass_hg(donor_mass, z, donor_age)
    x = giant_branch.mass_radius_exponent(z)
    return (1.67 - x + 2.0 * (core_mass / donor_mass) ** 5) / 2.13


def classify_rlof(
    donor_mass: float,
    companion_mass: float,
    separation: float,
    z: float,
    donor_age: float,
    q_crit_ms: float = Q_CRIT_MS,
) -> RLOFOutcome:
    """Classify the RLOF state/outcome for a binary at a given moment.

    `donor_mass`/`donor_age` describe the star being tested as the
    potential donor (HTP02's "primary", M1) -- call twice, once per
    star, to check both. Returns PHASE_NOT_MODELLED if the donor's
    phase isn't reachable at all (`main_sequence.phase()` raises), or
    if it's an HG donor whose radius/core mass this module can't
    compute (M >= M_FGB, or M >= CORE_MASS_BGB_MAX_MASS -- see
    `main_sequence.hg_radius`/`core_mass_hg`).

    MS donors (k=0,1): dynamically unstable RLOF (q1 > q_crit_ms)
    merges immediately -- MS donors are not in HTP02's CE-eligible
    donor list (Sec. 2.7.1). HG donors (k=2): dynamically unstable
    RLOF (q1 > hg_q_crit(...)) forms a COMMON_ENVELOPE instead -- HG
    donors *are* CE-eligible. The CE outcome (survive vs. merge) is
    not resolved here -- that's the energy-balance solve, not yet
    implemented (see the proposal doc).
    """
    try:
        donor_phase = main_sequence.phase(donor_mass, z, donor_age)
    except ValueError:
        return RLOFOutcome.PHASE_NOT_MODELLED

    q1 = donor_mass / companion_mass
    r_l1 = roche_lobe_radius(separation, q1)

    if donor_phase in (0, 1):
        donor_radius = main_sequence.ms_radius(donor_mass, z, donor_age)
        if donor_radius < r_l1:
            return RLOFOutcome.DETACHED
        if q1 > q_crit_ms:
            return RLOFOutcome.IMMEDIATE_MERGER
        return RLOFOutcome.STABLE_MASS_TRANSFER

    # donor_phase == 2 (HG) -- the only other value phase() returns.
    try:
        donor_radius = main_sequence.hg_radius(donor_mass, z, donor_age)
        q_crit = hg_q_crit(donor_mass, z, donor_age)
    except ValueError:
        return RLOFOutcome.PHASE_NOT_MODELLED

    if donor_radius < r_l1:
        return RLOFOutcome.DETACHED
    if q1 > q_crit:
        return RLOFOutcome.COMMON_ENVELOPE
    return RLOFOutcome.STABLE_MASS_TRANSFER


def find_rlof_onset(
    m1: float,
    m2: float,
    separation: float,
    z: float,
    epsilon: float = 1e-6,
    q_crit_ms: float = Q_CRIT_MS,
) -> tuple[float, RLOFOutcome, bool]:
    """Find the earliest RLOF onset for a binary across the MS and HG,
    checking both stars as the potential donor.

    Named `find_rlof_onset` (not `find_rlof_onset`, an earlier,
    now-inaccurate name) since it searches both phases -- see the
    "Update, 2026-08-24" note below for when HG search was added.

    Returns (t_rlof, outcome, donor_is_star1):
        t_rlof: time (Myr) the earlier-overflowing star first fills its
            Roche lobe, or np.inf if neither star does so while on the
            MS or HG.
        outcome: RLOFOutcome.DETACHED if t_rlof is inf; RLOFOutcome.
            PHASE_NOT_MODELLED if a crossing was found but this
            module can't classify it (an HG donor whose mass is
            outside `core_mass_hg`'s supported range -- radius is
            still computable there, via `M_FGB`, but core mass/q_crit
            is not, via `CORE_MASS_BGB_MAX_MASS`); otherwise the
            `classify_rlof()`-equivalent result at t_rlof.
        donor_is_star1: which star is the donor (True -> m1/star 1).

    Root-finding relies on R_donor(t) being monotonically increasing
    across the MS (confirmed by
    tests/test_hurley_main_sequence.py::test_ms_radius_and_luminosity_monotonically_increase)
    and, separately, across the HG (confirmed by
    tests/test_hertzsprung_gap.py), with continuity at the MS/HG
    boundary (R_HG(t_MS) == R_TMS exactly, by construction of eq. 27) --
    so at most one crossing exists per star across the whole MS+HG
    track, and a single bisection per phase actually searched is
    sufficient.

    Note on lifetime inconsistency: the search upper bound is Hurley
    et al. (2000)'s own `t_BGB(mass, z)`, which is NOT necessarily
    identical to Realta's separate, pre-existing `LifetimeTable`
    (Schaerer et al. 1993-based) used for the primary/secondary
    supernova timing elsewhere in `BinaryPopulation`. The two
    stellar-lifetime prescriptions are independently sourced and not
    reconciled -- see docs/provenance.md. In practice this means a
    predicted RLOF event can, for a given binary, fall after that
    binary's Schaerer-table-based supernova already occurred; the
    caller (`BinaryPopulation.evolve()`) naturally suppresses such
    stale predictions via its own `nturn == 0` gate rather than this
    function trying to account for it.

    Update, 2026-08-24: extended to also search the HG, using
    `hg_radius()`. GB (k=3) and later phases remain unreachable -- a
    donor that never overflows across MS+HG is treated as never
    overflowing at all within this module's scope, even though a real
    star would continue evolving past the GB.
    """

    def _radius_root(donor_mass: float, companion_mass: float) -> float | None:
        q1 = donor_mass / companion_mass
        r_l1 = roche_lobe_radius(separation, q1)

        r_at_start = main_sequence.ms_radius(donor_mass, z, epsilon)
        if r_at_start >= r_l1:
            return 0.0

        t_ms_donor = main_sequence.t_ms(donor_mass, z)
        t_ms_end = t_ms_donor * (1.0 - epsilon)
        r_at_ms_end = main_sequence.ms_radius(donor_mass, z, t_ms_end)
        if r_at_ms_end >= r_l1:

            def f_ms(t: float) -> float:
                return main_sequence.ms_radius(donor_mass, z, t) - r_l1

            return brentq(f_ms, epsilon, t_ms_end)

        # No MS crossing -- check the HG, if this module can compute a
        # radius for this donor there (raises for M >= M_FGB).
        try:
            t_bgb_donor = main_sequence.t_bgb(donor_mass, z)
            hg_duration = t_bgb_donor - t_ms_donor
            t_hg_start = t_ms_donor + epsilon * hg_duration
            t_hg_end = t_ms_donor + (1.0 - epsilon) * hg_duration
            r_at_hg_start = main_sequence.hg_radius(donor_mass, z, t_hg_start)
            r_at_hg_end = main_sequence.hg_radius(donor_mass, z, t_hg_end)
        except ValueError:
            return None  # HG radius not computable for this donor (M >= M_FGB)

        if r_at_hg_start >= r_l1:
            return t_hg_start  # crosses essentially right at HG onset
        if r_at_hg_end < r_l1:
            return None  # never overflows within MS+HG (GB+ out of scope)

        def f_hg(t: float) -> float:
            return main_sequence.hg_radius(donor_mass, z, t) - r_l1

        return brentq(f_hg, t_hg_start, t_hg_end)

    candidates: list[tuple[float, float, float, bool]] = []
    t1 = _radius_root(m1, m2)
    if t1 is not None:
        candidates.append((t1, m1, m2, True))
    t2 = _radius_root(m2, m1)
    if t2 is not None:
        candidates.append((t2, m2, m1, False))

    if not candidates:
        return np.inf, RLOFOutcome.DETACHED, True

    t_cross, donor_mass, companion_mass, donor_is_star1 = min(
        candidates, key=lambda c: c[0]
    )
    # Determine the outcome directly from q1/phase rather than
    # re-querying classify_rlof() at t_cross itself: at the exact root,
    # donor_radius == r_l1 to within floating-point noise, and
    # classify_rlof()'s strict "<" detached check can tip the wrong way
    # right at the boundary. The outcome only depends on the donor's
    # phase and mass ratio once RLOF has been confirmed to occur (see
    # classify_rlof()'s own logic), so this is equivalent and avoids
    # that edge case.
    q1 = donor_mass / companion_mass
    donor_phase = main_sequence.phase(donor_mass, z, t_cross)
    if donor_phase in (0, 1):
        outcome = (
            RLOFOutcome.IMMEDIATE_MERGER
            if q1 > q_crit_ms
            else RLOFOutcome.STABLE_MASS_TRANSFER
        )
    else:  # HG
        try:
            q_crit = hg_q_crit(donor_mass, z, t_cross)
        except ValueError:
            # Radius was computable (M < M_FGB) but core mass/q_crit
            # is not (M >= CORE_MASS_BGB_MAX_MASS) -- same gap
            # classify_rlof() itself reports as not modelled.
            outcome = RLOFOutcome.PHASE_NOT_MODELLED
        else:
            outcome = (
                RLOFOutcome.COMMON_ENVELOPE
                if q1 > q_crit
                else RLOFOutcome.STABLE_MASS_TRANSFER
            )
    return t_cross, outcome, donor_is_star1


def _widened_separation(
    separation_i: float, m1_i: float, m2_i: float, m1_f: float, m2_f: float
) -> float:
    """Orbital separation after conservative mass transfer.

    Conservative transfer conserves total mass and orbital angular
    momentum, L = M1*M2*sqrt(G*a/(M1+M2)) (standard two-body reduced
    formula; not itself an HTP02-specific result -- basic Keplerian
    mechanics). With M1+M2 constant, L=const gives
    a_f = a_i * (M1_i*M2_i / (M1_f*M2_f))^2.
    """
    return separation_i * (m1_i * m2_i / (m1_f * m2_f)) ** 2


def apply_stable_mass_transfer(
    donor_mass: float,
    companion_mass: float,
    separation: float,
    donor_radius: float,
) -> tuple[float, float, float]:
    """Instantaneous conservative stable mass transfer.

    Not itself an HTP02 prescription -- HTP02 rate-integrates mass
    transfer via Kelvin-Helmholtz/nuclear time-scales (eqs. 58-61),
    which does not fit Realta's instantaneous-event architecture (see
    docs/science/rlof-ce-classifier-proposal.md's "Decision" on
    instantaneous vs. rate-integrated treatment). This is a named
    simplification: mass moves from donor to companion, conservatively
    (total mass and orbital angular momentum both conserved -- see
    `_widened_separation`), until the widened orbit's Roche-lobe
    radius for the donor exactly equals its current radius --
    physically, the point at which the system would newly detach.
    `donor_radius` is treated as fixed during this instantaneous jump
    (no post-mass-loss stellar-structure response is modelled --
    that is the Brček, Hirai, Mandel & Lower (2026) concern named in
    the task brief but not available in this session; see the
    proposal doc).

    Only reachable for donor_mass < companion_mass (q1 < 1) -- the
    only case classify_rlof() ever labels STABLE_MASS_TRANSFER, since
    q1 > q_crit_ms (> ... in practice mostly > 1, see
    find_rlof_onset's module-level note) is IMMEDIATE_MERGER
    instead. For q1 < 1, conservative transfer widens the orbit (mass
    moves from the lighter to the heavier star), growing the Roche
    lobe until it exceeds the donor's (unchanged) radius -- this is
    exactly why q1 < q_crit_ms is the *stable* regime.

    Returns (new_donor_mass, new_companion_mass, new_separation).
    """
    if donor_mass >= companion_mass:
        raise ValueError(
            "apply_stable_mass_transfer expects donor_mass < companion_mass "
            f"(got donor={donor_mass}, companion={companion_mass}) -- this is "
            "the only regime classify_rlof() labels STABLE_MASS_TRANSFER."
        )

    def f(delta_m: float) -> float:
        m1_f = donor_mass - delta_m
        m2_f = companion_mass + delta_m
        a_f = _widened_separation(separation, donor_mass, companion_mass, m1_f, m2_f)
        r_l1_f = roche_lobe_radius(a_f, m1_f / m2_f)
        return r_l1_f - donor_radius

    epsilon = 1e-9 * donor_mass
    delta_m = brentq(f, epsilon, donor_mass - epsilon)

    new_donor_mass = donor_mass - delta_m
    new_companion_mass = companion_mass + delta_m
    new_separation = _widened_separation(
        separation, donor_mass, companion_mass, new_donor_mass, new_companion_mass
    )
    return new_donor_mass, new_companion_mass, new_separation


def rejuvenate_ms_gainer(
    mass_before: float, mass_after: float, age_before: float, z: float
) -> float:
    """Fractional main-sequence-lifetime remaining after an MS star
    gains mass via stable mass transfer (B3, docs/science/
    paper1-detailed-work-breakdown.md).

    Source: Tout, Aarseth, Pols & Eggleton (1997, MNRAS 291, 732),
    Sec. 5.1 "Rejuvenation", eq. (41), verified directly against the
    paper (not from memory) before implementing:

        t' = (mu/mu') * (tau'_MS / tau_MS) * t

    where `t`/`t'` are the star's age just before/after the mass gain,
    `tau_MS`/`tau'_MS` are its MS lifetime at the old/new mass, and
    `mu = M` (old mass), `mu' = M'` (new mass) -- EXCEPT `mu'` is
    redefined to equal `mu` (collapsing the mass-ratio factor to 1)
    when `0.3 < M/Msun < 1.3` (radiative core; the range-check uses
    the OLD mass). Outside that range (convective or fully-convective
    core), the surviving `M/M'` factor makes the star appear younger
    still than fractional-age preservation alone would give -- Tout et
    al. attribute this to unburnt hydrogen being mixed into the
    convective core by the incoming material (their Sec. 5.1 text,
    citing Sandage 1953's blue-straggler mechanism). This is the same
    formula HTP02 (2002) Sec. 2.6.6.1 cites (via Hurley, Pols & Tout
    2000 Sec. 7.1) for the radiative-core/HG case, extended here to
    the convective-core case using Tout et al.'s own eq. (41) directly
    (HTP02 itself does not reproduce that part, only cites it).

    `tau_MS` is HTP02's own `main_sequence.t_ms(mass, z)` (needed for
    self-consistency with eq. 41, which is itself part of the same
    Hurley/Tout fitting-formula family) -- NOT Realta's separate,
    Schaerer-based `LifetimeTable` used elsewhere for the star's actual
    supernova-triggering clock (see docs/provenance.md's "two
    independent lifetime prescriptions" note). Returns a dimensionless
    remaining-lifetime FRACTION (not an absolute time), for the caller
    to apply against whatever lifetime source it uses for absolute
    timing -- see `binaries/population.py::evolve`'s STABLE_MASS_TRANSFER
    branch, which applies this fraction to `LifetimeTable.get_lifetime`
    at the new mass, matching the existing post-interaction-reset
    convention (Hurley-sourced *fraction*, Schaerer-sourced *absolute*
    scale) rather than switching the star onto a third, inconsistent
    timing system.

    Assumes `age_before` is the star's true age since formation (i.e.
    it has not previously had its lifetime clock reset by an earlier
    interaction) -- true for Realta's current scope, since an RLOF
    event is only ever processed once per binary (`rlof_processed`).
    Clamped to a small positive floor rather than exactly 0, to avoid
    handing the caller a zero-duration remaining phase.
    """
    t_ms_before = main_sequence.t_ms(mass_before, z)
    t_ms_after = main_sequence.t_ms(mass_after, z)

    if 0.3 < mass_before < 1.3:
        mass_ratio_factor = 1.0
    else:
        mass_ratio_factor = mass_before / mass_after

    age_after = mass_ratio_factor * (t_ms_after / t_ms_before) * age_before
    remaining_fraction = 1.0 - age_after / t_ms_after
    return max(1e-6, min(1.0, remaining_fraction))


def apply_common_envelope(
    donor_mass: float,
    companion_mass: float,
    separation: float,
    z: float,
    donor_age: float,
    alpha_ce: float = ALPHA_CE,
    lambda_ce: float = LAMBDA_CE,
) -> tuple[bool, float, float, float | None]:
    """Resolve a common-envelope event for an HG donor with an MS
    companion (HTP02 Sec. 2.7.1, eqs. 69-73).

    Returns (survives, new_donor_mass, new_companion_mass,
    new_separation):
        survives=True: the envelope is fully ejected before the cores
            touch. The donor is stripped to its bare core
            (new_donor_mass = M_c1), the companion is unaffected
            (new_companion_mass = companion_mass), new_separation is
            the post-CE orbit (a_f, eq. 72) -- always tighter than the
            input `separation`, since CE is an inspiral.
        survives=False: the cores coalesce before the envelope is
            fully ejected. new_separation is None (there is no
            surviving binary); the caller should merge instead --
            new_donor_mass is still the donor's core mass (M_c1),
            since the merge should combine that with the companion,
            not the donor's pre-CE mass (the envelope has already
            been [partially] stripped by the time coalescence occurs).

    Scope: assumes the companion is an MS star (k2 in {0,1}) at
    `donor_age` -- HTP02's "effective core" treatment for an MS
    secondary (M'_c2 = companion_mass, R'_c2 = its actual radius,
    M'_env2 = 0, so it contributes no envelope-binding-energy term of
    its own, eq. 69). Does not check the companion's phase; callers
    must ensure this themselves (mirrors `apply_stable_mass_transfer`,
    which makes the same kind of assumption about its inputs).

    Named simplification for the merge case: unlike
    `apply_stable_mass_transfer` (which reaches an exact, self-
    consistent physical state), the merge branch here does NOT
    implement HTP02's eqs. 74-77 (the partial-envelope-retention
    Newton-Raphson solve for the merged star's final mass) -- that
    needs `R_i`, "the radius the system would have if it were to
    coalesce immediately," which HTP02 does not define operationally
    in a way this module could implement without further, separate
    study. Instead, the envelope is assumed fully lost at coalescence
    (a bare-core merger: `merge_stellar_masses(core_mass,
    companion_mass)`), the same "no partial retention" simplification
    already used for `merge_stellar_masses` itself in the
    IMMEDIATE_MERGER case. This is a real, documented gap, not an
    oversight -- see docs/science/rlof-ce-classifier-proposal.md.
    """
    core_mass = main_sequence.core_mass_hg(donor_mass, z, donor_age)
    core_radius_val = remnant.core_radius(core_mass, donor_mass, z)
    donor_radius = main_sequence.hg_radius(donor_mass, z, donor_age)
    companion_radius = main_sequence.ms_radius(companion_mass, z, donor_age)

    envelope_mass = donor_mass - core_mass
    # M'_env2 = 0 for an MS companion (its "effective core" is its
    # whole mass) -- eq. 69's second term vanishes; companion_radius
    # is not needed for E_bind,i itself, only for the coalescence
    # check below.
    e_bind_i_over_g = -(1.0 / lambda_ce) * (donor_mass * envelope_mass / donor_radius)

    e_orb_i_over_g = -0.5 * core_mass * companion_mass / separation
    e_orb_f_over_g = e_bind_i_over_g / alpha_ce + e_orb_i_over_g
    a_f = -0.5 * core_mass * companion_mass / e_orb_f_over_g

    # Coalescence check (HTP02, the paragraph introducing eq. 73):
    # would either the donor's bare core or the companion fill its
    # own Roche lobe before the orbit widens... narrows, rather, back
    # out to a_f? As the CE inspiral proceeds (separation decreasing
    # from the initial value toward a_f), the FIRST Roche-lobe-filling
    # condition reached -- i.e. the one at the LARGEST separation --
    # is what actually happens; a_f is only reached if it is smaller
    # than (i.e. the inspiral gets there before) that coalescence
    # separation.
    q_companion = companion_mass / core_mass
    q_core = core_mass / companion_mass
    a_l_companion = companion_radius / roche_lobe_radius(1.0, q_companion)
    a_l_core = core_radius_val / roche_lobe_radius(1.0, q_core)
    a_l = max(a_l_companion, a_l_core)

    if a_f > a_l:
        return True, core_mass, companion_mass, a_f
    return False, core_mass, companion_mass, None


def merge_stellar_masses(donor_mass: float, companion_mass: float) -> float:
    """Combined mass of an immediate/dynamical MS-MS merger.

    Conservative merger (no mass loss) -- HTP02 gives an explicit
    mass-loss prescription for common-envelope mergers (the envelope
    binding-energy balance, eqs. 69-77) but not for a direct/dynamical
    MS-MS collision merger, which is the only merger channel reachable
    by this module's current MS-only scope (see module docstring).
    Absent an explicit prescription, conservative merging is the
    natural default; this is a named simplification, not a citation.
    """
    return donor_mass + companion_mass
