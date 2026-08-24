from __future__ import annotations

import logging

import numpy as np
from scipy.stats import truncnorm

from realta.binaries.interaction import (
    RLOFOutcome,
    apply_common_envelope,
    apply_stable_mass_transfer,
    find_rlof_onset,
    merge_stellar_masses,
    rejuvenate_ms_gainer,
    roche_lobe_radius,
)
from realta.config import SimulationConfig
from realta.imf.factory import get_imf
from realta.io.tables import IonizingPhotonTable, LifetimeTable, RemnantTable
from realta.stellar import main_sequence
from realta.xray.luminosity import XRayLuminosity

logger = logging.getLogger("realta")

# Realta's imetal convention (config.py: "1=Z=0, 2=Z=0.008, 3=Z=0.02"),
# mapped to the numeric Z the Hurley/Tout stellar formulae need. imetal=1
# (Z=0) is deliberately absent -- those formulae are undefined at Z=0
# (they involve log(Z)); see the use_rlof_classifier handling below.
_IMETAL_TO_Z = {2: 0.008, 3: 0.02}


class BinaryPopulation:
    """Vectorized binary population manager.

    The core Monte Carlo engine, based on the coeval globular-cluster
    HMXB population model of Power et al. (2009) -- see
    docs/provenance.md for the full paper-equation -> implementation ->
    test traceability table this class's methods are part of.
    """

    PFAC = 365.229126
    AFAC = 0.0193852859

    # `self.a` (from AFAC/PFAC above) is in AU, not Rsun -- confirmed by
    # evaluating AFAC's own formula for an Earth-Sun-like case (M=1
    # Msun, P=365.25 days): it returns a~=0.99, matching 1 AU, not 1
    # Rsun (which would put the Earth inside the Sun). PFAC=365.229126
    # itself is essentially the number of days in a sidereal year --
    # i.e. these constants encode Kepler's third law in the standard
    # P(yr)^2 = a(AU)^3/M(Msun) astronomical convention, the same one
    # Power et al.'s original Fortran uses (see docs/provenance.md's
    # semi-major-axis row for the ~1% precision-convention note on
    # AFAC itself).
    #
    # binaries/interaction.py's RLOF/CE module (added this session) was
    # built and unit-tested entirely with hand-picked Rsun-scale
    # separations -- it compares `separation` directly against donor/
    # companion radii from the Hurley/Tout stellar-radius fits, which
    # are explicitly in Rsun. Passing `self.a` (AU) into it unconverted
    # made every donor look ~215x closer to its Roche lobe than it
    # really is -- a real bug found by running the full Paper 1
    # pipeline end-to-end (see docs/provenance.md's "known gaps"/RLOF
    # section): with the un-converted units, essentially every massive
    # binary classified as IMMEDIATE_MERGER, regardless of period.
    # RSUN_PER_AU converts at that boundary -- self.a itself STAYS in
    # AU everywhere else (the SN1 mass-loss orbit-widening code below,
    # and every existing pinned regression value, is untouched).
    #
    # 1 au = 1.495978707e11 m (IAU 2012 exact definition); R_sun =
    # 6.957e8 m (IAU 2015 nominal solar radius) -> ratio = 215.032.
    RSUN_PER_AU = 215.032

    # Converts total X-ray luminosity (in units of `lunit`) to an
    # ionising photon rate, using a model spectral shape (Power et al.
    # 2013): 6.2415e11 is erg->eV; 13.6 eV is the hydrogen ionisation
    # energy; the log ratio is a spectral-shape correction for photons
    # distributed between 13.6 eV and 1500 eV, referenced against an
    # upper bound of 1e6 eV.
    NPHOT_PER_LUMX = (6.2415e11 / 13.6) * np.log(1500.0 / 13.6) / np.log(1.0e6 / 13.6)

    def __init__(self, config: SimulationConfig):
        self.config = config
        self.m1: np.ndarray = np.array([])
        self.m2: np.ndarray = np.array([])
        self.period: np.ndarray = np.array([])
        self.a: np.ndarray = np.array([])
        self.turnoff_time: np.ndarray = np.array([])
        self.t2_lifetime: np.ndarray = np.array([])
        self.nturn: np.ndarray = np.array([])
        self.is_survived: np.ndarray = np.array([])
        self.lum_xray: np.ndarray = np.array([])
        # Paper 1 merger-channel bookkeeping (see
        # docs/science/paper1-binary-interaction-proposal.md) -- inert
        # (all False/nan) unless config.p_merge > 0. Minimal event
        # record kept now because a future Figure 3 will need it;
        # merges happen once, at formation, so merge_time is always 0.0
        # for merged systems.
        self.did_merge: np.ndarray = np.array([])
        self.merge_time: np.ndarray = np.array([])
        # MS Roche-lobe-overflow bookkeeping (see
        # docs/science/rlof-ce-classifier-proposal.md,
        # binaries/interaction.py) -- inert (rlof_time all inf) unless
        # config.use_rlof_classifier is True. Precomputed once per
        # binary in generate_population(), the same pattern as
        # turnoff_time, so evolve() only needs a cheap tnow comparison.
        self.rlof_time: np.ndarray = np.array([])
        self.rlof_outcome: np.ndarray = np.array([])
        self.rlof_donor_is_star1: np.ndarray = np.array([])
        self.rlof_processed: np.ndarray = np.array([])
        # Total mass formed across the WHOLE sampled IMF (mmin-mmax),
        # not just the M >= mcut binary progenitors -- set in
        # generate_population(). Needed to correctly normalize
        # population-luminosity comparisons such as MSLuminosityTable
        # against this specific run's ntot/mmin/mmax choice.
        self.total_mass_msun: float = 0.0

        self.imf = get_imf(config.imf_type, slope=config.imf_slope)
        self.np_rng = np.random.default_rng(config.iseed)

        self.lifetime_table = LifetimeTable(config.imetal, config.data_dir)
        self.remnant_table = RemnantTable(config.data_dir)
        self.ionizing_table = IonizingPhotonTable(config.data_dir)

        self.xray_calc = XRayLuminosity(
            lxmin=10.0 ** (config.lxmin - np.log10(config.lunit)),
            lxmax=10.0 ** (config.lxmax - np.log10(config.lunit)),
            lunit=config.lunit,
            distribution=config.xray_distribution,
        )

        self.generate_population()

    def generate_population(self):
        cfg = self.config
        logger.info(f"Generating {cfg.ntot} stars with IMF type {cfg.imf_type}")

        # Vectorized IMF mass generation
        raw_masses = self.imf.sample(cfg.ntot, cfg.mmin, cfg.mmax, self.np_rng)
        self.total_mass_msun = float(raw_masses.sum())

        # Filter primary mass cutoff upfront
        mask_massive = raw_masses >= cfg.mcut
        m1 = raw_masses[mask_massive]
        n_massive = len(m1)

        # Primordial binary fraction (A1, docs/science/paper1-detailed-
        # work-breakdown.md): binary_fraction=1.0 (default) reproduces
        # the pre-existing Power et al. (2009) baseline exactly (every
        # M >= mcut star paired, Sec. 2.1) -- the RNG draw below is
        # skipped entirely in that case, not drawn-and-discarded, so
        # the RNG stream is untouched by default (same pattern as
        # p_merge==0 below). `fsur` is NOT a formation-time filter: it
        # is applied once, later, as the HMXB activation probability at
        # primary supernova (see evolve()).
        #
        # "single" prescription (see
        # docs/science/paper1-binary-interaction-proposal.md): forces
        # no companion for any massive star (has_companion all False),
        # regardless of any separately-configured binary_fraction --
        # unlike every other prescription, m1 stays fully populated
        # (not emptied to n_massive=0, as this prescription used to do
        # before A3). Emptying the array was a reasonable shortcut
        # while the only things computed from m1 were L_X/HMXB-related
        # (correctly zero for single stars either way), but A3's
        # massive-star Q_H(t) also reads m1/nturn to know which
        # massive stars exist -- with an empty array it would (and, in
        # an intervening version of this code, briefly did) see none
        # at all and silently report Q_H=0 for single-star populations,
        # which is physically wrong: massive stars ionize regardless of
        # binarity. Migrating "single" onto the same has_companion
        # mechanism A1 already built fixes this for free -- L_X/HMXB
        # activation stay exactly zero (m2=0 blocks that unconditionally,
        # same as before), but L_bol/Q_H tracking now works correctly.
        if cfg.binary_prescription == "single":
            has_companion = np.zeros(n_massive, dtype=bool)
        elif cfg.binary_fraction >= 1.0:
            has_companion = np.ones(n_massive, dtype=bool)
        else:
            has_companion = self.np_rng.random(n_massive) < cfg.binary_fraction

        # Period distribution: "log_uniform" (default, unchanged
        # baseline) or "log_normal" (generic, pmin/pmax-derived
        # parameters -- see config.py::SimulationConfig.period_distribution
        # for the exact convention and why it's not literature-sourced).
        if cfg.period_distribution == "log_normal":
            log10_pmin, log10_pmax = np.log10(cfg.pmin), np.log10(cfg.pmax)
            mu = 0.5 * (log10_pmin + log10_pmax)
            sigma = (log10_pmax - log10_pmin) / 6.0
            a_trunc = (log10_pmin - mu) / sigma
            b_trunc = (log10_pmax - mu) / sigma
            log10_periods = truncnorm.rvs(
                a_trunc,
                b_trunc,
                loc=mu,
                scale=sigma,
                size=n_massive,
                random_state=self.np_rng,
            )
            periods = 10.0**log10_periods
        else:  # "log_uniform"
            log_pmin, log_pmax = np.log(cfg.pmin), np.log(cfg.pmax)
            periods = np.exp(
                log_pmin + (log_pmax - log_pmin) * self.np_rng.random(n_massive)
            )

        # Mass-ratio distribution: "uniform" (default, unchanged
        # baseline) or "flat_q" (flat in q = m2/m1 rather than in
        # absolute companion mass -- see
        # config.py::SimulationConfig.mass_ratio_distribution).
        if cfg.mass_ratio_distribution == "flat_q":
            m2 = m1 * self.np_rng.random(n_massive)
        else:  # "uniform"
            if cfg.mcomp < 0:
                m2 = cfg.mmin + (m1 - cfg.mmin) * self.np_rng.random(n_massive)
            else:
                m2 = cfg.mcomp + (m1 - cfg.mcomp) * self.np_rng.random(n_massive)
            m2 = np.minimum(m2, m1)

        # Stars without a companion (has_companion=False) have no real
        # orbit -- m2/period/a are placeholder zeros, not physical
        # values, and are excluded from the RLOF classifier and the
        # pre-SN merger channel below (both require a real companion).
        m2 = np.where(has_companion, m2, 0.0)
        periods = np.where(has_companion, periods, 0.0)

        # Pre-SN merger channel (enhanced_mergers prescription; inert
        # -- did_merge all False -- when config.p_merge == 0, which is
        # every other prescription's default). Merges happen once, here,
        # at formation: a merged system's companion is folded into the
        # primary (rejuvenating it -- its lifetime is recomputed below
        # from the merged mass) and permanently zeroed out, which is
        # sufficient on its own to keep it out of the HMXB channel in
        # evolve() (activation there requires m2 > |mcomp|). See the
        # proposal doc for why this is NOT paper-derived physics.
        # p_merge == 0 is every prescription's default except
        # enhanced_mergers -- skip the RNG draw entirely in that case
        # rather than drawing-and-discarding, so the RNG stream (and
        # therefore every pinned regression value) is untouched when
        # the merger channel isn't in use. `has_companion` additionally
        # excludes no-companion stars (period=0 would otherwise trivially
        # satisfy `periods < p_merge_max_period`).
        if cfg.p_merge > 0.0:
            did_merge = (
                has_companion
                & (periods < cfg.p_merge_max_period)
                & (self.np_rng.random(n_massive) < cfg.p_merge)
            )
        else:
            did_merge = np.zeros(n_massive, dtype=bool)
        m1 = np.where(did_merge, m1 + cfg.f_merge * m2, m1)
        m2 = np.where(did_merge, 0.0, m2)

        # Vectorized Semi-Major Axis calculation
        q = np.divide(m2, m1, out=np.zeros_like(m1), where=m1 > 0)
        a = self.AFAC * (m1 ** (1 / 3)) * ((1.0 + q) ** (1 / 3)) * (periods ** (2 / 3))

        # Vectorized lifetimes (merged primaries use their post-merger
        # mass, i.e. they are rejuvenated rather than retaining the
        # unmerged primary's lifetime)
        t_off1 = np.array([self.lifetime_table.get_lifetime(m) for m in m1])
        t_off2 = np.array(
            [self.lifetime_table.get_lifetime(m) if m > 0 else 0.0 for m in m2]
        )

        # Pre-sort population ONCE by primary turnoff time
        sort_idx = np.argsort(t_off1)

        self.m1 = m1[sort_idx]
        self.m2 = m2[sort_idx]
        self.period = periods[sort_idx]
        self.a = a[sort_idx]
        self.turnoff_time = t_off1[sort_idx]
        self.t2_lifetime = t_off2[sort_idx]
        self.did_merge = did_merge[sort_idx]
        self.merge_time = np.where(self.did_merge, 0.0, np.nan)

        # Tracking arrays
        self.nturn = np.zeros(n_massive, dtype=np.int8)
        self.is_survived = np.ones(n_massive, dtype=bool)
        self.lum_xray = np.zeros(n_massive, dtype=np.float64)

        # MS Roche-lobe-overflow precomputation (opt-in,
        # config.use_rlof_classifier -- see
        # docs/science/rlof-ce-classifier-proposal.md). Inert
        # (rlof_time all inf, never processed) when disabled, which is
        # every existing config's default -- this branch does not run
        # at all in that case, so it cannot perturb the pinned baseline.
        self.rlof_time = np.full(n_massive, np.inf)
        # NOT np.full(n_massive, RLOFOutcome.DETACHED, dtype=object) --
        # np.full() silently truncates/corrupts a str-Enum fill value
        # even with dtype=object (confirmed: the array ends up holding
        # a plain, truncated str that fails equality against the real
        # enum member). Per-element assignment (done below, in the
        # precompute loop) is unaffected -- only the bulk fill is
        # broken. List-construction avoids it entirely.
        self.rlof_outcome = np.array([RLOFOutcome.DETACHED] * n_massive, dtype=object)
        self.rlof_donor_is_star1 = np.ones(n_massive, dtype=bool)
        self.rlof_processed = np.zeros(n_massive, dtype=bool)

        if cfg.use_rlof_classifier:
            z = _IMETAL_TO_Z.get(cfg.imetal)
            if z is None:
                logger.warning(
                    "use_rlof_classifier=True but imetal=%s (Z=0): the "
                    "Hurley/Tout stellar formulae are undefined at Z=0 "
                    "-- skipping RLOF classification for this run.",
                    cfg.imetal,
                )
            else:
                for i in range(n_massive):
                    if self.did_merge[i] or self.m2[i] <= 0.0:
                        # Already merged at formation, or never had a
                        # companion at all (A1's binary_fraction<1
                        # placeholder m2=0 -- find_rlof_onset/
                        # classify_rlof divide by companion_mass, which
                        # would be a divide-by-zero/inf-q1 corruption
                        # here, not a real RLOF-eligible system).
                        continue
                    t_rlof, outcome, donor_is_star1 = find_rlof_onset(
                        self.m1[i],
                        self.m2[i],
                        self.a[i] * self.RSUN_PER_AU,
                        z,
                        q_crit_ms=cfg.q_crit_ms,
                    )
                    self.rlof_time[i] = t_rlof
                    self.rlof_outcome[i] = outcome
                    self.rlof_donor_is_star1[i] = donor_is_star1

        logger.info(f"Initialized {n_massive} massive binaries in vectorized memory.")

    def evolve(self, tnow: float, dt: float) -> tuple[float, float, int, int]:
        """Advance the population by one timestep.

        Returns (lumx_tot, nphot_tot, nactive, ndead):
            lumx_tot: summed active-HMXB X-ray luminosity, in erg/s.
            nphot_tot: effective ionising photon rate from HMXBs, in s^-1
                (see the Phase 3 comment below for its provenance).
            nactive: number of primary supernovae during *this* timestep
                (a formation-rate count, not a running census -- see the
                Phase 3 comment below).
            ndead: cumulative number of binaries with both stars now
                compact remnants.
        """
        mcomp_abs = abs(self.config.mcomp)

        # --- Phase 0: MS Roche-lobe overflow (opt-in, use_rlof_classifier) ---
        # Inert when disabled (rlof_time all inf -> mask never matches).
        # See docs/science/rlof-ce-classifier-proposal.md. Gated on
        # nturn == 0 (neither star has had its primary SN yet, per
        # Realta's own -- independently sourced -- LifetimeTable): a
        # predicted RLOF event whose time falls after this binary's own
        # supernova already occurred is naturally suppressed rather
        # than retroactively applied, since the two stellar-lifetime
        # prescriptions are not reconciled (see
        # binaries/interaction.py::find_rlof_onset's docstring).
        if self.config.use_rlof_classifier:
            rlof_mask = (
                (~self.rlof_processed) & (self.nturn == 0) & (tnow >= self.rlof_time)
            )
            if np.any(rlof_mask):
                idx0 = np.where(rlof_mask)[0]
                z = _IMETAL_TO_Z.get(self.config.imetal)
                for i in idx0:
                    self.rlof_processed[i] = True
                    if self.rlof_outcome[i] == RLOFOutcome.IMMEDIATE_MERGER:
                        # Same treatment as the existing formation-time
                        # merger channel: fold the companion in, zero
                        # it out (this alone keeps the system out of
                        # the HMXB channel below, since activation
                        # requires m2 > |mcomp|), and reset the
                        # lifetime clock for the merged mass starting
                        # now -- a full-reset simplification, not
                        # partial (Tout et al. 1997/Brček et al.
                        # 2026) rejuvenation; see the proposal doc.
                        merged_mass = merge_stellar_masses(self.m1[i], self.m2[i])
                        self.m1[i] = merged_mass
                        self.m2[i] = 0.0
                        self.did_merge[i] = True
                        self.merge_time[i] = tnow
                        self.turnoff_time[i] = tnow + self.lifetime_table.get_lifetime(
                            merged_mass
                        )
                        self.t2_lifetime[i] = 0.0
                    elif self.rlof_outcome[i] == RLOFOutcome.STABLE_MASS_TRANSFER:
                        # Instantaneous conservative transfer to the new
                        # detachment point -- see
                        # binaries/interaction.py::apply_stable_mass_transfer
                        # and the proposal doc's "Decision" on
                        # instantaneous vs. rate-integrated treatment.
                        if self.rlof_donor_is_star1[i]:
                            donor_mass, companion_mass = self.m1[i], self.m2[i]
                        else:
                            donor_mass, companion_mass = self.m2[i], self.m1[i]
                        # find_rlof_onset() can find a stable-MT
                        # crossing during either the MS or the HG
                        # (added when HG search was wired in) -- use
                        # whichever radius function actually applies
                        # at the donor's phase at rlof_time, not
                        # unconditionally ms_radius.
                        donor_phase = main_sequence.phase(
                            donor_mass, z, self.rlof_time[i]
                        )
                        if donor_phase in (0, 1):
                            donor_radius = main_sequence.ms_radius(
                                donor_mass, z, self.rlof_time[i]
                            )
                        else:
                            donor_radius = main_sequence.hg_radius(
                                donor_mass, z, self.rlof_time[i]
                            )
                        new_donor, new_companion, new_a_rsun = (
                            apply_stable_mass_transfer(
                                donor_mass,
                                companion_mass,
                                self.a[i] * self.RSUN_PER_AU,
                                donor_radius,
                            )
                        )
                        if self.rlof_donor_is_star1[i]:
                            self.m1[i], self.m2[i] = new_donor, new_companion
                        else:
                            self.m2[i], self.m1[i] = new_donor, new_companion
                        self.a[i] = new_a_rsun / self.RSUN_PER_AU
                        mtot = new_donor + new_companion
                        self.period[i] = self.PFAC * np.sqrt((self.a[i] ** 3) / mtot)
                        # Donor's lifetime clock: full-reset
                        # simplification (same as mergers) -- no
                        # verified response/rejuvenation prescription
                        # for a mass-LOSING MS/HG donor is applied here.
                        # Companion's (gainer's) lifetime clock:
                        # rejuvenated via Tout et al. (1997) eq. 41 (B3,
                        # docs/science/paper1-detailed-work-breakdown.md)
                        # instead of a full reset -- see
                        # binaries/interaction.py::rejuvenate_ms_gainer
                        # -- but only when the companion is genuinely
                        # MS-phase at tnow (that formula's own scope);
                        # falls back to the full-reset simplification
                        # otherwise (an HG+ companion, or one past this
                        # module's t_BGB scope -- classify_rlof() places
                        # no phase constraint on the companion for
                        # STABLE_MASS_TRANSFER, so this is a real,
                        # if presumably uncommon, case to guard).
                        if self.rlof_donor_is_star1[i]:
                            donor_new_mass, companion_new_mass = (
                                self.m1[i],
                                self.m2[i],
                            )
                        else:
                            donor_new_mass, companion_new_mass = (
                                self.m2[i],
                                self.m1[i],
                            )
                        try:
                            companion_is_ms = main_sequence.phase(
                                companion_mass, z, tnow
                            ) in (0, 1)
                        except ValueError:
                            companion_is_ms = False
                        if companion_is_ms:
                            remaining_fraction = rejuvenate_ms_gainer(
                                companion_mass, companion_new_mass, tnow, z
                            )
                            companion_new_time = (
                                tnow
                                + remaining_fraction
                                * self.lifetime_table.get_lifetime(companion_new_mass)
                            )
                        else:
                            companion_new_time = (
                                tnow
                                + self.lifetime_table.get_lifetime(companion_new_mass)
                            )
                        donor_new_time = tnow + self.lifetime_table.get_lifetime(
                            donor_new_mass
                        )
                        if self.rlof_donor_is_star1[i]:
                            self.turnoff_time[i] = donor_new_time
                            self.t2_lifetime[i] = companion_new_time
                        else:
                            self.turnoff_time[i] = companion_new_time
                            self.t2_lifetime[i] = donor_new_time
                    elif self.rlof_outcome[i] == RLOFOutcome.COMMON_ENVELOPE:
                        # HG donors dynamically unstable at RLOF (see
                        # binaries/interaction.py::hg_q_crit) resolve
                        # via the alpha-lambda energy-balance solve --
                        # see apply_common_envelope's docstring.
                        if self.rlof_donor_is_star1[i]:
                            donor_mass, companion_mass = self.m1[i], self.m2[i]
                        else:
                            donor_mass, companion_mass = self.m2[i], self.m1[i]
                        survives, new_donor, new_companion, new_a_rsun = (
                            apply_common_envelope(
                                donor_mass,
                                companion_mass,
                                self.a[i] * self.RSUN_PER_AU,
                                z,
                                self.rlof_time[i],
                                alpha_ce=self.config.alpha_ce,
                                lambda_ce=self.config.lambda_ce,
                            )
                        )
                        if survives:
                            # Donor stripped to its bare core; companion
                            # unaffected; orbit tightened to a_f -- same
                            # full-reset lifetime-clock simplification
                            # used for stable mass transfer above.
                            if self.rlof_donor_is_star1[i]:
                                self.m1[i], self.m2[i] = new_donor, new_companion
                            else:
                                self.m2[i], self.m1[i] = new_donor, new_companion
                            self.a[i] = new_a_rsun / self.RSUN_PER_AU
                            mtot = new_donor + new_companion
                            self.period[i] = self.PFAC * np.sqrt(
                                (self.a[i] ** 3) / mtot
                            )
                            # Donor stripped to its bare core -- full-
                            # reset lifetime-clock simplification
                            # (matching STABLE_MASS_TRANSFER's donor
                            # treatment above). The companion is
                            # genuinely mass-unaffected by a surviving
                            # CE (apply_common_envelope's own
                            # docstring: "the companion is unaffected")
                            # -- its lifetime clock is left entirely
                            # untouched here, not reset. A reset would
                            # incorrectly de-age it (pretend it just
                            # formed anew at tnow) despite nothing
                            # having physically happened to it -- found
                            # while implementing B3's rejuvenation fix
                            # for STABLE_MASS_TRANSFER (a related but
                            # distinct issue: that companion genuinely
                            # gains mass, this one doesn't change at
                            # all), see docs/provenance.md Section 12.
                            if self.rlof_donor_is_star1[i]:
                                self.turnoff_time[i] = (
                                    tnow + self.lifetime_table.get_lifetime(self.m1[i])
                                )
                            else:
                                self.t2_lifetime[i] = (
                                    tnow + self.lifetime_table.get_lifetime(self.m2[i])
                                )
                        else:
                            # Cores coalesce before the envelope is
                            # fully ejected -- merge the donor's core
                            # (new_donor, already reduced from its
                            # pre-CE mass) with the companion, same
                            # treatment as IMMEDIATE_MERGER above.
                            merged_mass = merge_stellar_masses(new_donor, new_companion)
                            self.m1[i] = merged_mass
                            self.m2[i] = 0.0
                            self.did_merge[i] = True
                            self.merge_time[i] = tnow
                            self.turnoff_time[i] = (
                                tnow + self.lifetime_table.get_lifetime(merged_mass)
                            )
                            self.t2_lifetime[i] = 0.0

        # --- Phase 1: Primary Supernova Transitions ---
        # Triggers when system is on MS (nturn == 0) and tnow exceeds primary lifetime
        sn1_mask = (
            (self.nturn == 0) & (tnow >= self.turnoff_time) & (self.turnoff_time > 0.0)
        )
        if np.any(sn1_mask):
            idx1 = np.where(sn1_mask)[0]
            for i in idx1:
                remnant1 = self.remnant_table.get_remnant_mass(self.m1[i])
                deltam = self.m1[i] - remnant1
                floss = (
                    deltam / (self.m1[i] + self.m2[i])
                    if (self.m1[i] + self.m2[i]) > 0
                    else 1.0
                )

                self.m1[i] = remnant1
                # Transition to HMXB phase: next turnoff is when secondary completes MS
                self.turnoff_time[i] = self.t2_lifetime[i]
                self.nturn[i] = 1

                if floss <= 0.5:
                    # Deterministic sudden-mass-loss survival criterion
                    # (Power et al. 2009, Sec. 2.1) -- no additional
                    # stochastic term in this criterion itself.
                    mtot = self.m1[i] + self.m2[i]
                    if mtot > 0:
                        self.a[i] *= deltam / mtot
                        self.period[i] = self.PFAC * np.sqrt((self.a[i] ** 3) / mtot)
                    self.is_survived[i] = True

                    # HMXB activation gate: f_sur, the probability that a
                    # surviving, sufficiently massive binary is observed
                    # as an active HMXB (Power et al. 2009, Sec. 2.1).
                    # The X-ray luminosity is drawn exactly once, here,
                    # and held fixed for the rest of the binary's active
                    # HMXB lifetime -- it is NOT redrawn every timestep.
                    #
                    # interaction_boost multiplies fsur, but -- since
                    # the reconciliation with the physics-based RLOF
                    # classifier (docs/science/rlof-ce-classifier-
                    # proposal.md "Decision 3") -- ONLY for binaries
                    # that the classifier actually found underwent
                    # stable mass transfer on the MS (a genuine
                    # interaction event), not unconditionally for every
                    # surviving binary as in the original placeholder
                    # version. interaction_boost is 1.0 (no effect) for
                    # every prescription except standard_interaction/
                    # enhanced_interaction, and rlof_outcome is always
                    # DETACHED when use_rlof_classifier=False (every
                    # other prescription's default), so this does not
                    # perturb the pre-existing baseline -- see
                    # docs/provenance.md Section 6 for the old-vs-new
                    # pinned values for the three affected prescriptions.
                    had_stable_mt = (
                        self.config.use_rlof_classifier
                        and self.rlof_processed[i]
                        and self.rlof_outcome[i] == RLOFOutcome.STABLE_MASS_TRANSFER
                    )
                    boost = self.config.interaction_boost if had_stable_mt else 1.0
                    fsur_eff = min(1.0, self.config.fsur * boost)
                    if self.m2[i] > mcomp_abs and self.np_rng.random() <= fsur_eff:
                        self.lum_xray[i] = self.xray_calc.get_lumx(
                            self.m1[i],
                            self.m2[i],
                            self.period[i],
                            self.a[i],
                            rng=self.np_rng,
                        )
                else:
                    self.is_survived[i] = False

        # --- Phase 1.5: Post-SN secondary Roche-lobe overflow (opt-in,
        # config.use_post_sn_rlof) ---
        # docs/science/paper1-followup-prompt.md. Inert when disabled
        # (mask never matches, since it's gated on the config flag
        # itself), so it cannot perturb the pre-existing baseline.
        # Checked every timestep (not precomputed like Phase 0's
        # rlof_time) because the secondary's radius grows continuously
        # and the trigger is a live comparison against its own Roche
        # lobe -- matching the existing SN1/SN2 phases' own live-check
        # style, not Phase 0's precompute-once pattern (which needs a
        # root-finder because Phase 0 must know the crossing time in
        # advance to place it correctly relative to nturn transitions;
        # here nturn==1 is already the gate, so no such ordering
        # problem exists).
        if self.config.use_post_sn_rlof:
            z2 = _IMETAL_TO_Z.get(self.config.imetal)
            if z2 is None:
                logger.warning(
                    "use_post_sn_rlof=True but imetal=%s (Z=0): the "
                    "Hurley/Tout stellar formulae are undefined at Z=0 "
                    "-- skipping post-SN secondary RLOF for this run.",
                    self.config.imetal,
                )
            else:
                # Not yet active (lum_xray==0) doubles as "hasn't
                # already triggered via this channel or via the
                # stochastic fsur draw above" -- once either sets
                # lum_xray>0, this mask naturally stops matching.
                psr_mask = (self.nturn == 1) & self.is_survived & (self.lum_xray == 0.0)
                if np.any(psr_mask):
                    idx_psr = np.where(psr_mask)[0]
                    for i in idx_psr:
                        if self.m2[i] <= mcomp_abs:
                            continue
                        try:
                            donor_phase = main_sequence.phase(self.m2[i], z2, tnow)
                            if donor_phase in (0, 1):
                                donor_radius = main_sequence.ms_radius(
                                    self.m2[i], z2, tnow
                                )
                            else:
                                donor_radius = main_sequence.hg_radius(
                                    self.m2[i], z2, tnow
                                )
                        except ValueError:
                            # Past this module's t_BGB scope, or mass
                            # outside its supported range -- not
                            # modelled, skip (same as classify_rlof's
                            # PHASE_NOT_MODELLED handling).
                            continue
                        q2 = self.m2[i] / self.m1[i]
                        r_l2 = roche_lobe_radius(self.a[i] * self.RSUN_PER_AU, q2)
                        if donor_radius >= r_l2:
                            self.lum_xray[i] = self.xray_calc.get_lumx(
                                self.m1[i],
                                self.m2[i],
                                self.period[i],
                                self.a[i],
                                rng=self.np_rng,
                            )

        # --- Phase 2: Secondary Supernova Transitions ---
        # Triggers ONLY when secondary star completes lifetime (tnow >= t2_lifetime)
        sn2_mask = (self.nturn == 1) & (tnow >= self.turnoff_time)
        if np.any(sn2_mask):
            idx2 = np.where(sn2_mask)[0]
            for i in idx2:
                self.m2[i] = self.remnant_table.get_remnant_mass(self.m2[i])
                self.turnoff_time[i] = 0.0
                self.nturn[i] = 2
                self.is_survived[i] = False
                self.lum_xray[i] = 0.0

        # --- Phase 3: Aggregate observables ---
        # self.lum_xray already holds each active HMXB's persistent
        # luminosity (drawn once, at activation, in Phase 1 above) --
        # sum as-is rather than redrawing every timestep.
        #
        # self.lum_xray is stored internally normalized by `lunit`
        # (matching xray_calc's internal lxmin/lxmax and Eddington-limit
        # bookkeeping). evolve()'s returned lumx_tot is scaled back to
        # actual erg/s here so it means what its name says.
        lumx_tot = float(np.sum(self.lum_xray)) * self.config.lunit

        # --- Counts ---
        # `nactive` is reset to zero every timestep and counts only the
        # primary supernovae that occur *during this step*, unconditional
        # on survival/activation -- it is a formation-rate column, not a
        # running census of currently-active HMXBs. sn1_mask, computed
        # above, is exactly that set of events.
        nactive = int(np.count_nonzero(sn1_mask))
        ndead = int(np.count_nonzero(self.nturn == 2))

        # --- Ionising photon rate ---
        # An empirical conversion of the *X-ray* luminosity itself into a
        # photoionising rate, assuming a spectrum between 13.6 eV and
        # 1500 eV referenced against an upper bound of 1e6 eV (Power et
        # al. 2013 -- see NPHOT_PER_LUMX above). The ionizing_table
        # machinery (io/tables.py::IonizingPhotonTable) is retained but
        # is not used to compute nphot_tot -- it estimates a different
        # quantity (the main-sequence ionising photon budget), which is
        # not currently wired into this evolution loop.
        nphot_tot = lumx_tot * self.NPHOT_PER_LUMX

        return lumx_tot, nphot_tot, nactive, ndead
