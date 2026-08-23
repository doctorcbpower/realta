from __future__ import annotations

import logging

import numpy as np

from realta.config import SimulationConfig
from realta.imf.factory import get_imf
from realta.io.tables import IonizingPhotonTable, LifetimeTable, RemnantTable
from realta.xray.luminosity import XRayLuminosity

logger = logging.getLogger("realta")


class BinaryPopulation:
    """Vectorized binary population manager.

    The core Monte Carlo engine, based on the coeval globular-cluster
    HMXB population model of Power et al. (2009) -- see
    docs/provenance.md for the full paper-equation -> implementation ->
    test traceability table this class's methods are part of.
    """

    PFAC = 365.229126
    AFAC = 0.0193852859

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
        # Total mass formed across the WHOLE sampled IMF (mmin-mmax),
        # not just the M >= mcut binary progenitors -- set in
        # generate_population(). Needed to correctly normalize
        # population-luminosity comparisons such as MSLuminosityTable
        # against this specific run's ntot/mmin/mmax choice.
        self.total_mass_msun: float = 0.0

        self.imf = get_imf(config.imf_type)
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

        # Primordial binary fraction is 100% for stars above mcut
        # (Power et al. 2009, Sec. 2.1) -- every massive star gets a
        # companion and a period. `fsur` is NOT a formation-time
        # filter: it is applied once, later, as the HMXB activation
        # probability at primary supernova (see evolve()).
        log_pmin, log_pmax = np.log(cfg.pmin), np.log(cfg.pmax)
        periods = np.exp(
            log_pmin + (log_pmax - log_pmin) * self.np_rng.random(n_massive)
        )

        # Vectorized companion mass assignment
        if cfg.mcomp < 0:
            m2 = cfg.mmin + (m1 - cfg.mmin) * self.np_rng.random(n_massive)
        else:
            m2 = cfg.mcomp + (m1 - cfg.mcomp) * self.np_rng.random(n_massive)

        m2 = np.minimum(m2, m1)

        # Vectorized Semi-Major Axis calculation
        q = np.divide(m2, m1, out=np.zeros_like(m1), where=m1 > 0)
        a = self.AFAC * (m1 ** (1 / 3)) * ((1.0 + q) ** (1 / 3)) * (periods ** (2 / 3))

        # Vectorized lifetimes
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

        # Tracking arrays
        self.nturn = np.zeros(n_massive, dtype=np.int8)
        self.is_survived = np.ones(n_massive, dtype=bool)
        self.lum_xray = np.zeros(n_massive, dtype=np.float64)

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
                    if (
                        self.m2[i] > mcomp_abs
                        and self.np_rng.random() <= self.config.fsur
                    ):
                        self.lum_xray[i] = self.xray_calc.get_lumx(
                            self.m1[i],
                            self.m2[i],
                            self.period[i],
                            self.a[i],
                            rng=self.np_rng,
                        )
                else:
                    self.is_survived[i] = False

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
