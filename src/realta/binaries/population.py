from __future__ import annotations

import logging

import numpy as np

from realta.config import SimulationConfig
from realta.imf.factory import get_imf
from realta.io.tables import IonizingPhotonTable, LifetimeTable, RemnantTable
from realta.xray.luminosity import XRayLuminosity

logger = logging.getLogger("realta")


class BinaryPopulation:
    """Vectorized binary population manager."""

    PFAC = 365.229126
    AFAC = 0.0193852859

    def __init__(self, config: SimulationConfig):
        self.config = config
        self.imf = get_imf(config.imf_type)
        self.np_rng = np.random.default_rng(config.iseed)

        self.lifetime_table = LifetimeTable(config.imetal, config.data_dir)
        self.remnant_table = RemnantTable(config.data_dir)
        self.ionizing_table = IonizingPhotonTable(config.data_dir)

        self.xray_calc = XRayLuminosity(
            lxmin=10.0 ** (config.lxmin - np.log10(config.lunit)),
            lxmax=10.0 ** (config.lxmax - np.log10(config.lunit)),
            lunit=config.lunit,
        )

        self.generate_population()

    def generate_population(self):
        cfg = self.config
        logger.info(f"Generating {cfg.ntot} stars with IMF type {cfg.imf_type}")

        # Vectorized IMF mass generation
        raw_masses = self.imf.sample(cfg.ntot, cfg.mmin, cfg.mmax, self.np_rng)

        # Filter primary mass cutoff upfront
        mask_massive = raw_masses >= cfg.mcut
        m1 = raw_masses[mask_massive]
        n_massive = len(m1)

        # Vectorized binary fraction & period sampling
        is_binary = self.np_rng.random(n_massive) <= cfg.fbin
        log_pmin, log_pmax = np.log(cfg.pmin), np.log(cfg.pmax)
        periods = (
            np.exp(log_pmin + (log_pmax - log_pmin) * self.np_rng.random(n_massive))
            * is_binary
        )

        # Vectorized companion mass assignment
        if cfg.mcomp < 0:
            m2 = cfg.mmin + (m1 - cfg.mmin) * self.np_rng.random(n_massive)
        else:
            m2 = cfg.mcomp + (m1 - cfg.mcomp) * self.np_rng.random(n_massive)

        m2 = np.minimum(m2, m1) * is_binary

        # Vectorized Semi-Major Axis calculation
        q = np.divide(m2, m1, out=np.zeros_like(m1), where=m1 > 0)
        a = np.where(
            is_binary,
            self.AFAC * (m1 ** (1 / 3)) * ((1.0 + q) ** (1 / 3)) * (periods ** (2 / 3)),
            0.0,
        )

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
        fsur_val = getattr(self.config, "fsur", getattr(self.config, "f_sur", 1.0))
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
                    mtot = self.m1[i] + self.m2[i]
                    if mtot > 0:
                        self.a[i] *= deltam / mtot
                        self.period[i] = self.PFAC * np.sqrt((self.a[i] ** 3) / mtot)
                    self.is_survived[i] = self.np_rng.random() <= fsur_val
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

        # --- Phase 3: Evaluate Active HMXBs ---
        # Systems currently in HMXB phase (nturn == 1) that have NOT reached secondary turnoff
        active_hmxb_mask = (self.nturn == 1) & self.is_survived & (self.m2 >= mcomp_abs)

        self.lum_xray.fill(0.0)
        if np.any(active_hmxb_mask):
            hmxb_idx = np.where(active_hmxb_mask)[0]
            for i in hmxb_idx:
                self.lum_xray[i] = self.xray_calc.get_lumx(
                    self.m1[i],
                    self.m2[i],
                    self.period[i],
                    self.a[i],
                    iseed=None,
                    use_weibull=True,
                )

        lumx_tot = float(np.sum(self.lum_xray))

        # --- Counts & Photon Totals ---
        nactive = int(np.count_nonzero(self.nturn == 1))
        ndead = int(np.count_nonzero(self.nturn == 2))

        nphot_tot = 0.0
        ms_mask = self.nturn == 0
        if np.any(ms_mask):
            for i in np.where(ms_mask)[0]:
                ng1 = self.ionizing_table.get_ngamma(self.m1[i])
                ng2 = self.ionizing_table.get_ngamma(self.m2[i])
                t1, t2 = self.turnoff_time[i], self.t2_lifetime[i]
                if t1 > 0 and t2 > 0:
                    nphot_tot += 10.0 ** (
                        ng1 + np.log10(dt / t1) + ng2 + np.log10(dt / t2) - 60
                    )

        return lumx_tot, nphot_tot, nactive, ndead
