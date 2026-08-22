"""
HMXRB Simulation - Python Refactoring
=====================================

Complete refactoring of the FORTRAN code from Power et al. 2009 for simulating
High-Mass X-ray Binaries (HMXRBs) in globular clusters.

Original FORTRAN files analyzed:
- main.f: Primary simulation driver
- make_stars.f: Stellar population generation
- get_lumx.f: X-ray luminosity calculations
- lifetime.f: Stellar lifetime estimates
- get_mremnant.f: Remnant mass estimates
- get_ngamma.f: Ionizing photon calculations
- salpeter.f, kroupa.f, log_normal_IMF.f: Initial Mass Functions
- choose_metallicity.f, read_mremnant.f, read_ionise.f: Data readers
- locate.f, indexx.f, ran3.f: Utilities

Author: Chris Power
Date: 2026
"""

import numpy as np
import yaml
import os
import logging
import importlib.resource
from typing import Tuple, List, Dict, Optional, Union
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('hmxrb')


# =============================================================================
# CONFIGURATION
# =============================================================================

@dataclass
class SimulationConfig:
    """Configuration for the HMXRB simulation."""
    # Simulation parameters
    ntot: int = 100000
    mmin: float = 0.01
    mmax: float = 100.0
    mcut: float = 8.0
    tmax: float = 100.0
    dt: float = 0.01
    
    # IMF type: 1=Salpeter, 2=Kroupa, 3=Chabrier
    imf_type: int = 2
    
    # Binary parameters
    pmin: float = 0.1
    pmax: float = 1000.0
    mcomp: float = 0.5
    fbin: float = 0.5
    
    # Metallicity: 1=Z=0, 2=Z=0.008, 3=Z=0.02
    imetal: int = 2
    
    # X-ray luminosity
    lxmin: float = 33.0
    lxmax: float = 39.0
    lunit: float = 1.0e33
    
    # Random seed
    iseed: int = -12345
    
    # Data directory
    data_dir: str = 'data'


def load_config(config_path: Optional[str] = None) -> SimulationConfig:
    """Load configuration from YAML file or use defaults."""
    if config_path and os.path.exists(config_path):
        with open(config_path, 'r') as f:
            config_dict = yaml.safe_load(f)
        
        # Convert to SimulationConfig
        config = SimulationConfig()
        for key, value in config_dict.items():
            if hasattr(config, key):
                setattr(config, key, value)
        return config
    return SimulationConfig()


# =============================================================================
# RANDOM NUMBER GENERATOR
# =============================================================================

class RandomGenerator:
    """
    Random number generator replacing FORTRAN's ran3.
    Uses NumPy's random number generator for better performance and reproducibility.
    """
    
    def __init__(self, seed: int):
        self.rng = np.random.default_rng(seed)
    
    def random(self) -> float:
        """Generate a random number between 0 and 1."""
        return self.rng.random()
    
    def uniform(self, low: float = 0.0, high: float = 1.0) -> float:
        """Generate a uniform random number in [low, high)."""
        return self.rng.uniform(low, high)
    
    def exponential(self, scale: float = 1.0) -> float:
        """Generate an exponential random number."""
        return self.rng.exponential(scale)


# =============================================================================
# INITIAL MASS FUNCTIONS
# =============================================================================

class IMF(ABC):
    """Abstract base class for Initial Mass Functions."""
    
    @abstractmethod
    def cdf(self, m: float, mmin: float, mmax: float) -> float:
        """Cumulative distribution function P(<m)."""
        pass
    
    def sample(self, n: int, mmin: float, mmax: float, rng: RandomGenerator) -> np.ndarray:
        """
        Sample n masses from the IMF using inverse transform sampling.
        
        This replaces the FORTRAN approach of tabulating the CDF and using locate().
        """
        # Generate uniform random numbers
        u = np.array([rng.random() for _ in range(n)])
        
        # Use numerical inversion of CDF
        # For efficiency, we use a simple approach: binary search on mass
        masses = []
        for ui in u:
            # Binary search for mass where cdf(mass) = ui
            low, high = mmin, mmax
            for _ in range(100):  # 100 iterations should be enough
                mid = (low + high) / 2
                cdf_mid = self.cdf(mid, mmin, mmax)
                if cdf_mid < ui:
                    low = mid
                else:
                    high = mid
                if abs(high - low) < 1e-10:
                    break
            masses.append((low + high) / 2)
        
        return np.array(masses)


class SalpeterIMF(IMF):
    """
    Salpeter IMF: dN/dm ∝ m^(-alpha)
    
    Original FORTRAN: salpeter.f
    Alpha = 2.35 (standard Salpeter)
    """
    
    def __init__(self, alpha: float = 2.35):
        self.alpha = alpha
        self.beta = 1.0 - alpha  # = -1.35 for alpha=2.35
    
    def cdf(self, m: float, mmin: float, mmax: float) -> float:
        """
        Cumulative distribution function for Salpeter IMF.
        
        P(<m) = (m^beta - mmin^beta) / (mmax^beta - mmin^beta)
        where beta = 1 - alpha
        """
        if m <= mmin:
            return 0.0
        if m >= mmax:
            return 1.0
        
        numerator = m**self.beta - mmin**self.beta
        denominator = mmax**self.beta - mmin**self.beta
        
        return numerator / denominator


class KroupaIMF(IMF):
    """
    Kroupa piecewise IMF.
    
    Original FORTRAN: kroupa.f
    
    The Kroupa IMF has 4 power-law segments:
    - m < 0.08: beta = 0.3 (shallow)
    - 0.08 <= m < 0.5: beta = 1.3
    - 0.5 <= m < 1.0: beta = 2.3
    - m >= 1.0: beta = 2.3
    """
    
    def __init__(self):
        self.m0 = 0.01
        self.m1 = 0.08
        self.m2 = 0.5
        self.m3 = 1.0
        
        self.beta0 = 0.3
        self.beta1 = 1.3
        self.beta2 = 2.3
        self.beta3 = 2.3
    
    def _integrated_density(self, m: float, mmin: float, mmax: float, beta: float) -> float:
        """Integrate m^(-beta) from mmin to m."""
        if beta == 1.0:
            return np.log(m) - np.log(mmin)
        else:
            return (m**(1.0 - beta) - mmin**(1.0 - beta)) / (1.0 - beta)
    
    def cdf(self, m: float, mmin: float, mmax: float) -> float:
        """
        Cumulative distribution for Kroupa IMF.
        
        The CDF is normalized across all segments.
        """
        if m <= mmin:
            return 0.0
        if m >= mmax:
            return 1.0
        
        # Calculate the total normalization
        # dm0 = integral from m0 to m1
        # dm1 = integral from m1 to m2
        # dm2 = integral from m2 to m3
        # dm3 = integral from m3 to mmax
        
        # But we need to respect the user's mmin and mmax
        # Adjust the break points to be within [mmin, mmax]
        breaks = [mmin, self.m1, self.m2, self.m3, mmax]
        betas = [self.beta0, self.beta1, self.beta2, self.beta3, self.beta3]
        
        # Filter breaks that are outside [mmin, mmax]
        valid_breaks = [b for b in breaks if mmin <= b <= mmax]
        if valid_breaks[0] != mmin:
            valid_breaks.insert(0, mmin)
        if valid_breaks[-1] != mmax:
            valid_breaks.append(mmax)
        
        # Recalculate betas for valid segments
        valid_betas = []
        for i in range(len(valid_breaks) - 1):
            # Find which original segment this corresponds to
            for j in range(len(breaks) - 1):
                if breaks[j] <= valid_breaks[i] < breaks[j+1]:
                    valid_betas.append(betas[j])
                    break
            else:
                valid_betas.append(betas[-1])
        
        # Calculate normalization
        dmtot = 0.0
        for i in range(len(valid_breaks) - 1):
            dmtot += self._integrated_density(
                valid_breaks[i+1], valid_breaks[i], mmax, valid_betas[i]
            )
        
        # Calculate CDF at m
        cdf_val = 0.0
        for i in range(len(valid_breaks) - 1):
            if m <= valid_breaks[i+1]:
                cdf_val += self._integrated_density(
                    min(m, valid_breaks[i+1]), valid_breaks[i], mmax, valid_betas[i]
                )
                break
            else:
                cdf_val += self._integrated_density(
                    valid_breaks[i+1], valid_breaks[i], mmax, valid_betas[i]
                )
        
        return cdf_val / dmtot


class ChabrierIMF(IMF):
    """
    Chabrier IMF: Log-normal for m <= M_tr, power-law for m > M_tr.
    
    Original FORTRAN: log_normal_IMF.f
    
    For primordial IMF (as in the original code):
    - Mc = 3.5 Msun
    - M_tr = 4.0 Msun
    - sigma = 0.2
    - x = 1.7 (power-law slope for m > M_tr)
    """
    
    def __init__(self, mc: float = 3.5, m_tr: float = 4.0, 
                 sigma: float = 0.2, x: float = 1.7):
        self.mc = mc
        self.m_tr = m_tr
        self.sigma = sigma
        self.x = x
        self.pi = np.pi
    
    def _erf(self, z: float) -> float:
        """Error function approximation."""
        from scipy.special import erf
        return erf(z)
    
    def cdf(self, m: float, mmin: float, mmax: float) -> float:
        """
        Cumulative distribution for Chabrier IMF.
        
        This is a direct translation of the FORTRAN code.
        """
        if m <= mmin:
            return 0.0
        if m >= mmax:
            return 1.0
        
        arg_erf1 = np.log10(self.m_tr / self.mc) / (np.log10(10.0) * self.sigma * np.sqrt(2.0))
        arg_erf2 = np.log10(mmin / self.mc) / (np.log10(10.0) * self.sigma * np.sqrt(2.0))
        arg_erf = np.log10(m / self.mc) / (np.log10(10.0) * self.sigma * np.sqrt(2.0))
        
        term1 = self.sigma * np.sqrt(self.pi / 2.0) * np.log10(10.0) * (self._erf(arg_erf1) - self._erf(arg_erf2))
        term2 = np.exp(-(np.log10(self.m_tr / self.mc))**2 / (2 * self.sigma**2)) / (self.x * self.m_tr**(-self.x))
        term2 *= (mmax**(-self.x) - self.m_tr**(-self.x))
        
        denom = term1 - term2
        
        A = 1.0 / denom
        B = A * np.exp(-((np.log10(self.m_tr / self.mc) / self.sigma)**2.0) / 2.0) / self.m_tr**(-self.x)
        
        if m <= self.m_tr:
            cdf_val = A * np.sqrt(self.pi / 2.0) * self.sigma * np.log10(10.0)
            cdf_val *= (self._erf(arg_erf) - self._erf(arg_erf2))
        else:
            cdf_val = A * np.sqrt(self.pi / 2.0) * self.sigma * np.log10(10.0)
            cdf_val *= (self._erf(arg_erf1) - self._erf(arg_erf2))
            cdf_val -= B * (m**(-self.x) - self.m_tr**(-self.x)) / self.x
        
        return cdf_val


def get_imf(imf_type: int) -> IMF:
    """Factory function to get IMF instance."""
    if imf_type == 1:
        return SalpeterIMF()
    elif imf_type == 2:
        return KroupaIMF()
    elif imf_type == 3:
        return ChabrierIMF()
    else:
        logger.warning(f"Unknown IMF type {imf_type}, defaulting to Salpeter")
        return SalpeterIMF()


# =============================================================================
# DATA LOADING AND INTERPOLATION
# =============================================================================

class DataTable:
    """
    Base class for loading and interpolating tabulated data.
    
    Replaces FORTRAN routines: choose_metallicity.f, read_mremnant.f, read_ionise.f
    """
    def __init__(self, data_dir: Optional[str] = None):
        if data_dir is None:
            # Resolves path inside the installed package bundle
            self.data_dir = (
                importlib.resources.files("realta") / "data"
            )  #[cite: 2]
        else:
            self.data_dir = Path(data_dir)
        self.loaded = False
    
    def load(self):
        """Load data from file."""
        raise NotImplementedError
    
    def interpolate(self, x: float) -> float:
        """Linear interpolation in log-log space (as in original FORTRAN)."""
        raise NotImplementedError


class LifetimeTable(DataTable):
    """
    Stellar lifetime data table.
    
    Original FORTRAN: choose_metallicity.f, lifetime.f (get_lifetime)
    
    Loads lifetime data for different metallicities and provides interpolation.
    """
    
    METAL_FILES = {
        1: 'lifetimes_z0.dat',
        2: 'lifetimes_z8e-3.dat',
        3: 'lifetimes_z2e-2.dat'
    }
    
    def __init__(self, imetal: int = 2, data_dir: str = 'data'):
        super().__init__(data_dir)
        self.imetal = imetal
        self.mass = None
        self.lifetime = None
        self.load()
    
    def load(self):
        """Load lifetime data from file."""
        filename = self.METAL_FILES.get(self.imetal, self.METAL_FILES[1])
        filepath = self.data_dir / filename
        
        if not filepath.exists():
            logger.warning(f"Lifetime data file {filepath} not found. "
                          "Using placeholder data.")
            self._create_placeholder_data()
            return
        
        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()
            
            # Skip header (6 lines in original FORTRAN)
            data_lines = lines[6:]
            
            masses = []
            lifetimes = []
            
            for line in data_lines:
                if line.strip() and not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 2:
                        masses.append(float(parts[0]))
                        lifetimes.append(float(parts[1]))
            
            # Convert to log10 as in original FORTRAN
            self.mass = np.log10(np.array(masses))
            self.lifetime = np.log10(np.array(lifetimes))
            self.loaded = True
            
        except Exception as e:
            logger.error(f"Error loading lifetime data: {e}")
            self._create_placeholder_data()
    
    def _create_placeholder_data(self):
        """Create placeholder data for testing."""
        # Simple approximation: lifetime ∝ mass^(-2.5)
        masses = np.logspace(-1, 2, 100)  # 0.1 to 100 Msun
        lifetimes = 10.0 * masses**(-2.5)  # in Myr
        
        self.mass = np.log10(masses)
        self.lifetime = np.log10(lifetimes)
        self.loaded = True
    
    def get_lifetime(self, star_mass: float) -> float:
        """
        Get stellar lifetime in Myr using linear interpolation in log-log space.
        
        Original FORTRAN: get_lifetime() function
        """
        if not self.loaded:
            self.load()
        
        lmass = np.log10(star_mass)
        
        # Find the interval (replaces the do-while loop in FORTRAN)
        idx = np.searchsorted(self.mass, lmass, side='right') - 1
        idx = max(0, min(idx, len(self.mass) - 2))
        
        # Linear interpolation: y = ax + b
        a = (self.lifetime[idx+1] - self.lifetime[idx]) / (self.mass[idx+1] - self.mass[idx])
        b = self.lifetime[idx] - a * self.mass[idx]
        
        log_lifetime = a * lmass + b
        return 10.0**log_lifetime


class RemnantTable(DataTable):
    """
    Remnant mass data table.
    
    Original FORTRAN: read_mremnant.f, get_mremnant.f
    """
    
    def __init__(self, data_dir: str = 'data'):
        super().__init__(data_dir)
        self.minit = None  # Initial mass
        self.mfin = None   # Final (remnant) mass
        self.load()
    
    def load(self):
        """Load remnant mass data from file."""
        filepath = self.data_dir / 'remnant_masses.dat'
        
        if not filepath.exists():
            logger.warning(f"Remnant data file {filepath} not found. "
                          "Using placeholder data.")
            self._create_placeholder_data()
            return
        
        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()
            
            # Skip header (5 lines in original FORTRAN)
            data_lines = lines[5:]
            
            minit = []
            mfin = []
            
            for line in data_lines:
                if line.strip() and not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 2:
                        minit.append(float(parts[0]))
                        mfin.append(float(parts[1]))
            
            # Convert to log10 as in original FORTRAN
            self.minit = np.log10(np.array(minit))
            self.mfin = np.log10(np.array(mfin))
            self.loaded = True
            
        except Exception as e:
            logger.error(f"Error loading remnant data: {e}")
            self._create_placeholder_data()
    
    def _create_placeholder_data(self):
        """Create placeholder remnant mass data."""
        # Simple approximation
        minit = np.logspace(0, 2, 100)  # 1 to 100 Msun
        # Remnant mass approximation
        mfin = np.zeros_like(minit)
        for i, m in enumerate(minit):
            if m < 8:
                mfin[i] = 0.0  # No remnant for low-mass stars
            elif m < 20:
                mfin[i] = 1.4  # Neutron star
            elif m < 40:
                mfin[i] = 5.0  # Low-mass black hole
            else:
                mfin[i] = 10.0  # High-mass black hole
        
        self.minit = np.log10(minit)
        self.mfin = np.log10(mfin)
        self.loaded = True
    
    def get_remnant_mass(self, star_mass: float) -> float:
        """
        Get remnant mass in solar masses.
        
        Original FORTRAN: get_mremnant() function
        """
        if not self.loaded:
            self.load()
        
        lmass = np.log10(star_mass)
        
        # Find the interval
        idx = np.searchsorted(self.minit, lmass, side='right') - 1
        idx = max(0, min(idx, len(self.minit) - 2))
        
        # Linear interpolation
        a = (self.mfin[idx+1] - self.mfin[idx]) / (self.minit[idx+1] - self.minit[idx])
        b = self.mfin[idx] - a * self.minit[idx]
        
        log_mfin = a * lmass + b
        return 10.0**log_mfin


class IonizingPhotonTable(DataTable):
    """
    Ionizing photon data table.
    
    Original FORTRAN: read_ionise.f, get_ngamma.f
    """
    
    # Constants from original FORTRAN
    MUNIT = 1.99e30  # Solar mass in grams
    MATOM = 1.67e-27  # Atomic mass unit in grams
    
    def __init__(self, data_dir: str = 'data'):
        super().__init__(data_dir)
        self.mstar = None  # Stellar mass
        self.ngamma = None # log10(number of ionizing photons per solar mass)
        self.load()
    
    def load(self):
        """Load ionizing photon data from file."""
        filepath = self.data_dir / 'ionise.dat'
        
        if not filepath.exists():
            logger.warning(f"Ionizing photon data file {filepath} not found. "
                          "Using placeholder data.")
            self._create_placeholder_data()
            return
        
        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()
            
            # Skip header (1 line in original FORTRAN)
            data_lines = lines[1:]
            
            mstar = []
            ngamma_raw = []
            
            for line in data_lines:
                if line.strip() and not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 2:
                        mstar.append(float(parts[0]))
                        ngamma_raw.append(float(parts[1]))
            
            # Convert as in original FORTRAN:
            # ngamma = log10(ngamma_raw) + log10(mstar) + log10(munit) - log10(matom)
            mstar_arr = np.array(mstar)
            ngamma_arr = np.array(ngamma_raw)
            
            self.mstar = np.log10(mstar_arr)
            self.ngamma = np.log10(ngamma_arr) + np.log10(mstar_arr) + \
                         np.log10(self.MUNIT) - np.log10(self.MATOM)
            self.loaded = True
            
        except Exception as e:
            logger.error(f"Error loading ionizing photon data: {e}")
            self._create_placeholder_data()
    
    def _create_placeholder_data(self):
        """Create placeholder ionizing photon data."""
        mstar = np.logspace(0, 2, 100)  # 1 to 100 Msun
        # Simple approximation: ngamma ∝ mass^2 for massive stars
        ngamma_raw = np.zeros_like(mstar)
        for i, m in enumerate(mstar):
            if m >= 8:
                ngamma_raw[i] = m**2 * 1e48  # photons per second
            else:
                ngamma_raw[i] = 0.0
        
        self.mstar = np.log10(mstar)
        self.ngamma = np.log10(ngamma_raw) + np.log10(mstar) + \
                     np.log10(self.MUNIT) - np.log10(self.MATOM)
        self.loaded = True
    
    def get_ngamma(self, star_mass: float) -> float:
        """
        Get number of ionizing photons per unit solar mass.
        
        Original FORTRAN: get_ngamma() function
        Returns -10.0 for stars with mass < 8 Msun
        """
        if not self.loaded:
            self.load()
        
        if star_mass < 8.0:
            return -10.0
        
        lmass = np.log10(star_mass)
        
        # Find the interval
        idx = np.searchsorted(self.mstar, lmass, side='right') - 1
        idx = max(0, min(idx, len(self.mstar) - 2))
        
        # Linear interpolation
        a = (self.ngamma[idx+1] - self.ngamma[idx]) / (self.mstar[idx+1] - self.mstar[idx])
        b = self.ngamma[idx] - a * self.mstar[idx]
        
        return a * lmass + b


# =============================================================================
# X-RAY LUMINOSITY
# =============================================================================

class XRayLuminosity:
    """
    X-ray luminosity calculator for binaries.
    
    Original FORTRAN: get_lumx.f
    """
    
    def __init__(self, lxmin: float, lxmax: float, lunit: float = 1.0e33):
        self.lxmin = lxmin
        self.lxmax = lxmax
        self.lunit = lunit
        self.eta = 0.1  # Accretion efficiency
        
        # Weibull distribution parameters (from FORTRAN)
        self.lambda_ = 0.5
        self.k = 1.9
    
    def eddington_luminosity(self, mass: float) -> float:
        """
        Calculate Eddington luminosity in erg/s.
        
        Original FORTRAN: ledd = 10**(alog10(1.3)+38-alog10(lunit))*(massp)
        """
        return 10.0**(np.log10(1.3) + 38 - np.log10(self.lunit)) * mass
    
    def get_lumx(self, massp: float, masss: float, period: float, a: float,
                  iseed: Optional[int] = None, use_weibull: bool = True) -> float:
        """
        Get X-ray luminosity for a binary.
        
        Original FORTRAN: get_lumx() function
        
        Parameters:
        -----------
        massp : float
            Primary mass in solar masses
        masss : float
            Secondary mass in solar masses
        period : float
            Orbital period in days
        a : float
            Semi-major axis in AU
        iseed : int, optional
            Random seed (if None, uses uniform distribution)
        use_weibull : bool
            If True, use Weibull distribution; if False, use uniform
            
        Returns:
        --------
        float
            X-ray luminosity in erg/s
        """
        ledd = self.eddington_luminosity(massp)
        
        if self.lxmin == self.lxmax:
            return self.lxmin
        
        if iseed is None or iseed >= 0:
            # Uniform distribution
            # Original: 10**(alog10(lxmin) + ((alog10(lxmax)-alog10(lxmin))*ran3(iseed)))
            # Using numpy random
            u = np.random.random()
            log_lx = np.log10(self.lxmin) + (np.log10(self.lxmax) - np.log10(self.lxmin)) * u
            lumx = 10.0**log_lx
        else:
            # Peaked (Weibull) distribution
            # This is a direct translation of the FORTRAN code
            # Note: The original FORTRAN has a loop (goto 232) to ensure lumx <= ledd
            while True:
                xmprob = self.lambda_ * ((self.k - 1.0) / self.k)**(1.0 / self.k)
                # Generate random number
                u = np.random.random()
                xprob = self.lambda_ * (-np.log(u))**(1.0 / self.k)
                
                get_lumx = (np.log10(self.lxmax) - np.log10(self.lxmin)) * xprob / xmprob
                get_lumx = np.log10(self.lxmin) + get_lumx
                lumx = 10.0**get_lumx
                
                if lumx <= ledd:
                    break
        
        return lumx


# =============================================================================
# BINARY POPULATION
# =============================================================================

@dataclass
class Binary:
    """
    Represents a binary star system.
    
    Attributes:
        primary_mass : float
            Primary star mass in solar masses
        secondary_mass : float
            Secondary star mass in solar masses
        period : float
            Orbital period in days
        a : float
            Semi-major axis in AU
        turnoff_time : float
            Time at which the primary turns off the main sequence (Myr)
        nturn : int
            Number of stars that have turned off MS (0, 1, or 2)
        lum_xray : float
            Current X-ray luminosity in erg/s
        index : int
            Index in the population (for sorting)
    """
    primary_mass: float
    secondary_mass: float
    period: float
    a: float
    turnoff_time: float = 0.0
    nturn: int = 0
    lum_xray: float = 0.0
    index: int = 0


class BinaryPopulation:
    """
    Manages a population of binary stars.
    
    Original FORTRAN: make_stars.f, plus evolution logic from main.f
    """
    
    # Constants from original FORTRAN
    PFAC = 365.229126  # Days to years
    AFAC = 0.0193852859  # Conversion factor for semi-major axis
    
    def __init__(self, config: SimulationConfig):
        self.config = config
        self.binaries: List[Binary] = []
        self.rng = RandomGenerator(config.iseed)
        self.imf = get_imf(config.imf_type)
        
        # Data tables
        self.lifetime_table = LifetimeTable(config.imetal, config.data_dir)
        self.remnant_table = RemnantTable(config.data_dir)
        self.ionizing_table = IonizingPhotonTable(config.data_dir)
        
        # X-ray calculator
        self.xray_calc = XRayLuminosity(
            lxmin=10.0**(config.lxmin - np.log10(config.lunit)),
            lxmax=10.0**(config.lxmax - np.log10(config.lunit)),
            lunit=config.lunit
        )
        
        # Generate population
        self.generate_population()
    
    def generate_population(self):
        """
        Generate initial binary population.
        
        Original FORTRAN: make_stars() subroutine
        """
        ntot = self.config.ntot
        mmin = self.config.mmin
        mmax = self.config.mmax
        mcut = self.config.mcut
        pmin = self.config.pmin
        pmax = self.config.pmax
        mcomp = self.config.mcomp
        fbin = self.config.fbin
        
        logger.info(f"Generating {ntot} stars with IMF type {self.config.imf_type}")
        
        # Sample masses from IMF
        masses = self.imf.sample(ntot, mmin, mmax, self.rng)
        
        # Initialize arrays
        primary_masses = []
        secondary_masses = []
        periods = []
        semi_major_axes = []
        
        nmass = 0
        nhmxb = 0
        
        for mass in masses:
            primary_masses.append(mass)
            
            # Determine if this is a binary
            if self.rng.random() <= fbin and mass >= mcut:
                nmass += 1
                
                # Draw binary period from logarithmic distribution
                # Original FORTRAN: pstar = exp(log(pmax/pmin)*ran3(iseed))
                log_period = np.log(pmax / pmin) * self.rng.random()
                period = np.exp(log_period)
                periods.append(period)
                
                # Draw companion mass
                # Original FORTRAN: qmass = mcmpct + (mass - mcmpct) * ran3(iseed)
                if mcomp < 0:
                    companion_mass = mmin + (mass - mmin) * self.rng.random()
                else:
                    companion_mass = mcomp + (mass - mcomp) * self.rng.random()
                
                # Ensure companion mass doesn't exceed primary
                companion_mass = min(companion_mass, mass)
                secondary_masses.append(companion_mass)
                
                if companion_mass >= abs(mcomp):
                    nhmxb += 1
                
                # Calculate semi-major axis
                # Original FORTRAN: a = afac * mass**(1/3) * (1 + companion/mass)**(1/3)
                # Then: a = a * period**(2/3)
                a_val = (self.AFAC * mass**(1.0/3.0) * 
                        (1.0 + companion_mass / mass)**(1.0/3.0))
                a_val *= period**(2.0/3.0)
                semi_major_axes.append(a_val)
                
            else:
                # Isolated star
                secondary_masses.append(0.0)
                periods.append(0.0)
                semi_major_axes.append(0.0)
        
        # Sort by primary mass (descending) as in original FORTRAN
        # This uses the indexx() approach from Numerical Recipes
        indices = np.argsort(-np.array(primary_masses))
        
        # Create binary objects
        self.binaries = []
        for i in indices:
            binary = Binary(
                primary_mass=primary_masses[i],
                secondary_mass=secondary_masses[i],
                period=periods[i],
                a=semi_major_axes[i],
                index=i
            )
            
            # Calculate turnoff time for primary
            if primary_masses[i] >= self.config.mcut:
                binary.turnoff_time = self.lifetime_table.get_lifetime(primary_masses[i])
            
            self.binaries.append(binary)
        
        # Filter to only massive stars (primary_mass >= mcut)
        self.binaries = [b for b in self.binaries if b.primary_mass >= mcut]
        
        logger.info(f"Generated {len(self.binaries)} massive binaries, "
                   f"{nhmxb} HMXB progenitors")
    
    def sort_by_turnoff_time(self):
        """Sort binaries by turnoff time (ascending)."""
        self.binaries.sort(key=lambda b: b.turnoff_time)
        # Update indices
        for i, b in enumerate(self.binaries):
            b.index = i
    
    def evolve(self, tnow: float, dt: float) -> Tuple[float, float, int, int]:
        """
        Evolve the binary population to time tnow.
        
        Returns:
            Tuple containing:
            - Total X-ray luminosity (erg/s)
            - Total number of ionizing photons
            - Number of active binaries
            - Number of dead binaries (both stars are remnants)
        """
        lumx_tot = 0.0
        nphot_tot = 0.0
        nactive = 0
        ndead = 0
        
        # Sort by turnoff time for efficient processing
        self.sort_by_turnoff_time()
        
        # Process each binary
        for binary in self.binaries:
            # Initialize X-ray luminosity
            binary.lum_xray = 0.0
            
            # Check if we've passed the turnoff time
            if tnow < binary.turnoff_time:
                continue
            
            # Both stars are remnants (binary is dead)
            if binary.turnoff_time == 0.0:
                ndead += 1
                continue
            
            # Determine which star is turning off
            if binary.nturn == 0:
                # Primary is turning off
                primary_remnant = self.remnant_table.get_remnant_mass(binary.primary_mass)
                deltam = binary.primary_mass - primary_remnant
                floss = deltam / (binary.primary_mass + binary.secondary_mass)
                
                # Update primary mass to remnant
                binary.primary_mass = primary_remnant
                
                # Update turnoff time to secondary's turnoff time
                binary.turnoff_time = self.lifetime_table.get_lifetime(binary.secondary_mass)
                binary.nturn = 1
                nactive += 1
                
            elif binary.nturn == 1:
                # Secondary is turning off
                secondary_remnant = self.remnant_table.get_remnant_mass(binary.secondary_mass)
                deltam = binary.secondary_mass - secondary_remnant
                floss = deltam / (binary.primary_mass + binary.secondary_mass)
                
                # Update secondary mass to remnant
                binary.secondary_mass = secondary_remnant
                
                # Both are now remnants
                binary.turnoff_time = 0.0
                binary.nturn = 2
                ndead += 1
                
                # Skip X-ray calculation for dead binary
                continue
            
            # Update semi-major axis if mass loss fraction is <= 0.5
            if floss <= 0.5:
                binary.a *= deltam / (binary.primary_mass + binary.secondary_mass)
                
                # Update period
                binary.period = (self.PFAC * np.sqrt(binary.a**3 / 
                           (binary.primary_mass + binary.secondary_mass)))
            
            # Calculate X-ray luminosity if binary is still active
            if binary.turnoff_time > 0.0 and floss <= 0.5:
                if binary.secondary_mass >= abs(self.config.mcomp):
                    if self.rng.random() <= self.config.fbin:
                        binary.lum_xray = self.xray_calc.get_lumx(
                            binary.primary_mass,
                            binary.secondary_mass,
                            binary.period,
                            binary.a,
                            iseed=None,  # Use uniform distribution
                            use_weibull=True
                        )
            
            # Add to totals
            lumx_tot += binary.lum_xray
            
            # Calculate ionizing photons from both stars
            if binary.nturn == 0:
                # Both stars on MS
                ng1 = self.ionizing_table.get_ngamma(binary.primary_mass)
                ng2 = self.ionizing_table.get_ngamma(binary.secondary_mass)
                t1 = self.lifetime_table.get_lifetime(binary.primary_mass)
                t2 = self.lifetime_table.get_lifetime(binary.secondary_mass)
                
                nphot_tot += 10.0**(ng1 + np.log10(dt/t1) + ng2 + np.log10(dt/t2) - 60)
            elif binary.nturn == 1:
                # Primary is remnant, secondary on MS
                ng2 = self.ionizing_table.get_ngamma(binary.secondary_mass)
                t2 = self.lifetime_table.get_lifetime(binary.secondary_mass)
                
                nphot_tot += 10.0**(ng2 + np.log10(dt/t2) - 60)
        
        return lumx_tot, nphot_tot, nactive, ndead


# =============================================================================
# MAIN SIMULATION
# =============================================================================

class ClusterSimulation:
    """
    Main simulation class.
    
    Original FORTRAN: main.f
    
    This class orchestrates the entire simulation, including:
    - Loading configuration
    - Generating initial population
    - Running time evolution
    - Writing output
    """
    
    def __init__(self, config: Optional[SimulationConfig] = None):
        if config is None:
            self.config = load_config()
        else:
            self.config = config
        
        self.population: Optional[BinaryPopulation] = None
    
    def initialize(self):
        """Initialize the simulation."""
        logger.info("Initializing simulation...")
        self.population = BinaryPopulation(self.config)
        logger.info("Simulation initialized.")
    
    def run(self, output_dir: str = 'output') -> Dict:
        """
        Run the simulation.
        
        Returns:
            Dictionary containing simulation results
        """
        if self.population is None:
            self.initialize()
        
        # Create output directory
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        # Write initial conditions
        self._write_initial_conditions(output_dir)
        
        # Time evolution
        tmax = self.config.tmax
        dt = self.config.dt
        
        tnow = 0.0
        results = []
        
        logger.info(f"Starting time evolution to {tmax} Myr with dt={dt} Myr")
        
        while tnow <= tmax:
            lumx_tot, nphot_tot, nactive, ndead = self.population.evolve(tnow, dt)
            
            results.append({
                'time': tnow,
                'lumx_tot': lumx_tot,
                'nphot_tot': nphot_tot,
                'nactive': nactive,
                'ndead': ndead
            })
            
            tnow += dt
        
        # Write results
        self._write_results(results, output_dir)
        
        logger.info("Simulation complete.")
        return results
    
    def _write_initial_conditions(self, output_dir: str):
        """Write initial conditions to file."""
        imf_name = {1: 'Salpeter', 2: 'Kroupa', 3: 'Chabrier'}
        metal_name = {1: 'Z=0', 2: 'Z=0.008', 3: 'Z=0.02'}
        
        filename = Path(output_dir) / f"{imf_name[self.config.imf_type]}.init.dat"
        
        with open(filename, 'w') as f:
            f.write(f"# {imf_name[self.config.imf_type]} IMF\n")
            f.write(f"# ntot (mmin,mmax,mcut)/Msol (pmin,pmax)/days (lxmin,lxmax)/ergs/s\n")
            f.write(f"{self.config.ntot} {self.config.mmin} {self.config.mmax} "
                   f"{self.config.mcut} {self.config.pmin} {self.config.pmax} "
                   f"{self.config.lxmin} {self.config.lxmax}\n")
            f.write(f"# n (m1,m2)/M* P/days a/AU (t1,t2)/Myrs (mr1,mr2)/M*\n")
            
            for i, binary in enumerate(self.population.binaries):
                t1 = self.population.lifetime_table.get_lifetime(binary.primary_mass)
                t2 = self.population.lifetime_table.get_lifetime(binary.secondary_mass)
                mr1 = self.population.remnant_table.get_remnant_mass(binary.primary_mass)
                mr2 = self.population.remnant_table.get_remnant_mass(binary.secondary_mass)
                
                f.write(f"{i+1} {binary.primary_mass:12.4f} {binary.secondary_mass:12.4f} "
                       f"{binary.period:12.4f} {binary.a:12.4f} {t1:12.4f} {t2:12.4f} "
                       f"{mr1:12.4f} {mr2:12.4f}\n")
        
        logger.info(f"Initial conditions written to {filename}")
    
    def _write_results(self, results: List[Dict], output_dir: str):
        """Write simulation results to file."""
        imf_name = {1: 'Salpeter', 2: 'Kroupa', 3: 'Chabrier'}
        
        filename = Path(output_dir) / f"{imf_name[self.config.imf_type]}.tevol.dat"
        
        with open(filename, 'w') as f:
            f.write(f"# {imf_name[self.config.imf_type]} IMF\n")
            f.write(f"# ntot (mmin,mmax,mcut)/Msol (pmin,pmax)/days\n")
            f.write(f"{self.config.ntot} {self.config.mmin} {self.config.mmax} "
                   f"{self.config.mcut} {self.config.pmin} {self.config.pmax}\n")
            f.write(f"# t/Myrs lx_tot/ergs nphot npop ndead\n")
            
            for r in results:
                f.write(f"{r['time']:18.8e} {r['lumx_tot']:18.8e} "
                       f"{r['nphot_tot']:18.8e} {r['nactive']:9d} {r['ndead']:9d}\n")
        
        logger.info(f"Results written to {filename}")


# =============================================================================
# COMMAND-LINE INTERFACE
# =============================================================================

def main():
    """Command-line interface for running the simulation."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='HMXRB Simulation - Power et al. 2009 Python Refactoring'
    )
    parser.add_argument(
        '--config', '-c',
        type=str,
        default=None,
        help='Path to configuration YAML file'
    )
    parser.add_argument(
        '--output', '-o',
        type=str,
        default='output',
        help='Output directory'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging'
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Load configuration
    config = load_config(args.config)
    
    # Run simulation
    sim = ClusterSimulation(config)
    results = sim.run(args.output)
    
    # Print summary
    print(f"\nSimulation complete!")
    print(f"Output directory: {args.output}")
    print(f"Total binaries: {len(sim.population.binaries)}")
    print(f"Final time: {results[-1]['time']} Myr")
    print(f"Final X-ray luminosity: {results[-1]['lumx_tot']:.2e} erg/s")


if __name__ == '__main__':
    main()
