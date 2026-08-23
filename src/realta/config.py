from __future__ import annotations

import os
from dataclasses import dataclass

import yaml


@dataclass
class SimulationConfig:
    """Configuration for the HMXRB simulation."""

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

    # Probability that a surviving primary-supernova binary is observed as
    # an active HMXB (Power et al. 2009: get_lumx.f / main.f, gated by
    # `ran3(iseed).le.fbin`). NOT a primordial binary fraction -- every
    # massive star above `mcut` is assigned a companion at formation
    # (fpbin=1.0 in the reference make_stars.f); fbin only gates whether a
    # post-supernova binary is counted as X-ray active.
    fbin: float = 0.5

    # Reserved for a future natal-kick / disruption survival prescription
    # (see brief Level 2, "improved natal kicks"). The Power et al. (2009)
    # reference has no such term: binary survival after the primary
    # supernova is governed purely by the deterministic sudden-mass-loss
    # criterion (floss <= 0.5). Defaults to 1.0 (always survive, i.e. a
    # no-op) so the baseline reproduces the reference model; set below 1.0
    # only when deliberately exploring an improved kick-survival model.
    fsur: float = 1.0

    # Metallicity: 1=Z=0, 2=Z=0.008, 3=Z=0.02
    imetal: int = 2

    # X-ray luminosity
    lxmin: float = 33.0
    lxmax: float = 39.0
    lunit: float = 1.0e33

    # Shape of the per-binary X-ray luminosity draw (xray/luminosity.py).
    # "weibull": peaked distribution rejection-sampled below the Eddington
    #   luminosity -- matches every real run of the Fortran reference
    #   (get_lumx.f only takes its "uniform" branch when iseed is exactly
    #   -1, a debug/test sentinel never used by main.f, which always sets
    #   iseed = -abs(iseed)).
    # "uniform": flat log-uniform draw between lxmin and lxmax.
    xray_distribution: str = "weibull"

    # Random seed
    iseed: int = 12345

    # Data directory
    data_dir: str | None = None

    def __post_init__(self):
        if not 0.0 <= self.fbin <= 1.0:
            raise ValueError(f"fbin must be in [0, 1], got {self.fbin}")
        if not 0.0 <= self.fsur <= 1.0:
            raise ValueError(f"fsur must be in [0, 1], got {self.fsur}")
        if self.xray_distribution not in ("weibull", "uniform"):
            raise ValueError(
                "xray_distribution must be 'weibull' or 'uniform', "
                f"got {self.xray_distribution!r}"
            )
        if self.dt <= 0:
            raise ValueError(f"dt must be positive, got {self.dt}")
        if self.pmin >= self.pmax:
            raise ValueError(
                f"pmin ({self.pmin}) must be strictly less than pmax ({self.pmax})"
            )
        if self.mmin >= self.mmax:
            raise ValueError(
                f"mmin ({self.mmin}) must be strictly less than mmax ({self.mmax})"
            )


def load_config(config_path: str | None = None) -> SimulationConfig:
    """Load configuration from YAML file or use defaults."""
    if config_path and os.path.exists(config_path):
        with open(config_path, "r") as f:
            config_dict = yaml.safe_load(f)

        config = SimulationConfig()
        for key, value in config_dict.items():
            if hasattr(config, key):
                if key == "iseed":
                    value = abs(int(value))
                setattr(config, key, value)
        return config
    return SimulationConfig()
