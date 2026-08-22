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
    fbin: float = 0.5
    fsur: float = 0.1

    # Metallicity: 1=Z=0, 2=Z=0.008, 3=Z=0.02
    imetal: int = 2

    # X-ray luminosity
    lxmin: float = 33.0
    lxmax: float = 39.0
    lunit: float = 1.0e33

    # Random seed
    iseed: int = 12345

    # Data directory
    data_dir: str | None = None

    def __post_init__(self):
        if not 0.0 <= self.fbin <= 1.0:
            raise ValueError(f"fbin must be in [0, 1], got {self.fbin}")
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
