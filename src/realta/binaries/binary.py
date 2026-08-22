from dataclasses import dataclass


@dataclass
class Binary:
    """Represents a binary star system."""

    primary_mass: float
    secondary_mass: float
    period: float
    a: float
    turnoff_time: float = 0.0
    nturn: int = 0
    lum_xray: float = 0.0
    index: int = 0
