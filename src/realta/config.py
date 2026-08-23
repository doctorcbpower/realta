from __future__ import annotations

import dataclasses
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

    # Probability that a binary which remains bound after the primary
    # supernova (floss <= 0.5, the deterministic sudden-mass-loss
    # criterion -- Power et al. 2009, Section 2.1) goes on to be counted
    # as an active HMXB. This is the paper's own f_sur: "if the binary
    # remains bound, then it has a probability of f_sur that it will
    # evolve into a HMXB" (Sec. 2.1). The reference Fortran (main.f)
    # implements exactly this with a variable it happens to call `fbin`
    # (`ran3(iseed).le.fbin`), which does NOT match the paper's own
    # notation and is not a primordial binary fraction -- every massive
    # star above `mcut` is assigned a companion at formation regardless
    # (fpbin=1.0 in the reference make_stars.f). Renamed here to `fsur`
    # to match the paper. fsur=1 (Fig. 1's baseline) means every bound
    # binary becomes an active HMXB; the paper's Figs 2-3 explore
    # fsur < 1 as a post-hoc linear rescaling of a single fsur=1 run's
    # aggregate output rather than by re-drawing the population, since
    # the two are equivalent in expectation -- Realta instead re-draws
    # per run (matching main.f exactly), so re-run with a different
    # `fsur` rather than rescaling `lumx_tot`/`nphot_tot` externally.
    fsur: float = 0.5

    # Metallicity: 1=Z=0, 2=Z=0.008, 3=Z=0.02
    imetal: int = 2

    # X-ray luminosity bounds for the per-HMXB draw, as log10(erg/s) --
    # e.g. lxmin=33.0 means 1e33 erg/s, NOT the linear value 1e33 itself.
    # Matches the reference Fortran parameter file's convention exactly
    # (main.f: `lxmin = 10**(lxmin-alog10(lunit))`) and config.yml's
    # defaults (33.0, 39.0 -> 1e33-1e39 erg/s).
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
        # YAML's default resolver only recognizes scientific notation with
        # an explicit sign after 'e' (1.0e+33), not 1.0e33 -- the latter
        # silently loads as a str, which then fails confusingly deep in
        # downstream numpy calls (e.g. np.log10) rather than here. Coerce
        # int/float fields eagerly so a config typo like this fails loudly
        # and immediately, with a message that points at the actual field.
        for field in dataclasses.fields(self):
            if field.type not in ("int", "float"):
                continue
            value = getattr(self, field.name)
            target_type = int if field.type == "int" else float
            if not isinstance(value, target_type):
                try:
                    setattr(self, field.name, target_type(value))
                except (TypeError, ValueError) as exc:
                    raise ValueError(
                        f"Config field {field.name!r} must be a {field.type}, "
                        f"got {value!r} ({type(value).__name__}). If this came "
                        "from a YAML file, check for scientific notation "
                        "without an explicit sign (use 1.0e+33, not 1.0e33)."
                    ) from exc

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
        # lxmin/lxmax are log10(erg/s) exponents, not linear luminosities
        # (see the field docstring above) -- catch the easy mistake of
        # passing e.g. 1e38 instead of 38.0, which otherwise silently
        # overflows 10**(lxmax - log10(lunit)) to inf downstream.
        if not (0.0 < self.lxmin <= self.lxmax <= 60.0):
            raise ValueError(
                "lxmin/lxmax must be log10(luminosity in erg/s) -- e.g. "
                "lxmin=33.0 for 1e33 erg/s, not the linear value 1e33 -- "
                f"with 0 < lxmin <= lxmax <= 60. Got lxmin={self.lxmin}, "
                f"lxmax={self.lxmax}."
            )


def load_config(config_path: str | None = None) -> SimulationConfig:
    """Load configuration from a flat YAML file, or use defaults.

    The YAML file must be flat -- one key per SimulationConfig field,
    e.g. `ntot: 100000`, not grouped under headings like `simulation:`.
    An earlier version of this function silently accepted (and ignored)
    a grouped/nested file, and separately silently ignored any
    unrecognized key, in both cases falling back to defaults without
    any warning. Both are now errors: a typo, an old-format config file,
    or a nonexistent path should never produce a simulation that looks
    like it used your settings but silently didn't.
    """
    if config_path is None:
        return SimulationConfig()

    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path, "r") as f:
        config_dict = yaml.safe_load(f) or {}

    valid_fields = {field.name for field in dataclasses.fields(SimulationConfig)}
    unknown_keys = sorted(set(config_dict) - valid_fields)
    if unknown_keys:
        raise ValueError(
            f"Unknown key(s) in {config_path}: {unknown_keys}. "
            f"Valid keys are: {sorted(valid_fields)}."
        )

    if "iseed" in config_dict:
        config_dict["iseed"] = abs(int(config_dict["iseed"]))

    return SimulationConfig(**config_dict)
