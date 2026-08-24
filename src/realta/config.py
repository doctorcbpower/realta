from __future__ import annotations

import dataclasses
import os
from dataclasses import dataclass

import yaml

# Paper 1 binary-interaction/merger prescriptions -- see
# docs/science/paper1-binary-interaction-proposal.md for the full
# rationale. None of this is derived from Power et al. (2009/2013);
# it is a new, Realta-specific parameterization layered on the existing
# floss/fsur gate, reviewed and accepted 2026-08-24. "non_interacting"
# is the default and reproduces the pre-existing baseline exactly.
BINARY_PRESCRIPTIONS = (
    "single",
    "non_interacting",
    "standard_interaction",
    "enhanced_interaction",
    "enhanced_mergers",
)

# Illustrative, not-paper-derived defaults implied by each prescription
# for the parameters below, when the parameter itself is left
# unspecified (None) in the config. A prescription not listed for a
# given parameter leaves that parameter at its neutral (no-effect)
# value.
#
# Reconciled with the physics-based RLOF classifier (Hurley/Tout et al.
# -- see binaries/interaction.py, docs/science/rlof-ce-classifier-
# proposal.md "Decision 3") on 2026-08-24: standard_interaction/
# enhanced_interaction/enhanced_mergers now turn on
# use_rlof_classifier=True and drive their behaviour through real
# Roche-lobe physics rather than purely through interaction_boost/
# p_merge alone. This is a genuine scientific-behaviour change for
# these three prescriptions specifically (NOT for "single" or
# "non_interacting", which are unaffected -- use_rlof_classifier still
# defaults to False for those, so the pinned baseline is untouched) --
# see docs/provenance.md Section 6 for the old-vs-new pinned values and
# the explicit rationale, per the brief's "never disguise a scientific
# change as a refactor" principle. `p_merge` is deliberately NOT set
# for enhanced_mergers any more (was 0.2/10.0/0.5) -- the real
# classifier's IMMEDIATE_MERGER outcome, driven by a lowered
# `q_crit_ms`, replaces it as that prescription's merger driver; the
# formation-time p_merge/p_merge_max_period/f_merge channel remains
# fully available as an independent, explicit override for anyone who
# wants to layer both.
_PRESCRIPTION_DEFAULTS: dict[str, dict[str, float | bool]] = {
    "standard_interaction": {"interaction_boost": 1.5, "use_rlof_classifier": True},
    "enhanced_interaction": {"interaction_boost": 3.0, "use_rlof_classifier": True},
    "enhanced_mergers": {
        # Lower than HTP02's own fiducial q_crit_ms=0.695 (see
        # binaries/interaction.py::Q_CRIT_MS) -- more RLOF-ing systems
        # are classified as dynamically-unstable mergers rather than
        # stable mass transfer. NOT paper-derived; an illustrative,
        # named, revisable choice (same status as Q_CRIT_MS itself,
        # which the task's own brief flags as uncertain/literature-
        # varying), not a new independent free parameter.
        "q_crit_ms": 0.4,
        "use_rlof_classifier": True,
    },
}


@dataclass
class SimulationConfig:
    """Configuration for the HMXRB simulation."""

    ntot: int = 100000
    # Lower stellar mass bound, Msun. 0.1 is the practical stellar
    # lower-mass cutoff (objects below ~0.08 Msun are substellar, not
    # true hydrogen-burning stars) -- see imf/kroupa.py's class
    # docstring for how this interacts with the Kroupa IMF's breakpoints.
    mmin: float = 0.1
    mmax: float = 100.0
    mcut: float = 8.0
    tmax: float = 100.0
    dt: float = 0.01

    # IMF type: 1=Salpeter, 2=Kroupa, 3=Chabrier
    imf_type: int = 2

    # Continuous single power-law IMF slope (dN/dm ~ m^-imf_slope),
    # overriding SalpeterIMF's own default alpha=2.35 -- only applies
    # when imf_type=1; ignored otherwise. If None, SalpeterIMF's own
    # default is used, reproducing the pre-existing baseline exactly.
    # See imf/factory.py::get_imf's docstring for why this is
    # Salpeter-only (docs/science/paper1-detailed-work-breakdown.md,
    # item A4).
    imf_slope: float | None = None

    # Binary parameters
    pmin: float = 0.1
    pmax: float = 1000.0
    mcomp: float = 0.5

    # Fraction of M >= mcut stars that get a companion at all (A1,
    # docs/science/paper1-detailed-work-breakdown.md -- roadmap item
    # 7). Default 1.0 reproduces the pre-existing Power et al. (2009)
    # baseline exactly (every massive star paired, Sec. 2.1) -- the
    # per-star Bernoulli draw this introduces is skipped entirely at
    # binary_fraction=1.0 (see generate_population), so the RNG stream
    # is untouched by default, the same pattern already used for
    # p_merge==0. A star that doesn't get a companion (m2=0, period=0,
    # a=0 -- placeholder values, not a real orbit) is still tracked
    # through SN/lifetime bookkeeping, unlike the "single"
    # binary_prescription (which empties the tracked array entirely) --
    # this generalizes that special case to any fraction in [0, 1].
    binary_fraction: float = 1.0

    # Mass-ratio distribution for assigned companions: "uniform"
    # (default, unchanged baseline -- m2 ~ Uniform(mcomp, m1), Power et
    # al. 2009 Sec. 2.1) or "flat_q" (m2 = m1 * Uniform(0, 1), flat in
    # mass ratio q rather than in absolute companion mass -- the
    # alternative development-roadmap.md item 7's own example names).
    mass_ratio_distribution: str = "uniform"

    # Period distribution: "log_uniform" (default, unchanged baseline
    # -- log P ~ Uniform(log pmin, log pmax)) or "log_normal" (the
    # other alternative development-roadmap.md item 7 names). The
    # log_normal parameters are NOT literature-sourced -- a named,
    # generic simplification (like Q_CRIT_MS elsewhere): mean and
    # width are derived directly from pmin/pmax (mu = midpoint of
    # log10(pmin)/log10(pmax), sigma = range/6, i.e. pmin/pmax sit at
    # roughly +-3 sigma), truncated to [pmin, pmax] so it can never
    # sample outside the configured bounds. See
    # binaries/population.py::generate_population for the exact
    # formula.
    period_distribution: str = "log_uniform"

    # Probability that a binary which remains bound after the primary
    # supernova (floss <= 0.5, the deterministic sudden-mass-loss
    # criterion -- Power et al. 2009, Section 2.1) goes on to be counted
    # as an active HMXB. This is the paper's own f_sur: "if the binary
    # remains bound, then it has a probability of f_sur that it will
    # evolve into a HMXB" (Sec. 2.1). Named `fsur` here to match the
    # paper's own notation -- note this is NOT a primordial binary
    # fraction -- every massive star above `mcut` is assigned a
    # companion at formation regardless. fsur=1 (Fig. 1's baseline)
    # means every bound binary becomes an active HMXB; the paper's
    # Figs 2-3 explore fsur < 1 as a post-hoc linear rescaling of a
    # single fsur=1 run's aggregate output rather than by re-drawing
    # the population, since the two are equivalent in expectation --
    # Realta instead re-draws per run, so re-run with a different
    # `fsur` rather than rescaling `lumx_tot`/`nphot_tot` externally.
    fsur: float = 0.5

    # Paper 1 binary-interaction/merger prescription -- see
    # docs/science/paper1-binary-interaction-proposal.md. One of
    # BINARY_PRESCRIPTIONS above. "non_interacting" (default)
    # reproduces the pre-existing floss/fsur-only baseline exactly.
    # "single" additionally suppresses binary formation entirely (no
    # companion, no period, no HMXB channel at all).
    binary_prescription: str = "non_interacting"

    # interaction_boost: multiplicative boost on fsur at HMXB
    # activation (fsur_eff = min(1, fsur * interaction_boost)). If
    # None, resolved from binary_prescription's illustrative default
    # (1.0 -- no boost -- if the prescription doesn't specify one).
    # NOT paper-derived; see the proposal doc.
    interaction_boost: float | None = None

    # Pre-SN merger channel (only reachable when binary_prescription
    # implies nonzero p_merge, i.e. "enhanced_mergers", or when set
    # explicitly). A binary with period < p_merge_max_period is merged
    # with probability p_merge at formation; f_merge is the fraction of
    # m2 retained by the merged star. If None, each resolves from
    # binary_prescription's illustrative default (0.0 -- no mergers --
    # if unspecified). NOT paper-derived; see the proposal doc.
    p_merge: float | None = None
    p_merge_max_period: float | None = None
    f_merge: float | None = None

    # Physics-based MS Roche-lobe-overflow classifier (Hurley et al.
    # 2000/Hurley, Tout & Pols 2002 -- see
    # docs/science/rlof-ce-classifier-proposal.md and
    # binaries/interaction.py). If None, resolved from
    # binary_prescription: True for standard_interaction/
    # enhanced_interaction/enhanced_mergers, False (the existing
    # floss/fsur/interaction_boost-only baseline, exactly preserved)
    # for every other prescription -- see _PRESCRIPTION_DEFAULTS above
    # for the full reconciliation rationale and docs/provenance.md
    # Section 6 for what changed. When True, this adds a new MS-RLOF
    # event (stable mass transfer or immediate merger) to evolve(),
    # AND -- since this reconciliation -- interaction_boost is only
    # applied to binaries whose RLOF outcome was stable mass transfer,
    # not unconditionally to every surviving binary as before. Requires
    # imetal=2 or 3 (Z=0.008 or 0.02) -- the underlying Hurley/Tout
    # formulae are undefined at Z=0 (imetal=1); a run with
    # use_rlof_classifier=True and imetal=1 logs a warning and skips
    # RLOF classification entirely rather than crashing.
    use_rlof_classifier: bool | None = None

    # Post-SN secondary Roche-lobe overflow (docs/science/paper1-
    # followup-prompt.md) -- independent of use_rlof_classifier above,
    # which only covers PRE-SN interaction between two still-live
    # stars. This is the physically dominant real HMXB-formation
    # channel: the secondary's own later RLOF onto the by-then-compact
    # primary (Case B/C mass transfer), checked every timestep once
    # nturn==1 (primary already compact, secondary not yet exploded).
    # Deliberately minimal scope: a single RLOF-only trigger (no wind-
    # accretion/RLOF-fed spectral distinction), no consequence model
    # (secondary mass/envelope and the compact primary's mass are left
    # unchanged -- Hovis-Afflerbach et al. 2025's stripped-donor
    # properties remain an explicit, separate extension point, not
    # implemented here). On trigger, HMXB activation becomes CERTAIN
    # (not a stochastic `fsur` draw) rather than the existing wind-fed
    # approximation `fsur` represents -- real Roche-lobe accretion onto
    # a compact object is qualitatively different from what `fsur` is
    # meant to capture. Default False reproduces the pre-existing
    # baseline exactly (no new code path runs at all). Requires
    # imetal=2 or 3, same as use_rlof_classifier -- skipped with a
    # logged warning at imetal=1 (Z=0).
    use_post_sn_rlof: bool = False

    # Critical mass ratio (q1 = M_donor/M_companion) above which a
    # Roche-lobe-overflowing MS donor merges dynamically rather than
    # transferring mass stably (Hurley, Tout & Pols 2002, Sec. 2.6.4).
    # See binaries/interaction.py::Q_CRIT_MS for the caveat that this
    # value is stated there specifically for deeply-convective donors
    # and is extended here to radiative MS donors as a named
    # simplification. If None, resolved from binary_prescription: 0.4
    # (lower, more mergers) for enhanced_mergers, 0.695 (HTP02's own
    # fiducial value) otherwise.
    q_crit_ms: float | None = None

    # Common-envelope efficiency (alpha_CE) and binding-energy
    # parameter (lambda_CE), HTP02 eqs. 69-73 -- resolves a
    # COMMON_ENVELOPE classification into survive-vs-merge and the new
    # orbit (binaries/interaction.py::apply_common_envelope). See that
    # module's ALPHA_CE/LAMBDA_CE constants for the Zuo & Li (2014)
    # literature basis; both are physically uncertain enough to warrant
    # being overridable, not buried, matching the q_crit_ms pattern.
    # If None, resolved to those module defaults (0.9, 0.5) regardless
    # of binary_prescription -- unlike q_crit_ms, no prescription
    # currently varies these.
    alpha_ce: float | None = None
    lambda_ce: float | None = None

    # Metallicity: 1=Z=0, 2=Z=0.008, 3=Z=0.02
    imetal: int = 2

    # X-ray luminosity bounds for the per-HMXB draw, as log10(erg/s) --
    # e.g. lxmin=33.0 means 1e33 erg/s, NOT the linear value 1e33 itself.
    # config.yml's defaults (33.0, 39.0 -> 1e33-1e39 erg/s).
    lxmin: float = 33.0
    lxmax: float = 39.0
    lunit: float = 1.0e33

    # Shape of the per-binary X-ray luminosity draw (xray/luminosity.py).
    # "weibull": peaked distribution rejection-sampled below the
    #   Eddington luminosity -- Realta's default, matching Power et al.
    #   (2009), Sec. 2.2.
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

        if self.binary_prescription not in BINARY_PRESCRIPTIONS:
            raise ValueError(
                f"binary_prescription must be one of {BINARY_PRESCRIPTIONS}, "
                f"got {self.binary_prescription!r}"
            )
        defaults = _PRESCRIPTION_DEFAULTS.get(self.binary_prescription, {})
        if self.interaction_boost is None:
            self.interaction_boost = defaults.get("interaction_boost", 1.0)
        if self.p_merge is None:
            self.p_merge = defaults.get("p_merge", 0.0)
        if self.p_merge_max_period is None:
            self.p_merge_max_period = defaults.get("p_merge_max_period", 0.0)
        if self.f_merge is None:
            self.f_merge = defaults.get("f_merge", 0.0)
        if self.use_rlof_classifier is None:
            self.use_rlof_classifier = defaults.get("use_rlof_classifier", False)
        if self.q_crit_ms is None:
            self.q_crit_ms = defaults.get("q_crit_ms", 0.695)
        # Local import: realta.binaries's __init__ imports
        # BinaryPopulation, which imports this module -- a circular
        # import at module-load time if imported at the top of this
        # file instead. Deferred to here (runtime, not import time),
        # by which point both modules are already fully initialized.
        from realta.binaries.interaction import ALPHA_CE, LAMBDA_CE

        if self.alpha_ce is None:
            self.alpha_ce = defaults.get("alpha_ce", ALPHA_CE)
        if self.lambda_ce is None:
            self.lambda_ce = defaults.get("lambda_ce", LAMBDA_CE)

        if self.interaction_boost < 0.0:
            raise ValueError(
                f"interaction_boost must be >= 0, got {self.interaction_boost}"
            )
        if not 0.0 <= self.p_merge <= 1.0:
            raise ValueError(f"p_merge must be in [0, 1], got {self.p_merge}")
        if self.q_crit_ms <= 0.0:
            raise ValueError(f"q_crit_ms must be > 0, got {self.q_crit_ms}")
        if self.alpha_ce <= 0.0:
            raise ValueError(f"alpha_ce must be > 0, got {self.alpha_ce}")
        if self.lambda_ce <= 0.0:
            raise ValueError(f"lambda_ce must be > 0, got {self.lambda_ce}")
        if self.imf_slope is not None:
            if self.imf_slope <= 0.0:
                raise ValueError(f"imf_slope must be > 0, got {self.imf_slope}")
            if self.imf_slope == 1.0:
                # SalpeterIMF.cdf's denominator (mmax**beta - mmin**beta,
                # beta = 1 - alpha) is identically 0 at alpha=1.0 --
                # this pre-existing singularity was unreachable before
                # imf_slope made alpha user-configurable (the 2.35
                # default never hits it); reject explicitly rather than
                # silently producing NaN masses.
                raise ValueError(
                    "imf_slope must not be exactly 1.0 -- SalpeterIMF's "
                    "CDF is singular there (a log-form CDF would be "
                    "needed for that case, not implemented)"
                )
        if not 0.0 <= self.binary_fraction <= 1.0:
            raise ValueError(
                f"binary_fraction must be in [0, 1], got {self.binary_fraction}"
            )
        if self.mass_ratio_distribution not in ("uniform", "flat_q"):
            raise ValueError(
                "mass_ratio_distribution must be 'uniform' or 'flat_q', "
                f"got {self.mass_ratio_distribution!r}"
            )
        if self.period_distribution not in ("log_uniform", "log_normal"):
            raise ValueError(
                "period_distribution must be 'log_uniform' or 'log_normal', "
                f"got {self.period_distribution!r}"
            )
        if self.p_merge_max_period < 0.0:
            raise ValueError(
                f"p_merge_max_period must be >= 0, got {self.p_merge_max_period}"
            )
        if not 0.0 <= self.f_merge <= 1.0:
            raise ValueError(f"f_merge must be in [0, 1], got {self.f_merge}")

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
