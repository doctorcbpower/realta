from __future__ import annotations

import numpy as np


class XRayLuminosity:
    """X-ray luminosity calculator for binaries.

    Reference: Power et al. (2009), MNRAS, 395, 1146 (get_lumx.f).

    The Fortran reference selects between two draws based on the sign/value
    of the random seed passed to it: a flat log-uniform draw when
    ``iseed == -1`` (a debug/test sentinel), and a peaked (Weibull-shaped)
    draw, rejection-sampled below the Eddington luminosity, for every other
    seed value. Since the driver (main.f) always forces the working seed to
    ``-abs(iseed)`` -- never exactly -1 for a real run -- every production
    run of the reference model uses the peaked draw. Realta exposes this as
    an explicit ``distribution`` choice rather than overloading the seed
    value, defaulting to "weibull" to match real Power et al. (2009) runs.
    """

    def __init__(
        self,
        lxmin: float,
        lxmax: float,
        lunit: float = 1.0e33,
        distribution: str = "weibull",
    ):
        self.lxmin = lxmin
        self.lxmax = lxmax
        self.lunit = lunit
        self.eta = 0.1
        self.lambda_ = 0.5
        self.k = 1.9
        if distribution not in ("weibull", "uniform"):
            raise ValueError(f"Unknown X-ray luminosity distribution: {distribution!r}")
        self.distribution = distribution

    def eddington_luminosity(self, mass: float) -> float:
        return 10.0 ** (np.log10(1.3) + 38 - np.log10(self.lunit)) * mass

    def get_lumx(
        self,
        massp: float,
        masss: float,
        period: float,
        a: float,
        rng: np.random.Generator,
    ) -> float:
        """Draw one X-ray luminosity sample for an active HMXB.

        `masss`, `period`, and `a` are accepted for interface parity with
        the reference implementation but, like get_lumx.f, are not used in
        the calculation -- only the primary (compact-object) mass sets the
        Eddington limit that bounds the draw.

        `rng` must be the same seeded `numpy.random.Generator` used for the
        rest of the population, so that a run is fully reproducible from
        `SimulationConfig.iseed`.
        """
        ledd = self.eddington_luminosity(massp)

        if self.lxmin == self.lxmax:
            return self.lxmin

        if self.distribution == "uniform":
            u = rng.random()
            log_lx = (
                np.log10(self.lxmin) + (np.log10(self.lxmax) - np.log10(self.lxmin)) * u
            )
            return 10.0**log_lx

        # Weibull-shaped ("peaked") draw, rejection-sampled below the
        # Eddington limit -- matches the reference implementation's
        # behaviour for every real (non-sentinel) seed.
        xmprob = self.lambda_ * ((self.k - 1.0) / self.k) ** (1.0 / self.k)
        while True:
            u = rng.random()
            xprob = self.lambda_ * (-np.log(u)) ** (1.0 / self.k)
            log_lx = np.log10(self.lxmin) + (
                (np.log10(self.lxmax) - np.log10(self.lxmin)) * xprob / xmprob
            )
            lumx = 10.0**log_lx
            if lumx <= ledd:
                return lumx
