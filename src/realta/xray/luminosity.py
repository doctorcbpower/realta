from __future__ import annotations

import numpy as np


class XRayLuminosity:
    """X-ray luminosity calculator for binaries."""

    def __init__(self, lxmin: float, lxmax: float, lunit: float = 1.0e33):
        self.lxmin = lxmin
        self.lxmax = lxmax
        self.lunit = lunit
        self.eta = 0.1
        self.lambda_ = 0.5
        self.k = 1.9

    def eddington_luminosity(self, mass: float) -> float:
        return 10.0 ** (np.log10(1.3) + 38 - np.log10(self.lunit)) * mass

    def get_lumx(
        self,
        massp: float,
        masss: float,
        period: float,
        a: float,
        iseed: int | None = None,
        use_weibull: bool = True,
    ) -> float:
        ledd = self.eddington_luminosity(massp)

        if self.lxmin == self.lxmax:
            return self.lxmin

        if iseed is None or iseed >= 0:
            u = np.random.random()
            log_lx = (
                np.log10(self.lxmin) + (np.log10(self.lxmax) - np.log10(self.lxmin)) * u
            )
            lumx = 10.0**log_lx
        else:
            while True:
                xmprob = self.lambda_ * ((self.k - 1.0) / self.k) ** (1.0 / self.k)
                u = np.random.random()
                xprob = self.lambda_ * (-np.log(u)) ** (1.0 / self.k)

                get_lumx = (
                    (np.log10(self.lxmax) - np.log10(self.lxmin)) * xprob / xmprob
                )
                get_lumx = np.log10(self.lxmin) + get_lumx
                lumx = 10.0**get_lumx

                if lumx <= ledd:
                    break

        return lumx
