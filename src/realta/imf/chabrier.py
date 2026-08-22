import numpy as np
from scipy.special import erf

from realta.imf.base import IMF


class ChabrierIMF(IMF):
    def __init__(
        self,
        mc: float = 3.5,
        m_tr: float = 4.0,
        sigma: float = 0.2,
        x: float = 1.7,
    ):
        self.mc = mc
        self.m_tr = m_tr
        self.sigma = sigma
        self.x = x
        self.pi = np.pi

    def cdf(self, m: float, mmin: float, mmax: float) -> float:
        if m <= mmin:
            return 0.0
        if m >= mmax:
            return 1.0

        arg_erf1 = np.log10(self.m_tr / self.mc) / (
            np.log10(10.0) * self.sigma * np.sqrt(2.0)
        )
        arg_erf2 = np.log10(mmin / self.mc) / (
            np.log10(10.0) * self.sigma * np.sqrt(2.0)
        )
        arg_erf = np.log10(m / self.mc) / (np.log10(10.0) * self.sigma * np.sqrt(2.0))

        term1 = (
            self.sigma
            * np.sqrt(self.pi / 2.0)
            * np.log10(10.0)
            * (erf(arg_erf1) - erf(arg_erf2))
        )
        term2 = np.exp(
            -((np.log10(self.m_tr / self.mc)) ** 2) / (2 * self.sigma**2)
        ) / (self.x * self.m_tr ** (-self.x))
        term2 *= mmax ** (-self.x) - self.m_tr ** (-self.x)

        denom = term1 - term2
        A = 1.0 / denom
        B = (
            A
            * np.exp(-(((np.log10(self.m_tr / self.mc) / self.sigma) ** 2.0) / 2.0))
            / self.m_tr ** (-self.x)
        )

        if m <= self.m_tr:
            cdf_val = A * np.sqrt(self.pi / 2.0) * self.sigma * np.log10(10.0)
            cdf_val *= erf(arg_erf) - erf(arg_erf2)
        else:
            cdf_val = A * np.sqrt(self.pi / 2.0) * self.sigma * np.log10(10.0)
            cdf_val *= erf(arg_erf1) - erf(arg_erf2)
            cdf_val -= B * (m ** (-self.x) - self.m_tr ** (-self.x)) / self.x

        return cdf_val
