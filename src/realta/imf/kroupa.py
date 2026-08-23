import numpy as np

from realta.imf.base import IMF


class KroupaIMF(IMF):
    """Broken power-law IMF (Kroupa 2001, MNRAS 322, 231).

    Four segments in dN/dm ~ m^-beta, breaking at m1=0.08, m2=0.5,
    m3=1.0 Msun with beta0=0.3, beta1=1.3, beta2=beta3=2.3 -- the
    canonical Kroupa (2001) values, and an exact match to the reference
    Fortran's kroupa.f (`parameter(beta0=0.3,beta1=1.3,beta2=2.3,
    beta3=2.3)` / `parameter(m0=0.01,m1=0.08,m2=0.5,m3=1)`). Since
    beta2 == beta3, m3=1.0 is not a physically distinct break -- it is
    reproduced here only because the reference Fortran has it as a
    separate (numerically inert) segment. config.imf_type=2 selects
    this IMF (see realta.imf.factory.get_imf) and is Realta's default.

    FLAGGED DISCREPANCY (not fixed -- see class docstring convention in
    the development brief, "flag ambiguity rather than resolve it"):
    kroupa.f hardcodes its lowest break at m0=0.01 and never actually
    uses its own `mmin` argument in the CDF integral -- the low-mass
    leg always integrates from 0.01, regardless of what mmin is passed
    to make_stars(). This Python port instead uses the caller-supplied
    `mmin` as the true lower bound of the first segment (see `breaks`
    in cdf() below). With Realta's own config.yml default (mmin=0.01,
    matching m0), the two are numerically identical. If mmin is ever
    set to anything else, this port's normalization will differ from
    the reference Fortran's -- worth a decision on which behaviour is
    intended before mmin is used as a free parameter in production.
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

    def _integrated_density(
        self, m: float, mmin: float, mmax: float, beta: float
    ) -> float:
        if beta == 1.0:
            return np.log(m) - np.log(mmin)
        else:
            return (m ** (1.0 - beta) - mmin ** (1.0 - beta)) / (1.0 - beta)

    def cdf(self, m: float, mmin: float, mmax: float) -> float:
        """P(<m) for the broken power law, normalized on [mmin, mmax].

        See the class docstring for the flagged mmin-handling
        discrepancy versus the reference kroupa.f.
        """
        if m <= mmin:
            return 0.0
        if m >= mmax:
            return 1.0

        breaks = [mmin, self.m1, self.m2, self.m3, mmax]
        betas = [self.beta0, self.beta1, self.beta2, self.beta3, self.beta3]

        valid_breaks = [b for b in breaks if mmin <= b <= mmax]
        if valid_breaks[0] != mmin:
            valid_breaks.insert(0, mmin)
        if valid_breaks[-1] != mmax:
            valid_breaks.append(mmax)

        valid_betas = []
        for i in range(len(valid_breaks) - 1):
            for j in range(len(breaks) - 1):
                if breaks[j] <= valid_breaks[i] < breaks[j + 1]:
                    valid_betas.append(betas[j])
                    break
            else:
                valid_betas.append(betas[-1])

        dmtot = 0.0
        for i in range(len(valid_breaks) - 1):
            dmtot += self._integrated_density(
                valid_breaks[i + 1], valid_breaks[i], mmax, valid_betas[i]
            )

        cdf_val = 0.0
        for i in range(len(valid_breaks) - 1):
            if m <= valid_breaks[i + 1]:
                cdf_val += self._integrated_density(
                    min(m, valid_breaks[i + 1]),
                    valid_breaks[i],
                    mmax,
                    valid_betas[i],
                )
                break
            else:
                cdf_val += self._integrated_density(
                    valid_breaks[i + 1], valid_breaks[i], mmax, valid_betas[i]
                )

        return cdf_val / dmtot
