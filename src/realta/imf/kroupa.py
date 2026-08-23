import numpy as np

from realta.imf.base import IMF


class KroupaIMF(IMF):
    """Broken power-law IMF (Kroupa 2001, MNRAS 322, 231).

    Four segments in dN/dm ~ m^-beta, breaking at m1=0.08, m2=0.5,
    m3=1.0 Msun with beta0=0.3, beta1=1.3, beta2=beta3=2.3 -- the
    canonical Kroupa (2001) values. Since beta2 == beta3, m3=1.0 is not
    a physically distinct break, but is kept as a separate (numerically
    inert) segment for consistency with the original model.
    config.imf_type=2 selects this IMF (see realta.imf.factory.get_imf)
    and is Realta's default.

    cdf() normalizes on the caller-supplied [mmin, mmax], i.e. `mmin` is
    the true lower bound of the first segment (see `breaks` in cdf()
    below) -- it is a free parameter of the population being simulated,
    not a fixed property of the IMF itself. Realta's own config.yml
    default is mmin=0.1 Msun, the practical stellar lower-mass cutoff
    (below ~0.08 Msun objects are substellar, not stars) -- this
    excludes the beta0=0.3 segment (0.01-0.08 Msun in the canonical
    form) from being sampled at all by default, since mmin=0.1 > m1.
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
        """P(<m) for the broken power law, normalized on [mmin, mmax]."""
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
