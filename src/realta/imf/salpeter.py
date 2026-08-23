from realta.imf.base import IMF


class SalpeterIMF(IMF):
    """Single power-law IMF (Salpeter 1955, ApJ 121, 161): dN/dm ~ m^-alpha.

    Reference Fortran: salpeter.f, `beta=2.35` hardcoded -- alpha=2.35 is
    the default here to match exactly. config.imf_type=1 selects this
    IMF (see realta.imf.factory.get_imf).

    cdf(m, mmin, mmax) = (m^(1-alpha) - mmin^(1-alpha))
                          / (mmax^(1-alpha) - mmin^(1-alpha))
    i.e. a single segment normalized on [mmin, mmax] -- both bounds are
    the caller-supplied config.mmin/config.mmax, matching salpeter.f
    exactly (mmin is used directly there, unlike kroupa.f -- see
    KroupaIMF's docstring for a flagged discrepancy in that IMF only).
    """

    def __init__(self, alpha: float = 2.35):
        self.alpha = alpha
        self.beta = 1.0 - alpha

    def cdf(self, m: float, mmin: float, mmax: float) -> float:
        if m <= mmin:
            return 0.0
        if m >= mmax:
            return 1.0

        numerator = m**self.beta - mmin**self.beta
        denominator = mmax**self.beta - mmin**self.beta
        return numerator / denominator
