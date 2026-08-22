from realta.imf.base import IMF


class SalpeterIMF(IMF):
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
