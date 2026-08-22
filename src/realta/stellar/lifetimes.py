from realta.stellar.base import LifetimeModel


class TabulatedLifetimeModel(LifetimeModel):

    def __init__(self, Z: float = 0.008):
        pass

    def lifetime(self, mass):
        return 10.0
