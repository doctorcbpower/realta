from realta.stellar.base import IonisationModel


class TabulatedIonisationModel(IonisationModel):
    def photon_rate(self, mass):
        return 48.0
