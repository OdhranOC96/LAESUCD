from air import Air


class IsenthalpicUnit:

    def __init__(self, air: Air):
        self._mAir = air

    # // Functions for the base class (and derived classes)
    def calcOutTemperature(self, inTemp: float, inPres: float,
                           outPres: float) -> float:  # //calculate outlet temperature
        h1 = self.calcInEnthalpy(inTemp, inPres)
        return self._mAir.isenthalpic(h1, outPres)

    def calcInEnthalpy(self, inTemp: float, inPres: float) -> float:  # //calculate inlet enthalpy
        return self._mAir.enthalpy(inTemp, inPres)

    def calcGasQuality(self, inTemp: float, inPres: float, outPres: float) -> float:  # //calculate gas fraction
        outTemp = self.calcOutTemperature(inTemp, inPres, outPres)
        return self._mAir.gasFraction(outTemp, outPres)

    def calcLiquidQuality(self, inTemp: float, inPres: float, outPres: float) -> float:  # //calculate liquid quality
        return 1 - self.calcGasQuality(inTemp, inPres, outPres)

    def calcOutVapourEnthalpy(self, inTemp: float, inPres: float,
                              outPres: float) -> float:  # //calculate outlet vapour enthalpy
        outTemp = self.calcOutTemperature(inTemp, inPres, outPres)
        return self._mAir.vH(outTemp, outPres)

    def calcOutLiquidEnthalpy(self, inTemp: float, inPres: float,
                              outPres: float) -> float:  # //calculate outlet liquid enthalpy
        outTemp = self.calcOutTemperature(inTemp, inPres, outPres)
        return self._mAir.lH(outTemp, outPres)
