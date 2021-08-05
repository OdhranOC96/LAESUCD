import os,sys 
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from FLUIDS.air import Air
from FLUIDS.compound import Compound, PRIsentropic, PRH, PRS
from FLUIDS.conversion import kjkgtojmol


class IsentropicUnit:
    # //These are the definitions of two Constructors of IsentropicUnit
    # // mAir(air), mCompound(nullptr) are known as initialization lists.
    def __init__(self, air: Air = None, compound: Compound = None):
        self._mAir = air
        self._mCompound = compound

    # // Functions for the base class and its derived classes
    def calcOutTemperature(self, inTemp: float, inPres: float,
                           outPres: float) -> float:  # //calculates outlet Temperature
        s1 = self.calcInEntropy(inTemp, inPres)
        if self._mAir is not None:
            return self._mAir.isentropic(s1, outPres)
        else:
            value = PRIsentropic(self._mCompound, s1, outPres)
            return value

    def calcInEnthalpy(self, inTemp: float, inPres: float) -> float:  # //calculates inlet enthalpy
        if self._mAir is not None:
            return self._mAir.enthalpy(inTemp, inPres)
        else:
            value = PRH(self._mCompound, inTemp, inPres)
            return value

    def calcOutEnthalpy(self, inTemp: float, inPres: float, outPres: float) -> float:  # //calculates outlet enthalpy
        outTemp = self.calcOutTemperature(inTemp, inPres, outPres)
        return self.calcInEnthalpy(outTemp, outPres)

    def calcInEntropy(self, inTemp: float, inPres: float) -> float:  # //calculates inlet entropy
        if self._mAir:
            return self._mAir.entropy(inTemp, inPres)
        else:
            value = PRS(self._mCompound, inTemp, inPres)
            return value

    def calcOutEntropy(self, inTemp: float, inPres: float, outPres: float) -> float:  # //calculates outlet entropy
        outTemp = self.calcOutTemperature(inTemp, inPres, outPres)
        return self.calcInEntropy(outTemp, outPres)

    def calcSpecificWork(self, inTemp: float, inPres: float, outPres: float) -> float:  # //calculates specific work
        inEnthalpy = self.calcInEnthalpy(inTemp, inPres)
        outTemp = self.calcOutTemperature(inTemp, inPres, outPres)
        outEnthalpy = self.calcOutEnthalpy(inTemp, inPres, outPres)
        return outEnthalpy - inEnthalpy

    def calcPowerOutIn(self, inTemp: float, inPres: float, outPres: float,
                       flow: float) -> float:  # //calculates power requirement
        work = self.calcSpecificWork(inTemp, inPres, outPres)
        if self._mAir:
            work = kjkgtojmol(work, self._mAir.getMolarMass())
            return work * flow
        else:
            work = kjkgtojmol(work, self._mCompound.getMolarMass())
            return work * flow
