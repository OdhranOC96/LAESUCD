import os,sys 
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from FLUIDS.compound import Compound


class Water(Compound):
    tCritical = 373.946  # static constexpr double tCritical = 373.946; // K
    pCritical = 22063912.8  # static constexpr double pCritical = 22063912.8; // Pa
    acentricFactor = 0.343  # static constexpr double acentricFactor = 0.343;
    ch3 = 0  # static constexpr double ch3 = 0;
    oh = 0  # static constexpr double oh = 0;
    ch2 = 0  # static constexpr double ch2 = 0;
    ch = 0  # static constexpr double ch = 0;
    O = 1  # static constexpr double O = 1;
    antoineA = 3.55959  # static constexpr double antoineA = 3.55959;
    antoineB = 643.748  # static constexpr double antoineB = 643.748;
    antoineC = -198.043  # static constexpr double antoineC = -198.043;
    molarMass = 44.1  # static constexpr double molarMass = 44.1; // g/mol

    def getTCritical(self) -> float:  # //Get critical temperature
        return Water.tCritical

    def getPCritical(self) -> float:  # //get critical pressure
        return Water.pCritical

    def getAcentricFactor(self) -> float:  # //get acentric factor
        return Water.acentricFactor

    def getCH3(self) -> float:  # //gets #CH3 groups in compound
        return Water.ch3

    def getOH(self) -> float:  # //gets #OH groups in compound
        return Water.oh

    def getCH2(self) -> float:  # //gets #CH2 groups in compound
        return Water.ch2

    def getCH(self) -> float:  # //gets #CH groups in compound
        return Water.ch

    def getO(self) -> float:  # //gets #O (non-ring) groups in compound
        return Water.O

    def getAntoineA(self) -> float:  # //gets Antoine constant A
        return Water.antoineA

    def getAntoineB(self) -> float:  # //gets Antoine constant B
        return Water.antoineB

    def getAntoineC(self) -> float:  # //gets Antoine constant C
        return Water.antoineC

    def getMolarMass(self) -> float:  # //gets the molar mass
        return Water.molarMass
