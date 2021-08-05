import os,sys 
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from FLUIDS.compound import Compound


class PropyleneGlycol(Compound):
    tCritical = 614.38  # static constexpr double tCritical = 614.38; // K
    pCritical = 5791740  # static constexpr double pCritical = 5791740; // Pa
    acentricFactor = 1.107  # static constexpr double acentricFactor = 1.107;
    ch3 = 1  # static constexpr double ch3 = 1;
    oh = 2  # static constexpr double oh = 2;
    ch2 = 1  # static constexpr double ch2 = 1;
    ch = 1  # static constexpr double ch = 1;
    O = 0  # static constexpr double O = 0;
    antoineA = 6.07936  # static constexpr double antoineA = 6.07936;
    antoineB = 2692.187  # static constexpr double antoineB = 2692.187;
    antoineC = -17.94  # static constexpr double antoineC = -17.94;
    molarMass = 76.09  # static constexpr double molarMass = 76.09; // g/mol

    def getTCritical(self) -> float:  # //Get critical temperature
        return PropyleneGlycol.tCritical

    def getPCritical(self) -> float:  # //get critical pressure
        return PropyleneGlycol.pCritical

    def getAcentricFactor(self) -> float:  # //get acentric factor
        return PropyleneGlycol.acentricFactor

    def getCH3(self) -> float:  # //gets #CH3 groups in compound
        return PropyleneGlycol.ch3

    def getOH(self) -> float:  # //gets #OH groups in compound
        return PropyleneGlycol.oh

    def getCH2(self) -> float:  # //gets #CH2 groups in compound
        return PropyleneGlycol.ch2

    def getCH(self) -> float:  # //gets #CH groups in compound
        return PropyleneGlycol.ch

    def getO(self) -> float:  # //gets #O (non-ring) groups in compound
        return PropyleneGlycol.O

    def getAntoineA(self) -> float:  # //gets Antoine constant A
        return PropyleneGlycol.antoineA

    def getAntoineB(self) -> float:  # //gets Antoine constant B
        return PropyleneGlycol.antoineB

    def getAntoineC(self) -> float:  # //gets Antoine constant C
        return PropyleneGlycol.antoineC

    def getMolarMass(self) -> float:  # //gets the molar mass
        return PropyleneGlycol.molarMass
