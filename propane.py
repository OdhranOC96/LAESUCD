from compound import Compound


class Propane(Compound):
    tCritical = 369.9  # static constexpr double tCritical = 369.9; // K
    pCritical = 4301000  # static constexpr double pCritical = 4301000; // Pa
    acentricFactor = 0.152  # static constexpr double acentricFactor = 0.152;
    ch3 = 2  # static constexpr double ch3 = 2;
    oh = 0  # static constexpr double oh = 0;
    ch2 = 1  # static constexpr double ch2 = 1;
    ch = 0  # static constexpr double ch = 0;
    O = 0  # static constexpr double O = 0;
    antoineA = 3.98292  # static constexpr double antoineA = 3.98292;
    antoineB = 819.296  # static constexpr double antoineB = 819.296;
    antoineC = -24.417  # static constexpr double antoineC = -24.417;
    molarMass = 44.1  # static constexpr double molarMass = 44.1; // g/mol

    def getTCritical(self) -> float:  # //Get critical temperature
        return Propane.tCritical

    def getPCritical(self) -> float:  # //get critical pressure
        return Propane.pCritical

    def getAcentricFactor(self) -> float:  # //get acentric factor
        return Propane.acentricFactor

    def getCH3(self) -> float:  # //gets #CH3 groups in compound
        return Propane.ch3

    def getOH(self) -> float:  # //gets #OH groups in compound
        return Propane.oh

    def getCH2(self) -> float:  # //gets #CH2 groups in compound
        return Propane.ch2

    def getCH(self) -> float:  # //gets #CH groups in compound
        return Propane.ch

    def getO(self) -> float:  # //gets #O (non-ring) groups in compound
        return Propane.O

    def getAntoineA(self) -> float:  # //gets Antoine constant A
        return Propane.antoineA

    def getAntoineB(self) -> float:  # //gets Antoine constant B
        return Propane.antoineB

    def getAntoineC(self) -> float:  # //gets Antoine constant C
        return Propane.antoineC

    def getMolarMass(self) -> float:  # //gets the molar mass
        return Propane.molarMass
