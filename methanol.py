from compound import Compound


class Methanol(Compound):
    tCritical = 512.6  # static constexpr double tCritical = 512.6; // K
    pCritical = 8101350  # static constexpr double pCritical = 8101350; // Pa
    acentricFactor = 0.5625  # static constexpr double acentricFactor = 0.5625;
    ch3 = 1  # static constexpr double ch3 = 1;
    oh = 1  # static constexpr double oh = 1;
    ch2 = 0  # static constexpr double ch2 = 0;
    ch = 0  # static constexpr double ch = 0;
    O = 0  # static constexpr double O = 0;
    antoineA = 5.15853  # static constexpr double antoineA = 5.15853;
    antoineB = 1569.613  # static constexpr double antoineB = 1569.613;
    antoineC = -34.846  # static constexpr double antoineC = -34.846;
    molarMass = 32.04  # static constexpr double molarMass = 32.04; // g/mol

    def getTCritical(self) -> float:  # //Get critical temperature
        return Methanol.tCritical

    def getPCritical(self) -> float:  # //get critical pressure
        return Methanol.pCritical

    def getAcentricFactor(self) -> float:  # //get acentric factor
        return Methanol.acentricFactor

    def getCH3(self) -> float:  # //gets #CH3 groups in compound
        return Methanol.ch3

    def getOH(self) -> float:  # //gets #OH groups in compound
        return Methanol.oh

    def getCH2(self) -> float:  # //gets #CH2 groups in compound
        return Methanol.ch2

    def getCH(self) -> float:  # //gets #CH groups in compound
        return Methanol.ch

    def getO(self) -> float:  # //gets #O (non-ring) groups in compound
        return Methanol.O

    def getAntoineA(self) -> float:  # //gets Antoine constant A
        return Methanol.antoineA

    def getAntoineB(self) -> float:  # //gets Antoine constant B
        return Methanol.antoineB

    def getAntoineC(self) -> float:  # //gets Antoine constant C
        return Methanol.antoineC

    def getMolarMass(self) -> float:  # //gets the molar mass
        return Methanol.molarMass
