from fluid import Fluid
import math

import constants
from cubicsolver import cubicsolver


class Compound(Fluid):
    def getTCritical(self) -> float:  # //Get critical temperature
        return 0

    def getPCritical(self) -> float:  # //get critical pressure
        return 0

    def getAcentricFactor(self) -> float:  # //get acentric factor
        return 0

    def getCH3(self) -> float:  # //gets #CH3 groups in compound
        return 0

    def getOH(self) -> float:  # //gets #OH groups in compound
        return 0

    def getCH2(self) -> float:  # //gets #CH2 groups in compound
        return 0

    def getCH(self) -> float:  # //gets #CH groups in compound
        return 0

    def getO(self) -> float:  # //gets #O (non-ring) groups in compound
        return 0

    def getAntoineA(self) -> float:  # //gets Antoine constant A
        return 0

    def getAntoineB(self) -> float:  # //gets Antoine constant B
        return 0

    def getAntoineC(self) -> float:  # //gets Antoine constant C
        return 0

    def getMolarMass(self) -> float:  # //gets the molar mass
        return 0

    def cp(self, temp, pres):
        # //override: Redefining the behaviour of function present in the base class
        return PRCp(self, temp, pres)

    def enthalpy(self, temp, pres):
        # //override: Redefining the behaviour of function present in the base class
        return PRH(self, temp, pres)


def CpIG(comp: Compound, temp: float) -> list:
    """
    //The parameter ‘comp’ is a reference to an instance of the class ‘Compound’.
    //‘const’ means this function cannot modify that instance
    """
    CpA = (comp.getCH3() * constants.CH3[0]) + (comp.getCH2() * constants.CH2[0]) + (
            comp.getOH() * constants.OH[0]) - 37.93
    CpB = (comp.getCH3() * constants.CH3[1]) + (comp.getCH2() * constants.CH2[1]) + (
            comp.getOH() * constants.OH[1]) + 0.21
    CpC = (comp.getCH3() * constants.CH3[2]) + (comp.getCH2() * constants.CH2[2]) + (comp.getOH() * constants.OH[2]) - (
            3.91 * 1E-4)
    CpD = (comp.getCH3() * constants.CH3[3]) + (comp.getCH2() * constants.CH2[2]) + (comp.getOH() * constants.OH[3]) + (
            2.06 * 1E-7)
    res = CpA + (CpB * temp) + (CpC * math.pow(temp, 2)) + (CpD * math.pow(temp, 3))
    return [res, CpA, CpB, CpC, CpD]


def HIG(comp: Compound, temp: float) -> float:
    # // 25C reference temperature
    tempRef = 298.15

    cpIdealGas = CpIG(comp, temp)
    # // cpIdealGas[0] is not used
    term1 = 1 * cpIdealGas[1] * (temp - tempRef)
    term2 = 0.5 * cpIdealGas[2] * (math.pow(temp, 2) - math.pow(tempRef, 2))
    term3 = 0.33 * cpIdealGas[3] * (math.pow(temp, 3) - math.pow(tempRef, 3))
    term4 = 0.25 * cpIdealGas[4] * (math.pow(temp, 4) - math.pow(tempRef, 4))
    return term1 + term2 + term3 + term4


def SIG(comp: Compound, temp: float, pres: float) -> float:
    pres = pres / 1000000
    tempRef = 298.15
    # // In MPa, the reference pressure
    presRef = 0.101325

    cpIdealGas = CpIG(comp, temp)
    # // cpIdealGas[0] is not used
    term1 = 1 * cpIdealGas[1] * math.log(temp / tempRef)
    term2 = 1 * cpIdealGas[2] * (temp - tempRef)
    term3 = 0.5 * cpIdealGas[3] * (math.pow(temp, 2) - math.pow(tempRef, 2))
    term4 = 0.33 * cpIdealGas[4] * (math.pow(temp, 3) - math.pow(tempRef, 3))
    term5 = -1 * constants.R * math.log(pres / presRef)
    return term1 + term2 + term3 + term4 + term5


def AntoineTSat(comp: Compound, pres: float) -> float:
    pres = pres / 100000
    res = comp.getAntoineB() / (comp.getAntoineA() - math.log10(pres)) - comp.getAntoineC()
    return res


def find_cubic_roots(comp: Compound, temp: float, pres: float, a: float, b: float, c: float) -> list:
    tempSat = AntoineTSat(comp, pres)
    total_root, x0, x1, x2 = cubicsolver(a, b, c)
    if total_root == 1:
        # // Find one real root
        if temp < tempSat:
            return [x0, 0, 0]
        else:
            return [0, 0, x0]
    else:
        # // Find three real roots
        if temp < tempSat:
            return [x0, 0, 0]
        elif temp > tempSat:
            return [0, 0, x2]
        else:
            """// we'll sometimes get 3 roots even if its unsaturated
            // so this needs to be adjusted to be a saturation check
            // i.e. not trusting the function"""
            return [x0, x1, x2]


def PRLiquidDFH(comp: Compound, temp: float, roots: list, coefA: float, coefB: float, coefK: float,
                alpha: float) -> float:
    tempCrit = comp.getTCritical()
    tempReduced = temp / tempCrit

    # // Relation for the liquid enthalpy departure function terms (H-HIG)
    term1 = constants.R * tempCrit
    term2 = tempReduced * (roots[0] - 1)
    term3 = 2.078 * (1 + coefK) * math.pow(alpha, 0.5)
    term4 = roots[0] + (2.414 * coefB)
    term5 = roots[0] - (0.414 * coefB)
    term6 = math.log(term4 / term5)

    return term1 * (term2 - term3 * term6)


def PRVapourDFH(comp: Compound, temp: float, roots: list, coefA: float, coefB: float, coefK: float,
                alpha: float) -> float:
    tempCrit = comp.getTCritical()
    tempReduced = temp / tempCrit

    # // Relation for the vapour enthalpy departure function terms (H-HIG)
    term1 = constants.R * tempCrit
    term2 = tempReduced * (roots[2] - 1)
    term3 = 2.078 * (1 + coefK) * math.pow(alpha, 0.5)
    term4 = math.fabs(roots[2] + (2.414 * coefB))
    term5 = math.fabs(roots[2] - (0.414 * coefB))
    term6 = math.log(term4 / term5)

    return term1 * (term2 - term3 * term6)


def PRH(comp: Compound, temp: float, pres: float) -> float:
    tempCrit = comp.getTCritical()
    presCrit = comp.getPCritical()
    acentric = comp.getAcentricFactor()

    suba = 0.45724 * math.pow(constants.R, 2) * math.pow(tempCrit, 2) / presCrit
    subb = 0.07780 * constants.R * tempCrit / presCrit
    coefK = 0.37464 + (1.54226 * acentric) - (0.26992 * math.pow(acentric, 2))
    tempReduced = temp / tempCrit
    alpha = math.pow(1 + coefK * (1 - math.pow(tempReduced, 0.5)), 2)

    coefA = (alpha * suba * pres) / (math.pow(constants.R, 2) * math.pow(temp, 2))
    coefB = (subb * pres) / (constants.R * temp)

    m1 = 1
    m2 = -1 * (1 - coefB)
    m3 = coefA - (2 * coefB) - (3 * math.pow(coefB, 2))
    m4 = -1 * ((coefA * coefB) - math.pow(coefB, 2) - math.pow(coefB, 3))

    roots = find_cubic_roots(comp, temp, pres, m2, m3, m4)

    liquidHDF = PRLiquidDFH(comp, temp, roots, coefA, coefB, coefK, alpha)
    vapourHDF = PRVapourDFH(comp, temp, roots, coefA, coefB, coefK, alpha)

    idealGasH = HIG(comp, temp)

    liquidH = liquidHDF + idealGasH
    vapourH = vapourHDF + idealGasH
    # // std::cout << "Liquid HDF: " << liquidHDF << std::endl
    # // std::cout << "Vapour HDF: " << vapourHDF << std::endl
    # // std::cout << "Liquid H: " << liquidH << std::endl
    # // std::cout << "Vapour H: " << vapourH << std::endl

    tempSat = AntoineTSat(comp, pres)

    if temp > tempSat:
        return vapourH
    else:
        return liquidH


def PRLiquidDFS(comp: Compound, temp: float, roots: list, coefA: float, coefB: float, coefK: float,
                alpha: float) -> float:
    tempCrit = comp.getTCritical()
    tempReduced = temp / tempCrit

    # // Relation for the liquid entropy departure function terms (S-SIG)
    term1 = math.log(roots[0] - coefB)
    term2 = 2.078 * coefK
    term3 = ((1 + coefK) / math.pow(tempReduced, 0.5)) - coefK
    term4 = roots[0] + (2.414 * coefB)
    term5 = roots[0] - (0.414 * coefB)
    term6 = math.log(term4 / term5)
    return constants.R * (term1 - (term2 * term3 * term6))


def PRVapourDFS(comp: Compound, temp: float, roots: list, coefA: float, coefB: float, coefK: float,
                alpha: float) -> float:
    tempCrit = comp.getTCritical()
    tempReduced = temp / tempCrit

    # // Relation for the vapour entropy departure function terms (S-SIG)
    term1 = math.log(roots[2] - coefB)
    term2 = 2.078 * coefK
    term3 = ((1 + coefK) / math.pow(tempReduced, 0.5)) - coefK
    term4 = roots[2] + (2.414 * coefB)
    term5 = roots[2] - (0.414 * coefB)
    term6 = math.log(term4 / term5)

    return constants.R * (term1 - (term2 * term3 * term6))


def PRS(comp: Compound, temp: float, pres: float) -> float:
    tempCrit = comp.getTCritical()
    presCrit = comp.getPCritical()
    acentric = comp.getAcentricFactor()

    suba = 0.45724 * math.pow(constants.R, 2) * math.pow(tempCrit, 2) / presCrit
    subb = 0.07780 * constants.R * tempCrit / presCrit
    coefK = 0.37464 + (1.54226 * acentric) - (0.26992 * math.pow(acentric, 2))
    tempReduced = temp / tempCrit
    alpha = math.pow(1 + coefK * (1 - math.pow(tempReduced, 0.5)), 2)

    coefA = (alpha * suba * pres) / (math.pow(constants.R, 2) * math.pow(temp, 2))
    coefB = (subb * pres) / (constants.R * temp)

    m1 = 1
    m2 = -1 * (1 - coefB)
    m3 = coefA - (2 * coefB) - (3 * math.pow(coefB, 2))
    m4 = -1 * ((coefA * coefB) - math.pow(coefB, 2) - math.pow(coefB, 3))

    roots = find_cubic_roots(comp, temp, pres, m2, m3, m4)

    liquidSDF = PRLiquidDFS(comp, temp, roots, coefA, coefB, coefK, alpha)
    vapourSDF = PRVapourDFS(comp, temp, roots, coefA, coefB, coefK, alpha)

    idealGasS = SIG(comp, temp, pres)

    liquidS = liquidSDF + idealGasS
    vapourS = vapourSDF + idealGasS

    tempSat = AntoineTSat(comp, pres)

    if temp > tempSat:
        return vapourS
    else:
        return liquidS


def PRIsentropic(comp: Compound, s1: float, pres: float) -> float:
    """
    //PRS is already set up to give the right input entropy so we can surely
    //iterate to get T now. Which we couldnt do before because we couldnt get
    //the right S
    """

    # // Temperature initial guess
    tempL = 273
    tempH = 400
    # // Evaluation at lower bound
    value = PRS(comp, tempL, pres)
    sLow = value - s1

    while True:
        tempM = (tempL + tempH) / 2
        values1 = PRS(comp, tempM, pres)
        sDiff = values1 - s1
        if sDiff * sLow > 0:
            tempL = tempM
        else:
            tempH = tempM
        if math.fabs(sDiff) <= 5:
            break
    return tempM


def PRCp(comp: Compound, temp: float, pres: float) -> float:
    delta = 1
    value1 = PRH(comp, temp, pres)
    value2 = PRH(comp, temp + delta, pres)
    return math.fabs(value2 - value1) / delta
