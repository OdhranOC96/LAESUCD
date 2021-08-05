import os,sys 
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import math

from FLUIDS.fluid import Fluid
from FLUIDS import constants
from FLUIDS import conversion
import numpy as np


class Air(Fluid):
    # // Table 13 Coefficients and exponents for the Equation of State for air
    K = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
    IK = [1, 1, 1, 2, 3, 3, 4, 4, 4, 6, 1, 3, 5, 6, 1, 3, 11, 1, 3]
    LK = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3]
    JK = [0, 0.33, 1.01, 0, 0, 0.15, 0, 0.2, 0.35, 1.35, 1.6, 0.8, 0.95, 1.25, 3.6, 6, 3.25, 3.5, 15]
    NK = [0.118160747229, 0.713116392079, -1.61824192067, 0.0714140178971, -0.0865421396646,
          0.134211176704, 0.0112626704218, -0.0420533228842, 0.0349008431982, 0.000164957183186,
          -0.101365037912, -0.173813690970, -0.0472103183731, -0.0122523554253, -0.146629609713,
          -0.0316055879821, -0.000233594806142, 0.0148287891978, -0.00938782884667]

    NBPP = [0.2260724, -7.080499, 5.700283, -12.44017, 17.81926, -10.81364, 0,
            0]  # ///< Table 10 Bubble-Point Pressure Eq.1
    NDPP = [-0.1567266, -5.539635, 0, 0, 0.7567212, 0, 0, -3.514322]  # ///< Dew-Point Pressure Eq.1
    NBP = [44.3413, -240.073, 285.139, -88.3366, -0.892181]  # ///< Table 10 Bubble-Point density Eq.2 coefficients
    NDP = [-2.0466, -4.7520, -13.259, -47.652]  # ///< Table 10 Dew-Point density Eq.3
    # // Table 12 Coefficients for the Equation of State of Air
    NI = [0.0000000605719400, -0.0000210274769, -0.000158860716, -13.841928076, 17.275266575,
          -0.000195363420, 2.490888032, 0.791309509, 0.212236768, -0.197938904,
          25.36365, 16.90741, 87.31279]

    rhoCritical = 10.4477  # < Critical density (mol/dm3)
    tCritical = 132.6312  # < Max Condentherm Temperature (K)
    pCritical = 3785020  # < Max condentherm pressure (Pa)
    molarMass = 28.9586  # < g/mol
    epsilon = 0.01  # < Tolerance, can be tuned
    epsilonT = 0.1  # < Temperature tolerance
    epsilonH = 1  # < Enthalpy tolerance
    tCriticalO2 = 154.58 #Critical temperature of O2
    tCriticalN2 = 126.2 #Critical temperature of N2

    def gasFraction(self, temp: float, pres: float) -> float:  # //calculates the gas fraction at a given T and P
        if temp > Air.tCritical:
            return 1
        tempBubble = self.Tsb(pres)
        tempDew = self.Tsd(pres)
        x = (temp - tempBubble) / (tempDew - tempBubble)
        return x

    def lH(self, temp: float, pres: float) -> float:  # //saturated liquid enthalpy
        pb = self.FnPb(temp)
        pd = self.FnPd(temp)
        # // out of range
        if pres < pd or pres > pb:
            return 0

        rhobp = self.FnRhobp(temp)
        deltaa = self._delta(rhobp)
        tau = self._taw(temp)
        da0_dtau = self._calcDadTau(tau)
        dar_dtau = self._calcDardTau(deltaa, tau)
        dar_ddel = self._darddel(deltaa, tau)
        res = constants.R * temp * (tau * (da0_dtau + dar_dtau) + (deltaa * dar_ddel) + 1)
        return res

    def vH(self, temp: float, pres: float) -> float:  # //saturated vapour enthalpy
        pb = self.FnPb(temp)
        pd = self.FnPd(temp)
        # // out of range
        if pres < pd or pres > pb:
            return 0

        rhodp = self.FnRhodp(temp)
        deltaa = self._delta(rhodp)
        tau = self._taw(temp)
        da0_dtau = self._calcDadTau(tau)
        dar_dtau = self._calcDardTau(deltaa, tau)
        dar_ddel = self._darddel(deltaa, tau)
        res = constants.R * temp * (tau * (da0_dtau + dar_dtau) + (deltaa * dar_ddel) + 1)
        return res

    def ls(self, temp: float, pres: float) -> float:  # //saturated liquid entropy
        pb = self.FnPb(temp)
        pd = self.FnPd(temp)
        # // out of range
        if pres < pd or pres > pb:
            return 0

        rhobp = self.FnRhobp(temp)
        deltaa = self._delta(rhobp)
        tau = self._taw(temp)
        da0_dtau = self._calcDadTau(tau)
        dar_dtau = self._calcDardTau(deltaa, tau)

        # // a0
        a0 = self._FnAtheta(deltaa, tau)

        # // ar
        ar = self._FnAr(deltaa, tau)

        res = constants.R * (tau * (da0_dtau + dar_dtau) - a0 - ar)
        return res

    def vs(self, temp: float, pres: float) -> float:  # //saturated vapour entropy
        pb = self.FnPb(temp)
        pd = self.FnPd(temp)
        # // out of range
        if pres < pd or pres > pb:
            return 0

        rhodp = self.FnRhodp(temp)
        deltaa = self._delta(rhodp)
        tau = self._taw(temp)
        da0_dtau = self._calcDadTau(tau)
        dar_dtau = self._calcDardTau(deltaa, tau)

        # // a0
        a0 = self._FnAtheta(deltaa, tau)

        # // ar
        ar = self._FnAr(deltaa, tau)

        res = constants.R * (tau * (da0_dtau + dar_dtau) - a0 - ar)
        return res

    def Tsb(self, pres: float) -> float:  # //bubble point temperature
        a = 20
        b = 300
        fa = self._getFaOrFb(a, pres)
        c = self._getFcBp(pres, fa, a, b)
        if c <= 60:
            c = 62
        return c

    def Tsd(self, pres: float) -> float:  # //dew point temperature
        a = 20
        b = 300
        fa = self._getFaOrFb(a, pres)
        c = self._getFcDp(pres, fa, a, b)
        return c

    def FnPb(self, temp: float) -> float:  # //bubble point pressure
        """
        // pb not applicable above critical T.
        // return a very large number required by the saturation check.
        """
        if temp > Air.tCritical:
            return 100000000000
        bp = 0
        for i, val in enumerate(Air.NBPP):
            bp += val * math.pow(1 - (temp / Air.tCritical), int((i + 1) / 2))
        pb = math.exp(bp * self._taw(temp)) * Air.pCritical
        return pb

    def FnPd(self, temp: float) -> float:  # //dew point pressure
        """
         // pd not applicable above critical T.
        // return a very large number required by the saturation check.
        """
        if temp > Air.tCritical:
            return 100000000000
        dp = 0.0
        for i, val in enumerate(Air.NDPP):
            dp += val * math.pow(1 - (temp / Air.tCritical), int((i + 1) / 2))
        pd = math.exp(dp * self._taw(temp)) * Air.pCritical
        return pd

    def FnRhobp(self, temp: float) -> float:  # //bubble point density
        """
        // Above critical temperature, return an approximate value in the valid range
        """
        if temp > Air.tCritical:
            return 12

        theta = self._calcTheta(temp)
        term = Air.NBP[0] * math.pow(theta, 0.65) + Air.NBP[1] * math.pow(theta, 0.85) + \
               Air.NBP[2] * math.pow(theta, 0.95) + Air.NBP[3] * math.pow(theta, 1.1) + \
               Air.NBP[4] * math.log(temp / Air.tCritical) + 1
        return term * Air.rhoCritical

    def FnRhodp(self, temp: float) -> float:  # //dew point density
        """
        // Above critical temperature, return an approximate value in the valid range
        """
        if temp > Air.tCritical:
            return 22

        theta = self._calcTheta(temp)
        term = Air.NDP[0] * math.pow(theta, 0.41) + Air.NDP[1] * theta + Air.NDP[2] * math.pow(theta, 2.8) + Air.NDP[
            3] * math.pow(theta, 6.5)

        return math.exp(term) * Air.rhoCritical

    def saturationcheck(self, temp: float, pres: float) -> float:  # //check for saturation
        pb = self.FnPb(temp)
        pd = self.FnPd(temp)
        if pd <= pres <= pb and temp < Air.tCritical:
            return True
        else:
            return False

    def enthalpyrho(self, rho: float,
                    temp: float) -> float:  # //calculate enthalpy (saturated and unsaturated) with density input
        tau = self._taw(temp)
        deltaa = self._delta(rho)
        da0_dtau = self._calcDadTau(tau)
        dar_dtau = self._calcDardTau(deltaa, tau)
        dar_ddel = self._darddel(deltaa, tau)
        res = (constants.R * temp) * (tau * (da0_dtau + dar_dtau) + (deltaa * dar_ddel) + 1)
        return res

    def density(self, temp: float, pres: float) -> float:  # //calculate density (saturated/unsaturtated)
        if self.saturationcheck(temp, pres):
            x = self.gasFraction(temp, pres)
            rhobp = self.FnRhobp(temp)
            rhodp = self.FnRhodp(temp)
            v = x / rhodp + (1 - x) / rhobp
            return 1 / v

        pb = self.FnPb(temp)
        pd = self.FnPd(temp)
        if pres > pd:
            a = 40
            b = self.FnRhobp(temp) - 10
        else:
            a = self.FnRhodp(temp)
            b = 0.001
        tau = self._taw(temp)

        # // fa
        deltaa = self._delta(a)
        dar_ddela = self._darddel(deltaa, tau)
        fa = (deltaa * dar_ddela) + 1 - (pres / (a * constants.R * temp * 1000))

        while True:
            c = (a + b) / 2
            deltac = self._delta(c)
            dar_ddelc = self._darddel(deltac, tau)
            fc = (deltac * dar_ddelc) + 1 - (pres / (c * constants.R * temp * 1000))
            if fa * fc > 0:
                a = c
            else:
                b = c

            if math.fabs(fc) <= Air.epsilon:
                break
        return c

    def entropy(self, temp: float, pres: float) -> float:  # //calculate entropy (saturated/unsaturated)
        if self.saturationcheck(temp, pres):
            x = self.gasFraction(temp, pres)
            vapourS = self.vs(temp, pres)
            liquidS = self.ls(temp, pres)
            res = x * vapourS + (1 - x) * liquidS
        else:
            rho = self.density(temp, pres)
            res = self.entropyrho(rho, temp)
        return res

    def entropyrho(self, rho: float, temp: float) -> float:  # //calculate entropy (saturated/unsaturated) density
        deltaa = self._delta(rho)
        tau = self._taw(temp)
        da0_dtau = self._calcDadTau(tau)
        dar_dtau = self._calcDardTau(deltaa, tau)

        # // a0
        a0 = self._FnAtheta(deltaa, tau)

        # // ar
        ar = self._FnAr(deltaa, tau)
        res = constants.R * (tau * (da0_dtau + dar_dtau) - a0 - ar)
        return res

    def exergy(self, temp: float, pres: float) -> float:  # //calculate exergy
        t0 = 298
        p0 = 100000
        h0 = self.enthalpy(self.density(t0, p0), t0)
        s0 = self.entropy(self.density(t0, p0), t0)
        h = self.enthalpy(self.density(temp, pres), temp)
        s = self.entropy(self.density(temp, pres), temp)

        # // j / mol
        res = ((h - h0) - t0 * (s - s0))
        return conversion.jmolktokjkgk(res, self.molarMass)

    def densityliquid(self, temp: float,
                      pres: float) -> float:  # //calculate forced liquid density (saturated/unsaturated)
        a = 40
        b = self.FnRhobp(temp) - 15
        tau = self._taw(temp)
        # // fa
        deltaa = self._delta(a)
        dar_ddela = self._darddel(deltaa, tau)
        fa = (deltaa * dar_ddela) + 1 - (pres / (a * constants.R * temp * 1000))

        while True:
            c = (a + b) / 2
            deltac = self._delta(c)
            dar_ddelc = self._darddel(deltac, tau)
            fc = (deltac * dar_ddelc) + 1 - (pres / (c * constants.R * temp * 1000))
            if fa * fc > 0:
                a = c
            else:
                b = c
            if math.fabs(fc) <= Air.epsilon:
                break
        return c

    def densitygas(self, temp: float,
                   pres: float) -> float:  # //calculate forced liquid density (saturated/unsaturated)
        # // Sometimes this range is too narrow
        a = self.FnRhodp(temp)
        b = 0.001
        tau = self._taw(temp)
        # // fa
        deltaa = self._delta(a)
        dar_ddela = self._darddel(deltaa, tau)
        fa = (deltaa * dar_ddela) + 1 - (pres / (a * constants.R * temp * 1000))

        while True:
            c = (a + b) / 2
            deltac = self._delta(c)
            dar_ddelc = self._darddel(deltac, tau)
            fc = (deltac * dar_ddelc) + 1 - (pres / (c * constants.R * temp * 1000))

            if fa * fc > 0:
                a = c
            else:
                b = c
            if math.fabs(fc) <= Air.epsilon:
                break
        return c

    def isenthalpic(self, h1: float, pres: float) -> float:  #
        tempL = 65
        tempH = 500
        hEval = 12000  # // Set high value (higher than possible root values)
        while True:
            tempM = (tempL + tempH) / 2
            if not self.saturationcheck(tempM, pres):
                h = self.enthalpy(tempM, pres)
                hDiff = (h1) - h
            else:
                liquidFraction = 1 - self.gasFraction(tempM, pres)
                h = liquidFraction * self.lH(tempM, pres) + (1 - liquidFraction) * self.vH(tempM, pres)
                hDiff = (h1) - h

            if hDiff * hEval > 0:
                tempL = tempM
            else:
                tempH = tempM

            if math.fabs(hDiff) <= Air.epsilonH:
                break
        return tempM

    def isentropic(self, s1: float,
                   pres: float) -> float:  # //calculates the outlet temperature of an isentropic engine
        # // Temperature initial guess
        tempL = 60
        tempH = 1700
        # // Evaluation at lower bound
        sLow = self.entropy(tempL, pres) - s1

        while True:
            tempM = (tempL + tempH) / 2
            if not self.saturationcheck(tempM, pres):
                s = self.entropy(tempM, pres)
                sDiff = s - s1
            else:
                x = self.gasFraction(tempM, pres)
                s = x * self.vs(tempM, pres) + (1 - x) * self.ls(tempM, pres)
                sDiff = s - s1

            if sDiff * sLow >= 0:
                tempL = tempM
            else:
                tempH = tempM

            if math.fabs(sDiff) <= Air.epsilonT:
                break
        return tempM

    def cv(self, temp: float, pres: float) -> float:  # //calculates constant volume heat capacity
        tau = self._taw(temp)
        rho = self.density(temp, pres)
        deltaa = self._delta(rho)
        term1 = -constants.R * (math.pow(tau, 2))
        term2 = (self._d2a0dt2(tau)) + (self._d2ardt2(deltaa, tau))
        return term1 * term2

    def getMolarMass(self):  # //calculates the gas fraction at a given T and P
        return Air.molarMass

    def enthalpy(self, temp, pres):
        if self.saturationcheck(temp, pres):
            x = self.gasFraction(temp, pres)
            vapourH = self.vH(temp, pres)
            liquidH = self.lH(temp, pres)
            res = x * vapourH + (1 - x) * liquidH
        else:
            rho = self.density(temp, pres)
            res = self.enthalpyrho(rho, temp)
        return res

    def cp(self, temp, pres):
        tau = self._taw(temp)
        rho = self.density(temp, pres)
        deltaa = self._delta(rho)

        term1 = deltaa * self._darddel(deltaa, tau)
        term2 = -1 * deltaa * tau * self._d2ardtdd(deltaa, tau)
        term3 = 2 * deltaa * self._darddel(deltaa, tau)
        term4 = math.pow(deltaa, 2) * self._d2ardd2(deltaa, tau)
      

        Cv = self.cv(temp, pres)

        Cp = Cv + constants.R * math.pow((1 + term1 + term2), 2) / (1 + term3 + term4)
        return Cp

    # //intermediate calculation functions
    def _delta(self, rho: float) -> float:
        """
        // Converts density to reduced density
        """
        return rho / Air.rhoCritical

    def _taw(self, temp: float) -> float:
        return Air.tCritical / temp

    def _darddel(self, delta: float, tau: float) -> float:
        """
        // FUNCTION: Calculate repetitive derivative term
        // INPUTS: Reduced Density, Reduced Temperature
        """
        term1 = 0.0
        for k in range(10):
            term1 += Air.IK[k] * Air.NK[k] * math.pow(delta, Air.IK[k] - 1) * math.pow(tau, Air.JK[k])
        term2 = 0.0
        for k in range(10, 19, 1):
            term2 += Air.NK[k] * math.pow(delta, Air.IK[k] - 1) * math.pow(tau, Air.JK[k]) * math.exp(
                -1 * math.pow(delta, Air.LK[k])) * (Air.IK[k] - Air.LK[k] * math.pow(delta, Air.LK[k]))
        return term1 + term2

    def _calcDadTau(self, tau: float) -> float:
        """
        // FUNCTION: Calculate repetitive derivative term
        // INPUTS: Reduced Temperature
        """
        term1 = 0
        for i in range(5):
            term1 += (i - 3) * Air.NI[i] * math.pow(tau, i - 4)
        term2 = 1.5 * Air.NI[5] * math.pow(tau, 0.5)
        term3 = Air.NI[6] / tau
        term4 = (Air.NI[7] * Air.NI[10]) / (math.exp(Air.NI[10] * tau) - 1)
        term5 = (Air.NI[8] * Air.NI[11]) / (math.exp(Air.NI[11] * tau) - 1)
        term6 = (Air.NI[9] * Air.NI[12]) / ((2 / 3) * math.exp(-1 * Air.NI[12] * tau) + 1)
        return term1 + term2 + term3 + term4 + term5 + term6

    def _calcDardTau(self, delta: float, tau: float) -> float:
        """
        // FUNCTION: Calculate repetitive derivative term
        // INPUTS: Reduced Density, Reduced Temperature
        """
        term1 = 0.0
        for k in range(10):
            term1 += Air.JK[k] * Air.NK[k] * math.pow(delta, Air.IK[k]) * math.pow(tau, Air.JK[k] - 1)
        term2 = 0.0
        for k in range(10, 19, 1):
            term2 += Air.JK[k] * Air.NK[k] * math.pow(delta, Air.IK[k]) * math.pow(tau, Air.JK[k] - 1) * math.exp(
                -1 * math.pow(delta, Air.LK[k]))
        return term1 + term2

    def _calcTheta(self, temp: float) -> float:
        return 1 - (temp / Air.tCritical)

    def _getFaOrFb(self, aOrb: float, pres: float) -> float:
        """
        // FUNCTION: Calculates the value of the summation term in the bubble point relation equation
        // INPUTS: Temperature, Pressure
        // OUTPUTS: Bubble Point Summation Term
        """
        bp = 0.0
        for i, val in enumerate(Air.NBPP):
            bp += val * math.pow(1 - (aOrb / Air.tCritical), int((i + 1) / 2))
        FaOrFb = (Air.tCritical * bp) / math.log(pres / Air.pCritical) - aOrb
        return FaOrFb

    def _getFcBp(self, pres: float, fa: float, a: float, b: float) -> float:
        """
        // FUNCTION: Calculate summation term in bubble point temperature equation
        // INPUTS: Pressure, fa, a, b (terms calculated within the bubble point temperature function
        // that do not have physical meaning)
        // OUTPUTS: Summation term
        """
        while True:
            c = (a + b) / 2
            fc = self._getFaOrFb(c, pres)
            if fa * fc > 0:
                a = c
            else:
                b = c
            if math.fabs(fc) <= Air.epsilon:
                break
        return c

    def _getFcDp(self, pres: float, fa: float, a: float, b: float) -> float:
        """
        // FUNCTION: Calculate summation term in dewbubble point temperature equation
        // INPUTS: Pressure, fa, a, b (terms calculated within the dew point temperature function
        // that do not have physical meaning)
        // OUTPUTS: Summation term
        """
        while True:
            c = (a + b) / 2
            dp = 0.0
            for k, val in enumerate(Air.NDPP):
                dp += val * math.pow(1 - (c / Air.tCritical), int((k + 1) / 2))
            fc = (Air.tCritical * dp) / math.log(pres / Air.pCritical) - c
            if fa * fc > 0:
                a = c
            else:
                b = c
            if math.fabs(fc) <= Air.epsilon:
                break
        return c

    def _FnAr(self, delta: float, tau: float) -> float:
        """
        // FUNCTION: Calculate repetitive derivative term
        // INPUTS: Reduced Density, Reduced Temperature
        // OUTPUTS: Derivative Term
        """
        term1 = 0
        for k in range(10):
            term1 += Air.NK[k] * math.pow(delta, Air.IK[k]) * math.pow(tau, Air.JK[k])
        term2 = 0
        for k in range(10, 19, 1):
            term2 += Air.NK[k] * math.pow(delta, Air.IK[k]) * math.pow(tau, Air.JK[k]) * math.exp(
                -1 * math.pow(delta, Air.LK[k]))
        return term1 + term2

    def _FnAtheta(self, delta: float, tau: float) -> float:
        """
        // FUNCTION: Calculate repetitive derivative term
        // INPUTS: Reduced Density, Reduced Temperature
        // OUTPUTS: Derivative Term
        """
        term1 = math.log(delta)
        term2 = 0;
        for i in range(5):
            term2 += Air.NI[i] * math.pow(tau, i - 3)

        term3 = Air.NI[5] * math.pow(tau, 1.5)
        term4 = Air.NI[6] * math.log(tau)
        term5 = Air.NI[7] * math.log(1 - math.exp(-1 * Air.NI[10] * tau))
        term6 = Air.NI[8] * math.log(1 - math.exp(-1 * Air.NI[11] * tau))
        term7 = Air.NI[9] * math.log((2 / 3) + math.exp(Air.NI[12] * tau))
        return term1 + term2 + term3 + term4 + term5 + term6 + term7

    def _d2a0dt2(self, tau: float) -> float:
        """
        // FUNCTION: Calculate repetitive derivative term
        // INPUTS: Reduced Temperature
        // OUTPUTS: Derivative Term
        """
        term1 = 0

        for i in range(1,5):
            term1 += (i - 4) * (i - 5) * Air.NI[i - 1] * math.pow(tau, i - 6)

        term2 = 0.75 * Air.NI[5] * math.pow(tau, -0.5)
        term3 = -1 * Air.NI[6] / math.pow(tau, 2)
        term4 = (-1 * Air.NI[7] * math.pow(Air.NI[10], 2) * math.exp(Air.NI[10] * tau)) / math.pow(
            (math.exp(Air.NI[10] * tau) - 1), 2)
        term5 = (-1 * Air.NI[8] * math.pow(Air.NI[11], 2) * math.exp(Air.NI[11] * tau)) / math.pow(
            (math.exp(Air.NI[11] * tau) - 1), 2)
        term6 = ((2 / 3) * Air.NI[9] * math.pow(Air.NI[12], 2) * (math.exp(-1 * Air.NI[12] * tau))) / (
                (2 / 3) * (math.exp(-1 * Air.NI[12] * tau)) + 1)

        d2a0dt2 = term1 + term2 + term3 + term4 + term5 + term6
        return d2a0dt2

    def _d2ardt2(self, delta: float, tau: float) -> float:
        """
        // FUNCTION: Calculate repetitive derivative term
        // INPUTS: Reduced Density, Reduced Temperature
        // OUTPUTS: Derivative Term
        """
        term1 = 0.0
        term2 = 0.0
        for k in range(1, 11, 1):
            term1 += Air.JK[k - 1] * (Air.JK[k - 1] - 1) * Air.NK[k - 1] * math.pow(delta, Air.IK[k - 1]) * math.pow(
                tau, Air.JK[k - 1] - 2)
        for k in range(11, 19, 1):
            term2 = Air.JK[k - 1] * (Air.JK[k - 1] - 1) * Air.NK[k - 1] * math.pow(delta, Air.IK[k - 1]) * math.pow(tau,
                                                                                                                    Air.JK[
                                                                                                                        k - 1] - 2) * math.exp(
                -1 * math.pow(delta, Air.LK[k - 1]))

        d2ardt2 = (term1 + term2)
        return d2ardt2

    def _d2ardtdd(self, delta: float, tau: float) -> float:
        """
        // FUNCTION: Calculate repetitive derivative term
        // INPUTS: Reduced Density, Reduced Temperature
        // OUTPUTS: Derivative Term
        """
        term1 = 0.0
        term2 = 0.0
        for k in range(1, 11, 1):
            term1 += Air.IK[k - 1] * Air.JK[k - 1] * Air.NK[k - 1] * math.pow(delta, Air.IK[k - 1] - 1) * math.pow(tau,
                                                                                                                   Air.JK[
                                                                                                                       k - 1] - 1)
        for k in range(11, 19, 1):
            term2 += Air.JK[k - 1] * Air.NK[k - 1] * math.pow(delta, Air.IK[k - 1] - 1) * math.pow(tau, Air.JK[
                k - 1] - 1) * math.exp(-1 * math.pow(delta, Air.LK[k - 1])) * (
                             Air.IK[k - 1] - (Air.LK[k - 1]) * math.pow(delta, Air.LK[k - 1]))
        d2ardtdd = (term1 + term2)
        return d2ardtdd

    def _d2ardd2(self, delta: float, tau: float) -> float:
        """
        // FUNCTION: Calculate repetitive derivative term
        // INPUTS: Reduced Density, Reduced Temperature
        // OUTPUTS: Derivative Term
        """
        term1 = 0.0
        term2 = 0.0
        for k in range(1, 11, 1):
            term1 += Air.IK[k - 1] * (Air.IK[k - 1] - 1) * Air.NK[k - 1] * (
                math.pow(delta, (Air.IK[k - 1] - 2))) * math.pow(tau, Air.JK[k - 1])

        for k in range(11, 19, 1):
            term2 += Air.NK[k - 1] * math.pow(delta, Air.IK[k - 1] - 2) * math.pow(tau, Air.JK[k - 1]) * math.exp(
                -1 * math.pow(delta, Air.LK[k - 1])) * (
                             (Air.IK[k - 1] - (Air.LK[k - 1]) * math.pow(delta, Air.LK[k - 1])) * (
                             Air.IK[k - 1] - 1 - ((Air.LK[k - 1]) * math.pow(delta, Air.LK[k - 1]))) - math.pow(
                         (Air.LK[k - 1]), 2) * math.pow(delta, Air.LK[k - 1]))
        return term1 + term2

    def CpO2Liq(self, T):
        #The following properties were calculated using Perry's Handbook (2018) formulas
        CpO2L=(175430.0-6152.3*T+113.92*T**2-0.92382*T**3+0.0027963*T**4)/1000 #J/mol*K oxygen liquid heat capacity
        return CpO2L
    
    def CpN2Liq(self, T):
    #The following properties were calculated using Perry's Handbook (2018) formulas
        CpN2L=(281970.0-12281.0*T+248.0*T**2-2.2182*T**3+0.0074902*T**4)/1000 #J/mol*K nitrogen liquid heat capacity
        return CpN2L

    def CpO2V(self, Tout):
    #The following properties were calculated using Perry's Handbook (2018) formulas
        CpO2Vap=(0.29103*1.0e5+0.10040*1.0e5*(2.5265*1.0e3/Tout/np.sinh(2.5265*1.0e3/Tout))**2+0.09356*1.0e5*(1153.8/Tout/np.cosh(1153.8/Tout))**2)/1000 #J/mol*K oxygen vapor heat capacity
        return CpO2Vap

    def CpN2V(self, Tout):
    #The following properties were calculated using Perry's Handbook (2018) formulas
        CpN2Vap=(0.29105*1.0e5+0.08615*1.0e5*(1.7016*1.0e3/Tout/np.sinh(1.7016*1.0e3/Tout))**2+0.00103*1.0e5*(909.79/Tout/np.cosh(909.79/Tout))**2)/1000 #J/mol*K nitrogen vapor heat capacity  
        return CpN2Vap
    
    def HvapO2(self, TrO2):
    #The following properties were calculated using Perry's Handbook (2018) formulas
        HvaporO2=(0.9008*1e7*(1.0-TrO2)**(0.4542-0.4096*TrO2+0.3183*TrO2**2.0))/1000 #J/mol oxygen heat of vaporization
        return HvaporO2
    
    def HvapN2(self, TrN2):
    #The following properties were calculated using Perry's Handbook (2018) formulas
        HvaporN2=(0.74905*1e7*(1.0-TrN2)**(0.40106-0.317*TrN2+0.27343*TrN2**2.0))/1000 #J/mol nitrogen heat of vaporization
        return HvaporN2
    
    def PsatO2(self, Tout):
        Psat_O2=np.exp(51.245-1200.2/Tout-6.4361*np.log(Tout)+0.028405*Tout)/1e5 #From Perry's Chemical Engineer Handbook (2018)
        return Psat_O2

    def PsatN2(self, Tout):
        Psat_N2=np.exp(58.282-1084.1/Tout-8.3144*np.log(Tout)+0.044127*Tout)/1e5 #From Perry's Chemical Engineer Handbook (2018)
        return Psat_N2
