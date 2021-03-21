from fluid import Fluid


class HeatExchanger:
    """
    /**
     * Constructs a HeatExchanger with two Fluids, taking ownership of unique pointers
     */
    """

    def __init__(self, f1: Fluid, f2: Fluid):
        self._mFluid1 = f1
        self._mFluid2 = f2

    """
    /**
     * Calculates the mass flow of Fluid 2.
     */
    """

    def calcMassFlow2(self, mf1: float, inPres1: float, inTemp1: float, outTemp1: float, inPres2: float, inTemp2: float,
                      outTemp2: float) -> float:
        enthalpyIn1 = self._mFluid1.enthalpy(inTemp1, inPres1)
        enthalpyOut1 = self._mFluid1.enthalpy(outTemp1, inPres1)
        enthalpyIn2 = self._mFluid2.enthalpy(inTemp2, inPres1)
        enthalpyOut2 = self._mFluid2.enthalpy(outTemp2, inPres1)
        m2 = (mf1 * (enthalpyIn1 - enthalpyOut1)) / (enthalpyOut2 - enthalpyIn2)
        return m2

    """
     /**
     * Calculates the outlet temperature of Fluid 2.
     */
     """

    def calcOutTemperature2(self, mf1: float, inPres1: float, inTemp1: float, outTemp1: float, mf2: float,
                            inPres2: float, inTemp2: float, pinchTemp: float = 5) -> float:
        cp1 = self._mFluid1.cp(inTemp1, inPres1)
        cp2 = self._mFluid2.cp(inTemp2, inPres2)
        # // thermal power
        tp1 = mf1 * cp1
        tp2 = mf2 * cp2

        if tp2 > tp1:
            if inTemp2 < inTemp1:
                outTemp2 = inTemp1 - pinchTemp
            else:
                outTemp2 = inTemp1 + pinchTemp
        else:
            outTemp2 = inTemp2 + (tp1 * (inTemp1 - outTemp1) / tp2)

        return outTemp2

    """
    /**
     * Calculates the outlet temperatures of Fluid 1 and Fluid 2
     */
    """

    def calcOutTemperature(self, mf1: float, inPres1: float, inTemp1: float, mf2: float, inPres2: float, inTemp2: float,
                           pinchTemp: float = 5) -> list:
        cp1 = self._mFluid1.cp(inTemp1, inPres1)
        cp2 = self._mFluid2.cp(inTemp2, inPres2)

        tp1 = mf1 * cp1
        tp2 = mf2 * cp2

        if tp1 < tp2:
            if inTemp1 < inTemp2:
                outTemp1 = inTemp2 - pinchTemp
            else:
                outTemp1 = inTemp2 + pinchTemp
            outTemp2 = inTemp2 + (tp1 * (inTemp1 - outTemp1) / tp2)
        else:
            if inTemp1 < inTemp2:
                outTemp2 = inTemp1 + pinchTemp
            else:
                outTemp2 = inTemp1 - pinchTemp
            outTemp1 = inTemp1 - (tp2 * (outTemp2 - inTemp2) / tp1)
        return [outTemp1, outTemp2]
