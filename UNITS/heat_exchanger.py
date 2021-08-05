import os,sys 
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from FLUIDS.fluid import Fluid
from FLUIDS.air import Air


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
        print("Tin: ", enthalpyIn2)
        enthalpyOut2 = self._mFluid2.enthalpy(outTemp2, inPres1)
        print("Tout: ", enthalpyOut2)
        m2 = (mf1 * abs((enthalpyIn1 - enthalpyOut1))) / (abs(enthalpyOut2 - enthalpyIn2))
        return m2

    """
     /**
     * Calculates the outlet temperature of Fluid 2.
     */
     """

    def calcOutTemperature(self, mf1: float, inPres1: float, inTemp1: float, outTemp1: float, mf2: float,
                            inPres2: float, inTemp2: float, pinchTemp: float = 5, sat1in: float = 0, sat2in: float = 0) -> float:
        cp1 = self._mFluid1.cp(inTemp1, inPres1)
        cp2 = self._mFluid2.cp(inTemp2, inPres2)
        # // thermal power
        tp1 = abs(mf1 * cp1)
        tp2 = abs(mf2 * cp2)

        if sat1in == 0:#these sat parameters are used if one of the streams is saturated (the program gets confused otherwise)
            enthalpyIn1 = self._mFluid1.enthalpy(inTemp1, inPres1)
        elif sat1in == 1:
            enthalpyIn1 = self._mFluid1.vH(inTemp1,inPres1)
        elif sat1in == 2: 
            enthalpyIn1 = self._mFluid1.lH(inTemp1,inPres1)
            
            
        
        if sat2in == 0:
            enthalpyIn2 = self._mFluid2.enthalpy(inTemp2, inPres2)
        elif sat2in == 1:
            enthalpyIn2 = self._mFluid2.vH(inTemp2,inPres2)
        elif sat2in == 2:
            enthalpyIn2 = self._mFluid2.lH(inTemp2,inPres2)

        
        enthalpyOut1 = self._mFluid1.enthalpy(outTemp1,inPres1)

        if tp2 < tp1:
            if inTemp2 < inTemp1:
                outTemp2 = inTemp1 - pinchTemp
            else:
                outTemp2 = inTemp1 + pinchTemp
        else: 
            
            enthalpyOut2 = enthalpyIn2 + (mf1 * (enthalpyIn1 - enthalpyOut1) / mf2)
            outTemp2 = self._mFluid2.isenthalpic(enthalpyOut2,inPres2)
        return outTemp2
       

    """
    /**
     * Calculates the outlet temperatures of Fluid 1 and Fluid 2
     */
    """

    def calcOutTemperature2(self, mf1: float, inPres1: float, inTemp1: float, mf2: float, inPres2: float, inTemp2: float,
                           pinchTemp: float = 5, sat1in: float = 0, sat2in: float = 0) -> list:
        cp1 = self._mFluid1.cp(inTemp1, inPres1)
        cp2 = self._mFluid2.cp(inTemp2, inPres2)

        tp1 = mf1 * cp1
        tp2 = mf2 * cp2

        if sat1in == 0:#these sat parameters are used if one of the streams is saturated (the program gets confused otherwise)
                enthalpyIn1 = self._mFluid1.enthalpy(inTemp1, inPres1)
        elif sat1in == 1:
                enthalpyIn1 = self._mFluid1.vH(inTemp1,inPres1)
        elif sat1in == 2: 
                enthalpyIn1 = self._mFluid1.lH(inTemp1,inPres1)
                
                
            
        if sat2in == 0:
                enthalpyIn2 = self._mFluid2.enthalpy(inTemp2, inPres2)
        elif sat2in == 1:
                enthalpyIn2 = self._mFluid2.vH(inTemp2,inPres2)
        elif sat2in == 2:
                enthalpyIn2 = self._mFluid2.lH(inTemp2,inPres2)
         
         #chance this code doesnt work below the same as func above as I didnt check it
        if tp1 < tp2:
            if inTemp1 < inTemp2:
                outTemp1 = inTemp2 - pinchTemp
            else:
                outTemp1 = inTemp2 + pinchTemp
            
            enthalpyOut1 = self._mFluid1.enthalpy(outTemp1, inPres1)
            enthalpyOut2 = (mf1/mf2)*(enthalpyIn1-enthalpyOut1) + enthalpyIn2
            outTemp2 = self._mFluid1.isenthalpic(enthalpyOut2,inPres1)
        else:
            if inTemp1 < inTemp2:
                outTemp2 = inTemp1 + pinchTemp
            else:
                outTemp2 = inTemp1 - pinchTemp
            
            enthalpyOut2 = self._mFluid2.enthalpy(outTemp2, inPres2)
            outEnthalpy1 = enthalpyIn1 - (mf2 * (enthalpyOut2 - enthalpyIn2) / mf1)
            outTemp1 = self._mFluid2.isenthalpic(outEnthalpy1,inPres2)

        return [outTemp1, outTemp2]
