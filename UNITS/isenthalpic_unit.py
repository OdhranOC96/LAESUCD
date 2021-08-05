import os,sys 
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from numpy import array
from FLUIDS.air import Air
import numpy as np
import scipy.optimize
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import plotly.graph_objects as go


class IsenthalpicUnit:

    def __init__(self, air: Air):
        self._mAir = air

    # // Functions for the base class (and derived classes)

    def ThrottleV(self, Tin: float, Po: float, P: float, L1: float, V1: float, y1O2: float, y1N2: float, x1O2: float, x1N2: float):
        Tref=75.0 #Reference Temperature (K)

        Z=1.0 #moles
        ZO2=0.21*Z #O2 moles (21% O2)
        ZN2=Z-ZO2 #N2 moles (79% N2)

        #Energy Balance in throttle valve
        #THout=120 #(K) initial T hot outlet value (HX), T inlet TV
        #t1=[]
        #while loop
        #while (abs(error2))>0.1:
        #Tin=THout #Temperature inlet of Throttle valve
  
        #TO DELETE ALL COMMENTS ABOVE (DOESNT MAKE SENSE). Or replace it with the right stuff
    
        def objective(X):
            L=X[0] #liquid moles
            xO2=X[1] #oxygen liquid molar fraction
            Tout=X[2] #Outlet T
            yO2=(ZO2-L*xO2)/(Z-L) #oxygen vapor molar fraction
            yN2=1.0-yO2 # nitrgogen vapor molar fraction
            TrO2=Tout/154.58 #Reduced Temperature O2
            TrN2=Tout/126.2 #Reduced temperature N2
            V=Z-L #Total vapour moles
            Psat=yO2*self._mAir.PsatO2(Tin)+yN2*self._mAir.PsatN2(Tin) #Saturation pressure of the mixture
    
            if Psat<=29.816: #?????
                E_B=V1*y1O2*self._mAir.CpO2V(Tin)*(Tin-Tref)+V1*y1N2*self._mAir.CpN2V(Tin)*(Tin-Tref)+L1*x1O2*self._mAir.CpO2Liq(Tin)*(Tin-Tref)+L1*x1N2*self._mAir.CpN2Liq(Tin)*(Tin-Tref)-L*xO2*(self._mAir.HvapO2(TrO2)+self._mAir.CpO2Liq(Tout)*(Tout-Tref))-L*(1-xO2)*(self._mAir.HvapN2(TrN2)+self._mAir.CpN2Liq(Tout)*(Tout-Tref))+60*28.9*V1-15*28.9*L#-(V)*yN2*(CpN2V(Tout)*(Tout-Tref))-(V)*yO2*(CpO2V(Tout)*(Tout-Tref))
            else:
                E_B=ZO2*self._mAir.CpO2V(Tin)*(Tin-Tref)+ZN2*self._mAir.CpN2V(Tin)*(Tin-Tref)-(V)*yN2*(self._mAir.CpN2V(Tout)*(Tout-Tref))-(V)*yO2*(self._mAir.CpO2V(Tout)*(Tout-Tref))
                return E_B
    
        def constraint1(X): #oxygen raoults law
            L=X[0] #liquid moles total
            xO2=X[1] #O2 moles
            Tout=X[2] #outlet valve temperature
            yO2=(ZO2-L*xO2)/(Z-L) #vapour moles O2
            yN2=1-yO2 #vapour moles N2
            return xO2*self._mAir.PsatO2(Tout)-yO2*P #Raoult's law

        def constraint2(X): #nitrogen raoults law
            L=X[0] #liquid moles total 
            xO2=X[1] #liquid moles oxygen
            Tout=X[2] #valve outlet temperature 
            yO2=(ZO2-L*xO2)/(Z-L) #oxygen vapour moles
            yN2=1-yO2 #nitrogen vapour moles
            xN2=1-xO2 #nitrogen liquid moles
            return xN2*self._mAir.PsatN2(Tout)-yN2*P #nitrogen raoults law

        def constraint3(X): #valve energy balance
            L=X[0] #liquid moles
            xO2=X[1] #oxygen liquid molar fraction
            Tout=X[2] #Outlet T
            yO2=(ZO2-L*xO2)/(Z-L) #oxygen vapor molar fraction
            yN2=1.0-yO2 # nitrgogen vapor molar fraction
            TrO2=Tout/154.58 #oxygen reduced temperature 
            TrN2=Tout/126.2 #nitrogen reduced temeprature
            V=Z-L #total vapour mples 
            Psat=yO2*self._mAir.PsatO2(Tin)+yN2*self._mAir.PsatN2(Tin) #mixture saturation pressure
            
            if Psat<=29.816: 
                E_B=V1*y1O2*self._mAir.CpO2V(Tin)*(Tin-Tref)+V1*y1N2*self._mAir.CpN2V(Tin)*(Tin-Tref)+L1*x1O2*self._mAir.CpO2Liq(Tin)*(Tin-Tref)+L1*x1N2*self._mAir.CpN2Liq(Tin)*(Tin-Tref)-L*xO2*(self._mAir.HvapO2(TrO2)+self._mAir.CpO2Liq(Tout)*(Tout-Tref))-L*(1-xO2)*(self._mAir.HvapN2(TrN2)+self._mAir.CpN2Liq(Tout)*(Tout-Tref))+60*28.9*V1-15*28.9*L##-(V)*yN2*(CpN2V(Tout)*(Tout-Tref))-(V)*yO2*(CpO2V(Tout)*(Tout-Tref))
            else:
                E_B=ZO2*self._mAir.CpO2V(Tin)*(Tin-Tref)+ZN2*self._mAir.CpN2V(Tin)*(Tin-Tref)-V*yN2*(self._mAir.CpN2V(Tout)*(Tout-Tref))-V*yO2*(self._mAir.CpO2V(Tout)*(Tout-Tref))
            return E_B
            
        #Initial values
        x0=[0.5,0.2,126.0] #liquid moles, oxygen liquid molar fraction, outlet temperature of valve
        #Boundaries
        bnds=((0,Z),(0,1),(54,154)) #setting up bounds for the max and min possible values
        #constraints
        con1={'type':'eq','fun':constraint1} #function 1, constraint 1
        con2={'type':'eq','fun':constraint2} #function 2, constraint 2
        con3={'type':'eq','fun':constraint3} #function 3, constraint 3
        con=[con1,con2,con3] #array of constraints
        solmin=minimize(objective,x0,method='SLSQP',bounds=bnds,constraints=con,options={'maxiter':10000}) #solving the objective function with the constraints that f1,f2,f3 must also simultaneously be solved (L,xO2,TOut)
        #print(solmin)
        #Results of the model
        L=solmin.x[0] 
        xO2=solmin.x[1]
        Tout=solmin.x[2]
        xN2=1-xO2
        V=Z-L
        TrO2=Tout/154.58
        TrN2=Tout/126.2
        yO2=(ZO2-L*xO2)/(Z-L) #oxygen vapor molar fraction
        yN2=1.0-yO2 # nitrgogen vapor molar fraction
        Psat=yO2*self._mAir.PsatO2(self, Tin)+yN2*self._mAir.PsatN2(self, Tin)
        if Psat<=29.816:
            #Enthalpy of liquid
            Enthalpy_L=L1*x1O2*self._mAir.CpO2Liq(self, Tin)*(Tin-Tref)+L1*x1N2*self._mAir.CpN2Liq(self, Tin)*(Tin-Tref)-L*xO2*(-1*self._mAir.HvapO2(self, TrO2)+self._mAir.CpO2Liq(self, Tout)*(Tout-Tref))-L*(1-xO2)*(-1*self._mAir.HvapN2(self, TrN2)+self._mAir.CpN2Liq(self, Tout)*(Tout-Tref))
            Enthalpy_V=V1*y1O2*self._mAir.CpO2V(self, Tin)*(Tin-Tref)+V1*y1N2*self._mAir.CpN2V(self, Tin)*(Tin-Tref)-(V)*yN2*(self._mAir.CpN2V(self, Tout)*(Tout-Tref))-(V)*yO2*(self._mAir.CpO2V(self, Tout)*(Tout-Tref))#
            
        else:
            
        #Enthalpy of Liquid
            Enthalpy_L=L*xO2*self._mAir.CpO2Liq(self, Tin)*(Tin-Tref)+L*xN2*self._mAir.CpN2Liq(self, Tin)*(Tin-Tref)-15*28.9*L-L*xO2*-1*self._mAir.HvapO2(self, TrO2)-L*xN2*-1*self._mAir.HvapN2(self, TrN2)#J
        #Enthalpy of vapor
            Enthalpy_V=ZO2*self._mAir.CpO2V(self, Tin)*(Tin-Tref)+ZN2*self._mAir.CpN2V(self, Tin)*(Tin-Tref)-V*yN2*(self._mAir.CpN2V(self, Tout)*(Tout-Tref))-V*yO2*(self._mAir.CpO2V(self, Tout)*(Tout-Tref))
        
        
        return [L,xO2,Tout,xN2,Enthalpy_L,Enthalpy_V]


    #THESE ARE ALL WRONG
    def calcOutTemperature(self, inTemp: float, inPres: float,
                           outPres: float) -> float:  # //calculate outlet temperature
        h1 = self.calcInEnthalpy(inTemp, inPres)
        print("h1: ", h1)
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
    
   
