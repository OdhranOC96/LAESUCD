from UNITS.expansion_valve import ExpansionValve
from FLUIDS.air import Air
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt 



def main():
    
    air = Air()
    ev1 = ExpansionValve(air)

    Tin = 126
    Po = 30
    P = 2
 
    [L,xO2,Tout,xN2,Enthalpy_L,Enthalpy_V] = ev1.ThrottleV(Tin,P,Po, 1,0,0.21,0.79,0,0)
    
    return Tin

if __name__ == '__main__':
    main()
