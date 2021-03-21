from air import Air
from compound import PRH
from isentropic_unit import IsentropicUnit
from methanol import Methanol
from water import Water
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt 



def main():
    # testing for air.py
    air = Air()
    print("enthalpy	298	101325\n", air.enthalpy(298, 101325))
    print("entropy	298	101325\n", air.entropy(298, 101325))
    print("density	298	101325\n", air.density(298, 101325))
    print("Cp	298	101325\n", air.cp(298, 101325))

    # testing for IsentropicUnit
    isentropic_unit = IsentropicUnit(air=air)
    print("calcOutTemperature	298	101325	500000\n", isentropic_unit.calcOutTemperature(298, 101325, 500000))
    print("WATER PRH 298	101325\n", PRH(Water(), 298, 101325))
    print("METHANOL PRH 298	101325\n", PRH(Methanol(), 298, 101325))

    x = np.arange(80, 293,1)
    y = [air.cv(_,101325) for _ in x]
    plt.plot(x,y)
    plt.show()
        
    


if __name__ == '__main__':
    main()
