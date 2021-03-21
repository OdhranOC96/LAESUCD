from air import Air
from compound import PRH
from isentropic_unit import IsentropicUnit
from methanol import Methanol
from water import Water
import numpy as np 
import matplotlib
import matplotlib.pyplot as plt 

air = Air()
x = np.arange(80, 293,1)
y = [air.cv(_,101325) for _ in x]
plt.plot(x,y)
plt.show()