import os,sys 
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from FLUIDS.air import Air
from UNITS.isenthalpic_unit import IsenthalpicUnit


class ExpansionValve(IsenthalpicUnit):

    def __init__(self, air: Air):
        super().__init__(air)
        pass
