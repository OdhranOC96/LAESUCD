import os,sys 
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from FLUIDS.air import Air
from FLUIDS.compound import Compound
from UNITS.isentropic_unit import IsentropicUnit


class Pump(IsentropicUnit):

    def __init__(self, air: Air = None, compound: Compound = None):
        super().__init__(air, compound)
        pass
