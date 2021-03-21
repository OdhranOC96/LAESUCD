from air import Air
from compound import Compound
from isentropic_unit import IsentropicUnit


class Turbine(IsentropicUnit):

    def __init__(self, air: Air = None, compound: Compound = None):
        super().__init__(air, compound)
        pass
