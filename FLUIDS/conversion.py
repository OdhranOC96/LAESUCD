def moldm3tokgm3(density, molarMass):
    return density * molarMass


def jmoltokjkg(enthalpy, molarMass):
    return enthalpy / molarMass


def kjkgtojmol(enthalpy, molarMass):
    return enthalpy * molarMass


def jmolktokjkgk(entropy, molarMass):
    return entropy / molarMass
