import pyneb as pn
import numpy as np
import matplotlib.pyplot as plt


def get_OoH_2(error = None):
    
    obs = pn.Observation('NGC300.dat', corrected=True, errIsRelative=False)
    if error == '-':
        for i in np.arange(obs.n_lines):
            obs.lines[i].corrIntens *= 1. - obs.lines[i].corrError
    if error == '+':
        for i in np.arange(obs.n_lines):
            obs.lines[i].corrIntens *= 1. + obs.lines[i].corrError

    O3 = pn.Atom('O', 3)
    r_O3 = O3.getEmissivity(1e4, 1e3, wave=4959)/O3.getEmissivity(1e4, 1e3, wave=5007)
    
    I = lambda label: obs.getLine(label=label).corrIntens
    R23 = (I('O2_3727A+') + I('O3_5007A') * (1 + r_O3)) / 100.
    x = np.log10(R23)
    OsH = 9.265 - 0.33 * x - 0.202 * x**2 - 0.207 * x**3 - 0.333 * x**4
    return OsH

def get_OoH_3():
    
    obs = pn.Observation('NGC300.dat', corrected=True, errIsRelative=False)

    I = lambda label: obs.getLine(label=label).corrIntens
    O3N2 = np.log10((I('O3_5007A') / 100.) / (I('N2_6584A') / (286.)))
    
    OsH = 8.73 - 0.32 * O3N2
    return OsH

def get_OoH_4(Tlow, Thigh, Ne):
    
    obs = pn.Observation('NGC300.dat', corrected=True, errIsRelative=False)
    O3 = pn.Atom('O', 3)
    r_O3 = O3.getEmissivity(1e4, 1e3, wave=4959)/O3.getEmissivity(1e4, 1e3, wave=5007)
    t2 = Tlow / 1e4
    t3 = Thigh / 1e4
    
    I = lambda label: obs.getLine(label=label).corrIntens
    R2 = I('O2_3727A+') / 100
    R3 = I('O3_5007A') * (1 + r_O3) / 100.
    R = I('O3_4363A') / 100
    R23 = R2 + R3
    X23 = np.log10(R23)
    x2 = 1e-4 * Ne * t2**(-0.5)
    Opp = np.log10(R3) + 6.174 + 1.251 / t3 - 0.55 * np.log10(t3)
    Op =  np.log10(R2) + 5.890 + 1.676 / t2 - 0.40 * np.log10(t2) + np.log10(1+1.35*x2)    
    
    OsH = np.log10(10**Op + 10**Opp)
    return OsH

