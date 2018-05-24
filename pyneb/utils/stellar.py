#/opt/local/bin/python
# Compute stellar Zanstra temperature and stellar luminosity.
# Last update: 2013-Oct-30 by RAS

import numpy as np
import scipy
import scipy.optimize as op
from scipy import constants as phy
from scipy.integrate import quad

WAV_H_IONIZ   = 911.7634
WAV_HE2_IONIZ = WAV_H_IONIZ/4.
ABMAG0     = 48.60
STMAG0     = 21.10
#FLUX_MAG0 = 3.68e-9    # from Pottasch (1984), outdated
FLUX_MAG0  = 3.63e-9    # or pow(10,-STMAG0/2.5)
MBOL_Sun   = 4.75

'''Determine the stellar Zanstra temperature from a nebular H or He emission 
   line flux and the stellar flux density. The method follows 
       Pottasch, S. 1984, "Planetary Nebulae" (Dordrecht: D. Reidel)
   with some modern updates. Also determine Bolometric correction and stellar 
   luminosity given temperature and distance. 

   Limitations: 
    - assumes Blackbody stellar SED
    - no correction for interstellar extinction.
'''

def logL(V, D_pc, T_eff):
    '''Calculate the stellar luminosity in Solar units, given an (unreddened) V 
       apparent magnitude, a distance (pc), and the stellar temperature (K).
    '''
    return -0.4*(V - 2.5*np.log10(D_pc) + BC(T_eff) - MBOL_Sun)


def BC(T_eff):
    '''Bolometric correction for very hot stars, from Vacca, Garmany & Shull,
       1996, ApJ, 460, 914.
    '''
    return 27.66 - 6.84*np.log10(T_eff)
#    return 28.46 - 7.08*np.log10(T_eff) + 0.08*log_g  # Fit incl. gravity

def integrand(x):
    return x**2/(np.exp(x)-1)

def G(ion,T):
    '''Compute integral over frequency for Planck function.
       ion - Reference Ion ("H | HE')
       T   - Stellar temperature (K)
       See Pottasch (1984), p. 169, eq. VII-8 for details. 
    '''

    if ion.upper()=='H':
        nu = phy.c / (WAV_H_IONIZ*1.e-10)
    elif ion.upper()=='HE2':
        nu = phy.c / (WAV_HE2_IONIZ*1.e-10)
    else:
        print('Unrecognized ion: %s' % ion)
        return None

    return quad(integrand, phy.h*nu/(phy.k*T), np.Inf)[0]


class Zanstra():
    def __init__(self, ion, logFlux=-10., magStar=10.):
        s = ion.upper()
        if s in ('H','HE2'):
            self.ion = s
        else:
            self.ion = 'H'
        self.logFlux = logFlux
        self.magStar = magStar

        if self.ion == 'H':
            self.ionScale = 3.95e-11
        elif self.ion == 'HE2':
            self.ionScale = 8.49e-11


    def z_ratio(self, T):
        '''Compute the ratio of the flux in a nebular emission line to the stellar 
           flux density for a given stellar (black-body) temperature T.
        '''

        return self.ionScale * T**3 * G(self.ion, T) * (np.exp(2.665e+4/T) - 1)


    def Func(self, T):
        '''Return the discriminant of the stellar Zanstra temperature function.
            '''

        r_obs = 10**self.logFlux / (FLUX_MAG0 * 10**(-self.magStar/2.5))
        return r_obs - self.z_ratio(T)


    def Solve(self):
        '''Determine the stellar Zanstra temperature from a nebular H or He 
           emission line flux and the stellar flux density.
        '''

        guess = 5.e+4
        return op.broyden1(self.Func, guess)
