#!/usr/bin/env python

import numpy as np

def angToInvMicron(wave):
    """Helper function to convert a numpy array of wavelengths in Angstroms to 
       inverse wavelengths in 1/um.

       @param wave   an input array (or list) of wavelengths in Angstroms
       @return       a numpy array of inverse wavelengths (microns)
    """

    return 1. / (wave[:] / 10000.)

class ExtClass(object):
    """Interstellar extinction class"""

    __EXTINCTION_LAWS__ = { \
        'GAL_CCM':    ('Avg Galactic: Cardelli, Clayton & Mathis (1989, ApJ, 345, 245)', None),\
        'GAL_GCC':    ('Avg Galactic: Gordon, Cartledge & Clayton (2009, ApJ, 705, 1320)', None),\
        'GAL_JBK':    ('Avg Galactic: Kaler (1976, ApJS, 31, 517)',          'Gal_Kaler.txt'),\
        'GAL_SM':     ('Avg Galactic: Savage & Mathis (1979, ARA&A, 17, 73)','Gal_SM79.txt'),\
        'LMC_HOWARTH':('LMC avg: Howarth (1983, MNRAS, 203, 301)',           'LMC_Howarth.txt'),\
        'LMC_GCMLW':  ('LMC avg: Gordon et al. (2003, ApJ, 594,279)',        'LMC_Gordon.txt'),\
        'SMC_GCMLW':  ('SMC bar: Gordon et al. (2003, ApJ, 594,279)',        'SMC_Gordon.txt'),\
        'SMC_PREVOT': ('SMC avg: Prevot et al. (1984, A&A, 132, 389)',       'SMC_Prevot.txt')\
        }

    def __init__(self):
        """Constructor

        Objects include in-memory interpolants for all supported extinction functions."""

        self.extFnDict = {}
        for extLaw in ExtClass.__EXTINCTION_LAWS__.keys():
            self.extFnDict.update({extLaw: self.__makeExtinctionFn__(extLaw)})


    def printExtinctionLaws(self):
        """Print supported extinction law acronyms and literature references."""

        for extLaw in ExtClass.__EXTINCTION_LAWS__.keys():
            print "%11s: %s" % (extLaw, self.__EXTINCTION_LAWS__[extLaw][0])


    def getExtLawNames(self):
        """Return a list of keyword names for the extinction laws."""

        return self.__EXTINCTION_LAWS__.keys()


    def f_lambda(self, extLaw, wave):
        """Compute extinction function at selected wavelengths for a named extinction law. 

           @param extLaw choice of extinction law
           @param wave   an input array (or list) of wavelengths in Angstroms
           @return       computed numpy array of reddening function values @wave
        """

        # Derive inverse wavelenths in microns from supplied wavelengths (in Angstroms)
        xx = angToInvMicron(np.asarray(wave))

        extFunc = self.extFnDict[extLaw.strip().upper()]
        return np.interp(xx, extFunc[:,0], extFunc[:,1])


    def deredden(self, extLaw, c, wave, flux):
        """Correct a numpy array of emission line fluxes for interstellar extinction.  

           @param extLaw choice of extinction law
           @param c      extinction constant
           @param wave   an input array (or list) of wavelengths in Angstroms
           @param flux   an input array (or list) of observed emission line fluxes
           @return       a numpy array of corrected emission line fluxes
        """

        assert (len(wave) == len(flux)), "Wavelength and flux arrays must have identical size"
        return flux * pow(10, c * self.f_lambda(extLaw, wave))


    def __gcc_ExtLaw__(self):
        """Generate a sampled extinction function A(lambda)/A(V) based on the Average 
           Galactic interstellar extinction function of Gordon, Cartledge & Clayton 
           (2009, ApJ, 705, 1320)

           This extinction law was originally defined as a set of piecewise continuous 
           polynomials in X (inverse microns). It is evaluated here on a fine grid in 
           wavenumber (steps of 0.05) in order to present the same interface as the 
           other tabulated extinction laws."""

        # For wavelengths longward of the L passband, linearly interpolate to 0.0 
        # at 1/lambda = 0.  (a=0.08, b=-0.0734 @ x=0.29)
        extFn_fir = np.asarray([[1.e-9, 0.], [0.29, 0.0563]])

        # IR region (0.30 x < 1.1)
        x_start=0.30; x_end=1.1; nSteps=int((x_end-x_start)*20)+1
        x_ir = np.linspace(x_start, x_end, num=nSteps)
        ext_ir = (0.574 - 0.527/3.1)*pow(x_ir, 1.61)
        extFn_ir = np.column_stack((x_ir, ext_ir))

        # Optical (1.1 < x <= 3.3)
        x_start=1.15; x_end=3.3; nSteps=int((x_end-x_start)*20)+1
        xx = np.linspace(x_start, x_end, num=nSteps)
        y = xx - 1.82
        a = 1 + y*(0.17699 + y*(-0.50447 + y*(-0.02427 + y*(0.72085 + \
                y*(0.01979 + y*(-0.77530 + y*0.32999))))))
        b = y*(1.41338 + y*(2.28305 + y*(1.07233 + y*(-5.38434 + \
            y*(-0.62251 + y*(5.30260 - y*2.09002)))))) 
        extFn_o = np.column_stack((xx, a + b/3.1))

        # UV (3.3 < x <= 5.9)
        # The coefficients are obviously not correct; 
        # the substitutions explore the possibility of a typo in the UV-bump term
        x_start=3.35; x_end=5.9; nSteps=int((x_end-x_start)*20)+1
        xx = np.linspace(x_start, x_end, num=nSteps)
#        a =  1.896 - 0.372*xx - 0.0108 / ((xx - 4.57)**2 + 0.0422)
        a =  1.896 - 0.372*xx - 0.108 / ((xx - 4.57)**2 + 0.0422)
        b = -3.503 + 2.057*xx + 0.7180 / ((xx - 4.59)**2 + 0.0530)
        extFn_uv = np.column_stack((xx, a + b/3.1))

        # Far-UV (5.9 < x <= 11.0)
        x_start=5.95; x_end=11.0; nSteps=int((x_end-x_start)*20)+2
        xx = np.linspace(x_start, x_end, num=nSteps)
#        a =  1.896 - 0.372*xx - 0.0108 / ((xx - 4.57)**2 + 0.0422)
        a =  1.896 - 0.372*xx - 0.108 / ((xx - 4.57)**2 + 0.0422)
        b = -3.503 + 2.057*xx + 0.7180 / ((xx - 4.59)**2 + 0.0530)
        y = xx - 5.9
        a += -(0.110 + 0.0099*y) * y**2
        b +=  (0.537 + 0.0530*y) * y**2
        extFn_fuv = np.column_stack((xx, a + b/3.1))

        # Contatenate piece-wise, sampled extinction function
        return np.concatenate((extFn_fir, extFn_ir, extFn_o, extFn_uv, extFn_fuv))


    def __ccm_ExtLaw__(self):
        """Generate a sampled extinction function A(lambda)/A(V) based on the 
           average Galactic interstellar extinction function of Cardelli, 
           Clayton & Mathis (1989, ApJ, 345, 245)

           This extinction law was originally defined as a set of piecewise continuous 
           polynomials in X (inverse microns). It is evaluated here on a fine grid in 
           wavenumber (steps of 0.05) in order to present the same interface as the 
           other tabulated extinction laws."""

        # For wavelengths longward of the L passband, linearly interpolate to 0.0 
        # at 1/lambda = 0.  (a=0.08, b=-0.0734 @ x=0.29)
        extFn_fir = np.asarray([[1.e-9, 0.], [0.29, 0.0563]])

        # IR region (0.30 x < 1.1)
        x_start=0.30; x_end=1.1; nSteps=int((x_end-x_start)*20)+1
        x_ir = np.linspace(x_start, x_end, num=nSteps)
        ext_ir = (0.574 - 0.527/3.1)*pow(x_ir, 1.61)
        extFn_ir = np.column_stack((x_ir, ext_ir))

        # Optical (1.1 < x <= 3.3)
        x_start=1.15; x_end=3.3; nSteps=int((x_end-x_start)*20)+1
        x_o = np.linspace(x_start, x_end, num=nSteps)
        y = x_o - 1.82
        a = 1 + y*(0.17699 + y*(-0.50447 + y*(-0.02427 + y*(0.72085 + \
                y*(0.01979 + y*(-0.77530 + y*0.32999))))))
        b = y*(1.41338 + y*(2.28305 + y*(1.07233 + y*(-5.38434 + \
            y*(-0.62251 + y*(5.30260 - y*2.09002)))))) 
        extFn_o = np.column_stack((x_o, a + b/3.1))

        # Near-UV (3.3 < x <= 5.9)
        x_start=3.35; x_end=5.9; nSteps=int((x_end-x_start)*20)+1
        x_nuv = np.linspace(x_start, x_end, num=nSteps)
        a =  1.752 - 0.316*x_nuv - 0.104 / ((x_nuv - 4.67)**2 + 0.341)
        b = -3.090 + 1.825*x_nuv + 1.206 / ((x_nuv - 4.62)**2 + 0.263)
        extFn_nuv = np.column_stack((x_nuv, a + b/3.1))

        # Mid-UV (5.9 < x <= 8.0)
        x_start=5.95; x_end=8.0; nSteps=int((x_end-x_start)*20)+1
        x_uv = np.linspace(x_start, x_end, num=nSteps)
        a =  1.752 - 0.316*x_uv - 0.104 / ((x_uv - 4.67)**2 + 0.341)
        b = -3.090 + 1.825*x_uv + 1.206 / ((x_uv - 4.62)**2 + 0.263)
        y = x_uv - 5.9
        a += -(0.04473 + 0.009779*y) * y**2
        b +=  (0.21300 + 0.1207*y) * y**2
        extFn_uv = np.column_stack((x_uv, a + b/3.1))

        # Far-UV (8.0 < x <= 10.0)
        x_start=8.05; x_end=10.0; nSteps=int((x_end-x_start)*20)+2
        x_fuv = np.linspace(x_start, x_end, num=nSteps)
        y = x_fuv - 8.0
        a = -1.073 - y*(0.628 - y*(0.137 - y*0.070))
        b = 13.670 + y*(4.257 - y*(0.420 - y*0.374))
        extFn_fuv = np.column_stack((x_fuv, a + b/3.1))

        # Contatenate piece-wise, sampled extinction function
        return np.concatenate((extFn_fir, extFn_ir, extFn_o, extFn_nuv, extFn_uv, extFn_fuv))


    def __makeExtinctionFn__(self, extLaw):
        """Return a normalized extinction function interpolant of choice.
 
           @param extLaw  choice of extinction law
           @return        a numpy array of the tabulated, normalized extinction function

           Normalization is: A(H-beta) = 0.; and A(Inf) = -1."""

        # Ensure that a valid extinction law was specified, or raise an exception. 
        # Note: for now we use an "assert" statement, but this needs to improve. 
        extLaw = extLaw.strip().upper()
        assert (extLaw in ExtClass.__EXTINCTION_LAWS__.keys()), \
            "Invalid choice of extinction law: %s" % extLaw

        # Get a tabulated extinction function from the named file (where applicable) 
        # and return a 2-D numpy array, with: 
        #     the first field populated with wavelengths, 
        #     second field populated with extinction @wavelength
        if (extLaw == 'GAL_CCM'):
            extFn = self.__ccm_ExtLaw__()
        elif (extLaw == 'GAL_GCC'):
            extFn = self.__gcc_ExtLaw__()
        else:
            extFn = np.loadtxt(ExtClass.__EXTINCTION_LAWS__[extLaw][1])

        # Each literature reference uses a different normalization of the extinction law,
        # but we want the extinction relative to H-beta (4861). So, determine the 
        # re-normalization and offset for each law: 
 
        # Cardelli, Clayton & Mathis (1989) for average Galactic
        # A(x)/E(B-V) for x = 1/wave(in microns)
        if (extLaw == 'GAL_CCM'):
            scale = 1.16420605; offset = 1.

        # Gordon, Cartledge & Clayton (2009) for average Galactic
        # A(x)/E(B-V) for x = 1/wave(in microns)
        if (extLaw == 'GAL_GCC'):
            scale = 1.164; offset = 1.

        # Savage & Mathis (1979) for average Galactic
        # A(x)/E(B-V) for x = 1/wave(in microns)
        if (extLaw == 'GAL_SM'):
            scale = 3.1 * 1.169930; offset = 1.

        # Kaler (1976) for average Galactic
        # f_lambda(wave) for wave (in Ang)
        if (extLaw == 'GAL_JBK'):
            # Reverse the sequence order & convert to inverse wavelength in microns
            extFn = extFn[::-1]
            extFn[:,0] = angToInvMicron(extFn[:,0])
            scale = 1.0; offset = 0.

        # Howarth (1983) for LMC average
        # A(x)/E(B-V) for x = 1/wave(in microns), sampling is 0.05 wavenumbers
        elif (extLaw == 'LMC_HOWARTH'):
            # Violet
            x_start=1.85; x_end=2.75; nSteps=int((x_end-x_start)*20)+1
            x_v = np.linspace(x_start, x_end, num=nSteps)
            delt = x_v - 1.82
            ext_v = 3.1 + (2.04 + 0.094 * delt) * delt
            extFn_v = np.column_stack((x_v, ext_v))

            # Ultraviolet out to ~1110 A
            x_start=2.8; x_end=9.0; nSteps=int((x_end-x_start)*20)+1
            x_uv = np.linspace(x_start, x_end, num=nSteps)
            delt = x_uv - 4.557
            ext_uv  = 3.1 - 0.236 + (0.462 + 0.105*x_uv)*x_uv + 0.454 / (delt**2 + 0.293)
            extFn_uv = np.column_stack((x_uv, ext_uv))

            extFn = np.concatenate((extFn, extFn_v, extFn_uv))
            scale = 3.1 * 1.15771599; offset = 1.

        # Gordon et al. (2003) for LMC average
        # A(x)/A(V) for x = 1/wave(in microns)
        elif (extLaw == 'LMC_GCMLW'):
            scale = 1.15394602; offset = 1.

        # Gordon et al. (2003) for SMC bar
        # A(x)/A(V) for x = 1/wave(in microns)
        elif (extLaw == 'SMC_GCMLW'):
            scale = 1.196504471; offset = 1.

        # Prevot et al. (1984) for average SMC
        # E(x)/E(B-V) for x = 1/wave(in microns)
        elif (extLaw == 'SMC_PREVOT'):
#            extFn[:,1] = (extFn[:,1] / 3.1 + 1.) / 1.144282 - 1.
            scale = 3.1*1.144282; offset = 0.144282/1.144282

        # Normalize extinction at 4861 Ang (H-beta) to 0.
        extFn[:,1] = extFn[:,1] / scale - offset
        return extFn


##########
def extinctionTest():
    """Tests for ExtClass"""

    # Reference wavelengths, fluxes, and F_lambda for each extinction law
    WV_EVAL  = [3726.0, 4101.8, 4340.5, 4861.3, 5875.7, 6562.8, 7319.4,  89900.]
    FLX_EVAL = [28.0,     25.1,   45.6,  100.0,   15.8,  308.7,   11.0,   125.6]
    EXT_REF = {\
        'GAL_CCM':     [0.322, 0.229, 0.156, 0.0, -0.203, -0.298, -0.399, -0.981], \
        'GAL_GCC':     [0.322, 0.229, 0.156, 0.0, -0.203, -0.298, -0.399, -0.981], \
        'GAL_SM':      [0.275, 0.191, 0.143, 0.0, -0.210, -0.309, -0.405, -0.983], \
        'GAL_JBK':     [0.299, 0.195, 0.133, 0.0, -0.215, -0.331, -0.431, -0.991], \
        'LMC_HOWARTH': [0.374, 0.225, 0.145, 0.0, -0.208, -0.307, -0.404, -0.984], \
        'LMC_GCMLW':   [0.307, 0.195, 0.135, 0.0, -0.207, -0.320, -0.419, -0.994], \
        'SMC_GCMLW':   [0.386, 0.244, 0.166, 0.0, -0.233, -0.340, -0.441, -0.997], \
        'SMC_PREVOT':  [0.336, 0.203, 0.131, 0.0, -0.193, -0.293, -0.392, -0.983] \
    }

    a = ExtClass()
    a.printExtinctionLaws()
    a.getExtLawNames()

    # Use list of reference test wavelengths (in Angstroms)
    xx = angToInvMicron(np.asarray(WV_EVAL))

    # Results for extinction laws should agree with reference calculations
    max_dev = 1.e-3
    for key in EXT_REF.keys():
        diff = a.f_lambda(key, WV_EVAL) - EXT_REF[key]
        dev = abs(np.max(diff)) + abs(np.min(diff))
        if dev < max_dev:
            grade = 'PASS'
        else:
            grade = 'FAIL'
        print "%10s F(lambda) max deviation: %6.4f -- %s" % (key, dev, grade)

    # Invalid extinction law should raise an exception
    try:
        dummy = a.f_lambda('BOGUS', WV_EVAL)
    except:
        print "BOGUS caught as invalid extinction law"

    # Try de-reddening some fluxes
    fluxCorr = a.deredden('GAL_JBK', 0.1, WV_EVAL, FLX_EVAL)

    print "Deredden using extinction law GAL_JBK, extinction constant c=0.1"
    print "     Wave    Input   Output"
    for i in range(len(WV_EVAL)):
        print "%9.1f %8.1f %8.1f" % (WV_EVAL[i], FLX_EVAL[i], fluxCorr[i])
