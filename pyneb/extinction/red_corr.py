import numpy as np
import pyneb as pn

if pn.config.INSTALLED['plt']:
    import matplotlib.pyplot as plt
if pn.config.INSTALLED['scipy']:
    from scipy import interpolate
from pyneb.utils.misc import execution_path

def my_X2(wave, params=[5000., 1., 2., 3.]):
    return params[1] * (wave / params[0]) + params[2] * (wave / params[0]) ** -1 + params[3] * (wave / params[0]) ** -2

def poly(x, coeffs):
    res = 0
    for i, coeff in enumerate(coeffs):
        res += coeff * x ** i
    return res

class RedCorr(object):
    """
    Reddening correction
    RC = RedCorr()
    
    """
    ##
    # @todo Manage error in extinction
    # @todo print extinction laws with references
    def __init__(self, E_BV=0., R_V=3.1, law='No correction', cHbeta=None,
                 UserFunction=None):
        """
        Reddening correction tool.

        Usage:
            RC = RedCorr(E_BV = 1.)
            RC.plot(laws = 'all')

            def my_X(wave, params = [5000., 1., 2., 3.]):
                \"""
                Description of the user defined law
                \"""
                return params[1] * (wave/params[0]) + params[2] * (wave/params[0])**-1 + params[3] * (wave/params[0])**-2

            RC.UserFunction = my_X
            RC.UserParams = [6000., 0., 0., 1.]
            RC.getCorr(5007)
        
        Parameters:
            - E_BV [float]: differential extinction between bands B and V
            - R_V = AV/E_BV. Default value is 3.1
            - law [str]: one of the defined laws (available with RedCorr.getLaws()) 
            - cHbeta: logarithmic extinction at Hbeta (prevalence on E_BV)
            - UserFunction X(wave, param): A user-defined function that accept 2 parameters: wavelength(s) in Angstrom 
            and an optional parameter and return X(lambda) = A(lambda)/E_BV = R.A(lambda)/AV. 
            The correction is: 10**(0.4*E_bv*X)
            
        """
        
        self.log_ = pn.log_ 
        self.calling = 'RedCorr'
        self._laws_dict = {} # dictionary pointing to a reddening function depending on the key
        self._laws_dict['No correction'] = self._zeros
        self._laws_dict['CCM89'] = self._CCM89
        self._laws_dict['CCM89 Bal07'] = self._CCM89_Bal07
        self._laws_dict['CCM89 oD94'] = self._CCM89_oD94
        self._laws_dict['S79 H83 CCM89'] = self._S79_H83_CCM89
        #self._laws_dict['GCC 09'] = self._GCC09 
        self._laws_dict['K76'] = self._K76
        self._laws_dict['SM79 Gal'] = self._SM79_Gal
        self._laws_dict['G03 LMC'] = self._G03_LMC
        self._laws_dict['MCC99 FM90 LMC'] = self._MCC99_FM90_LMC
        self._laws_dict['F99-like'] = self._F99_like
        #self._laws_dict['F99-like IDL'] = self._F99_like_IDL
        self._laws_dict['F99'] = self._F99
        self._laws_dict['F88 F99 LMC'] = self._FM88_F99_LMC

        self.UserFunction = UserFunction
        self.UserParams = None
        self.FitzParams = None
        
        self._given = None
        self.law = law
        self.R_V = R_V
        if cHbeta is not None:
            self._given = 'cHbeta'
            self.cHbeta = cHbeta
        else:
            self._given = 'E_BV'
            self.E_BV = E_BV
 
 
    def cHbetaFromEbv(self, ebv):
        """
        Return cHbeta from E(BV)
        Usage: 
            (1-f_lambda).cHbeta = 0.4.EBV.X_lambda applied to lambda = 4861, with f_beta = 0.:
            cHbeta = 0.4 . EBV . X_beta
        
        Parameter:
            - ebv    E(B-V)
            
        """
        Xbeta = self.X(pn.CST.HBETA)
        return 0.4 * ebv * Xbeta
        #return np.asarray((-0.61 + (0.61 ** 2 + 4 * 0.024 * ebv) ** 0.5) / (2 * 0.024))
    
 
    def EbvFromCHbeta(self, cHbeta):
        """
        Return E(B-V) from cHbeta
        Using: 
            (1-f_lambda).cHbeta = 0.4.EBV.X_lambda applied to lambda = 4861, with f_beta = 0.:
            cHbeta = 0.4 . EBV . X_beta
        
        Parameter:
            - cHbeta
            
        """
        Xbeta = self.X(pn.CST.HBETA)
        if Xbeta != 0.:
            return cHbeta * 2.5 / Xbeta
        else:
            return np.zeros_like(cHbeta)
        #return np.asarray(0.61 * cHbeta + 0.024 * cHbeta ** 2.)
    
    
    def getLaws(self):
        """
        Return the dictionary keys for the extinction laws
        
        """
        return self._laws_dict.keys()
    
    
    def printLaws(self):
        """
        Print out the extinction laws
        
        """
        for law in self._laws_dict.keys():
            try:
                doc = self._laws_dict[law].__doc__
            except:
                doc = ''
            print("'{0}': {1}".format(law, doc))
                         
    def _get_e_bv(self):
        return self.__E_BV
    def _get_r_v(self):
        return self.__R_V
    def _get_law(self):
        return self.__law
    def _get_cHbeta(self):
        return self.__cHbeta
    def _get_AV(self):
        return self.__AV
    def _get_uf(self):
        return self.__user_function
    
    def _set_e_bv(self, value):
        self.__E_BV = np.asarray(value)
        self.__cHbeta = self.cHbetaFromEbv(self.__E_BV)
        self.__AV = self.__E_BV * self.R_V
        self._given = 'E_BV'
    def _set_cHbeta(self, value):
        self.__cHbeta = np.asarray(value)
        self.__E_BV = self.EbvFromCHbeta(self.__cHbeta)
        self.__AV = self.__E_BV * self.R_V
        self._given = 'cHbeta'
    def _set_AV(self, value):
        self.__AV = np.asarray(value)
        self.__E_BV = self.__AV / self.R_V
        self.__cHbeta = self.cHbetaFromEbv(self.__E_BV)
        self._given = 'AV'
    def _set_r_v(self, value):       
        self.__R_V = np.asarray(value)
    def _set_law(self, value):
        if value not in self._laws_dict.keys():
            self.log_.error('Unknown extinction law reference: {0}'.format(value), calling=self.calling)
            self.__law = None
            self.X = None
        else:
            self.__law = value
            self.X = self._laws_dict[self.law]
            if self._given == 'E_BV':
                self.E_BV = self.E_BV
            elif self._given == 'cHbeta':  
                self.cHbeta = self.cHbeta
            elif self._given == 'AV':  
                self.AV = self.AV
    def _set_uf(self, value):
        self.__user_function = value
        if value is None:
            if 'user' in self._laws_dict:
                del self._laws_dict['user']
        else:
            def _uf2(wave):
                #This transform the user function with 2 parameters in a function of one single parameter
                return np.asarray(self.__user_function(wave, self.UserParams))
            _uf2.__doc__ = self.__user_function.__doc__
            self._laws_dict['user'] = _uf2
        
    
    E_BV = property(_get_e_bv, _set_e_bv, None, None)
    R_V = property(_get_r_v, _set_r_v, None, None)
    law = property(_get_law, _set_law, None, None)
    cHbeta = property(_get_cHbeta, _set_cHbeta, None, None)
    AV = property(_get_AV, _set_AV, None, None)
    UserFunction = property(_get_uf, _set_uf, None, None)
    
        
    def getCorr(self, wave, rel_wave=None): 
        """
        Return the extinction correction as:
            correction = 10**(0.4 * EBV * Xx) = 10**(A_lambda / 2.5)

        Usage:
            RC.getCorr(5007)
            RC.getCorr(5007, 4861)
            
        Parameters:
            - wave      wavelength (Angstrom)
            - rel_wave  wavelength (Angstrom) for a relative correction

        """
        if self.law is None:
            self.log_.warn('No extinction law defined.', calling=self.calling)
            return None
        if self._laws_dict[self.law] is None:
            self.log_.warn('No user defined extinction law.', calling=self.calling)
            return None
        else:
            if rel_wave is None:
                rel_corr = 1.
            else:
                rel_corr = self.getCorr(rel_wave) 
            X = self.X(wave)
            return np.squeeze(10. ** (0.4 * np.outer(X, self.E_BV).reshape(X.shape + self.E_BV.shape))) / rel_corr
    
    
    def getCorrHb(self, wave):
        """
        Return the extinction correction normalized to the correction at 4861AA.
            
        Parameter:
            - wave      wavelength (Angstrom)
        
        """
        return self.getCorr(wave, np.ones_like(wave) * pn.CST.HBETA)
    
    
    def getErrCorr(self, wave, err_E_BV, rel_wave=None):
        """
        Return the error on the correction for a given wavelength, given the error on E(B-V)
        
        Parameters:
            - wave         wavelength(s)
            - err_E_BV     error on E(B-V)
            - rel_wave     reference wavelength for the normalization (optional)

        """
        if rel_wave is None:
            rel_X = 0.
        else:
            rel_X = self.X(rel_wave)
        return  np.log(10) * abs(self.X(wave) - rel_X) * 0.4 * err_E_BV * self.E_BV


    def getErrCorrHb(self, wave, err_E_BV):
        """
        Return the the error on the correction relative to Hbeta for a given wavelength, 
            given the error on E(B-V)
        
        Parameters:
            - wave         wavelength(s)
            - err_E_BV     error on E(B-V)

        """
        return self.getErrCorr(self, wave, err_E_BV, rel_wave=pn.CST.HBETA)

    
    def setCorr(self, obs_over_theo, wave1, wave2):
        """
        Determination of the correction using the ratio of two observed line intensities 
            relative to the theoretical value.
        
        Usage:
            rc.setCorr(6.5/2.85, 6563., 4861.)
            
        Parameters:
            - obs_over_theo    ration of the observed ratio over the theoretical ratio
            - wave1, wave2     wavelengths at which the line rations are taken.

        """
        COR = RedCorr(E_BV= -2.5, R_V=self.R_V, law=self.law, UserFunction=self.UserFunction)
        f1 = np.log10(COR.getCorr(wave1))
        f2 = np.log10(COR.getCorr(wave2))
        if f1 != f2:
            self.E_BV = 2.5 * np.log10(obs_over_theo) / (f1 - f2)
        else:
            self.E_BV = 0.
        del COR

    def plot(self, w_inf=1000., w_sup=10000., laws=None, ax=None, **kwargs):
        """
        Plot extinction laws

        Parameters:
            - w_inf [float] lower limit of plot
            - w_sup [float] upper limit of plot
            - laws [list of strings] list of extinction law labels. If set to 'all', all the laws are plotted
            - ax : an axis object. If None, onr is created
            - **kwargs arguments to plot

        """
        if ax is None:
            f, ax = plt.subplots()
        colors = ['r', 'g', 'b', 'y', 'm', 'c']
        styles = ['-', '--', '-:', ':']
        if not pn.config.INSTALLED['plt']:
            pn.log_.error('matplotlib.pyplot not available for plotting', calling=self.calling)
        old_E_BV = self.E_BV
        old_law = self.law
        self.E_BV = 2.5
        w = np.linspace(w_inf, w_sup, 1000)
        if laws is None:
            laws = self.law
        elif laws == 'all':
            laws = self.getLaws()
        if type(laws) is str:
            laws = [laws]
        laws.sort()    
        for i, law in enumerate(laws):
            self.law = law
            corr = self.getCorrHb(w)
            if corr is not None:
                ax.plot(w, np.log10(corr), label=law, c=colors[i % 6], linestyle=styles[i // 6], **kwargs)
        ax.legend()
        ax.set_xlabel('Wavelength (A)')
        ax.set_ylabel('X')
        self.law = old_law
        self.E_BV = old_E_BV
    
                            
    def _CCM89(self, wave):
        """
        Cardelli, Clayton & Mathis 1989, ApJ 345, 245
        http://adsabs.harvard.edu/abs/1989ApJ...345..245C

        Comments: Depends on R_V, default value being 3.1

        Scope: Applicable to both dense and diffuse ISM
        Range: UV through IR
        
        """
        x = 1e4 / np.asarray([wave]) # inv microns
        a = np.zeros_like(x)
        b = np.zeros_like(x)
        
        tt = (x > 0.3) & (x <= 1.1)
        a[tt] = 0.574 * x[tt] ** 1.61 
        b[tt] = -0.527 * x[tt] ** 1.61
    
        tt = (x > 1.1) & (x <= 3.3)
        yg = x[tt] - 1.82
        a[tt] = (1. + 0.17699 * yg - 0.50447 * yg ** 2. - 0.02427 * yg ** 3. + 0.72085 * yg ** 4. + 
                 0.01979 * yg ** 5. - 0.7753 * yg ** 6. + 0.32999 * yg ** 7.)
        b[tt] = (0. + 1.41338 * yg + 2.28305 * yg ** 2. + 1.07233 * yg ** 3. - 5.38434 * yg ** 4. - 
                 0.622510 * yg ** 5. + 5.3026 * yg ** 6. - 2.09002 * yg ** 7.)
        
        tt = (x > 3.3) & (x <= 5.9)
        a[tt] = 1.752 - 0.316 * x[tt] - 0.104 / ((x[tt] - 4.67) ** 2. + 0.341)
        b[tt] = -3.090 + 1.825 * x[tt] + 1.206 / ((x[tt] - 4.62) ** 2 + 0.263)
        
        tt = (x > 5.9) & (x <= 8.0)
        a[tt] = (1.752 - 0.316 * x[tt] - 0.104 / ((x[tt] - 4.67) ** 2. + 0.341) - 
                 0.04473 * (x[tt] - 5.9) ** 2. - 0.009779 * (x[tt] - 5.9) ** 3.)
        b[tt] = (-3.090 + 1.825 * x[tt] + 1.206 / ((x[tt] - 4.62) ** 2. + 0.263) + 
                 0.2130 * (x[tt] - 5.9) ** 2. + 0.1207 * (x[tt] - 5.9) ** 3.)
        
        tt = (x > 8.0) & (x < 10.0)
        a[tt] = (-1.073 - 0.628 * (x[tt] - 8) + 0.137 * (x[tt] - 8) ** 2. - 
                 0.070 * (x[tt] - 8) ** 3.)
        b[tt] = (13.670 + 4.257 * (x[tt] - 8) - 0.420 * (x[tt] - 8) ** 2. + 
                 0.374 * (x[tt] - 8) ** 3.)
        
        Xx = self.R_V * a + b
        return np.squeeze(Xx)


    def _CCM89_Bal07(self, wave):
        """
        Galactic extinction law based on Cardelli et al 1989, modified by Blagrave et al 2007 
        for 3.3 < x < 8 (1250 < lambda < 3030)
        
        Blagrave et al 2007, ApJ, 655, 299 
        http://adsabs.harvard.edu/abs/2007ApJ...655..299B
        Cardelli, Clayton & Mathis 1989, ApJ 345, 245
        http://adsabs.harvard.edu/abs/1989ApJ...345..245C

        Comments:
        Same as CCM89 for x<3.3 and x>8
        Revised values for 3.3<x<8
        Based on observation of Orion stars
        Depends on R_V, default value being 3.1

        Range: UV through IR
        
        """
        x = 1e4 / np.asarray([wave]) # inv microns
        a = np.zeros_like(x)
        b = np.zeros_like(x)
        
        tt = (x > 0.3) & (x <= 1.1)
        a[tt] = 0.574 * x[tt] ** 1.61 
        b[tt] = -0.527 * x[tt] ** 1.61
    
        tt = (x > 1.1) & (x <= 3.3)
        yg = x[tt] - 1.82
        a[tt] = (1. + 0.17699 * yg - 0.50447 * yg ** 2. - 0.02427 * yg ** 3. + 0.72085 * yg ** 4. + 
                 0.01979 * yg ** 5. - 0.7753 * yg ** 6. + 0.32999 * yg ** 7.)
        b[tt] = (0. + 1.41338 * yg + 2.28305 * yg ** 2. + 1.07233 * yg ** 3. - 5.38434 * yg ** 4. - 
                 0.622510 * yg ** 5. + 5.3026 * yg ** 6. - 2.09002 * yg ** 7.)
        
        tt = (x > 3.3) & (x <= 5.9)
        a[tt] = 1.752 - 0.316 * x[tt] - 0.104 / ((x[tt] - 4.67) ** 2. + 0.341)
        b[tt] = -2.9 + 1.825 * x[tt] + 0.93 / ((x[tt] - 4.65) ** 2 + 0.263)
        
        
        tt = (x > 5.9) & (x <= 8.0)
        a[tt] = (1.752 - 0.316 * x[tt] - 0.104 / ((x[tt] - 4.67) ** 2. + 0.341) - 
                 0.04473 * (x[tt] - 5.9) ** 2. - 0.009779 * (x[tt] - 5.9) ** 3.)
        b[tt] = (-2.9 + 1.825 * x[tt] + 0.93 / ((x[tt] - 4.65) ** 2 + 0.263) + 
                 0.2130 * (x[tt] - 5.9) ** 2. + 0.1207 * (x[tt] - 5.9) ** 3.)
        
        tt = (x > 8.0) & (x < 10.0)
        a[tt] = (-1.073 - 0.628 * (x[tt] - 8) + 0.137 * (x[tt] - 8) ** 2. - 
                 0.070 * (x[tt] - 8) ** 3.)
        b[tt] = (13.670 + 4.257 * (x[tt] - 8) - 0.420 * (x[tt] - 8) ** 2. + 
                 0.374 * (x[tt] - 8) ** 3.)
        
        Xx = self.R_V * a + b
        return np.squeeze(Xx)
    

    def _CCM89_oD94(self, wave):
        """
        Galactic extinction law based on Cardelli et al 1989, modified by O'Donnell 1994
        for 1.1 < x < 3.3 (9100 < lambda < 3030)
        
        O'Donnell 1994, ApJ, 422, 1580
        http://adsabs.harvard.edu/abs/1994ApJ...422..158O
        Cardelli, Clayton & Mathis 1989, ApJ 345, 245
        http://adsabs.harvard.edu/abs/1989ApJ...345..245C

        Comments:
        Same as CCM89 for x<1.1 and x>3.3
        Revised values for 1.1<x<3.3
        Produces lower correction in the near UV at low R_V
        
        Scope: Galactic
        Range: UV through IR
        
        """
        x = 1e4 / np.asarray([wave]) # inv microns
        a = np.zeros_like(x)
        b = np.zeros_like(x)
        
        tt = (x > 0.3) & (x <= 1.1)
        a[tt] = 0.574 * x[tt] ** 1.61 
        b[tt] = -0.527 * x[tt] ** 1.61
    
        tt = (x > 1.1) & (x <= 3.3)
        yg = x[tt] - 1.82
        a[tt] = (1. + 0.104 * yg - 0.609 * yg ** 2. + 0.701 * yg ** 3. + 1.137 * yg ** 4. - 
                 1.718 * yg ** 5. - 0.827 * yg ** 6. + 1.647 * yg ** 7. - 0.505 * yg ** 8.)
        b[tt] = (0. + 1.952 * yg + 2.908 * yg ** 2. - 3.989 * yg ** 3. - 7.985 * yg ** 4. + 
                 11.102 * yg ** 5. + 5.491 * yg ** 6. - 10.805 * yg ** 7. + 3.347 * yg ** 8.)
        
        tt = (x > 3.3) & (x <= 5.9)
        a[tt] = 1.752 - 0.316 * x[tt] - 0.104 / ((x[tt] - 4.67) ** 2. + 0.341)
        b[tt] = -3.090 + 1.825 * x[tt] + 1.206 / ((x[tt] - 4.62) ** 2 + 0.263)
        
        tt = (x > 5.9) & (x <= 8.0)
        a[tt] = (1.752 - 0.316 * x[tt] - 0.104 / ((x[tt] - 4.67) ** 2. + 0.341) - 
                 0.04473 * (x[tt] - 5.9) ** 2. - 0.009779 * (x[tt] - 5.9) ** 3.)
        b[tt] = (-3.090 + 1.825 * x[tt] + 1.206 / ((x[tt] - 4.62) ** 2. + 0.263) + 
                 0.2130 * (x[tt] - 5.9) ** 2. + 0.1207 * (x[tt] - 5.9) ** 3.)
        
        tt = (x > 8.0) & (x < 10.0)
        a[tt] = (-1.073 - 0.628 * (x[tt] - 8) + 0.137 * (x[tt] - 8) ** 2. - 
                 0.070 * (x[tt] - 8) ** 3.)
        b[tt] = (13.670 + 4.257 * (x[tt] - 8) - 0.420 * (x[tt] - 8) ** 2. + 
                 0.374 * (x[tt] - 8) ** 3.)
        
        Xx = self.R_V * a + b
        return np.squeeze(Xx)


    def _S79_H83_CCM89(self, wave):
        """
        Galactic extinction law (0-33000 A range):
        - In the UV, from Seaton 1979
        - In the opt/NIR (3600-9100) Howarth 1983
        - In the FIR (9100-33000) Cardelly et al 1989
        
        Seaton 1979, MNRAS, 187, 73) and 
        http://adsabs.harvard.edu/abs/1979MNRAS.187P..73S
        Howarth 1983, MNRAS, 203, 301) Galactic law
        http://adsabs.harvard.edu/abs/1983MNRAS.204.1091H
        Cardelli, Clayton and Mathis 1989, ApJ, 345, 245
        http://adsabs.harvard.edu/abs/1989ApJ...345..245C
        
        Scope: Galactic
        Range: UV through IR
        
        """
        x = 1e4 / np.asarray([wave]) # inv microns
        Xx = np.zeros_like(x)

        # Cardelli, Clayton & Mathis 1989 33000-9000
        tt = (x > 0.3) & (x <= 1.1)
        Xx[tt] = self.R_V * (0.574 * x[tt] ** 1.61) - 0.527 * x[tt] ** 1.61
        
        # Howarth 1983, Galactic 9000-5000
        tt = (x > 1.1) & (x <= 1.83)
        Xx[tt] = (self.R_V - 3.1) + ((1.86 - 0.48 * x[tt]) * x[tt] - 0.1) * x[tt]
        
        # Howarth 1983, Galactic 5000-3600
        tt = (x > 1.83) & (x <= 2.75)
        Xx[tt] = self.R_V + 2.56 * (x[tt] - 1.83) - 0.993 * ((x[tt] - 1.83) ** 2)

        # Seaton 1979, Galactic
        tt = (x > 2.75) & (x <= 3.65)
        Xx[tt] = (self.R_V - 3.2) + 1.56 + 1.048 * x[tt] + 1.01 / (((x[tt] - 4.6) ** 2) + 0.280)
        
        # Seaton 1979, Galactic
        tt = (x > 3.65) & (x <= 7.14)
        Xx[tt] = (self.R_V - 3.2) + 2.29 + 0.848 * x[tt] + 1.01 / (((x[tt] - 4.6) ** 2) + 0.280)
        
        # Seaton 1979, Galactic
        tt = (x > 7.14) 
        Xx[tt] = (self.R_V - 3.2) + 16.17 - 3.20 * x[tt] + 0.2975 * x[tt] ** 2
 
        return np.squeeze(Xx)
    
    
    # Removed from the code because it is most probably wrong
    #def _GCC09(self, wave):
        """
        Gordon, Cartledge & Clayton 2009, ApJ, 705, 1320
        http://adsabs.harvard.edu/abs/2009ApJ...705.1320G
        
        Comments:
        Extinction function R.A(wave)/A(V) 
        R_V dependent
        WARNING: This law seems buggy in the 2200AA region.
        
        Scope:
        UV
           
        """
        """
        x = 1e4 / np.asarray([wave]) # inv microns
        Xx = np.zeros_like(x)

        tt = (x > 0.3) & (x <= 1.1)
        Xx[tt] = (self.R_V * 0.574 - 0.527) * x[tt] ** 1.61
        
        tt = (x > 1.1) & (x <= 3.3)
        y = x[tt] - 1.82
        a = 1 + y * (0.17699 + y * (-0.50447 + y * (-0.02427 + y * (0.72085 + \
                y * (0.01979 + y * (-0.77530 + y * 0.32999))))))
        b = y * (1.41338 + y * (2.28305 + y * (1.07233 + y * (-5.38434 + \
            y * (-0.62251 + y * (5.30260 - y * 2.09002)))))) 
        Xx[tt] = self.R_V * a + b

        ##
        # @bug The coefficients are obviously not correct; 
        tt = (x > 3.3) & (x <= 5.9)
        a = 1.896 - 0.372 * x[tt] - 0.0108 / ((x[tt] - 4.57) ** 2 + 0.0422)
        b = -3.503 + 2.057 * x[tt] + 0.7180 / ((x[tt] - 4.59) ** 2 + 0.0530)
        Xx[tt] = self.R_V * a + b

        tt = (x > 5.9) & (x <= 11.0)
        a = 1.896 - 0.372 * x[tt] - 0.0108 / ((x[tt] - 4.57) ** 2 + 0.0422)
        b = -3.503 + 2.057 * x[tt] + 0.7180 / ((x[tt] - 4.59) ** 2 + 0.0530)
        y = x[tt] - 5.9
        a += -(0.110 + 0.0099 * y) * y ** 2
        b += (0.537 + 0.0530 * y) * y ** 2
        Xx[tt] = self.R_V * a + b

        return np.squeeze(Xx)
        """

    def _K76(self, wave):
        """
        Kaler 1976, ApJS, 31, 517
        http://adsabs.harvard.edu/abs/1976ApJS...31..517K
        
        Comments:
        This function returns the correction relative to Hbeta (f_lambda) and not
        the extinction law (X(1/lambda)).
        It cannot be used for absolute correction.
        
        Range: 3000 to >20000
        
        """
        w = np.asarray([wave]) # inv microns
        if self.E_BV == 0.:
            return np.ones_like(w)
        f_tab = np.loadtxt(execution_path('Gal_Kaler.txt'))
        f = np.interp(w, f_tab[:, 0], f_tab[:, 1])
        return np.squeeze(f * self.cHbeta / 0.4 / self.E_BV)
    

    def _SM79_Gal(self, wave):
        """
        Galactic extinction law
        Savage & Mathis 1979, ARA&A, 17, 73
        http://adsabs.harvard.edu/abs/1979ARA%26A..17...73S
        
        Comments:
        Average of several extinction laws
        R_V=3.1
        
        Scope: Galactic
        Range: UV through IR
        
        """
        x = 1e4 / np.asarray([wave]) # inv microns

        X_tab = np.loadtxt(execution_path('Gal_SM79.txt'))
        Xx = np.interp(x, X_tab[:, 0], X_tab[:, 1])
        return np.squeeze(Xx)
    
    
    def _G03_LMC(self, wave):
        """
        Extinction curve for the LMC
        Gordon et al. (2003, ApJ, 594,279)
        http://adsabs.harvard.edu/abs/2003ApJ...594..279G
        
        Comments:
        Average curve for the LMC
        R_V = 3.41

        Scope: LMC
        Range: 1200 through fIR   

        """
        x = 1e4 / np.asarray([wave]) # inv microns

        X_tab = np.loadtxt(execution_path('LMC_Gordon.txt'))
        Xx = self.R_V * np.interp(x, X_tab[:, 0], X_tab[:, 1])
        return np.squeeze(Xx)
    
    
    def _F99_like(self, wave):
        """
        In the UV, it returns the Fitzpatrick & Massa 1990 law.
        In the opt/IR, it returns the Fitzpatrick & Massa 1990 law. 
        
        Fitzpatrick 1999, PASP, 11, 63
        http://adsabs.harvard.edu/abs/1999PASP..111...63F
        Fitzpatrick & Massa 1990, ApJS, 72, 163
        http://adsabs.harvard.edu/abs/1990ApJS...72..163F
        
        Comments:
        The FM90 depends on 6 parameters which must be set by the user and are stored in RedCorr.FitzParams.
        For the predefined set of parameters defined in FM99, use instead the F99 law.
        R_V must be provided, as the law depends on it. The dependence with R_V follows Table 4 in the F99 paper 

        Range: UV through IR
        
        """
        def fit_UV(x):

            Xx = c1 + c2 * x
            Xx += c3 * x ** 2 / ((x ** 2 - x0 ** 2) ** 2 + (x * gamma) ** 2)
            tt2 = (x > 5.9) 
            if tt2 is not False:
                Xx[tt2] += c4 * (0.5392 * (x[tt2] - 5.9) ** 2 + 0.05644 * (x[tt2] - 5.9) ** 3)
            Xx += self.R_V
            return Xx
            
        x = 1e4 / np.asarray([wave]) # inv microns
        Xx = np.zeros_like(x)
        if self.FitzParams is None:
            pn.log_.warn('Fitzpatrick law requires FitzParams', calling=self.calling)
            return None
        x0 = self.FitzParams[0]
        gamma = self.FitzParams[1]
        c1 = self.FitzParams[2]
        c2 = self.FitzParams[3]
        c3 = self.FitzParams[4]
        c4 = self.FitzParams[5]
        
        # UV from the 1988 paper:
        xcutuv = 10000.0 / 2700.0
        tt = (x >= xcutuv) 
        Xx[tt] = fit_UV(x[tt])
        
        l2x = lambda l: 1e4 / l
        x_opir = np.array([0, l2x(26500.0), l2x(12200.0), l2x(6000.0), l2x(5470.0), l2x(4670.0), l2x(4110.0),
                  l2x(2700.), l2x(2600.)])
        norm = self.R_V / 3.1

        # Opt and IR from the 1999 paper
        y_opir = np.array([0., 0.265 * norm, 0.829 * norm, -0.426 + 1.0044 * self.R_V,
                           - 0.050 + 1.0016 * self.R_V , 0.701 + 1.0016 * self.R_V,
                           1.208 + 1.0032 * self.R_V - 0.00033 * self.R_V ** 2, fit_UV(l2x(2700.)), fit_UV(l2x(2600.))])
        tt = x < xcutuv
        if tt.sum() > 0:
            tck = interpolate.splrep(x_opir, y_opir)
            Xx[tt] = interpolate.splev(x[tt], tck, der=0)
        return np.squeeze(Xx)


# Commented out because it duplicates _F99_like
#    def _F99_like_IDL(self, wave):
        """
        Same as F99_like, but with a different function in the opt/IR fitting, based on an IDL program
        provided by F99. The results should be identical.
    
        In the UV, it returns the Fitzpatrick & Massa 1990 law.
        In the opt/IR, it returns the Fitzpatrick 1990 law. 
        
        Fitzpatrick 1999, PASP, 11, 63
        http://adsabs.harvard.edu/abs/1999PASP..111...63F
        Fitzpatrick & Massa 1990, ApJS, 72, 163
        http://adsabs.harvard.edu/abs/1990ApJS...72..163F
        
        Comments:
        The FM90 depends on 6 parameters which must be set by the user and are stored in RedCorr.FitzParams.
        For the predefined set of parameters defined in FM99, use instead the F_99 method.
        R_V must be provided, as the law depends on it. The dependence with R_V follows Table 4 in the F99 paper 

        Scope:    
        Range: UV through IR
        
        
        """
        """
        def fit_UV(x):

            Xx = c1 + c2 * x
            Xx += c3 * x ** 2 / ((x ** 2 - x0 ** 2) ** 2 + (x * gamma) ** 2)
            tt2 = (x > 5.9) 
            if tt2 is not False:
                Xx[tt2] += c4 * (0.5392 * (x[tt2] - 5.9) ** 2 + 0.05644 * (x[tt2] - 5.9) ** 3)
            Xx += self.R_V
            return Xx
            
        x = 1e4 / np.asarray([wave]) # inv microns
        Xx = np.zeros_like(x)
        if self.FitzParams is None:
            pn.log_.error('Fitzpatrick law requires FitzParams to be set', calling=self.calling)
            return None
        x0 = self.FitzParams[0]
        gamma = self.FitzParams[1]
        c1 = self.FitzParams[2]
        c2 = self.FitzParams[3]
        c3 = self.FitzParams[4]
        c4 = self.FitzParams[5]
        
        xcutuv = 10000.0 / 2700.0
        tt = (x >= xcutuv) 
        Xx[tt] = fit_UV(x[tt])
        
        l2x = lambda l: 1e4 / l
        x_opir = np.array([0, l2x(26500.0), l2x(12200.0), l2x(6000.0), l2x(5470.0), l2x(4670.0), l2x(4110.0),
                  l2x(2700.), l2x(2600.)])
        norm = self.R_V / 3.1
        
        # Opt and IR computed with the IDL program provided by Fitzpatrick 1999:        
        y_opir = np.array([0.0, 0.26469 * norm, 0.82925 * norm, poly(self.R_V, [-4.22809e-01, 1.00270, 2.13572e-04]),
                    poly(self.R_V, [-5.13540e-02, 1.00216, -7.35778e-05]),
                    poly(self.R_V, [ 7.00127e-01, 1.00184, -3.32598e-05]),
                    poly(self.R_V, [ 1.19456, 1.01707, -5.46959e-03, 7.97809e-04, -4.45636e-05]),
                    fit_UV(l2x(2700.)), fit_UV(l2x(2600.))])
        
        tt = x < xcutuv
        if tt.sum() > 0:
            tck = interpolate.splrep(x_opir, y_opir)
            Xx[tt] = interpolate.splev(x[tt], tck, der=0)
        return np.squeeze(Xx)
        """    
            
    def _MCC99_FM90_LMC(self, wave):
        """
        In the UV, this method returns the extinction curve proposed for the LMC 
        by Misselt et al 1999 based on the 1990 variant of the Fitzpatrick & Massa law
        In the opt/IR, it returns the Fitzpatrick & Massa 1990 law.

        Misselt, Clayton & Gordon 1999 , ApJ, 515, 128 
        http://adsabs.harvard.edu/abs/1999ApJ...515..128M
        Fitzpatrick & Massa 1990, ApJS, 72, 163
        http://adsabs.harvard.edu/abs/1990ApJS...72..163F
        
        Comments:
        The Fitzpatrick & Massa 1990 law in the UV depends on 6 parameters, stored in RedCorr.FitzParams.
        The method sets RedCorr.FitzParams to the values of set in the Fitzpatrick 1999 paper, 
        which includes an explicit dependence on R_V.
        R_V must be provided, as the law depends on its value.
        We refer to FM90 and not to the original FM88 because the value of a constant in F(lambda) slightly changed (0.0564 -> 0.05644)
        The value of another constant of F(lambda) appears to change from FM90 to MCC99, but it is probably a typo (0.5392 -> 0.5329) 
        
        Scope: LMC

        """
    
        x0 = 4.596  
        gamma = 0.91
        c1 = -1.28
        c2 = 1.11    
        c3 = 2.73    
        c4 = 0.64    
        self.FitzParams = [x0, gamma, c1, c2, c3, c4]
        return self._F99_like(wave)


    def _F99(self, wave):
        """
        This method returns the R-dependent IR-through-UV extinction curve proposed by Fitzpatrick 1999.

        Fitzpatrick 1999, PASP, 11, 63
        http://adsabs.harvard.edu/abs/1999PASP..111...63F
        based on: 
        Fitzpatrick & Massa 1990, ApJS, 72, 163
        http://adsabs.harvard.edu/abs/1990ApJS...72..163F
        
        Comments:
        The Fitzpatrick & Massa 1990 law in the UV depends on 6 parameters, stored in RedCorr.FitzParams.
        The method sets RedCorr.FitzParams to the values of set in the Fitzpatrick 1999 paper, 
        which includes an explicit dependence on R_V.
        R_V must be provided, as the law depends on its value.
        
        Range: UV through IR

        """
        x0 = 4.596  
        gamma = 0.99    
        c3 = 3.23    
        c4 = 0.41    
        c2 = -0.824 + 4.717 / self.R_V # 0.7 if RV=3.1
        c1 = 2.030 - 3.007 * c2 # -0.0677 if RV = 3.1
        self.FitzParams = [x0, gamma, c1, c2, c3, c4]
        return self._F99_like(wave)


    def _FM88_F99_LMC(self, wave):
        """
        This method returns:
        - in the UV, the average LMC extinction curve derived by Fitzpatrick & Massa 1988
        - in the opt/IR, the R-dependent extinction curve proposed by Fitzpatrick 1999.
        
        Fitzpatrick 1999, PASP, 11, 63
        http://adsabs.harvard.edu/abs/1999PASP..111...63F
        Fitzpatrick & Massa 1988, ApJ, 328, 734
        http://adsabs.harvard.edu/abs/1988ApJ...328..734F
        
        Comments:
        The Fitzpatrick and Massa law in the UV depends on 6 parameters, stored in RedCorr.FitzParams and 
        here set to the LMC values derived in FM88
        R_V must be provided, as the law depends on it
        
        Scope: LMC
        Range: UV through IR

        """
        x0 = 4.608
        gamma = 0.994
        c1 = -0.687
        c2 = 0.891
        c3 = 2.55
        c4 = 0.504  
        self.FitzParams = [x0, gamma, c1, c2, c3, c4]
        return self._F99_like(wave)
    
    
    def _zeros(self, wave):
        """
        No correction, return 0.0
        
        """
        return np.zeros_like(wave)
    
