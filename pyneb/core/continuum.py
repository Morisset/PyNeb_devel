import pickle
import sys
import numpy as np
from scipy.interpolate import interp1d
from .pynebcore import RecAtom
from ..utils.physics import CST
from ..utils.misc import execution_path
    
class Continuum(object):
    
    def __init__(self, wl_min = 3500, wl_max = 4000):
        """
        
        """
        self.wl_min = wl_min
        self.wl_max = wl_max
        self.H1 = RecAtom('H',1)
        self.BE = None
        
    def __make_cont_Ercolano(self, tem, case, wl):
        """
        Adapted from http://adsabs.harvard.edu/abs/2006MNRAS.372.1875E
        """
        n_wl = len(wl)
        hnu =  CST.CLIGHT * 1e8 / wl * CST.HPLANCK  #!phy.c_ang_sec/wl*!phy.h
        if self.BE is None:
            with open(execution_path('../atomic_data/coeff_ercolano.pickle'), 'rb') as handle:
                if sys.version[0] == '2':
                    self.BE = pickle.load(handle)
                else:
                    self.BE = pickle.load(handle, encoding="latin-1")
                
        if case == 'H':
            tab_T = 10**self.BE['th']
            D = self.BE['dh']
        elif case == 'He1':
            tab_T = 10**self.BE['the1']
            D = self.BE['dhe1']
        elif case == 'He2':
            tab_T = 10**self.BE['the2']
            D = self.BE['dhe2']
        else:
            print('Invalid case {0}'.format(case))
            return None
        if (tem < np.min(tab_T)) or (tem > np.max(tab_T)):
            print('Invalid temperature {0}'.format(tem))
            return None
        
        BE_E_Ry = D[:,1]
        BE_E_erg = BE_E_Ry * CST.RYD_erg
        BE_E_Thr = BE_E_erg[D[:,0] == 1]
        Delta_E = np.zeros(n_wl)
        for i in np.arange(n_wl):
            DE = hnu[i] - BE_E_Thr
            Delta_E[i] = np.min(DE[DE > 0])
            
        n_T_sup = np.min(np.where(tab_T >= tem)[0])
        n_T_inf = n_T_sup - 1
        T_sup = tab_T[n_T_sup]
        T_inf = tab_T[n_T_inf]
        
        BE_coeff_sup = D[:, n_T_sup+2]
        BE_coeff_inf = D[:, n_T_inf+2]
        
        coeff_sup = interp1d(BE_E_erg, BE_coeff_sup)(hnu)
        coeff_inf = interp1d(BE_E_erg, BE_coeff_inf)(hnu)

        C_interp= (np.log10(tem) - np.log10(T_inf)) / (np.log10(T_sup) - np.log10(T_inf))
        
        coeff = coeff_sup * C_interp + coeff_inf*(1. - C_interp)
        
        cont = coeff * 1e-34 * tem**(-1.5) * np.exp(-Delta_E / tem / CST.BOLTZMANN) / wl**2. * CST.CLIGHT * 1e8 # erg/s.cm3/A
        return cont

    def __two_photon(self, tem, den, wl):
        y = 1215.7 / wl
        A = 202.0 * (y * (1. - y) * (1. -(4. * y * (1 - y))**0.8) + 0.88 * ( y * (1 - y))**1.53 * (4. * y * (1 - y))**0.8)
        alfa_eff = 0.838e-13 * (tem / 1e4)**(-0.728) # fit DP de Osterbrock
        q = 5.31e-4 * (tem / 1e4)**(-0.17) # fit DP de Osterbrock
        n_crit = 8.226 / q
        twophot_cont = CST.HPLANCK * CST.CLIGHT * 1e8 / wl**3. * 1215.7 * A / 8.226 * alfa_eff / \
                       (1. + den/n_crit)                
        twophot_cont[~np.isfinite(twophot_cont)] = 0.
        return twophot_cont
            
    def gen_continuum(self, tem, den, He1_H, He2_H, wl = None, 
                      cont_HI = True, cont_HeI = True, cont_HeII = True,cont_2p = True):
        """
        
        """
        try:
            _ = (e for e in tem)
            try:
                _ = (e for e in den)
            except:
                den = np.ones_like(tem) * den
        except TypeError:
            tem = [tem]
            den = [den]
            He1_H = [He1_H]
            He2_H = [He2_H]

        if wl is None:
            wl = np.arange(self.wl_min, self.wl_max+1, 1)
        cont = 0.
        if cont_HI:
            cont += np.array(list(map(lambda t:self.__make_cont_Ercolano(tem = t, case = 'H', wl = wl ), tem))).T
        if cont_HeI:
            cont += He1_H * np.array(list(map(lambda t: self.__make_cont_Ercolano(tem = t, case = 'He1', wl = wl ), tem))).T
        if cont_HeII:
            cont += He2_H * np.array(list(map(lambda t: self.__make_cont_Ercolano(tem = t, case = 'He2', wl = wl ), tem))).T
        if cont_2p:
            cont += np.array(list(map(lambda t, n: self.__two_photon(tem = t, den = n, wl = wl), tem, den))).T
        
        return cont.squeeze()
    
    def BJ(self, tem, den, He1_H, He2_H, wl_bbj = 3643, wl_abj = 3861, HI_label='11_2',cont_2p = True):
        
        
        
        fl_bbj, fl_abj = self.gen_continuum(tem = tem, den = den, He1_H = He1_H, 
                                                He2_H = He2_H, wl = np.array([wl_bbj, wl_abj]),
                                                cont_2p = cont_2p)
        
        H_1 = self.H1.getEmissivity(tem, den, label = HI_label)
        BJ_H1 = (fl_bbj - fl_abj) / H_1
        return BJ_H1
    
    