import numpy as np
from scipy import interpolate
import os
import re
import pyneb as pn
if pn.config.Chianti_version_main == '8':
    from . import _chianti_tools_8 as _chianti_tools
elif pn.config.Chianti_version_main == '7':
    from . import _chianti_tools
from . import _chianti_constants as const
from .physics import sym2name, vactoair
from .manage_atomic_data import getLevelsNIST, atom2chianti
from .misc import parseAtom

def Chianti_getA(ion_chianti, NLevels=None):
    """
    Return As from Chianti database.
    Usage:
        As_O3 = Chianti_getA('o_3')
    """
    wgfa = _chianti_tools.wgfaRead(ion_chianti)
    NLevels_ori = np.max(wgfa['lvl2'])
    As = np.zeros((NLevels_ori, NLevels_ori))
    for i, j, a in zip(wgfa['lvl2'], wgfa['lvl1'], wgfa['avalue']):
        As[i-1, j-1] = a
    if NLevels is not None:
        As = As[0:NLevels, 0:NLevels]
    return As

def Chianti_getE(ion_chianti, NLevels=None):
    """
    Return energies from the Chianti database
    Usage:
        E_O3 = Chianti_getE('o_3')
    """
    
    Es = _chianti_tools.elvlcRead(ion_chianti)
    if not 'j' in Es:
        return None
    dtype = [('ecmth', '<f8'), ('ionS', 'S10'), ('term', 'S20'),  ('pretty', 'S30'), 
             ('spd', 'S10'), ('ecm', '<f8'), ('j', '<f8'), ('l', '<i8'), ('erydth', '<f8'), 
             ('conf', '<i8'), ('lvl', '<i8'), ('spin', '<i8'), ('eryd', '<f8'), ('mult', '<i8')]
    arr = np.recarray((len(Es['j']),), dtype=dtype)
    for k in Es.keys():
        if k not in ('status', 'ref', 'label', 'filename'):
            arr[k] = Es[k]
    if NLevels is not None:
        arr = arr[0:NLevels]
    return arr

def get_levs_order(atom, NLevels=None):
    """
    Returns a dictionary giving the correspondence of indices between NIST and Chianti energy level tables.
    
    """
    E_chianti = Chianti_getE(atom2chianti(atom))
    if E_chianti is None:
        return None
    E_NIST = getLevelsNIST(atom)
    if E_NIST is None:
        return None
    # Here follow a list of filter to try to find coincidence between the configuration used in Chanti and in NIST
    remove_stars = lambda str: re.sub('\*','', str)
    remove_par = lambda str: re.sub(r'\([^)]*\)','', str)
    change_2pto1p = lambda str: re.sub('\.\.', '.', str)
    change_2pto0p = lambda str: re.sub('\.\.', '', str)
    remove_point_before_par = lambda str: re.sub('\.\(', '(', str)
    remove_before_firstp = lambda str: re.sub('^[^.]*.', '', str)
    change_ptpsp = lambda str: re.sub('\.', ' ', str)
    # Define the Dictionary
    Chianti2NIST = {}
    for i_N, E in enumerate(E_NIST):
        if i_N >= len(E_chianti):
            break
        if NLevels is not None and i_N > NLevels:
            break
        
        term = remove_stars(E['term'])

        c1 = remove_stars(E['conf'])
        
        c12 = remove_par(c1)
        c123 = change_2pto1p(c12)
        c124 = change_2pto0p(c12)
        
        c15 = remove_point_before_par(c1) 
        
        c16 = remove_before_firstp(c1)
        c126 = remove_before_firstp(c12)
        c1236 = remove_before_firstp(c123)
        c1246 = remove_before_firstp(c124)
        c156 = remove_before_firstp(c15)
        
        c17 = change_ptpsp(c1)
        c127 = change_ptpsp(c12)
        c1237 = change_ptpsp(c123)
        c1247 = change_ptpsp(c124)
        c157 = change_ptpsp(c15)
        c167 = change_ptpsp(c16)
        c1267 = change_ptpsp(c126)
        c12367 = change_ptpsp(c1236)
        c12467 = change_ptpsp(c1246)
        c1567 = change_ptpsp(c156)
       
        confs = (c1, c12, c123, c124, c15, c16, c126, c1236, c1246, c156,
                 c17, c127, c1237, c1247, c157, c167, c1267, c12367, c12467, c1567)
        for conf in confs:   
            pretty = '{0} {1}{2}'.format(conf, term, E['J'])
            i_Ch = np.where(E_chianti['pretty'] == pretty)[0]
            if len(i_Ch) == 1:
                if NLevels is not None and i_Ch[0] <= NLevels:
                    Chianti2NIST[i_N] = i_Ch[0]
    if len(Chianti2NIST) == 0:
        Chianti2NIST = None
        
        ### Il faut que ni i_N ni i_Ch depasse NLvels. Et il faut que NLevels soit diminue dansla cas contraire.
    return Chianti2NIST
        
def Chianti_getOmega(ion_chianti, tem, lev1=None, lev2=None, Splups=None, NLevels=None):
    """
    Return the values of Upsilon at a given temperature (may be a table)
    Usage:
        Chianti_getOmega('o_3', 10000)
        Chianti_getOmega('o_3', np.linspace(5000, 20000, 20))
    lower and upper levels can be given as lev1 and lev2. If not, full table is return
    Splups is the result of _chianti_tools.splupsRead(ion_chianti), can be passed to avoid reloading it
    NLevelsMax reduce the number of levels considered.
    """
    temp=np.asarray(tem)
    if Splups is None:
        if pn.config.Chianti_version_main == '8':
            Splups = _chianti_tools.scupsRead(ion_chianti)
            print(Splups.keys())
        elif pn.config.Chianti_version_main == '7':
            Splups = _chianti_tools.splupsRead(ion_chianti)
    if NLevels is None:
        nsplups = len(Splups["lvl1"])
        nlevels = np.max(Splups["lvl2"])
    else:
        nsplups = np.min((len(Splups["lvl1"]), NLevels))
        nlevels = np.min((np.max(Splups["lvl2"]), NLevels))
    ntemp=temp.size
    if (lev1 is None) and (lev2 is None):
        Omega = np.zeros((nlevels, nlevels, ntemp))
    else:
        Omega = np.zeros(ntemp)
    lvl1 = np.array(Splups["lvl1"])
    lvl2 = np.array(Splups["lvl2"])
    if lev1 is None:
        range1 = np.arange(nlevels)
    else:
        range1 = np.array([lev1-1])
    if lev2 is None:
        range2 = np.arange(nlevels)
    else:
        range2 = np.array([lev2-1])
    for l1 in range1:
        for l2 in range2:
            isplups = np.where(((lvl1 == (l1+1)) & (lvl2 == (l2+1))))
            if type(isplups) == type(()):
                isplups = isplups[0]                
            if len(isplups) == 0:
                ups = np.zeros(ntemp)
            else:
                if len(isplups) > 1 or np.ndim(isplups)==1:
                    isplups = isplups[0]
                ttype=Splups["ttype"][isplups]
                cups=Splups["cups"][isplups]
                if pn.config.Chianti_version_main == '8':
                    nspl=Splups["ntemp"][isplups]
                    dx=1./(float(nspl)-1.)
                    xs = Splups['btemp'][isplups]
                    splups=Splups["bscups"][isplups]
                elif pn.config.Chianti_version_main == '7':
                    nspl=Splups["nspl"][isplups]
                    dx=1./(float(nspl)-1.)
                    xs=dx*np.arange(nspl)
                    splups=Splups["splups"][isplups]
                de=Splups['de'][isplups]
                kte = const.boltzmann*temp/(de*const.ryd2erg)
                der=0
                if ttype == 1:
                    st=1.-np.log(cups)/np.log(kte+cups)
                    y2=interpolate.splrep(xs,splups,s=0)
                    sups=interpolate.splev(st,y2,der=der)
                    ups=sups*np.log(kte+np.exp(1.))
                elif ttype == 2:
                    st=kte/(kte+cups)
                    y2=interpolate.splrep(xs,splups,s=0)
                    sups=interpolate.splev(st,y2,der=der)
                    ups=sups
                elif ttype == 3:
                    st=kte/(kte+cups)
                    y2=interpolate.splrep(xs,splups,s=0)
                    sups=interpolate.splev(st,y2,der=der)
                    ups=sups/(kte+1.)
                elif ttype == 4:
                    st=1.-np.log(cups)/np.log(kte+cups)
                    y2=interpolate.splrep(xs,splups,s=0)
                    sups=interpolate.splev(st,y2,der=der)
                    ups=sups*np.log(kte+cups)
                elif ttype == 5:
                    # dielectronic rates
                    st=kte/(kte+cups)
                    y2=interpolate.splrep(xs,splups,s=0)
                    sups=interpolate.splev(st,y2,der=der)
                    ups=sups/(kte+0.)
                elif ttype == 6:
                    st=kte/(kte+cups)
                    y2=interpolate.splrep(xs,splups,s=0)
                    sups=interpolate.splev(st,y2,der=der)
                    ups=10.**sups
                elif ttype > 6:
                    pn.log_.warning(' t_type ne 1,2,3,4,5={} {} {}'.format(ttype,l1,l2), calling='pn.chianti.Chianti_getOmega')
            if (lev1 is None) and (lev2 is None):
                Omega[l2, l1, :] = ups
            else:
                Omega = ups

    Omega = np.where(Omega > 0.,Omega,0.)
    return Omega

class _AtomChianti(object):
    
    def __init__(self, elem=None, spec=None, atom=None, NLevels=None):
        """
        Object dealing with As values from the Chianti database.
        The directory where to find the data must be given through the environment variable XUVTOP.
        This object is not aimed to be used by itself, it is called by Atom.
        """
        self.log_ = pn.log_
        if atom is not None:
            self.atom = atom
            self.elem = parseAtom(atom)[0]
            self.spec = int(parseAtom(atom)[1])
        else:
            self.elem = elem
            self.spec = int(spec)
            self.atom = elem + str(self.spec)
        self.name = sym2name[self.elem]
        self.ion_chianti =  atom2chianti(self.atom)
        self.calling = 'Atom {}'.format(self.atom)
        self.Chianti2NIST = None
        self.NLevels = NLevels
        self._loadChianti()
        self.initWaves()
        
    def _loadChianti(self):
        self.fullFileName = _chianti_tools.ion2filename(self.ion_chianti) + '.wgfa'
        
        self.atomFile = self.fullFileName.split('/')[-1]
        self.atomPath = self.fullFileName[:-len(self.atomFile)]
        
        if not os.path.exists(self.fullFileName):
            self.log_.error('File {0} not found'.format(self.fullFileName), calling=self.calling)
        
        self.Chianti_version = self.atomPath.split('/')[-4]
        self.comments = {}
        
        Chianti_A = Chianti_getA(self.ion_chianti, NLevels=self.NLevels)
        if self.Chianti2NIST is None:
            self.Chianti2NIST = get_levs_order(self.atom)
        if self.Chianti2NIST is not None:
            Chianti_A_tmp = Chianti_A.copy()
            for i_chianti in self.Chianti2NIST:
                if self.Chianti2NIST[i_chianti] != i_chianti:
                    Chianti_A[i_chianti,:] = Chianti_A_tmp[self.Chianti2NIST[i_chianti],:]
                    Chianti_A[:,i_chianti] = Chianti_A_tmp[:,self.Chianti2NIST[i_chianti]]
        self.log_.message('Reading atom data from Chianti {}'.format(self.atomFile), calling = self.calling)
        self.Chianti_E = Chianti_getE(self.ion_chianti, NLevels=self.NLevels)
        if self.NLevels is None:
            self.NLevels = len(self.Chianti_E)
        self.NIST = getLevelsNIST(self.atom, self.NLevels)
        need_NIST = True
        
        if self.NIST is not None:
            self.NLevels = np.min((len(self.NIST), self.NLevels))    
            energy = self.NIST['energy'] / 1e8
            stat_weight = 1 + 2 * self.NIST['J']
            self.comments['VACUUM'] = '1'
            self.comments['NOTE'] = 'Energy levels'
            self.comments['SOURCE'] = 'NIST 2014'
        elif need_NIST:
            pn.log_.error('NIST data are needed for this format of atomic data', calling=self.calling) 

        self.E_in_vacuum = True
        self._Energy = energy
        self._StatWeight = stat_weight
        self._A = Chianti_A
        self.atomNLevels = self.NLevels
            
    def initWaves(self):
        """
        Initialization of wave_Ang
        
        """
        self.wave_Ang = np.zeros((self.NLevels, self.NLevels))
        
        for i in range(1, self.NLevels):
            for j in range(i):
                wave = 1. / abs(self._Energy[i] - self._Energy[j])
                if self.E_in_vacuum:
                    wave = vactoair(wave)
                self.wave_Ang[i, j] = self.wave_Ang[j, i] = wave
 
    def _test_lev(self, level):
        """
        Test whether selected level is legal

        Parameters:
            - level        selected atom level

        """       
        if level < -1 or level == 0 or level > self.NLevels:
            self.log_.error('Wrong value for level: {0}, maximum = {1}'.format(level, self.NLevels),
                            calling=self.calling)
    
    def getA(self, lev_i= -1, lev_j= -1):
        """
        Return the transition probability data. 
        If no arguments are given, the whole array of A is returned.
        A specific A value can be obtained by giving either the upper and lower levels.
            
        Usage:
            A_O3 = O3.getA()          # The whole A array is stored in A_O3
            O3.getA(4, 2)      # A(4, 2) of the O3 atom is printed
            O3.getA(2, 4)      # Returns 0

        Parameters:
            - lev_i  upper level of transition (default= -1, returns complete array)
            - lev_j  lower level of transition (default= -1, returns complete array)
            
        """
        
        self._test_lev(lev_i)
        self._test_lev(lev_j)
        if (lev_i == -1):
            if (lev_j == -1):
                return self._A
            else:
                return (self._A[lev_i - 1])
        else:
            return (self._A[lev_i - 1, lev_j - 1])

    def getStatWeight(self, level= -1):
        """
        Returns the array of statistical weights of the ion (if no arguments are given) 
            or the statistical weight of level i (if i is given as an argument, 
            with the optional keyword level).
            
        Usage:
            O3.getStatWeight()
            O3.getStatWeight(level=4)
            O3.getStatWeight(4)
                      
        Parameters:
            - level  selected atomic level (default= -1, returns complete array)
            
        """
        self._test_lev(level)
        
        if level == -1:
            return self._StatWeight
        else:
            return self._StatWeight[level-1]
        
    
    def getEnergy(self, level= -1, unit='1/Ang'):
        """
        Return energy level of selected level (or array of energy levels, if level not given) 
            in Angstrom^-1 (default) or another unit
        
        Usage:
            O3.getEnergy(4, unit='eV')
        Parameters:
            - level  selected atomic level (default= -1, returns complete array)
            - unit   [str] one of '1/Ang' (default), 'eV', or 'Ryd'    
            
        """
        self._test_lev(level)

        unit_dict = {'1/Ang': 1.,
                     'Ryd': pn.CST.RYD_ANG,
                     'eV': pn.CST.RYD_ANG * pn.CST.RYD_EV,
                     'cm-1': 1e8}
        if unit not in unit_dict:
            self.log_.warn('Unit {0} unknown, using 1/Ang'.format(unit), calling=self.calling + '.getEnergy')
            unit = '1/Ang'        
        
        if level == -1:
            return self._Energy * unit_dict[unit]
        else:
            return self._Energy[level-1] * unit_dict[unit]
        
    def getSources(self):
        """
        Return bibliographic sources for atomic data, as listed in the headers of the fits files
        
        """
        sources = []
        sources.append('A-values from {0}'.format(self.Chianti_version))
        for ref in  _chianti_tools.elvlcRead(self.ion_chianti)['ref']:
            sources.append(ref)
        return sources

    def printSources(self):
        """
        Print bibliographic sources for atomic data, as listed in the headers of the fits files
        
        """
        print('A-values from {0}'.format(self.Chianti_version))
        for ref in  _chianti_tools.elvlcRead(self.ion_chianti)['ref']:
            print(ref)
               
class _CollChianti(object):
    
    def __init__(self, elem=None, spec=None, atom=None, NLevels=None, TemArray=np.logspace(2, 5, 20)):
        """
        Object dealing with Upsilon values from the Chianti database.
        The directory where to find the data must be given through the environment variable XUVTOP.
        This object is not aimed to be used by itself, it is called by Atom.
        """

        self.log_ = pn.log_
        
        if atom is not None:
            self.atom = atom
            self.elem = parseAtom(atom)[0]
            self.spec = int(parseAtom(atom)[1])
        else:
            self.elem = elem
            self.spec = int(spec)
            self.atom = elem + str(self.spec)
        self.name = sym2name[self.elem]
        self.ion_chianti =  atom2chianti(self.atom)
        self.calling = 'Atom ' + self.atom
        self.NLevels = NLevels
        self.tem_units = 'K'
        self._TemArray = TemArray
        self._loadChianti()
        
    def _loadChianti(self):
        if pn.config.Chianti_version_main == '8':
            self.fullFileName = _chianti_tools.ion2filename(self.ion_chianti) + '.scups'
        elif pn.config.Chianti_version_main == '7':
            self.fullFileName = _chianti_tools.ion2filename(self.ion_chianti) + '.splups'
        self.collFile = self.fullFileName.split('/')[-1]
        self.collPath = self.fullFileName[:-len(self.collFile)]
        
        if not os.path.exists(self.fullFileName):
            self.log_.error('File {0} not found'.format(self.fullFileName), calling=self.calling)
        
        self.Chianti_version = self.collPath.split('/')[-4]
        self.comments = {}
        if pn.config.Chianti_version_main == '8':
            self.Splups = _chianti_tools.scupsRead(self.ion_chianti)
        elif pn.config.Chianti_version_main == '7':
            self.Splups = _chianti_tools.splupsRead(self.ion_chianti)
        self._CollArray = Chianti_getOmega(self.ion_chianti, tem=self._TemArray, 
                                   Splups=self.Splups, NLevels=self.NLevels)
        
        self.NLevels = np.min((self._CollArray.shape[0]-1, self._CollArray.shape[1]-1, self.NLevels))
        self.Chianti2NIST = get_levs_order(self.atom, NLevels=self.NLevels)
        
        if self.Chianti2NIST is not None:
            _CollArray_tmp = self._CollArray.copy()
            for i_chianti in self.Chianti2NIST:
                if self.Chianti2NIST[i_chianti] != i_chianti:
                    self._CollArray[i_chianti,:,:] = _CollArray_tmp[self.Chianti2NIST[i_chianti],:,:]
                    self._CollArray[:,i_chianti,:] = _CollArray_tmp[:,self.Chianti2NIST[i_chianti],:]
        
        self.log_.message('Reading coll data from Chianti {}'.format(self.collFile), calling = self.calling)
        self.NLevels = self._CollArray.shape[0]


    def _test_lev(self, level):
        """
        Test whether selected level is legal

        Parameters:
            - level        selected atom level

        """       
        if level < -1 or level == 0 or level > self.NLevels:
            self.log_.error('Wrong value for level: {0} (maximum: {1})'.format(level, self.NLevels),
                            calling=self.calling)

    def getOmegaArray(self, lev_i= -1, lev_j= -1):
        """
        Return array of original tabulated collision strengths for a given transition, 
            as a function of temperature.
        
        Usage:
            O3.getOmegaArray()
        
        Parameters:
            - lev_j  lower level (default= -1, returns complete array)
            - lev_i  upper level (default= -1, returns complete array)

        """
        self._test_lev(lev_i)
        self._test_lev(lev_j)
        if (lev_i <= lev_j) and (lev_i != -1):
            self.log_.warn("wrong levels given {0} <= {1}".format(lev_i, lev_j), calling=self.calling)
            return None
        elif lev_i == -1 or lev_j == -1:
            return self._CollArray
        else:
            return self._CollArray[lev_i-1, lev_j-1,:]

    def getOmega(self, tem, lev_i= -1, lev_j= -1):
        """
        Return interpolated value of the collision strength value at the given temperature 
            for the complete array or a specified transition.

        Usage:
            O3.getOmega(15000.)
            O3.getOmega([8e3, 1e4, 1.2e4])
            O3.getOmega([8e3, 1e4, 1.2e4], 5, 4)
        
        Parameters:
            - tem    electronic temperature in K. May be an array.
            - lev_i  upper level
            - lev_j  lower level

        """
        self._test_lev(lev_i)
        self._test_lev(lev_j)
        if (lev_i <= lev_j) and (lev_i != -1):
            self.log_.warn("wrong levels given {0} <= {1}".format(lev_i, lev_j), calling=self.calling)
            return None            
        if self.Chianti2NIST is None:
            self.Chianti2NIST = get_levs_order(self.atom)
        if (lev_i == -1) and (lev_j == -1):
            Omega = Chianti_getOmega(self.ion_chianti, tem, Splups=self.Splups, NLevels=self.NLevels)
            if self.Chianti2NIST is not None:
                Omega_tmp = Omega.copy()
                for i_chianti in self.Chianti2NIST:
                    if self.Chianti2NIST[i_chianti] != i_chianti:
                        Omega[i_chianti,:,:] = Omega_tmp[self.Chianti2NIST[i_chianti],:,:]
                        Omega[:,i_chianti,:] = Omega_tmp[:,self.Chianti2NIST[i_chianti],:]
        else:
            if self.Chianti2NIST is not None:
                try:
                    Omega = Chianti_getOmega(self.ion_chianti, tem, lev1=self.Chianti2NIST[lev_j-1]+1, lev2=self.Chianti2NIST[lev_i-1]+1, 
                                             Splups=self.Splups, NLevels=self.NLevels)
                except:
                    Omega = 0.
            else:
                Omega = Chianti_getOmega(self.ion_chianti, tem, lev1=lev_j, lev2=lev_i, 
                             Splups=self.Splups, NLevels=self.NLevels)
                
        return np.squeeze(Omega)

    def getTemArray(self, keep_unit=True):
        """
        Return array of tabulated original temperature points (as in fits file) 
            of collision strengths.
        
        Parameters:
            - keep_unit   return temperature in file units (default) or change it to Kelvin (False)
        """
        if keep_unit:
            return self._TemArray
        else:            
            if (self.tem_units == "log(K)"):
                return pow(10., self._TemArray)
            else:
                if (self.tem_units == "K/10000"):
                    return self._TemArray * 1.e4
                else: 
                    return self._TemArray
    
    def getSources(self):
        """
        Return bibliographic sources for atomic data, as listed in the headers of the fits files
        
        """
        sources = []
        sources.append('Omega-values from {0}'.format(self.Chianti_version))
        for ref in  _chianti_tools.splupsRead(self.ion_chianti)['ref']:
            sources.append(ref)
        return sources

    def printSources(self):
        
        for source in self.getSources():
            print(source)    
        
