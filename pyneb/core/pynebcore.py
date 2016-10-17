"""@module pynebcore
Main PyNeb module
Tools to manage atoms, emission lines and observational data.

@class Atom           atom object
@class EmissionLine   emission line object
@class Observation    observation object

"""
from __future__ import print_function
import numpy as np
import warnings
import os

from pyneb import config, log_, atomicData
from ..utils.misc import int_to_roman, strExtract, parseAtom, quiet_divide, _returnNone, solve
from ..utils import chebyshev
from ..utils.init import LINE_LABEL_LIST, BLEND_LIST, SPEC_LIST, label2levelDict
from ..utils.physics import sym2name, gsLevelDict, gsFromAtom, vactoair, CST, Z, IP
from ..utils.manage_atomic_data import getLevelsNIST
from ..utils.pn_chianti import _AtomChianti, _CollChianti
from ..extinction.red_corr import RedCorr
from fractions import Fraction
if config.INSTALLED['scipy']:
    from scipy import interpolate
    from scipy.special import gamma
if config.INSTALLED['plt']: 
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.collections import LineCollection
    from matplotlib import colors
if config.INSTALLED['mp']:
    from multiprocessing import Queue, Process
    from ..utils.multiprocs import getTemDen_helper
if config.INSTALLED['pyfits']:
    import pyfits
elif config.INSTALLED['pyfits from astropy']:
    import astropy.io.fits as pyfits
if config.INSTALLED['h5py']:
    import h5py


# Change the profiler to 'cpu', 'mem' or None to profile the execution of Atom.
profiler = None
#profiler = 'cpu'
if profiler == 'mem':
    try:
        from memory_profiler import profile
    except:
        def profile(f):
            return f
elif profiler is None:
    def profile(f):
        return f

class _AtomDataNone(object):
    
    def __init__(self):

        self.atomFile = _returnNone()
        self.atomPath = _returnNone()
        self.atomFitsFile = _returnNone()
        self.atomFitsPath = _returnNone()
        self.wave_Ang = np.nan
        self.getStatWeight = _returnNone
        self.getEnergy = _returnNone
        self.getA = _returnNone
        self.atomNLevels = 0
        self.NLevels = 0

class _CollDataNone(object):
    
    def __init__(self):
        self.getOmegaArray = _returnNone
        self.getTemArray = _returnNone
        self.collFile = _returnNone()
        self.collPath = _returnNone()
        self.collFitsFile = _returnNone()
        self.collFitsPath = _returnNone()
        self.collNLevels = 0
        self.tem_units = _returnNone()
        self.NLevels = 0

class _AtomDataFits(object):
    
    def __init__(self, elem=None, spec=None, atom=None, NLevels=None):
        self.log_ = log_
        if atom is not None:
            self.atom = atom
            self.elem = parseAtom(atom)[0]
            self.spec = int(parseAtom(atom)[1])
        else:
            self.elem = elem
            self.spec = int(spec)
            self.atom = elem + str(self.spec)
        self.name = sym2name[self.elem]
        self.calling = 'Atom ' + self.atom
        self.NLevels = NLevels
        self._loadFit()
        self._A = self.getA() # index = quantum number - 1
        self._Energy = self.getEnergy() #Angstrom^-1
        self._StatWeight = self.getStatWeight()
        self.initWaves()
        
        
    def _loadFit(self):
        """
        Load fits file of atomic data.
        The fits files names are defined in the atomicData Object

        """                 
        self.atomFile = atomicData.getDataFile(self.atom, data_type='atom')
        self.atomFitsFile = self.atomFile
        if self.atomFitsFile is None:
            self.log_.error('No atom data for ion {0}'.format(self.atom), calling=self.calling)
            return None
        if self.atom in config.DataFiles:
            if self.atomFitsFile not in config.DataFiles[self.atom]:
                config.DataFiles[self.atom].append(self.atomFitsFile)
        else:
            config.DataFiles[self.atom] = [self.atomFitsFile]
            
        self.atomPath = atomicData.getDirForFile(self.atomFile)
        self.atomFitsPath = self.atomPath
        file_to_open = '{0}/{1}'.format(self.atomPath, self.atomFile)
        if not os.path.exists(file_to_open):
            self.log_.error('File {0} not found'.format(file_to_open), calling=self.calling)
        atomFits = pyfits.open(file_to_open, ignore_missing_end=True)
        self.log_.message('Reading atom data from {0}'.format(self.atomFile), calling=self.calling)
        
        AtomExt = atomFits[1]


        #Read headers
        self.AtomHeader = AtomExt.header

        self.gs = 'unknown'
        try:
            self.gs = self.AtomHeader['GSCONFIG']
        except:
            pass
        if 'SPECTRUM' in self.AtomHeader:
            if int(self.AtomHeader['SPECTRUM']) != self.spec:
                log_.error('The spectrum I read in the file {0} is {1}, but you are requesting {2}'.format(self.atomFitsFile, self.AtomHeader['SPECTRUM'],
                                                                                                    self.spec), calling=self.calling)
        if 'ATOM' in self.AtomHeader:
            if self.AtomHeader['ATOM'] != sym2name[self.elem]:
                log_.error('The element name I read in the file {0} is {1}, but I was expecting {2}. Check the keyword ATOM'.format(self.atomFitsFile, self.AtomHeader['ATOM'],
                                                                                                    sym2name[self.elem]), calling=self.calling)
                
        #Read data
        self._AtomData = AtomExt.data

        try:
            self.atomNLevels = self.AtomHeader['N_LEVELS']
        except:
            log_.error('N_LEVELS is not set in {0}'.format(self.atomFitsFile))
            
        self.log_.message('NLevels of atomic data: {0}'.format(self.atomNLevels),
                          calling=self.calling)
        
        self.NLevels = self.atomNLevels
        
        atomFits.close()


    def initWaves(self):
        """
        Initialization of wave_Ang 
        
        """
        self.wave_Ang = np.zeros((self.NLevels, self.NLevels))
        for i in range(1, self.NLevels):
            for j in range(i):
                self.wave_Ang[i, j] = self.wave_Ang[j, i] = 1. / abs(self._Energy[i] - self._Energy[j])

        
    def _test_lev(self, level):
        """
        Test whether selected level is legal

        Parameters:
            - level        selected atom level

        """       
        if level < -1 or level == 0 or level > self.NLevels:
            self.log_.error('Wrong value for level: {0}, maximum = {1}'.format(level, self.NLevels),
                            calling=self.calling)
        
    def getSources(self):
        """
        Return bibliographic sources for atomic data, as listed in the headers of the fits files
        """
        sources = [] 
        header = self.AtomHeader
        for i in range(len(header.items())):
            if 'SOURCE' in header.items()[i][0]:
                number = header.items()[i][0].lstrip('SOURCE')
                try:
                    sources.append(self.atom + ': ' + header.get('NOTE' + str(number)) + ':', header.get('SOURCE' + str(number)))
                except:
                    sources.append(self.atom + ': ' + 'Atomic data:', header.get('SOURCE' + str(number)))
        return sources
    
    
    def printSources(self):
        
        for source in self.getSources():
            print(source)    

    
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
                # Line below commented out because introduces bug if Nlevels(atom) > Nlevels(coll)
                #resultArray = np.zeros([self.NLevels, self.NLevels])
                resultArray = np.zeros([self.atomNLevels, self.atomNLevels])
                for i in range(self.atomNLevels - 1):
                    lev_i = i + 1
                    j = i + 1
                    while (j < self.atomNLevels):
                        lev_j = j + 1
                        resultArray[j][i] = self.getA(lev_j, lev_i)
                        j += 1
                return resultArray

            return (self._AtomData.field('A(*->{0})'.format(lev_j)))
        else:
            return (self._AtomData.field('A(*->{0})'.format(lev_j))[lev_i - 1])


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
        if (level == -1):
            return ((self._AtomData.field('Stat_Weight')) * 1.0)
        elif (level > self.NLevels  or level < 1):
            self.log_.warn('level {0} outside the [1-{1}] range: aborting'.format(level, self.NLevels)
                               , calling=self.calling)
#                level = level % self.NLevels 
        else:
            return ((self._AtomData.field('Stat_Weight')[level - 1]) * 1.0)
    

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
                     'Ryd': CST.RYD_ANG,
                     'eV': CST.RYD_ANG * CST.RYD_EV,
                     'erg': CST.HPLANCK * CST.CLIGHT * 1e8}
        if unit not in unit_dict:
            self.log_.warn('Unit {0} unknown, using 1/Ang'.format(unit), calling=self.calling + '.getEnergy')
            unit = '1/Ang'
        if (level == -1):
            return self._AtomData.field('Energy') * unit_dict[unit]
        elif (level > self.NLevels  or level < 1):
            self.log_.warn('level {0} outside the [1-{1}] range: aborting'.format(level, self.NLevels)
                               , calling=self.calling)
#                level = level % self.NLevels 
        else:
            return self._AtomData.field('Energy')[level - 1] * unit_dict[unit]


class _AtomDataAscii(object):
    
    def __init__(self, elem=None, spec=None, atom=None, NLevels=None):
        self.log_ = log_
        if atom is not None:
            self.atom = atom
            self.elem = parseAtom(atom)[0]
            self.spec = int(parseAtom(atom)[1])
        else:
            self.elem = elem
            self.spec = int(spec)
            self.atom = elem + str(self.spec)
        self.name = sym2name[self.elem]
        self.calling = 'Atom ' + self.atom
        self.NLevels = NLevels
        self._loadAscii()
        self.initWaves()
        
    def _loadAscii(self):
        
        self.atomFile = atomicData.getDataFile(self.atom, data_type='atom')
        if self.atomFile is None:
            self.log_.error('No atom data for ion {0}'.format(self.atom), calling=self.calling)
            return None
        if self.atom in config.DataFiles:
            if self.atomFile not in config.DataFiles[self.atom]:
                config.DataFiles[self.atom].append(self.atomFile)
        else:
            config.DataFiles[self.atom] = [self.atomFile]
            
        self.atomPath = atomicData.getDirForFile(self.atomFile)
        file_to_open = '{0}/{1}'.format(self.atomPath, self.atomFile)
        if not os.path.exists(file_to_open):
            self.log_.error('File {0} not found'.format(file_to_open), calling=self.calling)
            
        self.log_.message('Reading atom data from {0}'.format(self.atomFile), calling=self.calling)

    
        self.E_in_vacuum = False
        # Read data from ascii file
        #
        # Read energies and stat weights 
        f = open(file_to_open)
        data = f.readlines()
        f.close()
        if data[0].strip() == 'Aij':
            # This means that we are dealing with new format (no more energies nor stats weights
            # so NIST data are needed
            need_NIST = True
            at_data = np.array([d.split() for d in data if d[0:3]!='***'][2::], dtype='float')
            comments_tab = [d for d in data if d[0:3]=='***']
            self.comments = {}
            for com in comments_tab:
                key = com.split()[1]
                self.comments[key] = com.split(key)[1].strip()
            A = at_data.copy()
            if A.shape[0] != A.shape[1]:
                log_.error('Atomic data must be a NxN matrix', calling=self.calling) 
            if self.NLevels is not None:
                A = A[0:self.NLevels, 0:self.NLevels]
            self.NLevels = A.shape[0]
        else:
            # This is the old format
            need_NIST = False
            units = data[1].split()[0]
            at_data = np.array([d.split() for d in data if d[0:3]!='***'][2::], dtype='float')
            comments_tab = [d for d in data if d[0:3]=='***']
            self.comments = {}
            for com in comments_tab:
                key = com.split()[1]
                self.comments[key] = com.split(key)[1].strip()
            A = at_data.copy()
            if self.NLevels is not None:
                A = at_data[0:self.NLevels, 2:self.NLevels+2]
            self.NLevels = A.shape[0]
            
            
            # Read Es
            energy = at_data[:,0]
            if units == 'eV': 
                energy /= CST.RYD_EV * CST.RYD_ANG
            elif units == 'Rydberg':
                energy /= CST.RYD_ANG
            elif units == 'cm-1':
                energy /= 1e8
            # Read statistical weights
            stat_weight = at_data[:,1]
            
            # Read As
            A = np.zeros([self.NLevels, self.NLevels])
            A[:,:] = at_data[:,2:]

        
        self.NIST = getLevelsNIST(self.atom, self.NLevels)
        
        web = 'Ref. {0} of NIST 2014 (try this: http://physics.nist.gov/cgi-bin/ASBib1/get_ASBib_ref.cgi?db=el&db_id={0}&comment_code=&element={1}&spectr_charge={2}&'
        if self.NIST is not None:
            energy = self.NIST['energy'] / 1e8
            stat_weight = 1 + 2 * self.NIST['J']
            self.comments['VACUUM'] = '1'
            self.comments['NOTE'] = 'Energy levels'
            source = '\n    '
            for ref in np.unique(self.NIST['ref']):
                source = source + web.format(ref[1:], self.elem, self.spec) + ')\n  + ' 
            self.comments['SOURCE'] = source[0:-5]
        elif need_NIST:
            log_.error('NIST data are needed for this format of atomic data', calling=self.calling) 
        
        self._Energy = energy
        self._StatWeight = stat_weight
        self._A = A
        self.atomNLevels = self.NLevels
        if 'GSCONFIG' in self.comments:
            self.gs = self.comments['GSCONFIG']
        else:
            self.gs = 'unknown'
        if 'SPECTRUM' in self.comments:
            if int(self.comments['SPECTRUM']) != self.spec:
                log_.error('The spectrum I read in the file {0} is {1}, but you are requesting {2}'.format(self.atomFile, self.comments['SPECTRUM'],
                                                                                                    self.spec), calling=self.calling)
        if 'ATOM' in self.comments:
            if self.comments['ATOM'] != sym2name[self.elem]:
                log_.error('The element name I read in the file {0} is {1}, but I was expecting {2}. Check the keyword ATOM'.format(self.atomFile, self.comments['ATOM'],
                                                                                                    sym2name[self.elem]), calling=self.calling)
        if 'VACUUM' in self.comments:
            if self.comments['VACUUM'] == '1':
                self.E_in_vacuum = True
    
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

    
    def getSources(self):
        """
        Return bibliographic sources for atomic data, as listed in the headers of the fits files
        
        """
        sources = []
        sources_str = {}
        notes = {}
        for com in self.comments:
            if 'SOURCE' in com:
                number = com.split('SOURCE')[1]
                sources_str[number] = self.comments[com]
            if 'NOTE' in com:
                number = com.split('NOTE')[1]
                notes[number] = self.comments[com]
        for key in sources_str:
            try:
                sources.append('{0}: {1}: {2}'.format(self.atom, notes[key], sources_str[key])) 
            except:
                sources.append('{0}: {1}'.format(self.atom, sources[key]))
        return sources

    def printSources(self):
        
        for source in self.getSources():
            print(source)    
    
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
                     'Ryd': CST.RYD_ANG,
                     'eV': CST.RYD_ANG * CST.RYD_EV,
                     'cm-1': 1e8,
                     'erg': CST.HPLANCK * CST.CLIGHT * 1e8}
        if unit not in unit_dict:
            self.log_.warn('Unit {0} unknown, using 1/Ang'.format(unit), calling=self.calling + '.getEnergy')
            unit = '1/Ang'        
        
        if level == -1:
            return self._Energy * unit_dict[unit]
        else:
            return self._Energy[level-1] * unit_dict[unit]
        

class _AtomDataStout(object):
    
    def __init__(self):
        pass
        
class _CollDataFits(object):
    
    def __init__(self, elem=None, spec=None, atom=None, OmegaInterp='Cheb', noExtrapol = False, NLevels=None):
        self.log_ = log_
        if atom is not None:
            self.atom = atom
            self.elem = parseAtom(atom)[0]
            self.spec = int(parseAtom(atom)[1])
        else:
            self.elem = elem
            self.spec = int(spec)
            self.atom = elem + str(self.spec)
        self.name = sym2name[self.elem]
        self.noExtrapol = noExtrapol
        self.calling = 'Atom ' + self.atom
        self.NLevels = NLevels
        self._loadFit()
        self.initOmegas(OmegaInterp=OmegaInterp)
        self.tem_units = self.CollHeader['TUNIT1']
        if (self.tem_units == "log(K)"):
            def transfo_tem(x):
                y = np.log10(x)
                if np.ndim(y) == 0:
                    y = np.array(y)
                return y
        elif (self.tem_units == "K/10000"):
                def transfo_tem(x):
                    y = x / 1e4
                    if np.ndim(y) == 0:
                        y = np.array(y)
                    return y
        else: #T in Kelvin in the fits file
            def transfo_tem(x):
                y = x
                return x
        self._transfo_tem = transfo_tem
        
        
    def _loadFit(self):
        """
        Load fits file of atomic data.
        The fits files names are defined in the atomicData Object

        """                 
        self.collFile = atomicData.getDataFile(self.atom, data_type='coll')
        self.collFitsFile = self.collFile 
        if self.collFile is None:
            self.log_.error('No coll data for ion {0}'.format(self.atom), calling=self.calling)
            return None
        
        if self.atom in config.DataFiles:
            if self.collFile not in config.DataFiles[self.atom]:
                config.DataFiles[self.atom].append(self.collFile)
        else:
            config.DataFiles[self.atom] = [self.collFitsFile]
            
        self.collPath = atomicData.getDirForFile(self.collFile)
        self.collFitsPath = self.collPath 
        file_to_open = '{0}/{1}'.format(self.collPath, self.collFile)
        if not os.path.exists(file_to_open):
            self.log_.error('File {0} not found'.format(file_to_open), calling=self.calling)        
        collFits = pyfits.open(file_to_open, ignore_missing_end=True)
        self.log_.message('Reading coll data from {0}'.format(self.collFile), calling=self.calling)
        
        CollExt = collFits[1]

        #Read headers
        self.CollHeader = CollExt.header
        self.CollExtNames = CollExt.columns.names

        if 'SPECTRUM' in self.CollHeader:
            if int(self.CollHeader['SPECTRUM']) != self.spec:
                log_.error('The spectrum I read in the file {0} is {1}, but you are requesting {2}'.format(self.collFitsFile, self.CollHeader['SPECTRUM'],
                                                                                                    self.spec), calling=self.calling)
        if 'ATOM' in self.CollHeader:
            if self.CollHeader['ATOM'] != sym2name[self.elem]:
                log_.error('The element name I read in the file {0} is {1}, but I was expecting {2}. Check the keyword ATOM'.format(self.collFitsFile, self.CollHeader['ATOM'],
                                                                                                    sym2name[self.elem]), calling=self.calling)
                
        #Read data
        self._CollData = CollExt.data
        if self.NLevels is None:
            try:
                self.NLevels = self.CollHeader['N_LEVELS']
            except:
                log_.error('N_LEVELS is not set in {0}'.format(self.collFitsFile))
            
        self.log_.message('NLevels of collisional data: {0}'.format(self.NLevels),
                          calling=self.calling)
        collFits.close()


    def initOmegas(self, OmegaInterp='Cheb'):
        """
        Initialization of the Chebishev orders and coefficients. This method is called by __init__ and when 
        OmegaInterp is changed.

        Parameter:
            - OmegaInterp   one of ('Cheb', 'Linear')
        
        """
        self.OmegaInterp = OmegaInterp
        self._ChebOrder = self.getChebOrder()
        self.ChebCoeffs = self._getChebCoeffArray()


    def _test_lev(self, level):
        """
        Test whether selected level is legal

        Parameters:
            - level        selected atom level

        """       
        if level < -1 or level == 0 or level > self.NLevels:
            self.log_.error('Wrong value for level: {0} (maximum: {1})'.format(level, self.NLevels),
                            calling=self.calling)
        
        
    def getSources(self):
        """
        Return bibliographic sources for atomic data, as listed in the headers of the fits files.

        """
        sources = []
        header = self.CollHeader
        for i in range(len(header.items())):
            if 'SOURCE' in header.items()[i][0]:
                number = header.items()[i][0].lstrip('SOURCE')
                try:
                    sources.append(self.atom + ': ' + header.get('NOTE' + str(number)) + ':', header.get('SOURCE' + str(number)))
                except:
                    sources.append(self.atom + ': ' + 'Collision strengths:', header.get('SOURCE' + str(number)))
        return sources
    
    def printSources(self):
        
        for source in self.getSources():
            print(source)    
    
    def getChebOrder(self, lev_i= -1, lev_j= -1):
        """
        Return order of Chebyshev polynomial fitting the collision strengths of selected transition.
        If transition not specified, return result for all transitions.
        
        Usage:
            O3.getChebOrder(4, 2)
        
        Parameters:
            - lev_i  upper level (default= -1, returns complete array)
            - lev_j  lower level (default= -1, returns complete array)

        """
        self._test_lev(lev_i)
        self._test_lev(lev_j)
        if (lev_i == -1):
            if (lev_j == -1):
                ChebOrder = np.zeros([self.NLevels, self.NLevels])
                for i in range(self.NLevels - 1):
                    lev_i = i + 1
                    j = i + 1
                    while (j < self.NLevels):
                        lev_j = j + 1
                        fieldName = 'Omega({0}->{1})'.format(lev_j, lev_i)
                        index = self.CollExtNames.index(fieldName)
                        if self.OmegaInterp == 'Linear':
                            ChebOrder[j][i] = -1
                        else:
                            ChebOrder[j][i] = int(self.CollHeader['TUNIT{0}'.format(index + 1)])
                        j += 1
        else:   
            if self.OmegaInterp == 'Linear':
                ChebOrder = -1
            else:
                fieldName = 'Omega({0}->{1})'.format(lev_i, lev_j)
                index = self.CollExtNames.index(fieldName)
                ChebOrder = int(self.CollHeader['TUNIT{0}'.format(index + 1)])
        return ChebOrder


    def _getChebCoeffArray(self):
        """
        Return coefficients of Chebyshev fit of all transitions in a given atom
            (3D array of Chebyshev fits; the first two dimensions are atomic levels,
            the third is the array of coefficients of the Chebyshev polynomial
            for the corresponding transition.

        """

        def _chebyshevFit(temArray, omegaDict):
            """
            Return coefficients of Chebyshev fit of a given transition
            
            Parameters:
                - temArray    array of temperature points in collisional strengths data
                - omegaDict   order of Chebyshev fit + array of collisional strengths 
                                of a given transition

            """
            warnings.filterwarnings('ignore', "The fit may be poorly conditioned")

            # -1 because (IRAF's order) = (polynomial degree + 1)
            order = omegaDict['order'] - 1
            data = omegaDict['data']
            coeffs = chebyshev.chebfit(temArray, data, order)
            return coeffs
        
        if np.max(self._ChebOrder) == -1:
            return -1
        else:
            ChebCoeffArray = np.zeros([self.NLevels, self.NLevels, np.max(self._ChebOrder)])
            for i in range(self.NLevels - 1):
                lev_i = i + 1
                j = i + 1
                while (j < self.NLevels):
                    lev_j = j + 1
                    if self._ChebOrder[lev_j - 1, lev_i - 1] == -1:
                        ChebCoeffArray[i, j, 0: self._ChebOrder[lev_j - 1, lev_i - 1]] = -1
                    else:
                        myOmega = {'order': self._ChebOrder[lev_j - 1, lev_i - 1],
                                   'data': self.getOmegaArray(lev_j, lev_i)}
                        if myOmega['order'] != 0.0:
                            ChebCoeffArray[i, j, 0: self._ChebOrder[lev_j - 1, lev_i - 1]] = \
                                np.array(_chebyshevFit(self.getTemArray(), myOmega))
                    j += 1
            return ChebCoeffArray



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
            result = []
            for name in self._CollData.names:
                if 'Omega' in name:
                    result.append(name + ": " + str(self._CollData.field(name)))
            return result
        else:
            return (self._CollData.field('Omega({0}->{1})'.format(lev_i, lev_j)))


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
            self.log_.warn("wrong levels given: {0} <= {1}".format(lev_i, lev_j), calling=self.calling)
            return None            
        tem_in_file_units = self._transfo_tem(tem)
        TemArray = self.getTemArray()
        if (lev_i == -1) and (lev_j == -1):
            tem = np.asarray(tem)
            res_shape = [self.NLevels, self.NLevels]
            for sh in tem.shape:
                res_shape.append(sh)
            Omega = np.zeros(res_shape)
    
            for i in range(self.NLevels - 1):
                j = i + 1
                while (j < self.NLevels):
                    Omega[j][i] = self.getOmega(tem, j + 1, i + 1)
                    j += 1
        elif self._ChebOrder[lev_i - 1, lev_j - 1] != -1:
            fit = self.ChebCoeffs[lev_j - 1, lev_i - 1, 0: self._ChebOrder[lev_i - 1, lev_j - 1]]
            tem_eval = tem_in_file_units
            if self.noExtrapol or config.get_noExtrapol():
                leftExtrapol = np.NAN
                rightExtrapol = np.NAN
            else:
                leftExtrapol = TemArray[0]
                rightExtrapol = TemArray[-1]
            
            if np.ndim(tem_eval) > 0:
                left = (tem_eval < TemArray[0])
                tem_eval[left] = leftExtrapol
                right = (tem_eval > TemArray[-1])
                tem_eval[right] = rightExtrapol
            else:
                if tem_eval < TemArray[0]:
                    tem_eval = leftExtrapol
                elif tem_eval > TemArray[-1]:
                    tem_eval = rightExtrapol
            Omega = chebyshev.chebval(tem_eval, fit)
        else:
            OmegaArray = self.getOmegaArray(lev_i, lev_j)
            if self.noExtrapol or config.get_noExtrapol():
                leftExtrapol = np.NAN
                rightExtrapol = np.NAN
            else:
                leftExtrapol = OmegaArray[0]
                rightExtrapol = OmegaArray[-1]
            Omega = np.interp(tem_in_file_units, self.getTemArray(), OmegaArray,
                              left=leftExtrapol, right=rightExtrapol)
        return np.squeeze(Omega)


    def getTemArray(self, keep_unit=True):
        """
        Return array of tabulated original temperature points (as in fits file) 
            of collision strengths.
        
        Parameters:
            - keep_unit   return temperature in file units (default) or change it to Kelvin (False)


        """
        tem_in_file = self._CollData.field(self.CollHeader['TTYPE1'])
        if keep_unit:
            return tem_in_file
        else:            
            if (self.tem_units == "log(K)"):
                return pow(10., tem_in_file)
            else:
                if (self.tem_units == "K/10000"):
                    return tem_in_file * 1.e4
                else: #T in Kelvin in the fits file
                    return tem_in_file

class _CollDataAscii(object):

    def __init__(self, elem=None, spec=None, atom=None, OmegaInterp='Linear', noExtrapol = False, 
                 NLevels=None):
        self.log_ = log_
        
        if atom is not None:
            self.atom = atom
            self.elem = parseAtom(atom)[0]
            self.spec = int(parseAtom(atom)[1])
        else:
            self.elem = elem
            self.spec = int(spec)
            self.atom = elem + str(self.spec)
        self.NLevels = NLevels
        self.name = sym2name[self.elem]
        self.noExtrapol = noExtrapol
        self.calling = 'Atom ' + self.atom
        if OmegaInterp != 'Linear':
            self.log_.error('Ascii files does not support other interpolation than Linear', 
                            calling = self.calling)
        
        self._loadAscii()
        #self.initOmegas(OmegaInterp=OmegaInterp)
        if 'T_UNIT' in self.comments:
            self.tem_units = self.comments['T_UNIT']
        else:
            self.tem_units = "K/10000" # Default value
        if (self.tem_units == "log(K)"):
            def transfo_tem(x):
                y = np.log10(x)
                if np.ndim(y) == 0:
                    y = np.array(y)
                return y
        elif (self.tem_units == "K/10000"):
                def transfo_tem(x):
                    y = x / 1e4
                    if np.ndim(y) == 0:
                        y = np.array(y)
                    return y
        else: #T in Kelvin in the fits file
                def transfo_tem(x):
                    y = x
                    return x
        self._transfo_tem = transfo_tem

    def _loadAscii(self):
        
        self.collFile = atomicData.getDataFile(self.atom, data_type='coll')
        if self.collFile is None:
            self.log_.error('No coll data for ion {0}'.format(self.atom), calling=self.calling)
            return None
        if self.atom in config.DataFiles:
            if self.collFile not in config.DataFiles[self.atom]:
                config.DataFiles[self.atom].append(self.collFile)
        else:
            config.DataFiles[self.atom] = [self.collFile]
            
        self.collPath = atomicData.getDirForFile(self.collFile)
        file_to_open = '{0}/{1}'.format(self.collPath, self.collFile)
        if not os.path.exists(file_to_open):
            self.log_.error('File {0} not found'.format(file_to_open), calling=self.calling)
            
        self.log_.message('Reading coll data from {0}'.format(self.collFile), calling=self.calling)

    
        # Read data from ascii file
        #
        # Read energies and stat weights 
        f = open(file_to_open)
        data = f.readlines()
        f.close()
    
        coll_data = np.array([d.split() for d in data if d[0:3]!='***'], dtype='float')
        comments_tab = [d for d in data if d[0:3]=='***']
        self.comments = {}
        for com in comments_tab:
            key = com.split()[1]
            self.comments[key] = com.split(key)[1].replace("'", "").replace('"', '').strip()
        self._lev_is = coll_data[:,0]
        self._lev_js = coll_data[:,1]
        self._TemArray = coll_data[0,2:]
        if self.NLevels is None:
            self.NLevels = int(np.max(self._lev_js))
        else:
            self.NLevels = np.min((self.NLevels, int(np.max(self._lev_js))))

        self._CollArray = np.zeros((self.NLevels, self.NLevels, len(self._TemArray)))
        for i in range(len(self._lev_is)):
            lev_i = int(self._lev_is[i])
            lev_j = int(self._lev_js[i])
            if (lev_i != 0) and (lev_j != 0) and (lev_i <= self.NLevels) and (lev_j <= self.NLevels):
                self._CollArray[lev_i-1, lev_j-1, :] = coll_data[i,2:]
           
        if 'SPECTRUM' in self.comments:
            if int(self.comments['SPECTRUM']) != self.spec:
                log_.error('The spectrum I read in the file {0} is {1}, but you are requesting {2}'.format(self.collFile, self.comments['SPECTRUM'],
                                                                                                    self.spec), calling=self.calling)
        if 'ATOM' in self.comments:
            if self.comments['ATOM'] != sym2name[self.elem]:
                log_.error('The element name I read in the file {0} is {1}, but I was expecting {2}. Check the keyword ATOM'.format(self.collFile, self.comments['ATOM'],
                                                                                                    sym2name[self.elem]), calling=self.calling)        
        
    def getSources(self):
        sources = []
        sources_dic = {}
        notes = {}
        for com in self.comments:
            if 'SOURCE' in com:
                number = com.split('SOURCE')[1]
                sources_dic[number] = self.comments[com]
            if 'NOTE' in com:
                number = com.split('NOTE')[1]
                notes[number] = self.comments[com]
        for key in sources_dic:
            try:
                sources.append('{0}: {1}: {2}'.format(self.atom, notes[key], sources_dic[key])) 
            except:
                sources.append('{0}: {1}'.format(self.atom, sources[key]))
        return sources
                
    def printSources(self):
        
        for source in self.getSources():
            print(source)

    def _test_lev(self, level):
        """
        Test whether selected level is legal

        Parameters:
            - level        selected atom level

        """       
        if level < -1 or level == 0 or level > self.NLevels:
            self.log_.error('Wrong value for level: {0} (maximum: {1})'.format(level, self.NLevels),
                            calling=self.calling)

    @profile
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
            return self._CollArray[lev_j-1, lev_i-1,:]
        
    @profile
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
        tem_in_file_units = self._transfo_tem(tem)
        if (lev_i == -1) and (lev_j == -1):
            tem = np.asarray(tem)
            res_shape = [self.NLevels, self.NLevels]
            for sh in tem.shape:
                res_shape.append(sh)
            Omega = np.zeros(res_shape)
    
            for i in range(self.NLevels - 1):
                j = i + 1
                while (j < self.NLevels):
                    Omega[j][i] = self.getOmega(tem, j + 1, i + 1)
                    j += 1
        else:
            OmegaArray = self.getOmegaArray(lev_i, lev_j)
            if self.noExtrapol or config.get_noExtrapol():
                leftExtrapol = np.NAN
                rightExtrapol = np.NAN
            else:
                leftExtrapol = OmegaArray[0]
                rightExtrapol = OmegaArray[-1]
            Omega = np.interp(tem_in_file_units, self.getTemArray(), OmegaArray,
                              left=leftExtrapol, right=rightExtrapol)
        
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
                else: #T in Kelvin in the fits file
                    return self._TemArray


class _CollDataStout(object):
    
    def __init__(self, elem=None, spec=None, atom=None, OmegaInterp='Linear', noExtrapol = False, 
                 NLevelsMax=None):
        self.log_ = log_
        self.calling = '_CollDataStout'
        if not config.INSTALLED['Stout']:
            log_.error('The STOUT_DIR environment variable is not defined', calling=self.calling)
            return None
        
        if atom is not None:
            self.atom = atom
            self.elem = parseAtom(atom)[0]
            self.spec = int(parseAtom(atom)[1])
        else:
            self.elem = elem
            self.spec = int(spec)
            self.atom = elem + str(self.spec)
        self.name = sym2name[self.elem]
        self.noExtrapol = noExtrapol
        
        if OmegaInterp != 'Linear':
            self.log_.error('Stout files does not support other interpolation than Linear', 
                            calling = self.calling)
        
        self._loadStout()
    
    def _loadStout(self):
        st_ion = '{}_{}'.format(self.atom.lower(), self.spec)
        self.fullFileName = '{0}/{1}/{1}.coll'.format(config.Stout_dir, st_ion)
        
        self.collFile = self.fullFileName.split('/')[-1]
        self.collPath = self.fullFileName[:-len(self.collFile)]
        
        

class Atom(object):
    """
    Define the atom object, fill it with data, explore the data, and 
    compute quantities such as level populations and line emissivities.

    """
    # This following Class variable will hold the (unique) references of every Atom instance created.
    # This allows to list all the references used in a project.
    
    
    @profile
    def __init__(self, elem=None, spec=None, atom=None, OmegaInterp='Linear', noExtrapol = False, NLevels=None):
        """
        Atom constructor
        
        Parameters:
            - elem          symbol of the selected element
            - spec          ionization stage in spectroscopic notation (I = 1, II = 2, etc.)
            - atom          ion (e.g. 'O3').
            - OmegaInterp   one of ('Cheb', 'Linear'). Default is Linear. "Cheb" works only for fits files. 
            - noExtrapol    if set to False (default), Omega will be extrapolated above and below
                            the highest and lowest temperatures where it is defined. If set to True
                            a NaN will be return.
            
        Usage:
            O3 = pn.Atom('O',3)
            N2 = pn.Atom(atom='N2')
            S2 = pn.Atom(atom='S2', OmegaInterp='Linear')
        """        
        self.log_ = log_
        self.type = 'coll'
        if atom is not None:
            self.atom = str.capitalize(atom)
            self.elem = parseAtom(self.atom)[0]
            self.spec = int(parseAtom(self.atom)[1])
        else:
            if elem[0].isalpha():
                self.elem = str.capitalize(elem)
            else:
                self.elem = elem
            self.spec = int(spec)
            self.atom = self.elem + str(self.spec)
        self.name = sym2name[self.elem]
        try:
            self.Z = Z[self.elem]
        except:
            self.Z = -1
        if self.elem in IP:
            if self.spec == 0:
                self.IP = -1
                self.IP_up = -1
            elif self.spec == 1:
                self.IP = 0
                try:
                    self.IP_up = IP[self.elem][self.spec-1]
                except:
                    self.IP = -1                    
            else:
                try:
                    self.IP = IP[self.elem][self.spec-2]
                except:
                    self.IP = -1
                try:
                    self.IP_up = IP[self.elem][self.spec-1]
                except:
                    self.IP = -1                    
        else:
            self.IP = -1
            self.IP_up = -1
        self.calling = 'Atom ' + self.atom
        self.log_.message('Making atom object for {0} {1}'.format(self.elem, self.spec), calling=self.calling)
        self.NLevels = NLevels
        dataFile = atomicData.getDataFile(self.atom, data_type='atom')
        if dataFile is None:
            self.atomFileType = None
        else:
            self.atomFileType = dataFile.split('.')[-1]
        if self.atomFileType == 'fits':
            self.AtomData = _AtomDataFits(elem=self.elem, spec=self.spec, atom=self.atom, NLevels=self.NLevels)
        elif self.atomFileType == 'dat':
            self.AtomData = _AtomDataAscii(elem=self.elem, spec=self.spec, atom=self.atom, NLevels=self.NLevels)
        elif self.atomFileType == 'chianti':
            self.AtomData = _AtomChianti(elem=self.elem, spec=self.spec, atom=self.atom, NLevels=self.NLevels)
        elif self.atomFileType == 'stout':
            self.AtomData = _AtomDataStout(elem=self.elem, spec=self.spec, atom=self.atom, NLevels=self.NLevels)
        elif self.atomFileType is None:
            self.AtomData = _AtomDataNone()
        else:
            log_.error('Atom file extensions must be fits, dat or chianti')
                    
        self.atomFile = self.AtomData.atomFile
        self.atomPath = self.AtomData.atomPath
        self.atomFitsFile = self.atomFile # Obsolete
        self.atomFitsPath = self.atomPath # Obsolete
        self.wave_Ang = self.AtomData.wave_Ang
        self.getStatWeight = self.AtomData.getStatWeight
        self.getEnergy = self.AtomData.getEnergy
        self.atomNLevels = self.AtomData.NLevels


        dataFile = atomicData.getDataFile(self.atom, data_type='coll')
        if dataFile is None:
            self.collFileType = None
        else:
            self.collFileType = dataFile.split('.')[-1]
        if self.collFileType == 'fits':
            self.CollData = _CollDataFits(elem=self.elem, spec=self.spec, atom=self.atom, 
                                         OmegaInterp=OmegaInterp, noExtrapol = noExtrapol, NLevels=self.NLevels)
        elif self.collFileType == 'dat':
            self.CollData = _CollDataAscii(elem=self.elem, spec=self.spec, atom=self.atom, 
                                          OmegaInterp='Linear', noExtrapol = noExtrapol, NLevels=self.NLevels)
        elif self.collFileType == 'chianti':
            self.CollData = _CollChianti(elem=self.elem, spec=self.spec, atom=self.atom, NLevels=self.NLevels)
        elif self.collFileType == 'stout':
            self.CollData = _CollDataStout(elem=self.elem, spec=self.spec, atom=self.atom, NLevels=self.NLevels)            
        elif self.collFileType is None:
            self.CollData = _CollDataNone()
        try:
            self.CollHeader = self.CollData.CollHeader
        except:
            pass
        if "comments" not in self.CollData.__dict__.keys():
            self.CollData.comments = []
        self.getOmegaArray = self.CollData.getOmegaArray
        self.getTemArray = self.CollData.getTemArray
        self.collFile = self.CollData.collFile
        self.collPath = self.CollData.collPath
        self.collFitsFile = self.collFile # Obsolete
        self.collFitsPath = self.collPath # Obsolete
        self.collNLevels = self.CollData.NLevels
        self.tem_units = self.CollData.tem_units

        self.NLevels = np.min((self.atomNLevels, self.collNLevels))
            
        try:
            self.gs = self.AtomData.gs
        except:
            self.gs = gsFromAtom(self.atom)
        try:
            self.AtomHeader = self.AtomData.AtomHeader
        except:
            self.AtomHeader = None
        try:
            self.NIST = self.AtomData.NIST
        except:
            try:
                self.NIST = getLevelsNIST(self.atom, self.NLevels)
            except:
                self.NIST = None
        
        self.lineList = []
        for i in np.arange(self.NLevels):
            for j in np.arange(i):
                self.lineList.append(self.wave_Ang[i][j])
        self.lineList = np.array(self.lineList)
        
        self.energy_Ryd = quiet_divide(CST.RYD_ANG, self.wave_Ang)
        self.energy_eV = CST.RYD_EV * self.energy_Ryd
        
        self._A = self.getA() # index = quantum number - 1
        self._Energy = self.getEnergy() # Angstrom^-1
        self._StatWeight = self.getStatWeight()
        if self.NLevels > 0:
            self.EnergyNLevels = len(self._Energy)
        else:
            self.EnergyNLevels = None
            
    def getOmega(self, tem, lev_i= -1, lev_j= -1, wave= -1):
        """
        Return interpolated value of the collision strength value at the given temperature 
            for the complete array or a specified transition.
        If kappa is not None (non-maxwellian distribution of e-velocities), the collision 
            strength is corrected as in Mendoza & Bautista, 2014 ApJ 785, 91.

        Usage:
            O3.getOmega(15000.)
            O3.getOmega([8e3, 1e4, 1.2e4])
            O3.getOmega([8e3, 1e4, 1.2e4], 5, 4)
        
        Parameters:
            - tem    electronic temperature in K. May be an array.
            - lev_i  upper level
            - lev_j  lower level

        """
        """
        """
        
        
        
        if wave != -1:
            lev_i, lev_j = self.getTransition(wave)
        kappa = config.kappa 
        if kappa is None:
            to_return = self.CollData.getOmega(tem, lev_i, lev_j)
        else:
            #ToDo The Kappa correction should come AFTER the transformation into CS unit
            if (lev_i == -1) and (lev_j == -1):
                tem = np.asarray(tem)
                res_shape = [self.collNLevels, self.collNLevels]
                for sh in tem.shape:
                    res_shape.append(sh)
                Omega = np.zeros(res_shape)
        
                for i in range(self.collNLevels - 1):
                    j = i + 1
                    while (j < self.collNLevels):
                        Omega[j][i] = self.getOmega(tem, j + 1, i + 1)
                        j += 1
            else:
                OmegaMB = self.CollData.getOmega(tem, lev_i, lev_j)
                delta_E = self.getEnergy(lev_i, unit='eV') - self.getEnergy(lev_j, unit='eV')
                correc = ((kappa - 3./2.)**(-0.5) / kappa * gamma(kappa+1) / gamma(kappa-0.5) * 
                          (1 + delta_E/((kappa-1.5)*CST.BOLTZMANN_eVK*tem))**(-kappa)) * np.exp(delta_E/CST.BOLTZMANN_eVK/tem)
    
                Omega = correc * OmegaMB
                log_.message('Correcting for Kappa={0} by {1}'.format(kappa, correc), self.calling)

            to_return = np.squeeze(Omega)
        if 'COEFF' in self.CollData.comments:
            to_return *= float(self.CollData.comments['COEFF'])
        if 'O_UNIT' in self.CollData.comments:
            if self.CollData.comments['O_UNIT'] == 'DEEX RATE COEFF':
                to_return /= CST.KCOLLRATE / tem ** 0.5 / self.getStatWeight(lev_i)
            elif self.CollData.comments['O_UNIT'] == 'RATE COEFF':
                deltaE = self.getEnergy(lev_i, unit='erg') - self.getEnergy(lev_j, unit='erg')
                to_return *= (self.getStatWeight(lev_j) / self.getStatWeight(lev_i) * np.exp(deltaE /(CST.BOLTZMANN * tem))) #q21
                to_return /= CST.KCOLLRATE / tem ** 0.5 / self.getStatWeight(lev_i)
            elif self.CollData.comments['O_UNIT'] == 'COOLING':
                deltaE = self.getEnergy(lev_i, unit='erg') - self.getEnergy(lev_j, unit='erg')
                to_return /= deltaE # Loss to q12                
                to_return *= (self.getStatWeight(lev_j) / self.getStatWeight(lev_i) * np.exp(deltaE /(CST.BOLTZMANN * tem))) #q21
                to_return /= (CST.KCOLLRATE / np.sqrt(tem) / self.getStatWeight(lev_i)) # Omega
                
        return to_return
    
    @profile
    def getCollRates(self, tem, NLevels=None):
        """
        Return (n_levels x n_levels) array of collision rates at given temperature. 
        
        Usage:
            O3.getCollRates(tem=10000)
            O3.getCollRates([8e3, 1e4, 1.2e4])
        
        Parameters:
            - tem     electronic temperature in K. May be an array.


        """
        tem = np.asarray(tem)
        if NLevels is None:
            NLevels = np.min((self.collNLevels, self.EnergyNLevels))
        res_shape = [NLevels, NLevels]
        for sh in tem.shape:
            res_shape.append(sh)
        resultArray = np.zeros(res_shape)
        Omegas = self.getOmega(tem)
        for i in range(NLevels - 1):
            lev_i = i + 1
            j = i + 1
            energy_i = self._Energy[i]
            stat_weight_i = self._StatWeight[i]
            while (j < NLevels):
                lev_j = j + 1 
                energy_j = self._Energy[j]
                stat_weight_j = self._StatWeight[j]
                resultArray[j][i] = CST.KCOLLRATE / tem ** 0.5 / stat_weight_j * Omegas[lev_j-1, lev_i-1]
                resultArray[i][j] = ((stat_weight_j) / (stat_weight_i) * 
                                      np.exp((energy_i - energy_j) / (CST.BOLTZMANN_ANGK * tem)) * 
                                      resultArray[j][i])
                j += 1
        
        return np.squeeze(resultArray)

    
    def _Transition(self, wave, maxErrorA = 5.e-3, maxErrorm = 5.e-2):
        """
        Return an array with computed upper level, computed lower level, computed wavelength, 
            input wavelength
        
        Parameters:
            - wave       wavelength either in Angstrom (a float or a label: e.g., 5007, '5007A') 
                            or in micron (a label: '51.5m')
            - maxErrorA: tolerance if the input wavelength is in Angstrom
            - maxErrorm: tolerance if the input wavelength is in micron
                            
        """
        if str(wave)[-1] == 'A':
            inputWave = float(wave[:-1])
            label = '{0}_{1}'.format(self.atom, wave)
            maxError = maxErrorA
        elif str(wave)[-1] == 'm':
            inputWave = float(wave[:-1]) * 1e4
            label = '{0}_{1}'.format(self.atom, wave)
            maxError = maxErrorm
        else:
            inputWave = wave
            label = '{0}_{1}A'.format(self.atom, int(wave))
            maxError = maxErrorA
            
        if label in label2levelDict:
            result = [label2levelDict[label][0], label2levelDict[label][1], inputWave, inputWave]
            return(result)
        
        j, i = np.unravel_index(np.argmin(abs(self.wave_Ang - inputWave)), self.wave_Ang.shape)
        bestWave = self.wave_Ang[i, j]
        error = np.abs(bestWave - inputWave) / inputWave
        result = [i + 1, j + 1, bestWave, inputWave]
        if error > maxError:
            self.log_.warn('_Transition: wavelengths differ by more than {0:.2f}%: input = {1:.2f}, output = {2:.2f}'\
                           .format(100 * maxError, inputWave, bestWave), calling=self.calling)
        return(result)


    def getTransition(self, wave, maxErrorA = 5.e-3, maxErrorm = 5.e-2):
        """
        Return the indexes (upper level, lower level) of a transition for a given atom 
            from the wavelength.
        
        Usage:
            O3.getTransition(4959)
            
        Parameters:
            - wave      wavelength in Angstrom (a float or a label: e.g., 5007, '5007A') 
                or in micron (a label: '51.5m')
            - maxErrorA: tolerance if the input wavelength is in Angstrom
            - maxErrorm: tolerance if the input wavelength is in micron
                
        """ 
        res = self._Transition(wave, maxErrorA = maxErrorA, maxErrorm = maxErrorm)
        return(res[0], res[1])
        

    def printTransition(self, wave):
        """
        Print info on transition associated to input wavelength.
        
        Usage:
            O3.printTransition(4959)
            
        Parameters:
            - wave      wavelength in Angstrom (a float or a label: e.g., 5007, '5007A') 
                or in micron (a label: '51.5m')
                
        """
        closestTransition = self._Transition(wave)
        relativeError = closestTransition[3] / closestTransition[2] - 1
        print('Input wave: {0:.1F}'.format(closestTransition[3]))
        print('Closest wave found: {0:.1F}'.format(closestTransition[2]))
        print('Relative error: {0:.0E} '.format(relativeError))
        print('Transition: {0[0]} -> {0[1]}'.format(closestTransition))
        return
    
    def printSources(self):
        
        for source in self.getSources():
            print(source)    
    
    def getSources(self):
        sources = []
        if self.AtomData is not self.CollData:
            sources.extend(self.AtomData.getSources())
            sources.extend(self.CollData.getSources())
        else:
            sources.extend(self.AtomData.getSources())
        return sources
    
    def _test_lev(self, level):
        """
        Test whether selected level is legal

        Parameters:
            - level        selected atom level

        """       
        if level < -1 or level == 0 or level > self.NLevels:
            self.log_.error('Wrong value for level: {0}, maximum = {1}'.format(level, self.NLevels),
                            calling=self.calling)

    def getA(self, lev_i= -1, lev_j= -1, wave= -1):
        """
        Return the transition probability data. 
        If no arguments are given, the whole array of A is returned.
        A specific A value can be obtained by giving either the upper and lower levels or 
            the wavelength of the transition (keyword wave).
            
        Usage:
            A_O3 = O3.getA()          # The whole A array is stored in A_O3
            O3.getA(4, 2)      # A(4, 2) of the O3 atom is printed
            O3.getA(2, 4)      # Returns 0
            O3.getA(wave=4959)      

        Parameters:
            - lev_i  upper level of transition (default= -1, returns complete array)
            - lev_j  lower level of transition (default= -1, returns complete array)
            - wave   wavelength of transition. Takes precedence on lev_i and lev_j. Ignored if not set.
            
        """
        if wave != -1:
            lev_i, lev_j = self.getTransition(wave)
        
        return self.AtomData.getA(lev_i= lev_i, lev_j= lev_j)
        
    @profile
    def getPopulations(self, tem, den, product=True, NLevels=None):
        """
        Return array of populations at given temperature and density.
        The method returns a 1-, 2- or 3-D array containing the population of each level 
            for all temperatures and densities specified in the input vectors tem and den 
            (which can be n-element or 1-element vectors).
        If either quantity (tem or den) is a 1-element vector -that is, a single value-, 
            the resulting population array is collapsed along that dimension; 
            as a result, the result population array can be a 1-D, 2-D or 3-D array 
            (the three cases corresponding to situations in which both tem and den are single values; 
            one of them is a single value and the other an n-element vector; or both are multielement 
            vectors, respectively). In the general case, the level index is the first 
            [WARNING! It is not in physical unit, i.e. ground level = 0; to be normalized], 
            followed by the temperature index (if it exists) and the density index. 

        Usage:
            O3.getPopulations(1e4, 1e2)
            tem=np.array([10000., 12000., 15000., 20000]) # An array of four temperatures
            den=np.array([600., 800., 1000])      # An array of three densities
            O3.getPopulations(tem, den)           # is a (6, 4, 3) array
            O3.getPopulations(tem, den)[0,2,1]    # Returns the population of level 1 for T = 15000 
                                                    and Ne = 800
            tem = 20000                           # tem is no longer an array
            O3.getPopulations(tem, den)[0,2,1]  # Crashes: one index too much
            O3.getPopulations(tem, den)[0,1]    # Returns the population of level 1 for T = 20000 
                                                    and Ne = 800 [see warning]
            tem=np.array([10000., 15000., 20000]) # An array of three temperatures
            O3.getPopulations(tem, den, product = False)# is a (6, 3) array, tem and den beeing 
                                                            taken 2 by 2.
        
        Parameters:
            - tem       electronic temperature in K
            - den       electronic density in cm^-3
            - product   operate on all possible combinations of temperature and density 
                      (product = True, default case) or on those resulting from combining 
                      the i-th value of tem with the i-th value of den (product = False).
                      If product = False, then tem and den must be the same size.

        """
        tem = np.asarray(tem)
        den = np.asarray(den)
        if NLevels is None:
            n_level = self.NLevels
        else:
            n_level = NLevels
        if product:
            n_tem = tem.size
            n_den = den.size
            tem_ones = np.ones(n_tem)
            den_ones = np.ones(n_den)
            # q is vector-indexed (q(0, 1) = rate between levels 1 and 2)
            q = self.getCollRates(tem, n_level)
            Atem = np.outer(self._A[:n_level, :n_level], tem_ones).reshape(n_level, n_level, n_tem)
            pop_result = np.zeros((n_level, n_tem, n_den))
            sum_q_up = np.zeros((n_level, n_tem))
            sum_q_down = np.zeros((n_level, n_tem))
            sum_A = np.squeeze(Atem.sum(axis=1))
            self._critDensity = sum_A / q.sum(axis=1)
            for i in range(1, n_level):
                for j in range(i + 1, n_level):
                    sum_q_up[i] = sum_q_up[i] + q[i, j]
                for j in range(0, i):
                    sum_q_down[i] = sum_q_down[i] + q[i, j]
            coeff_matrix = ((np.outer(np.swapaxes(q, 0, 1), den) + 
                             np.outer(np.swapaxes(Atem, 0, 1), den_ones)).reshape(n_level, n_level, n_tem, n_den))
            coeff_matrix[0, :] = 1.
            for i in range(1, n_level):
                coeff_matrix[i, i] = (-(np.outer((sum_q_up[i] + sum_q_down[i]), den) + 
                                        np.outer(sum_A[i], den_ones)).reshape(1, 1, n_tem, n_den))
            vect = np.zeros(n_level)
            vect[0] = 1.
    
            for i_tem in range(n_tem):
                for i_den in range(n_den):
                    try:
                        pop_result[:, i_tem, i_den] = solve(np.squeeze(coeff_matrix[:, :, i_tem, i_den]), vect)
                    #except np.linalg.LinAlgError:
                    #    pop_result[:, i_tem, i_den] = np.nan
                    except:
                        self.log_.error('Error solving population matrix', calling=self.calling)
            pop = np.squeeze(pop_result)
        else:
            if tem.shape != den.shape:
                self.log_.error('tem and den must have the same shape', calling=self.calling)
                return None
            res_shape1 = [n_level]
            res_shape_rav1 = [n_level, tem.size]
            res_shape_rav2 = [n_level, n_level, tem.size]
            for sh in tem.shape:
                res_shape1.append(sh)
            tem_rav = tem.ravel()
            den_rav = den.ravel()
            q = self.getCollRates(tem_rav, n_level)
            A = self._A[:n_level, :n_level]
            pop_result = np.zeros(res_shape_rav1)
            coeff_matrix = np.ones(res_shape_rav2)
            sum_q_up = np.zeros(res_shape_rav1)
            sum_q_down = np.zeros(res_shape_rav1)
            sum_A = A.sum(axis=1)
            n_tem = tem_rav.size
            # Following line changed 29/11/2012. It made the code crash when atom_nlevels diff coll_nlevels
            #Atem = np.outer(self._A, np.ones(n_tem)).reshape(n_level, n_level, n_tem)
            Atem = np.outer(self._A[:n_level, :n_level], np.ones(n_tem)).reshape(n_level, n_level, n_tem)
            self._critDensity = Atem.sum(axis=1) / q.sum(axis=1)

            for i in range(1, n_level):
                for j in range(i + 1, n_level):
                    sum_q_up[i] = sum_q_up[i] + q[i, j]
                for j in range(0, i):
                    sum_q_down[i] = sum_q_down[i] + q[i, j]
            for row in range(1, n_level):
                # upper right half            
                for col in range(row + 1, n_level):
                    coeff_matrix[row, col] = den_rav * q[col, row] + A[col, row]
                # lower left half
                for col in range(0, row):
                    coeff_matrix[row, col] = den_rav * q[col, row]
                # diagonal
                coeff_matrix[row, row] = -(den_rav * (sum_q_up[row] + sum_q_down[row]) + sum_A[row])

            vect = np.zeros(n_level)
            vect[0] = 1.
            
            for i in range(tem.size):
                try:
                    pop_result[:, i] = solve(np.squeeze(coeff_matrix[:, :, i]), vect)
                except np.linalg.LinAlgError:
                    pop_result[:, i] = np.nan
                except:
                    self.log_.error('Error solving population matrix', calling=self.calling)
            
            pop = np.squeeze(pop_result.reshape(res_shape1))
            
        return pop

    
    def getLowDensRatio(self, lev_i1=-1, lev_i2=-1, wave1=-1, wave2=-1, to_eval=None):
        
        """
        Return the value of a diagostic ratio at the low density limit
        S2.getLowDensRatio(lev_i1 = 3, lev_i2 = 2)
        S2.getLowDensRatio(wave1 = 6716, wave2 = 6731)
        S2.getLowDensRatio(to_eval = 'L(6716)/L(6731)')
        """
        
        if wave1 != -1:
            lev_i1, lev_j1 = self.getTransition(wave1)
        if wave2 != -1:
            lev_i2, lev_j2 = self.getTransition(wave2)
            
        if to_eval is not None:
            L = lambda wave: self.getStatWeight(self.getTransition(wave)[0])
            return eval(to_eval)
            
        return self.getStatWeight(lev_i1) / self.getStatWeight(lev_i2)
        
        
    def getHighDensRatio(self, lev_i1=-1, lev_i2=-1, lev_j1=-1, lev_j2=-1, wave1=-1, wave2=-1, to_eval=None):
        
        """
        Return the value of a diagostic ratio at the high density limit
        S2.getHighDensRatio(lev_i1 = 3, lev_i2 = 2)
        S2.getHighDensRatio(wave1 = 6716, wave2 = 6731)
        S2.getHighDensRatio(to_eval = 'L(6716)/L(6731)')
        """
        
        if wave1 != -1:
            lev_i1, lev_j1 = self.getTransition(wave1)
        if wave2 != -1:
            lev_i2, lev_j2 = self.getTransition(wave2)
            
        if to_eval is not None:
            L = lambda wave: (self.getStatWeight(self.getTransition(wave)[0]) * 
                              self.getA(self.getTransition(wave)[0], self.getTransition(wave)[1]))
            return eval(to_eval)
            
        return (self.getStatWeight(lev_i1) / self.getStatWeight(lev_i2) *
                self.getA(lev_i1, lev_j1) / self.getA(lev_i2, lev_j2))
           
    def getDensityRange(self, lev_i1=-1, lev_i2=-1, lev_j1=-1, lev_j2=-1, wave1=-1, wave2=-1, 
                        to_eval=None, tol=0.1, tem=1e4):
        """
        Return the range of density where a given line ratio is between 10% and 90% of the low and high density limits
        """
        LowLim = self.getLowDensRatio(lev_i1, lev_i2, wave1, wave2, to_eval)
        HighLim = self.getHighDensRatio(lev_i1, lev_i2, lev_j1, lev_j2, wave1, wave2, to_eval)
        
        delta = abs(LowLim - HighLim)
        minRatio = min((LowLim, HighLim)) + tol * delta
        maxRatio = max((LowLim, HighLim)) - tol * delta
        dens1 = self.getTemDen(minRatio, tem=tem, lev_i1= lev_i1, lev_j1= lev_j1, lev_i2= lev_i2, lev_j2= lev_j2,
                  wave1= wave1, wave2= wave2, to_eval=to_eval)
        dens2 = self.getTemDen(maxRatio, tem=tem, lev_i1= lev_i1, lev_j1= lev_j1, lev_i2= lev_i2, lev_j2= lev_j2,
                  wave1= wave1, wave2= wave2, to_eval=to_eval)
        return(np.sort((dens1, dens2)))
               
    @profile
    def getCritDensity(self, tem, level= -1):
        """
        Return the critical density of selected level at given temperature. 
        If no transition is selected, return complete array.
        
        Usage:
            O3.getCritDensity(12000)
            O3.getCritDensity(12000, 4)
        
        Parameters:
            - tem    electronic temperature in K. May be an array.
            - level  selected atomic level (default= -1)

        """
        self._test_lev(level)
        self.getPopulations(tem, den=100.) # Any density would do
        if level != -1:
            return self._critDensity[level - 1]
        else:
            return self._critDensity
        
        
    @profile
    def getEmissivity(self, tem, den, lev_i= -1, lev_j= -1, wave= -1, product=True):
        """
        Return the line emissivity of selected transition or complete array of emissivities
        The transition is selected by the argument wave (if given); 
        if wave is not supplied, it is selected by the upper and lower levels (lev_i and lev_j); 
        if neither is given, the whole array is computed
        
        Usage:      
            O3.getEmissivity(12000, 100, 4, 2)         # (4, 2) transition
            O3.getEmissivity(10000, 10000, wave=5007)  # (4, 2) transition
            O3.getEmissivity(12000, 100)               # all transitions
            O3.getEmissivity([10000, 12000], [100, 500], 4, 2)
            O3.getEmissivity([10000, 12000], [100, 500])
            
        Parameters:
            - tem      electronic temperature in K. May be an array.
            - den      electronic density in cm^-3. May be an array.
            - lev_i    upper level (default= -1)
            - lev_j    lower level (default= -1)
            - wave     wavelength of transition. Takes precedence on lev_i and lev_j if set, 
                        ignored otherwise. It can also be a blend label.
            - product  Boolean. If True (default), all the combination of (tem, den) are used. 
                         If False, tem and den must have the same size and are joined.

        """
        if '{0}_{1}'.format(self.atom, wave) in BLEND_LIST:
            L = lambda wave: self.getEmissivity(tem, den, wave=wave, product=product)
            I = lambda lev_i, lev_j: self.getEmissivity(tem, den, lev_i=lev_i, lev_j=lev_j, product=product)
            try:
                res = eval(BLEND_LIST['{0}_{1}'.format(self.atom, wave)])
            except:
                log_.warn('{0} is not understood'.format(wave), calling=self.calling + 'getEmissivity')
                res = None
            return res
        self._test_lev(lev_i)
        self._test_lev(lev_j)
        tem = np.asarray(tem)
        den = np.asarray(den)
        if wave != -1:
            lev_i, lev_j = self.getTransition(wave)
        NLevels = self.NLevels
        if lev_i > NLevels or lev_j > NLevels:
            log_.error('The number of levels {} does not allow getting this emissivity. Consider changing the atomic data'.format(NLevels),
                          calling=self.calling) 
        if product:
            n_tem = tem.size
            n_den = den.size
            tem_ones = np.ones(n_tem)
            populations = self.getPopulations(tem, den, product=True)
            if ((lev_i == -1) and (lev_j == -1)):
                resultArray = np.zeros((NLevels, NLevels, n_tem, n_den))
                for i in range(NLevels):
                    lev_i = i + 1
                    j = i - 1 
                    while (j >= 0):
                        lev_j = j + 1
                        deltaE = (self._Energy[i] - self._Energy[j]) * CST.HPLANCK * CST.CLIGHT * 1.e8 
                        resultArray[i][j] = (deltaE * self._A[i, j] * populations[i].reshape(1, 1, n_tem, n_den) / 
                                             np.outer(tem_ones, den).reshape(1, 1, n_tem, n_den))
                        j -= 1
                return np.squeeze(resultArray)
            else:
                if (lev_i <= lev_j):
                    return 0.
                else:
                    i = lev_i - 1
                    j = lev_j - 1
                    deltaE = (self._Energy[i] - self._Energy[j]) * CST.HPLANCK * CST.CLIGHT * 1.e8 
                    return np.squeeze((populations[i] * deltaE * self._A[i, j]).reshape(1, 1, n_tem, n_den) / 
                                      np.outer(tem_ones, den).reshape(1, 1, n_tem, n_den))
        else:
            if tem.shape != den.shape:
                self.log_.error('tem and den must have the same shape', calling=self.calling)
                return None
            populations = self.getPopulations(tem, den, product=False)
            if (lev_i <= lev_j):
                return None
            else:
                i = lev_i - 1
                j = lev_j - 1
                deltaE = (self._Energy[i] - self._Energy[j]) * CST.HPLANCK * CST.CLIGHT * 1.e8 
                return populations[i] * deltaE * self._A[i, j] / den

        
    @profile
    def _getTemDen_1(self, int_ratio, tem= -1, den= -1, lev_i1= -1, lev_j1= -1, lev_i2= -1, lev_j2= -1,
                  wave1= -1, wave2= -1, maxError=1.e-3, method='nsect_recur', log=True, start_x= -1, end_x= -1,
                  to_eval=None, nCut=30, maxIter=20):

        ##
        # @todo manage blends
        self._test_lev(lev_i1)
        self._test_lev(lev_j1)
        self._test_lev(lev_i2)
        self._test_lev(lev_j2)
        tem = np.asarray(tem)
        den = np.asarray(den)
        if np.asarray(int_ratio).size != 1:
            shape = np.asarray(int_ratio).shape
            size = np.asarray(int_ratio).size
            result = np.zeros(size) 
            if np.asarray(tem).size != 1:
                if np.asarray(tem).shape != shape:
                    self.log_.error('getTemDen: int_ratio and tem/den must have the same size', calling=self.calling)
                    return None
                tem_ravel = np.asarray(tem).ravel()
                den_ravel = np.zeros(size) + den
            elif np.asarray(den).size != 1:
                if np.asarray(den).shape != shape:
                    self.log_.error('getTemDen: int_ratio and tem/den must have the same size', calling=self.calling)
                    return None
                den_ravel = np.asarray(den).ravel()
                tem_ravel = np.zeros(size) + tem
            else:
                den_ravel = np.zeros(size) + den
                tem_ravel = np.zeros(size) + tem
            int_ratio_ravel = np.asarray(int_ratio).ravel()
            for i in np.arange(size):
                result[i] = self._getTemDen_1(int_ratio_ravel[i], tem=tem_ravel[i], den=den_ravel[i],
                                           lev_i1=lev_i1, lev_j1=lev_j1, lev_i2=lev_i2, lev_j2=lev_j2,
                                           wave1=wave1, wave2=wave2, maxError=maxError, method=method,
                                           log=log, start_x=start_x, end_x=end_x,
                                           to_eval=to_eval)
            return result.reshape(shape)
            
        if wave1 != -1:
            lev_i1, lev_j1 = self.getTransition(wave1)
        if wave2 != -1:
            lev_i2, lev_j2 = self.getTransition(wave2)

        if to_eval is None:
            to_eval = 'I(' + str(lev_i1) + ',' + str(lev_j1) + ')/I(' + str(lev_i2) + ',' + str(lev_j2) + ')'


        if tem == -1:
            if start_x == -1:
                start_x = min(self.getTemArray(keep_unit=False))
                if log:
                    start_x = np.log10(start_x)
            if end_x == -1:
                end_x = max(self.getTemArray(keep_unit=False))
                if log:
                    end_x = np.log10(end_x)
            
            @profile
            def _func(x):
                """
                The function for which a root is looked for.
                It must return an array if x is an array, if an nsect-like method is used.
                The returned value is already normalized to the intensity ratio.

                """
                if log:
                    populations = self.getPopulations(10.**x, den)
                else:
                    populations = self.getPopulations(x, den)
                I = lambda lev_i, lev_j: (populations[lev_i - 1] * 
                                          (self._A[lev_i - 1, lev_j - 1] * 
                                           (self._Energy[lev_i - 1] - self._Energy[lev_j - 1])))
                L = lambda wave: (populations[self.getTransition(wave)[0] - 1] * 
                                  (self._A[self.getTransition(wave)[0] - 1, self.getTransition(wave)[1] - 1] * 
                                   (self._Energy[self.getTransition(wave)[0] - 1] - 
                                    self._Energy[self.getTransition(wave)[1] - 1])))
                result = eval(to_eval)
                return quiet_divide((result - int_ratio), int_ratio)
            
        elif den == -1:
            if start_x == -1:
                start_x = 1.e0
                if log: start_x = np.log10(start_x)
            if end_x == -1:
                end_x = 1e8
                if log: end_x = np.log10(end_x)

            @profile
            def _func(x):
                """
                The function for which a root is looked for.
                It must return an array if x is an array, if an nsect-like method is used.
                The returned value is already normalized to the intensity ratio.
                
                """
                if log:
                    populations = self.getPopulations(tem, pow(10., x))
                else:
                    populations = self.getPopulations(tem, x)
                I = lambda lev_i, lev_j: (populations[lev_i - 1] * 
                                          (self._A[lev_i - 1, lev_j - 1] * 
                                           (self._Energy[lev_i - 1] - self._Energy[lev_j - 1])))
                L = lambda wave: (populations[self.getTransition(wave)[0] - 1] * 
                                  (self._A[self.getTransition(wave)[0] - 1, self.getTransition(wave)[1] - 1] * 
                                   (self._Energy[self.getTransition(wave)[0] - 1] - 
                                    self._Energy[self.getTransition(wave)[1] - 1])))
                result = eval(to_eval)
                return quiet_divide((result - int_ratio), int_ratio)
        # improve exception handling (we must include cases where both tem, den = -1) 
        else:
            self.log_.error('ERROR in getTemDen: temperature and density cannot be simultaneously given',
                            calling=self.calling)

        if method == 'nsect_recur':
            """
            Recursive n-section method for finding a root
            It works by looking for the minumum of abs(f) using the fact that f(array) 
                returns an array

            """
            @profile
            def nsect_recur(f, x1, x2, nCut, maxIter, _iter=0):
                if _iter > maxIter:
                    return np.nan
                x = np.linspace(x1, x2, nCut)
                y = abs(f(x))
                x_min = np.argmin(y)
                if y[x_min] < maxError:
                    return x[x_min]
                else:
                    if x_min == 0:
                        x_min = 1
                    if x_min == nCut - 1:
                        x_min = nCut - 2
                    x1 = x[x_min - 1]
                    x2 = x[x_min + 1]
                    _iter = _iter + 1
                    return nsect_recur(f, x1, x2, nCut=nCut, maxIter=maxIter, _iter=_iter)

            result = nsect_recur(_func, start_x, end_x, nCut, maxIter)

        elif method == 'nsect_iter':
            """
            Iterative n-section method for finding a root (analogous to the bisection method)
            It works by looking for the minimum of abs(f) using the fact that f(array) returns an array

            """
            
            @profile
            def nsect_iter(f, x1, x2, nCut, maxIter):
                for i in range(maxIter):
                    x = np.linspace(x1, x2, nCut)
                    y = abs(f(x))
                    x_min = np.argmin(y)
                    if y[x_min] < maxError:
                        return x[x_min]
                    else:
                        if x_min == 0: x_min = 1
                        if x_min == nCut - 1: x_min = nCut - 2
                        x1 = x[x_min - 1]
                        x2 = x[x_min + 1]
                return np.nan

            result = nsect_iter(_func, start_x, end_x, nCut, maxIter)
            
        else:
            self.log_.error('ERROR in getTemDen: no valid method given', calling=self.calling)
            result = None

        if (log is True) and (result is not None):
            return pow(10., result)
        else:
            return result


    def _getTemDen_MP(self, int_ratio, tem= -1, den= -1, lev_i1= -1, lev_j1= -1, lev_i2= -1, lev_j2= -1,
                  wave1= -1, wave2= -1, maxError=1.e-3, method='nsect_recur', log=True, start_x= -1, end_x= -1,
                  to_eval=None, nCut=30, maxIter=20):
        
        if not config.INSTALLED['mp']:
            log_.error('_getTemDen_MP cannot be used in absence of multiprocessing package',
                          calling=self.calling)
            return None
        self._test_lev(lev_i1)
        self._test_lev(lev_j1)
        self._test_lev(lev_i2)
        self._test_lev(lev_j2)
        tem = np.asarray(tem)
        den = np.asarray(den)
        if np.asarray(int_ratio).size != 1:
            shape = np.asarray(int_ratio).shape
            size = np.asarray(int_ratio).size
            result = np.zeros(size) 
            if np.asarray(tem).size != 1:
                if np.asarray(tem).shape != shape:
                    self.log_.error('getTemDen: int_ratio and tem/den must have the same size', calling=self.calling)
                    return None
                tem_ravel = np.asarray(tem).ravel()
                den_ravel = np.zeros(size) + den
            elif np.asarray(den).size != 1:
                if np.asarray(den).shape != shape:
                    self.log_.error('getTemDen: int_ratio and tem/den must have the same size', calling=self.calling)
                    return None
                den_ravel = np.asarray(den).ravel()
                tem_ravel = np.zeros(size) + tem
            else:
                den_ravel = np.zeros(size) + den
                tem_ravel = np.zeros(size) + tem
            int_ratio_ravel = np.asarray(int_ratio).ravel()
        else:
            return self._getTemDen_1(int_ratio=int_ratio, tem=tem, den=den, lev_i1=lev_i1, lev_j1=lev_j1, lev_i2=lev_i2, lev_j2=lev_j2,
                  wave1=wave1, wave2=wave2, maxError=maxError, method=method, log=log, start_x=start_x,
                  end_x=end_x, to_eval=to_eval, nCut=nCut, maxIter=maxIter)
        
        Nprocs = config.Nprocs
        log_.message('number of CPUs = {0}'.format(Nprocs), calling=self.calling + '.getTemDenMP')
        gTDWorkerQ = Queue()
        gTDDoneQ = Queue()
        log_.message('Queues initialized', calling=self.calling + '.getTemDenMP')
        jobid = 0
        for int_rat1, tem1, den1 in zip(int_ratio_ravel, tem_ravel, den_ravel):
            gTDWorkerQ.put((jobid, int_rat1, tem1, den1))
            jobid += 1
        log_.message('put done', calling=self.calling + '.getTemDenMP')
        #gTDWorkerQSize = gTDWorkerQ.qsize() # this crash on OSX
        gTDWorkerQSize = size
        log_.message('Queue size {0}'.format(gTDWorkerQSize), calling=self.calling + '.getTemDenMP')

        
        gTDProcesses = []
        for i in range(Nprocs):
            p = Process(target=getTemDen_helper, args=(gTDWorkerQ, gTDDoneQ, self.atom, lev_i1, lev_j1, lev_i2, lev_j2,
                  wave1, wave2, maxError, method, log, start_x, end_x, to_eval, nCut, maxIter))
            p.start()
            gTDProcesses.append(p)
        log_.message('processes started', calling=self.calling + '.getTemDenMP')

        #
        result = []
        for i in range(gTDWorkerQSize):
            result.append(gTDDoneQ.get())
        log_.message('Result obtained', calling=self.calling + '.getTemDenMP')

        for p in gTDProcesses:
            p.join(timeout=0.1)
        log_.message('Joined', calling=self.calling + '.getTemDenMP')

        for i in range(Nprocs):
            gTDWorkerQ.put('STOP')
        result.sort()
        to_return = np.array(result)[:, 1]
        return to_return.reshape(int_ratio.shape)

    
    @profile
    def getTemDen(self, int_ratio, tem= -1, den= -1, lev_i1= -1, lev_j1= -1, lev_i2= -1, lev_j2= -1,
                  wave1= -1, wave2= -1, maxError=1.e-3, method='nsect_recur', log=True, start_x= -1, end_x= -1,
                  to_eval=None, nCut=30, maxIter=20):
        """
        Return either the temperature or the density given the other variable for a selected line ratio 
            of known intensity.
        The line ratio can involve two or more than two lines. 
        In the first case (only two lines), it can be specified giving either two transitions 
            (four atomic levels, i.e. two for each transition), or two wavelengths.
        In the general case (any number of lines), it can be specified as an algebraic expression 
            to be evaluated, involving either atomic levels or wavelengths.
        An array of values, rather than a single value, can also be given, in which case the result 
            will also be an array.
            
        Usage: 
            O3.getTemDen(150., den=100., wave1=5007, wave2=4363, maxError=1.e-2)
            O3.getTemDen(150., den=100., to_eval = '(I(4,3) + I(4,2) + I(4,1)) / I(5,4)')
            N2.getTemDen(150., den=100., to_eval = '(L(6584) + L(6548)) / L(5755)')
            O3.getTemDen([0.02, 0.04], den=[1.e4, 1.1e4], to_eval="I(5, 4) / (I(4, 3) + I(4, 2))")

        Parameters:
            - int_ratio    intensity ratio of the selected transition
            - tem          electronic temperature
            - den          electronic density
            - lev_i1       upper level of 1st transition
            - lev_j1       lower level of 1st transition
            - lev_i2       upper level of 2nd transition
            - lev_j2       lower level of 2nd transition
            - wave1        wavelength of 1st transition
            - wave2        wavelength of 2nd transition
            - maxError     tolerance on difference between input and computed ratio 
            - method       numerical method for finding the root (nsect_recur, nsect_iter)
            - log          switch of log (default = True). start_x and end_x are using this parameter.
            - start_x      lower end of the interval to explore. (default: lower end of collision 
                            strength temperature array for temperature, 1 if density)
            - end_x        higher end of the interval to explore. (default: higher end of collision 
                            strength temperature array for temperature, 1e8 if density)
            - to_eval      expression to be evaluated, using either I (for transitions identified 
                            through atomic levels) or L (for transitions identified through wavelengths)
            - nCut        number of sections in which each step is cut. 2 would be dichotomy.
            - maxIter     maximum number of iterations

        """        
        if config._use_mp:
            return self._getTemDen_MP(int_ratio=int_ratio, tem=tem, den=den, lev_i1=lev_i1, lev_j1=lev_j1, lev_i2=lev_i2, lev_j2=lev_j2,
                  wave1=wave1, wave2=wave2, maxError=maxError, method=method, log=log, start_x=start_x,
                  end_x=end_x, to_eval=to_eval, nCut=nCut, maxIter=maxIter)
        else:
            return self._getTemDen_1(int_ratio=int_ratio, tem=tem, den=den, lev_i1=lev_i1, lev_j1=lev_j1, lev_i2=lev_i2, lev_j2=lev_j2,
                  wave1=wave1, wave2=wave2, maxError=maxError, method=method, log=log, start_x=start_x,
                  end_x=end_x, to_eval=to_eval, nCut=nCut, maxIter=maxIter)

    def getIonAbundance(self, int_ratio, tem, den, lev_i= -1, lev_j= -1, wave= -1, to_eval=None, Hbeta=100.):
        """
        Compute the ionic abundance relative to H+ given the intensity of a line or sum of lines, 
        the temperature, and the density. 
        The line can be specified as a transition (i.e., giving the two atomic level involved), 
        as a wavelength, or as an algebraic expression. In the last case, a sum of lines 
        can also be supplied.

        Usage:
            O3.getIonAbundance(100, 1.5e4, 100., wave=5007)
            O3.getIonAbundance(130, 1.5e4, 100., to_eval='I(4,3) + I(4,2)')
            O3.getIonAbundance(np.array([100, 150]), np.array([1.5e4, 1.2e4]), np.array([100., 120]), 
                wave=5007)
            O2.getIonAbundance(int_ratio=100, tem=1.5e4, den=100., to_eval='L(3726)+L(3729)')

            
        Parameters:
            - int_ratio    relative line intensity (default normalization: Hbeta = 100). 
                            May be an array.
            - tem          electronic temperature in K. May be an array.
            - den          electronic density in cm^-3. May be an array.
            - lev_i        upper level of transition
            - lev_j        lower level of transition
            - wave         wavelength of transition. Takes precedence on lev_i and lev_j if set, 
                            ignored otherwise 
            - to_eval      expression to be evaluated. Takes precedence on wave if set, 
                            ignored otherwise.
            - Hbeta        line intensity normalization at Hbeta (default Hbeta = 100)
        
        """
        self._test_lev(lev_i)
        self._test_lev(lev_j)
        if np.ndim(tem) != np.ndim(den):
            self.log_.error('ten and den must have the same shape', calling=self.calling)
            return None
        if ((np.squeeze(np.asarray(int_ratio)).shape != np.squeeze(np.asarray(tem)).shape) | 
            (np.squeeze(np.asarray(den)).shape != np.squeeze(np.asarray(tem)).shape)):
            self.log_.warn('int_ratio, tem and den must does not have the same shape', calling=self.calling)
        if (lev_i == -1) & (lev_j == -1) & (wave == -1) & (to_eval is None):
            self.log_.error('At least one of lev_i, lev_j, wave or to_eval must be supplied', calling=self.calling)
            return None
        if to_eval == None:
            if wave != -1:
                lev_i, lev_j = self.getTransition(wave)     
            to_eval = 'I(' + str(lev_i) + ',' + str(lev_j) + ')' 
        I = lambda lev_i, lev_j: self.getEmissivity(tem, den, lev_i, lev_j, product=False)
        L = lambda wave: self.getEmissivity(tem, den, wave=wave, product=False)
        try:
            emis = eval(to_eval)
        except:
            self.log_.error('Unable to eval {0}'.format(to_eval), calling=self.calling)
            return None
        #int_ratio is in units of Hb = Hbeta keyword
        ionAbundance = ((int_ratio / Hbeta) * (getRecEmissivity(tem, den, 4, 2, atom='H1', product=False) / emis))
        return ionAbundance
    

    def printIonic(self, tem=None, den=None, printA=False, printPop=True, printCrit=True):
        """ 
        Print miscellaneous information (wavelengths, level populations, emissivities, 
            critical densities) for given physical conditions.
        If an electron temperature is given, level critical densities can be printed 
            (using printCrit=True)
        If an electron density is also given, line emissivities (also for Hbeta) are printed and level 
            populations can be printed (using printPop=True)
        
        Usage:
            O3.printIonic()
            O3.printIonic(printA=True)
            O3.printIonic(tem=10000., printCrit=True)
            O3.printIonic(tem=10000., den=1e3, printA=True, printPop=True, printCrit=True)

        Parameters:
            - tem          temperature
            - den          density
            - printA       also print transition probabilities (default=False)
            - printPop     also print level populations (needs tem and den)
            - printCrit    also print critical densities (needs tem)

        """
        print('elem = %s' % self.elem)
        print('spec = %i' % self.spec)
        if tem is not None:
            print('temperature = %6.1f K' % tem)
        if den is not None:
            print('density = %6.1f cm-3' % den)
        print("")
        if printPop and ((tem is None) or (den is None)):
            log_.warn('Cannot print populations as tem or den is missing', calling=self.calling)
            printPop = False
        if printCrit and (tem is None):
            log_.warn('Cannot print critical densities as tem is missing', calling=self.calling)
            printCrit = False
        to_print = ''
        if printPop:
            to_print += 'Level   Populations  '
        if printCrit:
            to_print += 'Critical densities'
        if to_print != '':
            print(to_print)
        if printCrit:
            critdens = self.getCritDensity(tem)
        if printPop:
            pop = self.getPopulations(tem, den)
        for i in range(0, self.NLevels):
            lev_i = i + 1
            to_print = 'Level %1i:  ' % (lev_i)
            if printPop:
                to_print += "%.3E  " % (pop[i])
            if printCrit:
                to_print += "%.3E" % critdens[i]
            if printPop or printCrit:
                print(to_print)
        if printPop or printCrit:
            print('')
            
        if (tem is not None) and (den is not None):
            emis = self.getEmissivity(tem, den)
        for i in range(1, self.NLevels):
            if printA:
                for j in range(i):
                    to_print = "{0:.3E}   ".format(np.float(self.getA(i + 1, j + 1)))
                    print(to_print, end="")
                print("")
            for j in range(i):
                if self.wave_Ang[i, j] > 10000.:
                    to_print = "%10.2fm " % (self.wave_Ang[i, j] / 1e4)
                else:
                    to_print = "%10.2fA " % self.wave_Ang[i, j]
                print(to_print, end="")
            print("")
            for j in range(i):
                print("    (%1i-->%1i) " % (i + 1, j + 1), end="")
            print("")
            if (tem is not None) and (den is not None):
                for j in range(i):
                    print("  %.3E " % emis[i,j], end="")
            print("\n")
        if (tem is not None) and (den is not None):
            try:
                H1 = RecAtom('H', 1)
                print("# H-beta volume emissivity:")
                print("%.3E N(H+) * N(e-)  (erg/s)" % H1.getEmissivity(tem, den, 4, 2))
            except:
                pass

    def printTemDen(self, int_ratio, tem= -1, den= -1, lev_i1= -1, lev_j1= -1, lev_i2= -1, lev_j2= -1, wave1= -1, wave2= -1,
                    maxError=1.e-3, method='nsect_recur', log=True, start_x= -1, end_x= -1, to_eval=None,
                    nCut=30, maxIter=20):
        """ 
        Print result of getTemDen function. See getTemDen for more details.
        
        Usage:
            O3.printTemDen(100, tem=10000, wave1=5007, wave2=4363)
        
        Parameters:
            - int_ratio    intensity ratio of the selected transition
            - tem          electronic temperature
            - den          electronic density
            - lev_i1       upper level of 1st transition
            - lev_j1       lower level of 1st transition
            - lev_i2       upper level of 2nd transition
            - lev_j2       lower level of 2nd transition
            - wave1        wavelength of 1st transition
            - wave2        wavelength of 2nd transition
            - maxError     tolerance on difference between input and computed ratio 
            - method       numerical method for finding the root (nsect_recur, nsect_iter)
            - log          log switch (default = True)
            - start_x      lower end of the interval to explore (default: lower end of collision 
                            strength temperature array)
            - end_x        higher end of the interval to explore (default: higher end of collision 
                            strength temperature array)
            - to_eval      expression to be evaluated, using either I (for transitions identified through 
                            atomic levels) or L (for transitions identified through wavelengths)
            - nCut        number of sections in which each step is cut. 2 would be dichotomy.
            - maxIter     maximum number of iterations

        """
        if tem == -1:
            option = 'temperature'
            assume = den
            assume_str = 'density'
            assume_unit = 'cm-3'
            unit = 'K'
        if den == -1:
            option = 'density'
            assume = tem
            assume_str = 'temperature'
            assume_unit = 'K'
            unit = 'cm-3'

        result = self.getTemDen(int_ratio, tem=tem, den=den, lev_i1=lev_i1, lev_j1=lev_j1, lev_i2=lev_i2, lev_j2=lev_j2,
                                wave1=wave1, wave2=wave2, maxError=maxError, method=method, log=log, start_x=start_x,
                                end_x=end_x, to_eval=to_eval, nCut=nCut, maxIter=maxIter)
        print('Ion = %s' % self.elem + " " + int_to_roman(self.spec))
        print('Option = %s' % option)
        print('Assumed %s = %.0f %s' % (assume_str, assume, assume_unit))
        if to_eval is None:
            print('Assumed I(%i)/I(%i) ratio = %.0f' % (wave1, wave2, int_ratio))
        else:
            print('Assumed %s = %.0f' % (to_eval, int_ratio))
#        print 'method = %s' % method
#        print 'maxError on ratio %f' % maxError
        print('Calculated %s = %.0f %s' % (option, result, unit))


    def plotEmiss(self, tem_min=1000, tem_max=30000, ionic_abund=1.0, den=1e3, style='-',
                  legend_loc=4, temLog=False, plot_total=False, plot_only_total=False, legend=True,
                  total_color='black', total_label='TOTAL', ax=None):
        """ 
        Plot the emissivity as a function of temperature of all the lines of the selected atom.  

        Usage: 
            O3.plotEmiss(tem_min=10000, tem_max=20000)
        
        Parameters:
            - tem_min         minimum value of the temperature range to span (default=1000)
            - tem_max         maximum value of the temperature range to span (default=30000)
            - ionic_abund     relative ionic abundance (default = 1.0)
            - den             electron density
            - style           line style of the plot (default: '-' [solid line])
            - legend_loc      localization of the legend (default: 4 = lower right; see plt.legend 
                                for more details)
            - temLog          linear (False) or logarithmic temperature axis (default = False)
            - plot_total      flag to also plot total emissivity (default = False)
            - plot_only_total flag to only plot total emissivity (default = False)
            - legend          flag to place legend (default = True)
            - total_color     color of the total emissivity (default = 'black')
            - total_label     label of the total emissivity (default = 'TOTAL')
            - ax              axis where to send the plot. If None, a new axis is done

        """
        if ax is None:
            f, ax = plt.subplots()
        if not config.INSTALLED['plt']: 
            log_.error('Matplotlib not available, no plot', calling=self.calling + '.plot')
            return None
        tem = np.logspace(np.log10(tem_min), np.log10(tem_max), 1000)
        total_emis = np.zeros_like(tem)
        for wave in self.lineList:
            lev_i, lev_j = self.getTransition(wave)
            if (lev_i <= self.NLevels) and (lev_j <= self.NLevels):
                color = ((np.log10(wave) - np.log10(np.min(self.lineList))) / 
                         (np.log10(np.max(self.lineList)) - np.log10(np.min(self.lineList)))) ** 0.4
                c = cm.jet(color, 1)
                if temLog:
                    x_to_plot = np.log10(tem)
                else:
                    x_to_plot = tem
                y_to_plot = ionic_abund * self.getEmissivity(tem, den, wave=wave)
                if not plot_only_total:
                    if (y_to_plot[0] > 0.):
                        ax.plot(x_to_plot, np.log10(y_to_plot),
                             label='{0:.0f}'.format(wave), color=c, linestyle=style)
                total_emis += y_to_plot
        if plot_total:
            ax.plot(x_to_plot, np.log10(total_emis),
                     label=total_label, color=total_color, linewidth=3, linestyle=style)
        if legend:
            ax.legend(loc=legend_loc)
        if temLog:
            ax.set_xlabel('log(T[K])')
        else:
            ax.set_xlabel('T[K]')
        ax.set_ylabel('log [erg.cm3/s]')
        ax.set_title('Line emissivities')

    def plotGrotrian(self, tem=1e4, den=1e2, thresh_int=1e-3, unit='eV', detailed=False, ax=None):
        """
        Draw a Grotrian plot of the selected atom, labelling only lines above a
        pecified intensity threshold (relative to the most intense line). 
        For ground state levels, the Russell-Saunders term symbol is also given.
        
        Usage:
            O3.plotGrotrian()
        
        Parameters:
            - tem          temperature at which the intensity threshold is to be computed 
            - den          density at which the intensity threshold is to be computed 
            - thresh_int   intensity threshold (relative to the most intense line, default: 1.e-3)
            - unit         one of 'eV' (default), '1/Ang' or 'Ryd'
            - ax           axis where to plot the result

        """
        if ax is None:
            f, ax = plt.subplots()
        if unit not in ['eV', 'Ryd', '1/Ang']:
            log_.error('Unit {0} not available'.format(unit))
            return None
        if not config.INSTALLED['plt']: 
            log_.error('Matplotlib not available, no plot', calling=self.calling + '.plot')
            return None
        color_list = ['b', 'r', 'y', 'c', 'm', 'g']
        energies = self.getEnergy(unit=unit)
         
        # VL 16 Jul 2013 - Detect level inversion, just for warning        
        level_multiplet = [0]
        multiplets = []
        if self.NIST is not None:
            term = []; j = []; en = []
            for item in self.NIST:
                term.append(item[1])
                j.append(item[2])
                en.append(item[3])
            sorted_indx = np.argsort(en) 
            term = np.array(term)[sorted_indx]
            j = np.array(j)[sorted_indx]
            en = np.array(en)[sorted_indx]
            for i in np.arange(1, len(en)):
                if (term[i] == term[i-1]):
                    level_multiplet.append(i)
                else:
                    multiplets.append(level_multiplet) 
                    level_multiplet = [i]
            multiplets.append(level_multiplet) # append last multiplet
# Mistakenly rounds off fractional J values. Corrected 26/12/2014            
#            levelLabels = ['$^{{{0}}}${1}$_{{{2:.0f}}}$'.format(l['term'][0],l['term'][1],l['J']) for l in self.NIST]
            levelLabels = ['$^{{{0}}}${1}$_{{{2}/{3}}}$'.format(l['term'][0], l['term'][1], Fraction(l['J']).numerator, Fraction(l['J']).denominator) for l in self.NIST]
        else:
            delta_e_max = 1.e-5  # Arbitrary limit to define hyperfine structure and identify multiplets 
            for i in np.arange(1, len(energies)):
                if (self.getEnergy()[i] - self.getEnergy()[i - 1] < delta_e_max):
                    # if levels close, append level to multiplet
                    level_multiplet.append(i) 
                else:
                    # if levels separated, add current multiplet to array and start a new multiplet
                    multiplets.append(level_multiplet) 
                    level_multiplet = [i]
            multiplets.append(level_multiplet) # append last multiplet
            warn_label = []
            levelLabels = gsLevelDict[self.gs]
            for i_multi in np.arange(len(multiplets)):
                stat_weights = []
                # check each multiplet for inversion separately
                for i_level in multiplets[i_multi]:
                    latex_sw = levelLabels[i_level]
                    sw = eval(latex_sw.split('$')[-2].replace('_', '').replace('{', '').replace('}', '.') + '*2.+1.')
                    stat_weights.append(sw)
                if (stat_weights != [self.getStatWeight()[k] for k in multiplets[i_multi]]):
                    warn_label.append(i_multi)
                    to_print = '\nLevel inversion in multiplet {0} of {1}'\
                             '\nThe labels of levels {2} are displayed in the wrong order' + \
                             '\nAssumed statistical weights: {3}' + \
                             '\nStatistical weights of data: {4}\n'
                    log_.warn(to_print.format(i_multi + 1, self.atom, multiplets[i_multi], 
                                                 stat_weights, [self.getStatWeight()[k] for k in multiplets[i_multi]]), 
                                 calling=self.calling + '.plot')

        # where level starts in the plot
       
        x_span = 0.75
        x_0 = (1 - x_span) / 4.    
        dx = 0.005
        max_en = np.max(energies)
        blow = 0.02 * max_en # blowup scale for multiplets
        ax.set_xlim((0, 1.))
        ax.set_ylim((np.max(energies) * -0.05, max_en * 1.05))
        ax.set_title('[{0} {1}]'.format(self.elem, int_to_roman(self.spec)))

        # Blow up of multiplets in the plot
        i_tot = 0
        for i_multi in np.arange(len(multiplets)):
            shift = (1-len(multiplets[i_multi])) / 2.
            for i in np.arange(len(multiplets[i_multi])):
                y = energies[multiplets[i_multi][i]] + shift * blow
                ax.plot([x_0, 1.5 * x_0, 1.6 * x_0, x_0 + x_span], 
                         [y, y, energies[multiplets[i_multi][i]], 
                          energies[multiplets[i_multi][i]]], lw = 1.2, 
                         color = color_list[np.mod(multiplets[i_multi][i], len(color_list))])
                ax.text(x_span + x_0 + dx, y, '{0:7.4f}'.format(energies[multiplets[i_multi][i]]), size='small', 
                         horizontalalignment='left', verticalalignment='center')
                try:
                    if (detailed is False):
                        ax.text(dx, y, levelLabels[i_tot], size='medium', verticalalignment='center')
                    else:
                        ax.text(dx, y, str(i_tot+1) + ": " + levelLabels[i_tot] , size='medium', verticalalignment='center', 
                                color = color_list[np.mod(multiplets[i_multi][i], len(color_list))])
                    i_tot += 1
                except:
                    pass
                shift += 1
        ax.text(1 - dx / 2., max(energies) * 0.5, 'E [{0:s}] '.format(unit), horizontalalignment='right', color='blue')
        ax.set_xlabel('Ground-state configuration: {0}'.format(self.gs), color="#004400")
               
        ax.xaxis.set_ticks_position("none")
        ax.yaxis.set_ticks_position("none")
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        
        all_emis = self.getEmissivity(tem, den)
        emis_max = np.max(all_emis)
        N_lines = (self.getEmissivity(1e4, 1e3) > (emis_max * thresh_int)).sum() # number of lines plotted
        x_pad = 0.1 * (x_span - 0.6 * x_0)     
        if N_lines > 1:
            delta_x = (x_span - 2 * x_pad) / (N_lines - 1)
            x = 1.6 * x_0 + 0.5 * x_pad
        else:
            delta_x = 0
            x = (1.6 * x_0 + x_span) / 2.
        cc = colors.ColorConverter()
        for j in np.arange(self.NLevels - 1) + 2:
            for i in np.arange(j - 1) + 1:
                emis = all_emis[j-1, i-1]
                if emis > (emis_max * thresh_int):
                    N_seg = 1000. #number of segments to draw emission line, to make it look smooth
                    xx = np.ones(N_seg+1) * x  # x coord of segmente making up one line
                    yy = np.linspace(self.getEnergy(i, unit=unit), self.getEnergy(j, unit=unit), N_seg+1)
                    scale = np.linspace(0, 1, N_seg+1)
                    alpha =  np.tanh((2*scale-1)+1)  # to make alpha vary as a smooth step function
                    points = np.array([xx, yy]).T.reshape(-1, 1, 2) 
                    segments = np.concatenate([points[:-1], points[1:]], axis=1) 
                    cmap1 = []
                    cmap2 = []
                    for seg in segments:
                            yy0 = seg.mean(0)[1] 
                            alpha0 = np.interp(yy0, yy, alpha)
                            rgb1 = cc.to_rgb(color_list[np.mod(j-1, len(color_list))]) 
                            rgb2 = cc.to_rgb(color_list[np.mod(i-1, len(color_list))]) 
                            cmap1.append([rgb1[0], rgb1[1], rgb1[2], alpha0]) 
                            cmap2.append([rgb2[0], rgb2[1], rgb2[2], 1-alpha0])
                    lc1 = LineCollection(segments, linewidths=3.5)
                    lc1.set_color(cmap1) 
                    lc2 = LineCollection(segments, linewidths=3.5)
                    lc2.set_color(cmap2) 
                    ax.add_collection(lc1) 
                    ax.add_collection(lc2) 
                    if self.wave_Ang[j - 1, i - 1] > 1e4:
                        to_print = '{0:.1f}m'.format(self.wave_Ang[j - 1, i - 1] / 1e4) + ' '
                    else:
                        to_print = '{0:.1f}A'.format(self.wave_Ang[j - 1, i - 1]) + ' '
                    ax.text(x + 0.05 * delta_x, self.getEnergy(j, unit=unit)-blow/2., to_print, size='small', rotation=90, verticalalignment='top')
                    x = x + delta_x

        # issue warning on plot if level inversion 
        try:
            for i in warn_label:
                ax.text(0, self.getEnergy(multiplets[i][0]+1, unit=unit), 'Warning ', ha='right', color='red')
        except:
            pass

    
    def __repr__(self):
        return 'Atom {0}{1} from {2} and {3}'.format(self.elem, self.spec, self.atomFile, self.collFile)


class RecAtom(object):
    
    def __init__(self, elem, spec, atom=None):
        """
        RecAtom class. Used to manage recombination data and compute emissivities.
        
        Usage:
            H1 = pn.RecAtom('H', 1)
        
        Parameters:
            - elem          symbol of the selected element
            - spec          ionization stage in spectroscopic notation (I = 1, II = 2, etc.)
         
        """
        self.log_ = log_
        self.type = 'rec'
        self.gs = None
        if atom is not None:
            self.atom = str.capitalize(atom)
            self.elem = parseAtom(self.atom)[0]
            self.spec = int(parseAtom(self.atom)[1])
        else:
            self.elem = str.capitalize(elem)
            self.spec = int(spec)
            self.atom = self.elem + str(self.spec)
        self.name = sym2name[self.elem]
        self.calling = 'Atom ' + self.atom
        self.log_.message('Making rec-atom object for {0} {1:d}'.format(self.elem, self.spec), calling=self.calling)
        try:
            self.Z = Z[self.elem]
        except:
            self.Z = -1
        if self.elem in IP:
            if self.spec == 1:
                self.IP = 0
            elif self.spec < len(IP[self.elem])+2:
                self.IP = IP[self.elem][self.spec-2]
            else:
                self.IP = -1
        else:
            self.IP = -1
        
        self.recFitsFile = atomicData.getDataFile(self.atom, 'rec')
        file_type = self.recFitsFile.split('.')[-1]
        if file_type == 'fits':
            self._loadFit()
        elif file_type == 'hdf5':
            self._loadHDF5()
        
        if 'trc' in atomicData.getDataFile()[self.atom].keys():
            self._loadTotRecombination()
        
        self.E_in_vacuum = True
        self.NIST = getLevelsNIST(self.atom)
        self.comments = {}
        web = 'Ref. {0} of NIST 2014 (try this: http://physics.nist.gov/cgi-bin/ASBib1/get_ASBib_ref.cgi?db=el&db_id={0}&comment_code=&element={1}&spectr_charge={2}&'
        if self.NIST is not None:
            self.NLevels = len(self.NIST)
            self._Energy = self.NIST['energy'] / 1e8
            self.comments['VACUUM'] = '1'
            self.comments['NOTE'] = 'Energy levels'
            source = '\n    '
            for ref in np.unique(self.NIST['ref']):
                source = source + web.format(ref[1:], self.elem, self.spec) + ')\n  + ' 
            self.comments['SOURCE'] = source[0:-5]
        self.initWaves()


    def _test_lev(self, level):
        """
        Test whether selected level is legal

        Parameters:
            - level        selected atom level

        """       
        if level < -1 or level == 0 or level > self.NLevels:
            self.log_.error('Wrong value for level: {0}, maximum = {1}'.format(level, self.NLevels),
                            calling=self.calling)

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


    def getWave(self, lev_i=None, lev_j=None):
        """
        Return the wavelength of a transition given the levels
        
        Usage:
            He2.getWave(4, 3) 
            
        Parameters:
            - lev_i, lev_j: upper and lower levels of the transition
                
        """ 
        self._test_lev(lev_i)
        self._test_lev(lev_j)
        return(self.wave_Ang[lev_i-1, lev_j-1])
        
    def _loadHDF5(self):
        """
        Method to read the atomic data hdf5 file and store the data
        Called by __init__
        """
        
        if not config.INSTALLED['h5py']:
            log_.error('You need to install h5py', calling=self.calling)
        self.recFitsFile = atomicData.getDataFile(self.atom, 'rec')
        if self.recFitsFile is None:
            log_.error('No hdf5 data for atom: {0}'.format(self.atom), calling=self.calling)
            return None
        self.recFitsFullPath = atomicData.getDataFullPath(self.atom, 'rec')
        try:
            hf5 = h5py.File(self.recFitsFullPath, 'r')
        except:
            log_.error('{0} recombination file not read'.format(self.recFitsFile), calling=self.calling)
        self._RecombData = hf5['updated_data'].value
        hf5.close()
        if self.atom in config.DataFiles:
            if self.recFitsFile not in config.DataFiles[self.atom]:
                config.DataFiles[self.atom].append(self.recFitsFile)
        else:
            config.DataFiles[self.atom] = [self.recFitsFile]
        try:
            self.temp = self._RecombData['TEMP']
        except:
            log_.error('No TEMP field in {0}'.format(self.recFitsFile))
        try:
            self.log_dens = self._RecombData['DENS']
        except:
            log_.error('No DENS field in {0}'.format(self.recFitsFile))
        self.labels = self._RecombData.dtype.names
        self.labels = tuple([l for l in self.labels if l not in ('TEMP', 'DENS')])
        if '_' in self._RecombData.dtype.names[0]:
            self.label_type = 'transitions'
        else:
            self.label_type = 'wavelengths'
        log_.message('{0} recombination data read from {1}'.format(self.atom, self.recFitsFile), calling=self.calling)
       
    def _loadFit(self):
        """
        Method to read the atomic data fits file and store the data
        Called by __init__

        """
        self.recFitsFile = atomicData.getDataFile(self.atom, 'rec')
        if self.recFitsFile is None:
            log_.error('No fits data for atom: {0}'.format(self.atom), calling=self.calling)
            return None
        self.recFitsFullPath = atomicData.getDataFullPath(self.atom, 'rec')
        try:
            hdu = pyfits.open(self.recFitsFullPath)
        except:
            log_.error('{0} recombination file not read'.format(self.recFitsFile), calling=self.calling)
        self._RecombData = hdu[1].data
        hdu.close()
        if self.atom in config.DataFiles:
            if self.recFitsFile not in config.DataFiles[self.atom]:
                config.DataFiles[self.atom].append(self.recFitsFile)
        else:
            config.DataFiles[self.atom] = [self.recFitsFile]
        try:
            self.temp = self._RecombData['TEMP']
        except:
            log_.error('No TEMP field in {0}'.format(self.recFitsFile))
        try:
            self.log_dens = self._RecombData['DENS']
        except:
            log_.error('No DENS field in {0}'.format(self.recFitsFile))
        self.labels = self._RecombData.names
        del self.labels[self.labels.index('TEMP')]
        del self.labels[self.labels.index('DENS')]
        if '_' in self._RecombData.names[0]:
            self.label_type = 'transitions'
        else:
            self.label_type = 'wavelengths'
        log_.message('{0} recombination data read from {1}'.format(self.atom, self.recFitsFile), calling=self.calling)


    def _loadTotRecombination(self):
        """
        Load the total recombination coefficient table. The case (A or B) is set by selecting the corresponding trc file.
            
        Parameters:
            None
                
        """ 
        self.TotRecFile = atomicData.getDataFullPath(self.atom, 'trc')
        f = open(self.TotRecFile)
        data = f.readlines()
        f.close()
        den_points = [np.float(d) for d in data[0].split()]
        tem_points = [np.float(d) for d in data[1].split()]
        self.alpha_grid = np.array([d.split() for d in data if d[0:3]!='***'][2::], dtype='float')
        self.lg_den_grid, self.lg_tem_grid = np.meshgrid(np.log10(den_points), np.log10(tem_points))

        
    def _Transition(self, wave, maxErrorA = 5.e-3, maxErrorm = 5.e-2):
        """
        Return an array with computed upper level, computed lower level, computed wavelength, 
            input wavelength
        
        Parameters:
            - wave       wavelength either in Angstrom (a float or a label: e.g., 5007, '5007A') 
                            or in micron (a label: '51.5m')
            - maxErrorA: tolerance if the input wavelength is in Angstrom
            - maxErrorm: tolerance if the input wavelength is in micron
                            
        """
        if str(wave)[-1] == 'A':
            inputWave = float(wave[:-1])
            label = '{0}_{1}'.format(self.atom, wave)
            maxError = maxErrorA
        elif str(wave)[-1] == 'm':
            inputWave = float(wave[:-1]) * 1e4
            label = '{0}_{1}'.format(self.atom, wave)
            maxError = maxErrorm
        else:
            inputWave = wave
            label = '{0}_{1}A'.format(self.atom, int(wave))
            maxError = maxErrorA
            
        if label in label2levelDict:
            result = [label2levelDict[label][0], label2levelDict[label][1], inputWave, inputWave]
            return(result)
        
        j, i = np.unravel_index(np.argmin(abs(self.wave_Ang - inputWave)), self.wave_Ang.shape)
        bestWave = self.wave_Ang[i, j]
        error = np.abs(bestWave - inputWave) / inputWave
        result = [i + 1, j + 1, bestWave, inputWave]
        if error > maxError:
            self.log_.warn('_Transition: wavelengths differ by more than {0:.2f}%: input = {1:.2f}, output = {2:.2f}'\
                           .format(100 * maxError, inputWave, bestWave), calling=self.calling)
        return(result)


    def getTotRecombination(self, tem, den, method='linear'):
        """
        Return the total recombination coefficient. The case (A or B) is set by selecting the corresponding trc file.
        
        Usage:
            atomicData.setDataFile('h_i_trc_SH95-caseA.dat')
            h1.getTotRecombination(tem=10000, den=5.e3)
            
        Parameters:
            - tem, den  temperature and density
            - method    interpolation method in the grid ('linear' = default, 'nearest', 'cubic')    
                
        """ 
        self.calling = 'getTotRecombination'
        if 'trc' in atomicData.getDataFile()[self.atom].keys():
            return interpolate.griddata((self.lg_den_grid.ravel(), self.lg_tem_grid.ravel()), self.alpha_grid.ravel(), (np.log10(den), np.log10(tem)), method=method)
        else:
            self.log_.warn('No recombination data available for {0} in the adopted dictionary (but data may exist: please query atomicData.getAllAvailableFiles("<atom>") for trc files)'.format(self.atom), calling=self.calling)
            return None
        
        
    def getTransition(self, wave, maxErrorA = 5.e-3, maxErrorm = 5.e-2):
        """
        Return the indexes (upper level, lower level) of a transition for a given atom 
            from the wavelength.
        
        Usage:
            O3.getTransition(4959)
            
        Parameters:
            - wave      wavelength in Angstrom (a float or a label: e.g., 5007, '5007A') 
                or in micron (a label: '51.5m')
            - maxErrorA: tolerance if the input wavelength is in Angstrom
            - maxErrorm: tolerance if the input wavelength is in micron
                
        """ 
        res = self._Transition(wave, maxErrorA = maxErrorA, maxErrorm = maxErrorm)
        return(res[0], res[1])
        
        
    def printTransition(self, wave):
        """
        Print info on transition associated to input wavelength.
        
        Usage:
            O3.printTransition(4959)
            
        Parameters:
            - wave      wavelength in Angstrom (a float or a label: e.g., 5007, '5007A') 
                or in micron (a label: '51.5m')
                
        """
        closestTransition = self._Transition(wave)
        relativeError = closestTransition[3] / closestTransition[2] - 1
        print('Input wave: {0:.1F}'.format(closestTransition[3]))
        print('Closest wave found: {0:.1F}'.format(closestTransition[2]))
        print('Relative error: {0:.0E} '.format(relativeError))
        print('Transition: {0[0]} -> {0[1]}'.format(closestTransition))
        return
    
              
    def printSources(self):
        """
        Print bibliographic sources for atomic data, as listed in the headers of the fits files.
        
        """       
        header = pyfits.open(self.recFitsFullPath, ignore_missing_end=True)[1].header
        for record in header.items():
            if 'SOURCE' in record[0]:
                number = record[0].lstrip('SOURCE')
                try:
                    print(self.atom + ': ' + header.get('NOTE' + str(number)) + ':', header.get('SOURCE' + str(number)))
                except:
                    print(self.atom + ': ' + 'Atomic data:', header.get('SOURCE' + str(number)))


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
                     'Ryd': CST.RYD_ANG,
                     'eV': CST.RYD_ANG * CST.RYD_EV,
                     'cm-1': 1e8}
        if unit not in unit_dict:
            self.log_.warn('Unit {0} unknown, using 1/Ang'.format(unit), calling=self.calling + '.getEnergy')
            unit = '1/Ang'        
        
        if level == -1:
            return self._Energy * unit_dict[unit]
        else:
            return self._Energy[level-1] * unit_dict[unit]

    
    def _checkLabel(self, label):
        """
        Return True if the label is compatible with the type of labels read form the data file
        'transitions' labels are of the form "J_I"
        'wavelengths' labels are of the form "1234.5" (dot is mandatory)

        """
        if ('_' in label) and (self.label_type == 'transitions'):
            return True
        elif ('.' in label) and (self.label_type == 'wavelengths'):
            return True
        else:
            return False
   
            
    def _getLabelStr(self, label, warn=True):
        """
        Returns a string containing the label. 
        If label is a float, it is transormed into a string .
        If label is not in the self.labels list, None is returned

        """
        if np.isreal(label):
            label_str = '{0:.1f}'.format(label)
        else:
            label_str = str(label)
        if label_str not in self.labels:
            if warn:
                log_.warn('Label {0} not in {1}.'.format(label_str, self.recFitsFile), calling=self.calling)
            return None
        else:
            return label_str
        
        
    def getEmissivity(self, tem, den, lev_i=None, lev_j=None, wave=None, label=None,
                      method='linear', product=True):
        """
        Return the emissivity of a recombination line. The arguments used to 
        define the line depend on whether the atom is an hydrogenoid or not. 
        In the first case, the transition can be specified either as a pair 
        of levels lev_i, lev_j or as a label. 
        In the second case, the transition can be specified either as a wavelength 
        or as a label.
        In either case, enter <atom>.labels to display the valid labels.
        
        Usage:
            H1 = pn.RecAtom('H', 1)
            H1.getEmissivity([1e4, 1.2e4], [1e3, 1e2], lev_i = 4, lev_j = 2)
            H1.getEmissivity([1e4, 1.2e4], [1e3, 1e2], label='4_2', product=False)
            tem = np.linspace(5000, 20000, 100)
            den = np.logspace(2, 6, 100)
            imHab = H1.getEmissivity(tem, den, label='3_2') / H1.getEmissivity(tem, den, label='4_2')

            He1 = pn.RecAtom('He', 1)
            He1.getEmissivity(1e4, 1e2, wave=4471.0)
            He1.getEmissivity(1e4, 1e2, label='4471.0')
            
        Parameters:
            - tem            temperature (K)
            - den            density (cm-3)
            - lev_i, lev_j   levels of the transition
            - wave           wavelength of the transition
            - label          label of the transition (e.g. "50_3", "1234.5")
            - method         interpolation method ('linear', 'nearest', 'cubic'), 
                             sent to scipy.interpolate.griddata    
            - product        Boolean. If True (default), all the combination of (tem, den) are used. 
                             If False, tem and den must have the same size and are joined.
        """
        
        if not config.INSTALLED['scipy']:
            log_.error('Scipy not installed, no RecAtom emissivities available',
                          calling=self.calling)
            return None
        tem = np.asarray(tem)
        den = np.asarray(den)
        if product:
            if tem.size == 1 or den.size == 1:
                temg = tem
                deng = den
            else:
                temg, deng = np.meshgrid(tem, den)
        else:
            if tem.size != den.size:
                log_.error('tem and den must have the same size', calling=self.calling)
                return None
            else:
                temg = tem
                deng = den
        if (lev_i is not None) and (lev_j is not None):
            label = '{0}_{1}'.format(lev_i, lev_j)
        if wave is not None:
            label = '{0:.1f}'.format(wave)
            label_str = self._getLabelStr(label, warn=False)
            if label_str is None:
                ij = self.getTransition(wave)
                label = '{}_{}'.format(ij[0], ij[1])
        label_str = self._getLabelStr(label, warn=False)
        if label_str is None:
            log_.warn('Wrong label {0}'.format(label), calling=self.calling)
            return None
        if not self._checkLabel(label_str):
            log_.warn('Wrong label {0}'.format(label_str), calling=self.calling)
            return None
        enu = self._RecombData[label_str]
            
        logd = np.log10(deng)
        temp_min = np.min(self.temp)
        temp_max = np.max(self.temp)
        log_dens_min = np.min(self.log_dens)
        log_dens_max = np.max(self.log_dens)
        tt = (logd < log_dens_min)
        if np.ndim(logd) == 0: 
            if tt == True:
                logd = log_dens_min
        else:
            logd[tt] = log_dens_min
        tt = (logd > log_dens_max)
        if np.ndim(logd) == 0:
            if tt == True:
                logd = log_dens_max
        else:
            logd[tt] = log_dens_max
        res = interpolate.griddata((self.temp.ravel(), self.log_dens.ravel()), enu.ravel(),
                                   (temg, logd), method=method)
        #return res.T
        return res


    def getIonAbundance(self, int_ratio, tem, den, lev_i= -1, lev_j= -1, wave= -1, label=None,
                        to_eval=None, Hbeta=100.):
        """
        Compute the ionic abundance relative to H+ given the temperature, the density and the 
            intensity of a line or sum of lines.
        The arguments used to define the line depend on whether the atom is an hydrogenoid or not. 
        For hydrogenoids, the transition can be specified either as a pair of levels 
            lev_i, lev_j, as a label, or as an I-type expression as the argument of to_eval 
            (e.g. to_eval='I(50, 19)' for the 50->19 transition). Wavelengths or L-type expressions 
            will not work in this case. 
        For non-hydrogenoids, the transition can be specified either as an integer wavelength, 
            as a label, or as an A-type expression as the argument of to_eval (e.g. to_eval='A(4471)' 
            for the lambda=4471 transition). Note that all these alternatives imply that the wavelength
            is known. Pairs of levels and I or L-type expressions would not work.

        The preferred method are the label and the I-type expressions, as the remaining parameters 
            are inherently fragile.

        Usage:
            He2.getIonAbundance(130, 1.5e4, 100., lev_i=5, lev_j=4)
            He2.getIonAbundance(130, 1.5e4, 100., label="5_4")
            He2.getIonAbundance(130, 1.5e4, 100., to_eval='I(4,3) + I(4,2)')
            He1.getIonAbundance(100, 1.5e4, 100., wave=5016)
            He1.getIonAbundance(100, 1.5e4, 100., label="5016.0")
            He1.getIonAbundance(100, 1.5e4, 100., to_eval='A(5016)')
            He1.getIonAbundance(np.array([100, 150]), np.array([1.5e4, 1.2e4]), np.array([100., 120]), 
                label="10830.0")
            
        Parameters:
            - int_ratio    relative line intensity (default normalization: Hbeta = 100). 
                            May be an array.
            - tem          electronic temperature in K. May be an array.
            - den          electronic density in cm^-3. May be an array.
            - lev_i        upper level of transition
            - lev_j        lower level of transition
            - wave         wavelength of transition. Takes precedence on lev_i and lev_j if set, 
                            ignored otherwise 
            - to_eval      expression to be evaluated. Takes precedence on wave if set, 
                            ignored otherwise.
            - Hbeta        line intensity normalization at Hbeta (default Hbeta = 100)

        
        """
        if np.ndim(tem) != np.ndim(den):
            self.log_.error('ten and den must have the same shape', calling=self.calling)
            return None
        if ((np.squeeze(np.asarray(int_ratio)).shape != np.squeeze(np.asarray(tem)).shape) | 
            (np.squeeze(np.asarray(den)).shape != np.squeeze(np.asarray(tem)).shape)):
            self.log_.warn('int_ratio, tem and den must does not have the same shape', calling=self.calling)
        if (lev_i == -1) & (lev_j == -1) & (wave == -1) & (to_eval is None) & (label is None):
            self.log_.error('At least one of lev_i, lev_j, or wave, or to_eval must be supplied', calling=self.calling)
            return None
        if to_eval == None:
            if wave != -1:
                to_eval = 'L({0})'.format(wave)
            elif label is not None:
                to_eval = 'A("{0}")'.format(label)
            else:
                to_eval = 'I({0}, {1})'.format(lev_i, lev_j)
        I = lambda lev_i, lev_j: self.getEmissivity(tem, den, lev_i, lev_j, product=False)
        L = lambda wave: self.getEmissivity(tem, den, wave=wave, product=False)
        A = lambda label: self.getEmissivity(tem, den, label=label, product=False)
        try:
            emis = eval(to_eval)
        except:
            self.log_.error('Unable to eval {0}'.format(to_eval), calling=self.calling)
            return None
        #int_ratio is in units of Hb = Hbeta keyword
        ionAbundance = ((int_ratio / Hbeta) * (getRecEmissivity(tem, den, 4, 2, atom='H1', product=False) / emis))
        return ionAbundance
    
    
    def __repr__(self):
        return 'Atom {0}{1} from {2}'.format(self.elem, self.spec, self.recFitsFile)
        

def getRecEmissivity(tem, den, lev_i=None, lev_j=None, atom='H1', method='linear', wave=None, product=True):
    """
    The function instantiates a RecAtom and store it into atomicData._RecombData for a further use.
    More possibilities are obtained using the RecAtom class.

    Usage:
        print(getRecEmissivity(1e4, 1e3, 3, 2, atom='H1') / getRecEmissivity(1e4, 1e3, 4, 2, atom='H1')) 
            # return Ha/Hb 

    Parameters:
        - tem           temperature in K
        - den           density in cm-3
        - lev_i, lev_j  levels (lev_i>lev_j, i_min=1)
        - atom          atom (e.g. 'H1', 'He1')
        - method        interpolation method ('linear', 'nearest', 'cubic'), sent to scipy.interpolate.griddata
        - wave          alternative way of identifying emision line.
        
    """
    calling = 'getRecEmissivity'
    
    if config.INSTALLED['scipy']:
        elem, spec = parseAtom(atom)
    
        if atom not in atomicData._RecombData:
            atomicData._RecombData[atom] = RecAtom(elem, spec)
        return atomicData._RecombData[atom].getEmissivity(tem=tem, den=den, lev_i=lev_i, lev_j=lev_j,
                                                             method=method, wave=wave, product=product)
    else:
        if (atom == 'H1') and (lev_i == 4) and (lev_j == 2):
            log_.warn('Scipy is missing, {0} returning Hbeta'.format(calling), calling)
            return getHbEmissivity(tem)
        else:
            log_.error('Only Hbeta emissivity available, as scipy not installed', calling)


def getHbEmissivity(tem= -1):
    """ 
    Compute Hbeta emissivity in erg/s/(N(H+)*N(e-)) for a given temperature with the formula 
        by Aller (1984)

    Usage:
        getHbemissivity(tem=10000)

    Parameters:
        - tem     electronic temperature in K

    """ 
    tem4 = np.asarray(tem) * 1.0e-4
    j_hb = 1.387 / pow(tem4, 0.983) / pow(10., 0.0424 / tem4) * 1.e-25

    # Remove jhb for tem4 < 0 or tem4 > 1e2
    ((tem4 < 0.) | (tem4 > 1e2)).choose(j_hb, -1)

    return j_hb


def getAtomDict(atom_list=None, elem_list=None, spec_list=None, **kwargs):
    """ 
    Initializes all atoms, according to the atomic files available.
    The elem objects are given conventional names elem+spec (e.g., O III is O3)

    Usage:
        all = pn.getAtomDict()
        print(all['S2'].name)
        some = pn.getAtomDict(elem_list=['C', 'N', 'O'])
        some = pn.getAtomDict(atom_list=['O2', 'O3', 'Ar3'])
    
    Parameters:
        - atom_list     a list of the ions for which the elem is to be computed 
                        Takes precedence on elem_list and spec_list
        - elem_list     a list of all the elements for which the elem is to be computed (all by default)
        - spec_list     a list of the spectra for which the elem is to be computed (all by default)
        _ **kwargs      argumentas passed to Atom, e.g. OmegaInterp

    """ 
    all_atoms = {}

    if atom_list is None:       
        if spec_list is None:
            spec_list = SPEC_LIST
        if elem_list is None:
            atom_list = atomicData.getAllAtoms()
        else:
            atom_list = []
            for elem in elem_list:
                for spec in spec_list:
                    atom_list.append(elem + str(spec)) 

    for atom in atom_list:
        elem, spec = parseAtom(atom)        
        try:
            all_atoms[atom] = Atom(elem, spec, **kwargs)
            log_.message('Including ' + atom, calling='getAtomDict')
        except:
            log_.message(atom + ' not found', calling='getAtomDict')
    return all_atoms


def getLineLabel(elem, spec, wave, blend=False):
    """
    Build a line label in the standard PyNeb format. Return atom_label, wave_label, and line_label 
            (strings representing the atom and wave fragments of the line label and the 
                complete line label)
    
    Parameters:
        - elem      symbol of the selected element
        - spec      ionization stage in spectroscopic notation (I = 1, II = 2, etc.)
        - wave      wavelength of the line
        - blend     blend flag (default = False)

    @see parseLineLabel
    
    """    
    if isinstance(wave, str):
        wave_label = wave
    else:
        wave4 = wave / 1.e4
        if wave < 10000.:
            wave_label = "{0:.0f}A".format(wave)
        else:
            wave_label = "{0:.1f}m".format(wave4)
    atom_label = elem + str(spec)
    if blend:
        blend_flag = '+'
    else:
        blend_flag = ''
    line_label = atom_label + '_' + wave_label + blend_flag

    return atom_label, wave_label, line_label


def parseLineLabel(lineLabel):
    """
    Parse the line label to extract the substrings referring to the atom (elem, spec and atom)
    and the numerical value of the wave
    
    Parameters:
        - label    line label in the standard PyNeb format
    
    Returns:
        - elem, spec, atom_label, wave, wave_label, blend  strings containing the elem, spec, atom and 
                                                              wave; numerical value of the wave, 
                                                              in Angstrom; and blend flag 
    """
    ##
    # @todo maybe rearrange the order so 1) it is compatible with getLineLabel, or 2) it lists all the strings first 

    # determine if the line is a blend or not
    if lineLabel[-1] == '+'  or lineLabel[-2] == '+':
        blend = True
    else:
        blend = False
    # extract information on the atom
    elem_spec = strExtract(lineLabel, ' ', '_')
    elem, spec = parseAtom(elem_spec)
    atom_label = elem + str(spec)
    if len(elem_spec) > 0:
        if elem_spec[-1] == 'r':
            atom_label += 'r'
    # extract information on the wave
    if blend:
        try:
            wave = float(strExtract(lineLabel[:-2], '_', ' '))
        except:
            wave = 0.
    else:
        try:
            wave = float(strExtract(lineLabel[:-1], '_', ' '))
        except:
            wave = 0.
    wave_label = strExtract(lineLabel, '_', ' ')
    if lineLabel[-1] == 'm':
        wave = wave * 1.e4
    
    return elem, spec, atom_label, wave, wave_label, blend

    
def isValid(line_label):
    """
    Return True if the line label correspond to a valid line from LINE_LABEL_LIST or BLEND_LIST
    
    Usage:
        isValid('O3_5007A')

    Parameters:
        - line_label    label to be tested

    """
    elem, spec, atom, wave, waveLabel, blend = parseLineLabel(line_label)
    if atom in LINE_LABEL_LIST:
            if (waveLabel in LINE_LABEL_LIST[atom]) or (line_label in BLEND_LIST):
                is_valid = True
            else:
                is_valid = False
    else:
        is_valid = False
    return is_valid


class EmissionLine(object):
    """
    Define the emission line object, which is defined by the line parameters and the intensity 
        parameters.
    The line parameters define the emitting ion, the wavelength or the transition, the label 
        in PyNeb format, and a flag defining whether the line is a blend or a single transition. 
    The intensity parameters describe the observed intensity, the observed uncertainty and 
        the corrected uncertainty
        
    Usage:
        line = pn.EmissionLine('O', 3, 5007, obsIntens=[1.4, 1.3])
        line = pn.EmissionLine(label = 'O3_5007A', obsIntens=320, corrected = True)
    
    Parameters:
        - elem        symbol of the selected element
        - spec        ionization stage in spectroscopic notation (I = 1, II = 2, etc.)
        - wave        wavelength of the line
        - blend       blend flag (boolean)
        - to_eval     algebraic expression describing the emission line in terms of single transitions
        - label       line label in the standard PyNeb format
        - obsIntens   observed intensity
        - obsError    uncertainty on the observed intensity
    
    """ 
    def __init__(self, elem=None, spec=None, wave=None, blend=False, to_eval=None, label=None,
                 obsIntens=None, obsError=None, corrected=False, _unit=None):
        
        self.log_ = log_ 
        self.calling = 'EmissionLine'
        self.corrected = corrected
        
        if label is None:
            self.elem = elem
            self.spec = int(spec)
            self.atom = self.elem + str(self.spec)
            self.wave = wave
            self.blend = blend
            self.atom, self.waveLabel, self.label = getLineLabel(elem, spec, wave, blend)
        else:
            self.label = label
            self.elem, self.spec, self.atom, self.wave, self.waveLabel, self.blend = parseLineLabel(label)
            if wave is not None:
                self.wave = wave

        if self.atom in LINE_LABEL_LIST:
            if (self.waveLabel in LINE_LABEL_LIST[self.atom]) or (self.label in BLEND_LIST):
                self.is_valid = True
                if to_eval is None:
                    if self.blend:
                        self.to_eval = BLEND_LIST[self.label]
                    else: 
                        self.to_eval = 'L(' + str(self.wave) + ')'
                else:
                    self.to_eval = to_eval
            else:
                self.is_valid = False
                self.to_eval = None
                self.log_.warn('line {0} for atom {1} not valid'.format(self.waveLabel, self.atom), calling=self.calling)
        else:
                self.is_valid = False
                self.to_eval = None
                self.log_.warn('Atom {0} not valid'.format(self.atom), calling=self.calling)
            
            
        self.obsIntens = np.asarray(obsIntens, dtype=float)
        if self.corrected:
            self.corrIntens = np.asarray(obsIntens, dtype=float)
        else:
            self.corrIntens = np.zeros_like(obsIntens)
            
        # the following is not public, as we still have to think about it. 
        # Don't une it...
        self._obsIntens_n = np.zeros_like(obsIntens)
        self._corrIntens_n = np.zeros_like(obsIntens)
        self._unit = _unit

        if obsError is None:
            self.obsError = np.zeros_like(self.obsIntens)
        else:
            self.obsError = np.asarray(obsError, dtype=float)

        if self.corrected:
            self.corrError = np.asarray(self.obsError, dtype=float)
        else:
            self.corrError = np.zeros_like(self.obsError)
    ##            
    # @var elem
    # Symbol of the element

    
    def correctIntens(self, RC, normWave=None):
        """
        Correct from extinction. The corrIntens and corrError values of the line 
        are updated according to the RedCorr object


        Parameters:
            - RC        an instantiation of the pn.RedCorr class 
            - normWave  a wavelength for the normalisation of the correction, e.g. 4861.  

        """
        if not isinstance(RC, RedCorr):
            log_.error('Trying to correct with something that is not a RedCor object',
                          calling=self.calling)
            return None
        if self.wave > 0.0:
            self.corrIntens = self.obsIntens * RC.getCorr(self.wave, normWave)     
        else:
            self.corrIntens = self.obsIntens
        self.corrError = self.obsError
        
        
    def addObs(self, newObsIntens, newObsError=None):
        """
        Add observed values to an existing line. 
        
        Parameters:
            - newObsIntens    observed intensity of the line
            - newObsError     error on the observed intensity (optional)
            
        """
        self.obsIntens = np.append(self.obsIntens, newObsIntens)
        if newObsError is None:
            self.obsError = np.append(self.obsError, np.zeros_like(newObsIntens))
        else:
            self.obsError = np.append(self.obsError, newObsError)
        if self.corrected:
            self.corrIntens = np.append(self.corrIntens, newObsIntens)
            if newObsError is None:
                self.corrError = np.append(self.corrError, np.zeros_like(newObsIntens))
            else:
                self.corrError = np.append(self.corrError, newObsError) 
        else:
            self.corrIntens = np.append(self.corrIntens, np.zeros_like(newObsIntens))
            self.corrError = np.append(self.corrError, np.zeros_like(newObsIntens))
    
    
    def printLine(self):
        """
        Provide information on the line: atom, label, to_eval, as well as the intensities and errors
        
        Usage:
            line.printLine()
            
        """
        print("""Line {0.atom} {0.label} evaluated as {0.to_eval}
Observed intensity: {0.obsIntens}
Observed error: {0.obsError}
Corrected intensity: {0.corrIntens}
Corrected error: {0.corrError}""".format(self))
     
     
    def __repr__(self):
        return 'Line {0.atom} {0.label}'.format(self)


class Observation(object):
    def __init__(self, obsFile=None, fileFormat='lines_in_cols', delimiter=None, err_default=0.10,
                 corrected=False, errIsRelative=True):
        """
        Define the observation object, which is a collection of observated intensities of one or more
        emission lines for one or more objects, with the corresponding errors.
        The observed intensities are read from a file or filled in by the addLine method.
        Includes an extinction correction object (pyneb.RedCorr) as Observation.extinction.
        
        Parameters:
            - obsFile       name of the file containing the observations. May be a file object or a file name 
            - fileFormat    lines_in_cols or lines_in_row depending on how the wavelengths are ordered 
                                (default: lines_in_cols)
            - delimiter     character separating entries 
            - err_default   default uncertainty assumed on intensities. Will overwrite the error from the file.
            - corrected     Boolean. True if the observed intensities are already corrected from extinction
                                (default: False)
            - errIsRelative Boolean. True if the errors are relative to the intensities, False if they
                                are in the same unit as the intensity (default: True)

        Example:
            Read a file containing corrected intensities:
                obs = pn.Observation('obs.dat', corrected = True)
            to obtain a dictionary with the observed  corrected intensities:
                i_cor = {label: obs.getLine(label = label).corrIntens for label in obs.lineLabels}
        
        """        
        self.log_ = log_ 
        self.calling = 'Observation'
        self.lines = []
        self.names = []
        self.extinction = RedCorr()
        self.corrected = corrected
        if self.corrected:
            self.extinction.law = 'No correction'
        if obsFile is not None:
            self.readData(obsFile=obsFile, fileFormat=fileFormat, delimiter=delimiter,
                          err_default=err_default, corrected=corrected, errIsRelative=errIsRelative)
    ##            
    # @var log_
    # myloggin object
    # @var extinction
    # RedCor object


    def addLine(self, line):
        """
        Add a line to an existing observation

        Parameters:
            - line    the selected emission line (an instance of EmissionLine)
            
        """
        if not isinstance(line, EmissionLine):            
            self.log_.error('Trying to add an inappropriate record to observations', calling=self.calling)
            return None
        if self.corrected:
            line.corrected = True
            self.correctData(line)
        self.lines.append(line)

    def fillObs(self, lineLabel, default=np.nan):
        """
        Create a fake observation of a given line, filled with a given value.
        Parameters:
            - lineLabel: the label of the new line. If the label corresponds to an already 
                defined observation, nothing is done and a warning is issued.
            - default: the value of the fake observations. Default is np.nan
        """

        if lineLabel not in self.lineLabels:
            newLine = EmissionLine(label=lineLabel, obsIntens=default*np.ones(self.n_obs))
            self.addLine(newLine)
        else:
            log_.warn('Line {0} already in obs'.format(lineLabel), calling = self.calling)

    def addObs(self, name, newObsIntens, newObsError=None):
        """
        Add an observation (i.e. a list of intensities corresponding to a new object) to the existing set.

        Parameters:
            - name            name of the new observation/object
            - newObsIntens    value(s) of the line intensities. Length must match Observation.n_lines
                
        """
        if np.ndim(newObsIntens) == 0:
            newObsIntens = [newObsIntens]
        if len(newObsIntens) != self.n_lines:
            self.log_.error('Length of observations to be added does not match n_lines = {0}'.format(self.n_lines),
                            calling=self.calling)
            return
        if name in self.names:
            self.log_.error('Name {0} already exists'.format(name))
            return
        for i, line in enumerate(self.lines):
            if newObsError is None:
                line.addObs(newObsIntens[i])
            else:
                line.addObs(newObsIntens[i], newObsError[i])
        self.names.append(name)

    @property
    def lineLabels(self):
        """
        Property
        Array of labels of the lines 
        
        """
        return np.asarray([line.label for line in self.lines])   

    @property
    def n_lines(self):
        """
        Property
        Number of lines
        
        """
        return len(self.lines)

    @property
    def n_valid_lines(self):
        """
        Property
        Number of valid lines (i.e., lines with labels recognized by PyNeb)
        
        """
        return len([line for line in self.lines if line.is_valid])

    @property
    def n_obs(self):
        """
        Property
        Number of observations. If the number of observations varies from one line to the other,
            returns the number of observations for each line as an array.
            
        """
        n_obs = np.asarray([l.obsIntens.size for l in self.lines])
        for n in n_obs:
            if n != n_obs[0]:
                return n_obs
                
        return n_obs[0]
    
    
    @property
    def n_obs_origin(self):
        """
        Property
        Number of observations which are not from MonteCarlo (i.e. without -MC- in the name)
        """
        return len([n for n in self.names if '-MC-' not in n])

    def getLine(self, elem=None, spec=None, wave=None, label=None, blend=False, i=None, j=None):
        """
        Return the lines corresponding to elem-spec-wave or to the label.
        
        """
        if label is None:
            label = getLineLabel(elem, spec, wave, blend)[2]
        lines = [line for line in self.lines if line.label == label.strip()]
        n_lines = len(lines)
        if n_lines == 0:
            self.log_.warn('No line for {0} from {1}{2} at wavelength {3} (blend={4})'.format(label, elem, spec, wave, blend),
                           calling=self.calling)
            return None
        elif n_lines == 1:
            return lines[0]
        else:
            return np.asarray(lines)
 

    def getSortedLines(self, crit='atom'):
        """
        Return a list of lines sorted by atoms or wavelengths.
        
        Parameters:
            - crit   criterion to sort the line list ('atom' [default] or 'wave')

        """
        if crit == 'atom':
            return sorted(self.lines, key=lambda line: line.atom + str(line.wave))
        elif crit == 'wave':
            return sorted(self.lines, key=lambda line: line.wave)
        else:
            self.log_.error('crit = {0} is not valid'.format(crit), calling=self.calling + '.getSortedLines')

    
    def getUniqueAtoms(self):
        """
        Return a numpy.ndarray of the atoms of the observed lines. If an atom emits 
        more than one line, it is returned only once (numpy.unique is applied 
        to the list before returning).

        """
        return np.unique([l.atom for l in self.lines])

    
    def readData(self, obsFile, fileFormat='lines_in_cols', delimiter=None, err_default=0.10, corrected=False,
                 errIsRelative=True):
        """
        Read observational data from an ascii file. The lines can be listed either in columns or in rows
        and the observed objects vary in the other direction. The uncertainty on the line intensities
        can be read from the file, or a constant relative value can be assumed.
        The lines must be identified by a label in PyNeb's format ion_wave (e.g., 'O3_5007'); the list of ions and
        corresponding wavelengths can also be found in pn.LINE_LABEL_LIST.
        The following optional fields may also be included (without quotes): 'NAME' (object's name, in one string), 
        'E(B-V)', 'cHbeta', and 'e' (observational error).
        
 
        Parameters:
            - obsFile        file containing the observations. May be a file object or a string 
            - fileFormat     emission lines vary across columns ('lines_in_cols', default) or 
                                across rows ('lines_in_rows'), or across rows with errors in columns 
                                ('lines_in_rows_err_cols')m in which case the column label must start with "err"
            - delimiter      field delimiter (default: None)  
            - err_default    default uncertainty on the line intensities
            - corrected      Boolean. True if the observed intensities are already corrected from extinction
                                 (default: False)
            - errIsRelative  Boolean. True if the errors are relative to the intensities, False if they
                                 are in the same unit as the intensity (default: False)

        """    
        format_list = ['lines_in_cols', 'lines_in_cols2', 'lines_in_rows', 'lines_in_rows_err_cols']
        if fileFormat not in format_list:
            self.log_.error('unknown format {0}'.format(fileFormat), calling='Observation.readData')

        if type(obsFile) is file:
            f = obsFile
            closeAfterUse = False
        else:
            f = open(obsFile, 'r')
            closeAfterUse = True
            
        if fileFormat == 'lines_in_cols':
            hdr = f.readline()
            labels = hdr.split(delimiter)
            labels = [l.strip() for l in labels]
            data = f.readlines()
            if closeAfterUse:
                f.close()
            self.names = [dd.split(delimiter)[0].strip() for dd in data]
            data_tab = np.asarray([[dd.split(delimiter)[i] for dd in data] for i in np.arange(len(labels))])
            
            for i, label in enumerate(labels):
                if label == 'NAME':
                    pass
                elif label == 'cHbeta':
                    self.extinction.cHbeta = data_tab[i].astype(np.float32)
                elif label == 'E(B-V)':
                    self.extinction.E_BV = data_tab[i].astype(np.float32)         
                elif label[-1] != 'e':
                    intens = data_tab[i].astype(np.float32)
                    try:
                        i_error = labels.index(label + 'e')
                        error = data_tab[i_error].astype(np.float32)
                        if not errIsRelative:
                            error = quiet_divide(error, intens)
                    except:
                        self.log_.message('No error found for line {0}'.format(label), calling=self.calling)
                        error = data_tab[1].astype(np.float32) * 0. + err_default
                    try:
                        line2add = EmissionLine(label=label, obsIntens=intens, obsError=error)
                    except:
                        self.log_.warn('Unknown line label {0}'.format(label), calling=self.calling)
                    try:
                        self.addLine(line2add)
                        self.log_.message('adding line {0}'.format(label), calling=self.calling)
                    except:
                        self.log_.warn('Impossible to add line'.format(label), calling=self.calling)
                        
        elif fileFormat == 'lines_in_cols2':
            if closeAfterUse:
                f.close()
            data_tab = np.genfromtxt(obsFile, dtype=None, delimiter=delimiter, names=True, deletechars='')
            for label in data_tab.dtype.names:
                if label == 'cHbeta':
                    self.extinction.cHbeta = data_tab[label]
                elif label == 'E(B-V)':
                    self.extinction.E_BV = data_tab[label]
                elif label == 'NAME':
                    self.names = list(data_tab[label])
                elif label == 'NAME2':
                    try:
                        names2 = list(data_tab[label])
                        self.names = [n1 + '_' + n2 for n1, n2 in zip(self.names, names2)]
                    except:
                        pass
                elif label[-1] != 'e':
                    if data_tab[label].dtype.type != np.string_:
                        intens = data_tab[label]
                        try:
                            error = data_tab[label + 'e']
                            if not errIsRelative:
                                error = error / intens
                        except:
                            self.log_.message('No error found for line {0}'.format(label), calling=self.calling)
                            error = np.ones_like(data_tab[label]) * err_default
                        try:
                            line2add = EmissionLine(label=label, obsIntens=intens, obsError=error)
                        except:
                            self.log_.warn('unkown line label {0}'.format(label), calling=self.calling)
                            print(label, intens, error)
                        try:
                            self.addLine(line2add)
                            self.log_.message('adding line {0}'.format(label), calling=self.calling)
                        except:
                            self.log_.warn('Impossible to add line'.format(label), calling=self.calling)
                            print(label, intens, error)
                    else:
                        self.log_.warn('Skipped {0}'.format(label), calling=self.calling)
            
        elif fileFormat == 'lines_in_rows':
            hdr = f.readline()
            self.names = hdr.split(delimiter)[1:]
            self.names =  [l.strip() for l in self.names]
            data = f.readlines()
            if closeAfterUse:
                f.close()

            labels = [dd.split(delimiter)[0].strip() for dd in data if len(dd.strip()) > 0]
            data_tab = np.asarray([[dd.split(delimiter)[i + 1] for dd in data if len(dd.strip()) > 0] for i in np.arange(len(self.names))])
            data_tab = data_tab.astype(np.float32)
            for i, label in enumerate(labels):
                if label == 'cHbeta':
                    self.extinction.cHbeta = data_tab[:, i]
                elif label == 'E(B-V)':
                    self.extinction.E_BV = data_tab[:, i]
                elif label[-1] != 'e':
                    intens = data_tab[:, i]
                    try:
                        i_error = labels.index(label + 'e')
                        error = data_tab[:, i_error]
                        if not errIsRelative:
                            error = error / intens
                    except:
                        self.log_.message('No error found for line {0}'.format(label), calling=self.calling)
                        error = np.ones_like(data_tab[:, 1]) * err_default
                    try:
                        line2add = EmissionLine(label=label, obsIntens=intens, obsError=error)
                    except:
                        self.log_.warn('unkown line label {0}'.format(label), calling=self.calling)
                        print(label, intens, error)
                    try:
                        self.addLine(line2add)
                        self.log_.message('adding line {0}'.format(label), calling=self.calling)
                    except:
                        self.log_.warn('Impossible to add line'.format(label), calling=self.calling)
                        print(label, intens, error)
                        
        elif fileFormat == 'lines_in_rows_err_cols':
            
            if closeAfterUse:
                f.close()
                
            data_tab = np.genfromtxt(obsFile, dtype=None, delimiter=delimiter, names=True)
            self.names = [name for name in data_tab.dtype.names[1::] if name[0:3] != 'err']
            error_names = [name for name in data_tab.dtype.names if name[0:3] == 'err']
            if len(self.names) != len(error_names):
                log_.error('Number of columns for intensities <> number of columns for errors',
                              calling=self.calling)
                return None
            #names_locations = [name in self.names for name in data_tab.dtype.names]
            #errors_locations = [name[0:3] == 'err' for name in data_tab.dtype.names]
            for i, label in enumerate(data_tab['LINE']):
                label = label.strip()
                if label == 'cHbeta':
                    self.extinction.cHbeta = np.array([data_tab[i][name] for name in self.names])
                elif label == 'E(B-V)':
                    self.extinction.E_BV = np.array([data_tab[i][name] for name in self.names])
                else:
                    intens = np.array([data_tab[i][name] for name in self.names])
                    error = np.array([data_tab[i][name] for name in error_names])
                    if not errIsRelative:
                        error = error / intens
                    try:
                        line2add = EmissionLine(label=label, obsIntens=intens, obsError=error)
                    except:
                        self.log_.warn('unkown line label {0}'.format(label), calling=self.calling)
                        print(label, intens, error)
                    try:
                        self.addLine(line2add)
                        self.log_.message('adding line {0}'.format(label), calling=self.calling)
                    except:
                        self.log_.warn('Impossible to add line'.format(label), calling=self.calling)
                        print(label, intens, error)
        
            
        if corrected:
            self.correctData()
            
            
    def getIntens(self, returnObs=False, obsName=None):
        """
        Return the line intensities in form of a dictionary with line labels as keys.
        
        Parameters:
            - returnObs  Boolean. If False (default), prints the corrected values. 
                            If True, prints the observed value. 
            - obsName    name of an observation. If not set or None, all the observations are printed

        """
        if obsName is not None:
            if obsName in self.names:
                obsIndex = self.names.index(obsName)
            else:
                self.log_.error('Name {} is not an Observation name'.format(obsName))
                return None
        else:
            obsIndex = np.arange(self.n_obs)
        to_return = {}
        for line in self.lines:
            if returnObs:
                to_return[line.label] = line.obsIntens[obsIndex]
            else:
                to_return[line.label] = line.corrIntens[obsIndex]
        return to_return
    
    
    def getError(self, returnObs=False, obsName=None):
        """
        Return the line intensity error in form of a dictionary with line labels as keys.
        
        Parameters:
            - returnObs  if False (default), prints the corrected values. 
                            If True, prints the observed value. 
            - obsName    name of an observation. If not set or None, all the observations are printed
            
        """
        if obsName is not None:
            if obsName in self.names:
                obsIndex = self.names.index(obsName)
            else:
                self.log_.error('Name {} is not an Observation name'.format(obsName))
                return None
        else:
            obsIndex = np.arange(self.n_obs)
        to_return = {}
        for line in self.lines:
            if returnObs:
                to_return[line.label] = line.obsError[obsIndex]
            else:
                to_return[line.label] = line.corrError[obsIndex]
        return to_return
    
    
    def printIntens(self, returnObs=False, obsName=None):
        """
        Print the line intensities.
        
        Parameters:
            - returnObs   if False (default), prints the corrected values. 
                            If True, prints the observed value. 
            - obsName     name of an observation. Is unset or None, all the observations are printed

        """    
        if obsName is not None:
            if obsName in self.names:
                obsIndex = self.names.index(obsName)
            else:
                self.log_.error('Name {} is not an Observation name'.format(obsName))
                return None
        else:
            obsIndex = np.arange(self.n_obs)
        for line in self.lines:
            if returnObs:
                print(line.label, line.obsIntens[obsIndex])
            else:
                print(line.label, line.corrIntens[obsIndex])
    
            
    def def_EBV(self, label1="H1_6563A", label2="H1_4861A", r_theo=2.85):
        """
        Define the extinction parameter using the ratio of 2 lines.
        Calls extinction.setCorr to set the EBV and cHbeta according to the parameters.
        Once this is done, one may call correctData to compute the EmissionLine.corrIntens
        
        Parameters:
            - label1 and label2 [EmissionLine.label] observed line whose intensities are used
            - r_theo [float] theoretical line ratio

        """
        line1 = self.getLine(label=label1)
        line2 = self.getLine(label=label2)
        if line1 is None:
            log_.error('{0} is not a valid label or is not observed'.format(line1), calling=self.calling)
            return None
        if line2 is None:
            log_.error('{0} is not a valid label or is not observed'.format(line2), calling=self.calling)
            return None

        obs_over_theo = (line1.obsIntens / line2.obsIntens) / r_theo 
        self.extinction.setCorr(obs_over_theo, line1.wave, line2.wave)
        
 
    def __normalize(self, label='H1_4861A'):
        """
        Normalize the line intensities to a reference line (Hbeta by default).
        Not yet implemented
        
        """
        if "=" in label:
            line_label, factor = label.split('=')
            factor = np.float(factor)
        else:
            line_label = label
            factor = 1.
        line_norm = self.getLine(label=line_label) 
        if line_norm is None:
            log_.warn('No normalization possible as {0} not found'.format(line_label), calling=self.calling)
            return None
        for line in self.lines:
            line._obsIntens_n = line.obsIntens / (line_norm.obsIntens * factor)
            line._unit = line_label.strip()
    
    
    def correctData(self, line=None, normWave=None):
        """
        Correct the line intensities with the correction computed with the RedCorr class (extinction.py)
        The result is stored in line.corrIntens (corrected intensity in absolute units).

        """
        if line is None:
            line = self.lines
        if np.ndim(line) != 0:
            for l in line:
                self.correctData(l, normWave=normWave)
        else:
            if not isinstance(line, EmissionLine):
                log_.error('Trying to correct something that is not a line', calling=self.calling)
                return None  
            line.correctIntens(self.extinction, normWave=normWave)

            """            
            # the following is commented for now, we'll see latter how to implement normalized intensities.
            if line._unit is not None:
                try:
                    line_norm = self.getLine(label=line._unit)
                    line._corrIntens_n = line.obsIntens_n * self.extinction.getCorr(line.wave, line_norm.wave)
                except:
                    log_.warn('No normalized correction for {0}'.format(line), calling='Observation.correctData')
                    line._corrIntens_n = None
            """
            
    def setAllErrors(self, err_default):
        """
        Set the relative uncertainty of all emission lines to a common constant value
        
        Parameters:
            - err_default     default value of the relative uncertainty
        
        """
        for line in self.lines:
            line.obs_err = np.ones_like(line.obs_err) * err_default

    
    
    def addMonteCarloObs(self, N=0, i_obs=None):
        """
        Adding MonteCarlo random-gauss values of fake observations to an obs object.
        The names of the fake observations will be OriginalName-MC-n, n ranging from 0 to N-1
        
        Parameters:
        N: number of new observations to be added for each original observation.
        i_obs: used in case only a given observations needs to be treated
        """
        n_lines = self.n_lines
        n_obs = self.n_obs
        if i_obs is None:
            log_.message('Entering', calling='addMonteCarloObs')
            for i in xrange(n_obs):
                self.addMonteCarloObs(i_obs=i, N=N)
            log_.message('Leaving', calling='addMonteCarloObs')
        else:
            intens = np.array([self.getIntens()[label] for label in self.lineLabels])[:,i_obs] # n_lines
            error = np.array([self.getError()[label] for label in self.lineLabels])[:,i_obs]
            all_new_obs = np.random.standard_normal((N, n_lines))
            
            for i in xrange(N):
                new_obs = intens * (all_new_obs[i,:] * error + 1)
                new_obs[new_obs < 0.] = 0.
                self.addObs('{0}-MC-{1}'.format(self.names[i_obs], i), new_obs, error)
    
