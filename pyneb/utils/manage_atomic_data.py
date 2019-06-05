import os
import pyneb as pn
import numpy as np
import re
from .misc import execution_path, parseAtom, roman_to_int, multi_split, int_to_roman, strExtract
from .init import ELEM_LIST, SPEC_LIST
from .physics import _predefinedDataFileDict

def atom2chianti(atom):
    elem = parseAtom(atom)[0]
    spec = int(parseAtom(atom)[1])
    return '{0}_{1}'.format(elem.lower(), spec)

class _ManageAtomicData(object):
    
    def __init__(self):
        self.log_ = pn.log_
        self.calling = '_ManageAtomicData'
        self._predefinedDataFileDict = _predefinedDataFileDict
        self.setDataFileDict()
        self.addDataFilePath()
        self.addDataFilePath('../atomic_data_fits/', inPyNeb=True)
        self.addDataFilePath('../atomic_data/', inPyNeb=True)
        self.addDataFilePath('./')
        self._RecombData = {}
        self._initChianti()
        self.read_gsconf()

    def includeFitsPath(self):
        self.addDataFilePath('../atomic_data_fits/old_fits/', inPyNeb=True)

    def removeFitsPath(self):
        self.removeDataFilePath('../atomic_data_fits/old_fits/')
        
    def getPredefinedDataFileDict(self, data_dict=None):
        """
        Retrieve all predefined dictionaries of atomic data (or a specified one if data_dict is specified)
        
        Parameters:
           - data_dict    name of the dictionary (default: None)

        """
#        Retrieve a predefined dictionary of atomic data (either the default one [data_dict unspecified]
#        or another one [identified by data_dict])
        if data_dict is None:
            return self._predefinedDataFileDict
        elif data_dict in self._predefinedDataFileDict:
            return self._predefinedDataFileDict[data_dict]
        else:
            self.log_.error('unknown DataFile dictionary {0}'.format(data_dict), calling=self.calling)
            return None


    def getAllPredefinedDict(self):
        """
        Return the labels of the available predefined dictionaries of atomic data

        """
        return pn.atomicData.getPredefinedDataFileDict().keys()

    
    def resetDataFileDict(self):
        """
        Reset the DataFileDict to the default value. This value is in defaultDict

        """
        self.setDataFileDict(self.defaultDict)


    def setDataFileDict(self, data_dict=None):
        """
        Use a dictionary for the atomic data (either predefined [string] or custom [dictionary]).
        
        Parameters:
           - data_dict    name of the dictionary 

        """
        if data_dict is None:
            self._DataFileDict = {}
            return None
        if type(data_dict) == type({}):
            for atom in data_dict:
                for data_type in data_dict[atom]:
                    self.setDataFile(data_dict[atom][data_type], atom, data_type)
            self.predefined = None
        elif data_dict in self.getPredefinedDataFileDict():
            self.setDataFileDict()
            self.setDataFileDict(self.getPredefinedDataFileDict(data_dict))
            self.predefined = data_dict
            self.log_.message('DataFileDict set to {0}.'.format(data_dict), calling=self.calling)
        else:
            self.log_.error('{0} is not a valid DataFileDict.'.format(data_dict), calling=self.calling)
 
            
    def addDataFilePath(self, fits_dir=None, inPyNeb=False):
        """
        Add a directory to the list of directories where atomic data files are searched for.
        
        Parameters:
           - fits_dir    directory
           - inPyNeb     Boolean.

        """
        if fits_dir is None:
            self._DataFilePaths = []
        else:
            try:
                if inPyNeb:
                    self._DataFilePaths.append(execution_path(fits_dir))
                else:
                    self._DataFilePaths.append(os.path.abspath(fits_dir))
            except:
                self.log_.warn('{0} could not be added to the path list', calling=self.calling)
 
 
    def removeDataFilePath(self, fits_dir=None):
        """
        Remove a directory from the list of directories where atomic data files are searched for.

        Parameters:
           - fits_dir    name of the directory to be removed

        """
        if fits_dir is None:
            pass
        else:
            try:
                self._DataFilePaths.remove(os.path.abspath(fits_dir))
            except:
                pass
            try:
                self._DataFilePaths.remove(execution_path(fits_dir))
            except:
                pass
 
                
    def getAllDataFilePaths(self):
        """
        Return the list of directories where atomic data files are searched for.

        """
        return self._DataFilePaths
 
            
    def getDirForFile(self, data_file):
        """
        Return the first directory from getDataFilePaths() where a file is found.
        If nothing is found, return None.

        Parameters:
           - data_file    name of the file

        """
        
        for dir in self._DataFilePaths:
            if data_file in os.listdir(dir):
                return dir
        if data_file.split('.')[-1] == 'chianti':
            strs = data_file.split('_')
            elem = strs[0].capitalize()
            spec = roman_to_int(strs[1])
            atom = elem + str(spec)
            strs = re.split(r'[_.]+',data_file)
            data_type = strs[2]
            if (data_type == 'atom') and (atom2chianti(atom) in self.ChiantiIONS['atom']):
                return('{}/{}/{}'.format(self.Chianti_path, elem.lower(), atom2chianti(atom)))
            if (data_type == 'coll') and (atom2chianti(atom) in self.ChiantiIONS['coll']):
                return('{}/{}/{}'.format(self.Chianti_path, elem.lower(), atom2chianti(atom)))
        else:
            return None
 
    
    def printDirForAllFiles(self):
        """
        Print the directories where all the files defined in the dictionary of 
            atomic data files are located.

        """
        for atom in self._DataFileDict:
            atom_file = self.getDataFile(atom, 'atom', warn=False)
            if atom_file is not None:
                print('{0} atom from {1} in {2}'.format(atom, atom_file , self.getDirForFile(atom_file)))
            coll_file = self.getDataFile(atom, 'coll', warn=False)
            if coll_file is not None:
                print('{0} coll from {1} in {2}'.format(atom, coll_file , self.getDirForFile(coll_file)))
            rec_file = self.getDataFile(atom, 'rec', warn=False)
            if rec_file is not None:
                print('{0} rec from {1} in {2}'.format(atom, rec_file , self.getDirForFile(rec_file)))
 
 
    def getAllAvailableFiles(self, atom=None, data_type=None):
        """
        Scan every directory in the list of paths, printing all the *atom*.fits, *atom*.dat,
            *coll*.fits and *rec*.fits files
            
        Parameters:
           - atom        atom name
           - data_type   either 'rec', 'atom', 'coll' or 'trc'

        """
        file_list = []
        if data_type is None:
            data_types = ['atom', 'coll', 'rec', 'trc']
        else:
            data_types = [data_type]
        if atom is not None:
            elem, spec = parseAtom(atom)
            elem = elem.lower()
            spec = int_to_roman(int(spec)).lower()
            atom_str = '{0}_{1}_'.format(elem, spec)
            if atom2chianti(atom) in self.ChiantiIONS['atom']:
                if 'atom' in data_types:
                    file_list.append('{0}_{1}_atom.chianti'.format(elem, spec))
            if atom2chianti(atom) in self.ChiantiIONS['coll']:
                if 'coll' in data_types:
                    file_list.append('{0}_{1}_coll.chianti'.format(elem, spec))
        else:
            atom_str = '.'
            for atom in self.ChiantiIONS['atom']:
                elem, spec = parseAtom(atom)
                elem = elem.lower()
                spec = int_to_roman(int(spec)).lower()
                if 'atom' in data_types:
                    file_list.append('{0}_{1}_atom.chianti'.format(elem, spec))
            for atom in self.ChiantiIONS['coll']:
                elem, spec = parseAtom(atom)
                elem = elem.lower()
                spec = int_to_roman(int(spec)).lower()
                if 'coll' in data_types:
                    file_list.append('{0}_{1}_coll.chianti'.format(elem, spec))
        for dir in self._DataFilePaths:
            files = os.listdir(dir)
            for file in files:
                if (('.fits' in file) or ('.func' in file) or ('.hdf5' in file) or ('.dat' in file) ) and (atom_str in file):
                    for dt in data_types:
                        if dt in file:
                            file_list.append(file)
                            
        file_list.sort()
        return file_list
    
    
    def getAllPossibleAtoms(self):
        """
        Return the list of all the possible collisional and recombination 
        atoms that can be built with the available datafiles (not all of 
        them necessarily included in the current atomic file dictionary) 
                    
        Parameters:
           - atom        atom name
           - data_type   either 'rec', 'atom' or 'coll'

        """
        atom_dict = {'atom':[], 'coll': [], 'rec': []}
        for path in self._DataFilePaths:
            files = os.listdir(path)
            for elem in ELEM_LIST:
                for spec in SPEC_LIST:
                    atom = elem + spec
                    atom_str = '{0}_{1}_'.format(elem.lower(), int_to_roman(int(spec)).lower())
                    for entry in files:
                        if (('.fits' in entry) or ('.dat' in entry) ) and (atom_str in entry):
                            if ('atom' in entry):
                                atom_dict['atom'].append(atom)
                            elif ('coll' in entry):
                                atom_dict['coll'].append(atom)
                            elif ('rec' in entry):
                                atom_dict['rec'].append(atom)
                            else:
                                pass
        coll_atom_list = list(set(atom_dict['atom']) & set(atom_dict['coll']))
        rec_atom_list = list(set(atom_dict['rec']))
        print('Data for the following collisional atoms exist:', coll_atom_list)
        print('Data for the following recombination atoms exist:', rec_atom_list)
        return [coll_atom_list, rec_atom_list] 
        
    
    def getAllAtoms(self, coll=True, rec=False):
        """
        Return a list of all the atoms (e.g. 'O3') included in the adopted dictionary 
        and for which dataFiles are available (a subset of all the ions returned by
        the getAllPossibleIons() command)

        Parameters:
            - coll     if True (default) includes the Atom objects
            - rec     if True (not the default) includes the RecAtoms

        """
        to_return = []
        for atom in self._DataFileDict:
            if coll and ('atom' in self._DataFileDict[atom]):
                to_return.append(atom)
            if rec and ('rec' in self._DataFileDict[atom]):
                to_return.append(atom)
        return to_return
            

    def getDataFile(self, atom=None, data_type=None, warn=True):
        """
        Return the name of the atomic data file associated to an atom and a type of data, 
        which is one of ('atom', 'coll', 'rec').
            
        Parameters:
            - atom        selected atom
            - data_type   type of atomic data ('atom', 'coll', 'rec', 'trc')
            - warn        warn if no associated data file is found
        
        """
        if atom is None:
            if data_type is None:
                return self._DataFileDict
            else:
                DFD = self._DataFileDict
                to_return = {}
                for at in DFD:
                    if data_type in DFD[at]:
                        to_return[at] = DFD[at][data_type]
                return to_return
        if data_type is None:
            return(self.getDataFile(atom, 'atom', warn=warn),
                   self.getDataFile(atom, 'coll', warn=warn),
                   self.getDataFile(atom, 'rec', warn=warn),
                   self.getDataFile(atom, 'trc', warn=warn))
        if data_type not in ('atom', 'coll', 'rec', 'trc'):
            if warn:
                self.log_.warn('{0} type unknown'.format(data_type), calling=self.calling)
            return None
        if atom not in self._DataFileDict:
            if warn:
                self.log_.warn('data for {0} not available'.format(atom), calling=self.calling)
            return None
        if data_type in self._DataFileDict[atom]:
            return self._DataFileDict[atom][data_type]
        else:
            if warn:
                self.log_.warn('{0} data not available for {1}'.format(data_type, atom), calling=self.calling)
            return None

        
    def getDataFullPath(self, atom, data_type=None, warn=True):
        """
        Return the full path of the atomic data file associated to an atom and a type of data, 
        which is one of ('atom', 'coll', 'rec').
        
        Parameters:
            - atom        selected atom
            - data_type   type of atomic data ('atom', 'coll', 'rec', 'trc')
            - warn        warn if no associated data file is found
        
        """
        data_file = self.getDataFile(atom, data_type=data_type, warn=warn)
        data_path = self.getDirForFile(data_file)
        if data_file is None:
            return None
        else:
            return '{0}/{1}'.format(data_path, data_file)


    def scanDirForDataFiles(self, fits_dir):
        """
        Scan a directory and associate any file named x_i_type_ref.fits to the atom Xi (i in roman), 
        type being one of ('atom', 'coll', 'rec').
        The directory must have been registered before using addDataFilePath(fits_dir)

        Parameters:
            - fits_dir        directory name

        """
        for file_ in os.listdir(fits_dir):
            names = multi_split(file_, ['_', '.'])
            elem = names[0].title()            
            if elem in ELEM_LIST:
                spec = str(roman_to_int(names[1]))
                if spec in SPEC_LIST:
                    atom = elem + spec
                    if atom not in self._DataFileDict:
                        self._DataFileDict[atom] = {}
                    if names[2] in ('atom', 'coll', 'rec', 'trc'):
                        self.setDataFile(atom, names[2], file_)
 
    def setDataFile(self, data_file=None, atom=None, data_type=None):
        """
        Associate an atomic data file to an atom and a type of data, which is one of 
            ('atom', 'coll', 'rec', 'func').

        Usage:
            pn.atomicData.setDataFile('cl_iii_atom_M83-KS86.fits')
            pn.atomicData.setDataFile('TEST.fits', 'O3', 'atom') # but the previous way 
                to name the files is preferred
            pn.atomicData.setDataFile('o_iii_atom.chianti')

        Parameters:
            - data_file    atomic data file. 
            - [atom]         selected atom
            - [data_type]    'atom', 'coll', 'rec'
            If atom and data_type not set, file format is assumed to be e.g. o_iii_coll_REF.fits

        """          
        if atom is None:
            strs = data_file.split('_')
            elem = strs[0].capitalize()
            spec = roman_to_int(strs[1])
            atom = elem + str(spec)
        else:
            elem, spec = parseAtom(atom)
        if data_type is None:
            strs = re.split(r'[_.]+',data_file)
            data_type = strs[2]
        if data_type not in ('atom', 'coll', 'rec', 'trc', 'func'):
            self.log_.error('{0} is not a valid type'.format(data_type))
            return None
        if elem not in ELEM_LIST:
            self.log_.error('{0} is not a valid element'.format(elem))
            return None
        if atom not in self._DataFileDict:
            self._DataFileDict[atom] = {}
        if self.getDirForFile(data_file) is not None:
            self._DataFileDict[atom][data_type] = data_file
            # as something changed, predefined is not valid anymore
            self.predefined = None
            if data_type == 'rec':
                if atom in self._RecombData:
                    del self._RecombData[atom]
            self.log_.message('Adding {0} {1} data for {2}'.format(data_file, data_type, atom), calling=self.calling)
        else:
            av_data = self.getAllAvailableFiles(atom, data_type)
            self.log_.error("""Unknown file {0} or corresponding File path not included in list. 
You may mean one of these files: {1}""".format(data_file, av_data),
                           calling=self.calling)
        if (atom == 'H1') and ('H1' in pn.atomicData._RecombData):
            del pn.atomicData._RecombData['H1']


    def printAllSources(self, at_set=None, predef=None):
        """
        Print bibliographic sources of the adopted data sets of a list of atoms.
        
        Usage:
        pn.atomicData.printAllSources(['O3', 'Ar4', 'S2'])
        pn.atomicData.printAllSources([O, S])
        pn.atomicData.printAllSources()

        Parameters: 
            - at_set   a list or tuple containing the atoms. If not set, 
                       print bibliographic sources for all the atoms in PyNeb's default dictionary
                       
        """        
        self.calling = 'printAllSources'
        if (type(at_set) == list) or (type(at_set) == tuple):
            at_dict = {}
            for item in at_set:
                atom = parseAtom(item)[0]
                spec = parseAtom(item)[1]
                if spec is '':
                    for ispec in SPEC_LIST:
                        try:
                            at_dict[atom+ispec] = pn.Atom(atom, ispec)
                        except:
                            pass
                else:
                    at_dict[item] = pn.Atom(atom, spec)
        elif at_set is not None:
            pn.log_.error('The argument at_set must be a list or a tuple', calling=self.calling)
            return None
        elif (at_set is None) and (predef is None):
            at_dict = pn.getAtomDict()
        elif predef in pn.atomicData.getAllPredefinedDict():
            current = pn.atomicData.predefined
            pn.atomicData.setDataFileDict(predef)
            at_dict = pn.getAtomDict()
            pn.atomicData.setDataFileDict(current)
        else: 
            pn.log_.error('The argument predef must be the label of a predefined dictionary', 
                          calling=self.calling)
        
        for item in sorted(at_dict):                                                              
            at_dict[item].printSources()
    
    def _initChianti(self):
        self.ChiantiIONS = {'atom':[], 'coll':[]}
        self.Chianti_path = None
        
        if pn.config.INSTALLED['Chianti']:
    
            self.Chianti_path = os.environ['XUVTOP']
            masterlist = '{0}/masterlist/masterlist.ions'.format(self.Chianti_path)
            try:
                with open(masterlist) as master:
                    for line in master.readlines():
                        atom = line.split()[0]
                        elem = atom.split('_')[0]
                        if pn.config.Chianti_version_main == '7':
                            coll_file = '{0}/{1}/{2}/{2}.splups'.format(self.Chianti_path, elem, atom)
                        elif pn.config.Chianti_version_main == '8':
                            coll_file = '{0}/{1}/{2}/{2}.scups'.format(self.Chianti_path, elem, atom)
                        elif pn.config.Chianti_version_main == '9':
                            coll_file = '{0}/{1}/{2}/{2}.scups'.format(self.Chianti_path, elem, atom)
                        else:
                            pn.log_.error('Unknown version of Chianti {}'.format(pn.config.Chianti_version),
                                          calling="_ManageAtomicData/_initChianti")
                        atom_file = '{0}/{1}/{2}/{2}.wgfa'.format(self.Chianti_path, elem, atom)
                        if os.path.exists(coll_file): 
                            self.ChiantiIONS['coll'].append(atom)
                        if os.path.exists(atom_file):
                            self.ChiantiIONS['atom'].append(atom)
            except:
                pn.log_.warn('File not found {}, no Chianti data available'.format(masterlist), 
                             calling='_initChianti')
    def addAllChianti(self):
        atoms = self.getAllAtoms()
        masterlist = '{0}/masterlist/masterlist.ions'.format(self.Chianti_path)
        try:
            with open(masterlist, 'r') as f:
                lines = f.readlines()
        except:
            pn.config.error('Masterlist file not read', calling='addAllChianti')
            return None
        chianti_ions = [l[0:7].strip() for l in lines]
        chianti_ions = [i.capitalize().replace('_','') for i in chianti_ions if i[-1] != 'd']
        for cion in chianti_ions:
            if cion not in atoms and cion not in ('H1', 'He1', 'He2'):
                try:
                    atfiles = self.getAllAvailableFiles(cion)
                    for atfile in atfiles:
                        self.setDataFile(atfile)
                except:
                    pass
    
    def read_gsconf(self):
        try:
            gsconf = np.genfromtxt(execution_path('../atomic_data/levels/gsconfs.dat'), 
                                   names=['atom', 'gsconf'], dtype='U5, U5')
            self.gsconf = {gs['atom']:gs['gsconf'] for gs in gsconf}
        except:
            self.gsconf = {}
        
    def printPoem(self, yr=0):
        """
        Print one of Wiese et al. spectroscopic poems
        
        Parameters:
            - yr     publication year. If no valid year is given, random print
 
        """
        if (yr == 96):
            print("""********************************************
            
We have used a near limitless data source,
It is the Opacity Project, of course.
Its results are certainly fine
For most any LS-coupled line.


(Wiese, Fuhr & Deters 1996)""")            
        elif (yr == 66):
            print("""********************************************

If there is no other data source            
Use the Coulomb approximation, of course.
The results should be certainly fine
For any moderately or highly excited line.

(Wiese, Smith & Glennon 1966)""")
        else:
            rnd = np.random.rand()
            if (rnd < 0.5):
                self.printPoem(yr=66)
            else:
                self.printPoem(yr=96)


def extract_flt(str_):
    """
    Input is a bytes (py3) or a string (py2)
    Return the first float in str_, removing the no digit (and no dot) leading part of it.
    Ex:
    extract_flt('(123.00?') -> 123.00
    """
    res = ''
    if len(str_) > 0:
        if str_.decode()[0] in ('(', '['):
            str_ = str_[1:]
    for l in str_.decode():
        if l.isdigit() or l == '.':
            res += l
        else:
            break
    if res == '':
        return np.nan
    else:
        return float(res)

def readNIST(NISTfile,NLevels=None):
    """
    The NIST levels file must be obtained from this page:
        http://physics.nist.gov/PhysRefData/ASD/levels_form.html
    with the following options:
        Level Units: cm-1
        Format output: ASCII
        Display output: in its entirely
        Energy ordered
        Level Information: Principal configuration, Principal term, Level, J
        Bibliographic references
    """
    data = np.genfromtxt(NISTfile, comments = '-', 
                    delimiter = '|', names = 'conf, term, J, energy, ref', 
                    dtype=('U23', 'U9', 'U4', 'float', 'U10'), autostrip=True,
                    converters = {'energy':extract_flt})
    mask = data['J'] != ''
    data = data[mask]
    previous_ref = ''
    previous_conf = ''
    previous_term = ''
    for d in data:
        if d['conf'] == '':
            d['conf'] = previous_conf
        else:
            previous_conf = d['conf']
        if d['term'] == '':
            d['term'] = previous_term
        else:
            previous_term = d['term']
        if d['ref'] == '':
            d['ref'] = previous_ref
        else:
            previous_ref = d['ref']
            
        if '?' in d['J']:
            d['J'] = d['J'].split('?')[0]
        if d['J'] == '':
            d['J'] = 0.
        elif '/' in d['J']:
            d['J'] = '{0:.1f}'.format(float(d['J'].split('/')[0]) / float(d['J'].split('/')[1]))
        elif ',' in d['J']:
            d['J'] = '{0:.1f}'.format(float(d['J'].split(',')[0]))
        else:
            d['J'] = '{0:.1f}'.format(float(d['J']))
    data = data[data['energy'].argsort()]
    data = data.astype([('conf', 'U23'), ('term', 'U9'), ('J', 'float'), ('energy', 'float'), ('ref', 'U10')])
    if NLevels is not None:
        data = data[0:NLevels]
    return data


def getLevelsNIST(atom, NLevels = None):
    """
    Return a numpy.array containing the NIST data related to the levels for the given atom
    The keys of the recarry are: conf, term, J, energy, ref
    Example:
        data = getLevelsNIST('O3')
        energies = data['energy']
        
    """
    elem, spec = parseAtom(atom)
    elem = elem.lower()
    spec = int_to_roman(int(spec)).lower()
    atom_str = '{0}_{1}'.format(elem, spec)
    level_file = '{0}_levels.dat'.format(atom_str)
    file_ = execution_path('../atomic_data/levels/{0}'.format(level_file))
    if os.path.isfile(file_):
        data = readNIST(file_, NLevels=NLevels)
        pn.log_.message('Reading energies and stat weights from {0}'.format(level_file), calling='getLevelsNIST')
        return data
    else:
        return None
    
def getNIST(elem, ion, fileout=None):
    """
    This will download the enegy levels from NIST and write them to a file
    example: getNIST('O',3) write the o_iii_levels.dat file
    """
    
    try:
        from selenium import webdriver
        from selenium.webdriver.support.ui import Select
    except:
        pn.log_.warn('selenium not installed', calling='get_NIST')
        return None
    import time
    atom_NIST = elem + ' ' + int_to_roman(int(ion))
    browser = webdriver.Firefox()
    try:
        browser.get('http://physics.nist.gov/PhysRefData/ASD/levels_form.html')
        Select(browser.find_element_by_name('format')).select_by_value('1')
        browser.find_element_by_css_selector("input[name='multiplet_ordered'][value='1']").click()
        browser.find_element_by_name('spectrum').send_keys(atom_NIST)
        browser.find_element_by_name('lande_out').click()
        browser.find_element_by_name('perc_out').click()
        browser.find_element_by_name('submit').click()
        time.sleep(1)
        tab = browser.find_element_by_tag_name('pre').text
    except:
        tab = None
    browser.close()
    
    if tab is not None:
        if fileout is None:
            fileout = '{}_{}_levels.dat'.format(elem.lower(), int_to_roman(int(ion)).lower())
        with open(fileout, 'w') as f:
            res = tab.split('\n')
            print_it = True
            i_line = 0
            for line in res[4:]:
                if '---' in line or '[' in line:
                    print_it = False
                if print_it:
                    f.write(line+'\n')
                    i_line += 1
        return i_line
    else:
        return None
            
    
def _atom_fits2ascii(filename):
    """
    Transform a *atom*.fits file into an *atom*.dat file.
    
    """
    strs = filename.split('/')[-1].split('_')
    elem = strs[0].capitalize()
    spec = roman_to_int(strs[1])
#    type = strs[2]
    pn.atomicData.setDataFile(filename.split('/')[-1])
    atom = pn.Atom(elem, spec)
    atom.printSources()
    print('')
    fileout = raw_input('Enter the name of the output file (between {0[0]}_{0[1]}_{0[2]}_ and .dat:'.format(strs))
    
    fout = open('{0[0]}_{0[1]}_{0[2]}_{1}.dat'.format(strs,fileout), 'w')
    fout.write('Aij\n')
#    str_ = ''
#    for i in np.arange(atom.NLevels):
#        str_ += '1/s '
#    str_ += '\n'
    str_ = '1/s ' * atom.NLevels + '\n'
    fout.write(str_)
    for a1 in atom.getA():
        str_ = ''
        for a2 in a1:
            str_ += '{0:10.7e}'.format(a2)
            str_ += ' '
        str_ += '\n'
        fout.write(str_)
    for item in atom.AtomData.AtomHeader.iteritems():
        if 'SPECTRUM' in item[0] or 'ATOM' in item[0] or 'GSCONFIG' in item[0]:
            fout.write('*** ' + item[0] + ' ' + str(item[1]) + '\n')
        if 'NOTE' in item[0]:
            if 'Energy' not in item[1]:
                N_Note = item[0].split('NOTE')[1]
                str_N_Note = str(N_Note.split('_')[0])
                fout.write('*** SOURCE' + str_N_Note + ' ' + atom.AtomData.AtomHeader.get('SOURCE'+str_N_Note)  + '\n')
                fout.write('*** ' + item[0] + ' ' + str(item[1]) + '\n')
    
    fout.close()
    
def _coll_fits2ascii(filename, overwrite=None):
    """
    Transform a *coll*.fits file into an *coll*.dat file.
    
    """
    strs = filename.split('/')[-1].split('_')
    elem = strs[0].capitalize()
    spec = roman_to_int(strs[1])
    pn.atomicData.setDataFile(filename.split('/')[-1])
    atom = pn.Atom(elem, spec)
    #atom.printSources()
    fileout = '{}.dat'.format(filename.split('.')[0])
    
    if overwrite is not True:
        if os.path.exists(fileout):
            erase = raw_input('{} exists, overwrite?'.format(fileout))
            if erase != 'y':
                return
    fout = open(fileout, 'w')
    fout.write('*** {0}{1} collision strengths data\n'.format(elem, strs[1].upper()))
    str_ = '0 0'
    for t in atom.getTemArray():
        str_ += ' {0:9.3e}'.format(t)
    str_ += '\n'
    fout.write(str_)
    for omega in atom.getOmegaArray():
        j = strExtract(omega, 'Omega(', '->')
        i = strExtract(omega, '->', ')')
        omegs = strExtract(omega, '[', ']')
        str_ = '{0} {1}'.format(i,j)
        for omeg in omegs.split():
            str_ += ' {0:9.3e}'.format(float(omeg))
        str_ += '\n'
        fout.write(str_)
        
    for item in atom.CollData.CollHeader.iteritems():
        if 'SPECTRUM' in item[0] or 'ATOM' in item[0] or 'GSCONFIG' in item[0]:
            fout.write('*** ' + item[0] + ' ' + str(item[1]) + '\n')
        if 'TUNIT1' == item[0]:
            fout.write('*** T_UNIT ' + str(item[1]) + '\n') 
        if 'SOURCE' in item[0]:
            N_Source = item[0].split('SOURCE')[1]
            str_N_Source = str(N_Source.split('_')[0])
            fout.write('*** ' + item[0] + ' ' + str(item[1]) + '\n')
            try:
                fout.write('*** NOTE' + str_N_Source + ' ' + atom.CollData.CollHeader.get('NOTE'+str_N_Source)  + '\n')
            except:
                fout.write('*** NOTE' + str_N_Source + ' collision strengths \n')
    print('{} done'.format(fileout))
    fout.close()
    
def print_stout_coll(Atom, file_):
    
    if type(Atom) is not pn.Atom:
        pn.log_.error('Atom must be a PyNeb.Atom object', calling='print_stout')
    
    temps = Atom.getTemArray(keep_unit=False)
    n_temp = len(temps)
    omegas = Atom.getOmegaArray()
    n_levels = omegas.shape[0]
    with open(file_, 'w') as f:
        f.write('11 10 14\n')
        f.write('TEMP {} \n'.format(' '.join(['{:.0f}'.format(t) for t in temps])))
        for i in range(n_levels):
            for j in range(n_levels):
                if omegas[i,j].sum() > 0.0:
                    f.write('CS ELECTRON {} {} {} \n'.format(i+1, j+1, ' '.join([str(o) for o in omegas[i,j]])))
        f.write('***************\n')
        f.write('Refs:\n')
        for s in Atom.CollData.getSources():
            f.write('{} \n'.format(s))


