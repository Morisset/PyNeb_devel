import os
from .logging import my_logging

class _Config(object):
    """
    This is the place where to put stuff that any module may need, a kind of COMMON.
    An instantiation is done in the main __init__ file using the "config" name.

    The CHIANTI and Stout database directories are read at import time from the
    XUVTOP and STOUT_DIR environment variables, and can be changed at runtime with
    set_chianti_path() and set_stout_path() (which also update the corresponding
    environment variable for this process, as the low-level readers rely on it).

    """


    def __init__(self):
        """
        Define variables that will be known from everywhere:
        - INSTALLED: a dictionary whose keys are the libraries that PyNeb may need
            e.g. 'plt' for matplotlib.pyplot, 'scipy' for scipy.interpolate. and whose values
            are Boolean
        
        - Datafile: a dictionary holding the HI atom, which can be used intensively in Atom.getIonAbundance
        
        - pypic_path: This is only set when needed (basically by getEmisGridDict).
            The default value for the pypic_path parameter of getEmisGridDict().
            The first try is ./.pypics: if it cannot be created, or it exists but it is not writable, 
            /tmp/pypics is tried; if it cannot be created, or it exists but it is not writable, the path
            is set to None and a pypic_path value shall be provided to getEmisGridDict().
        
        Parameters: none
        
        """
        self.log_ = my_logging()
        self.calling = '_Config'
                
        self.INSTALLED = {}
        try:
            import matplotlib.pyplot as plt
            self.INSTALLED['plt'] = True
        except:
            self.INSTALLED['plt'] = False
            self.log_.message('matplotlib not available', calling=self.calling)
        try:     
            from scipy import interpolate
            self.INSTALLED['scipy'] = True
        except:
            self.INSTALLED['scipy'] = False
            self.log_.message('scipy not available or interpolate not in scipy', calling=self.calling)
        try:
            import multiprocessing as mp
            self.INSTALLED['mp'] = True
            self.Nprocs = mp.cpu_count()
        except:
            self.INSTALLED['mp'] = False
            self.log_.message('multiprocessing not available', calling=self.calling)
            self.Nprocs = 1
        if 'XUVTOP' in os.environ:
            self.INSTALLED['Chianti'] = True
            self.Chianti_path = os.environ['XUVTOP']
            self.Chianti_version, self.Chianti_version_main = self._read_chianti_version(self.Chianti_path)
        else:
            self.INSTALLED['Chianti'] = False
            self.Chianti_path = None
            self.Chianti_version = None
            self.Chianti_version_main = None
        if 'STOUT_DIR' in os.environ:
            self.INSTALLED['Stout'] = True
            self.Stout_dir = os.environ['STOUT_DIR']
        else:
            self.INSTALLED['Stout'] = False
            self.Stout_dir = None
        try:
            import astropy.io.fits as pyfits
            self.INSTALLED['pyfits from astropy'] = True
            self.INSTALLED['pyfits'] = True
        except:
            self.INSTALLED['pyfits from astropy'] = False 
            try:
                import pyfits
                self.INSTALLED['pyfits'] = True
            except:
                self.INSTALLED['pyfits'] = False
        try:
            import h5py
            self.INSTALLED['h5py'] = True
        except:
            self.INSTALLED['h5py'] = False
        try:
            from astropy.table import Table
            self.INSTALLED['astropy Table'] = True
        except:
            self.INSTALLED['astropy Table'] = False
        try:
            import cvxopt
            self.INSTALLED['cvxopt'] = True
        except:
            self.INSTALLED['cvxopt'] = False   

        # self.INSTALLED['ai4neb'] = False

        self.DataFiles = {}
            
        self.unuse_multiprocs()
        
        self.kappa = None
        self.set_noExtrapol(False)
                    
        self.__pypic_path = None
        
        self.vactoair_low_wl = 2000. # UV in vacuum
        self.vactoair_high_wl = 1e30 # no upper limit, IR in air!!!
                        
    # def import_AI4Neb(self):
    #     try:
    #         from ai4neb import manage_RM
    #         self.INSTALLED['ai4neb'] = True
    #     except:
    #         self.INSTALLED['ai4neb'] = False

    def set_noExtrapol(self, value):
        self._noExtrapol = bool(value)

    def get_noExtrapol(self):
        return self._noExtrapol

    @staticmethod
    def _read_chianti_version(chianti_path):
        """
        Return (version, version_main) read from {chianti_path}/VERSION, or (None, None).
        """
        try:
            with open(os.path.join(chianti_path, 'VERSION')) as vFile:
                version = vFile.readline().strip()
            return version, version.split('.')[0]
        except Exception:
            return None, None

    def set_chianti_path(self, chianti_path=None):
        """
        Set (or unset, with None) the CHIANTI database directory at runtime,
        as an alternative to the XUVTOP environment variable.
        Side effect: sets/removes the XUVTOP environment variable for this process,
        since the low-level CHIANTI readers use it. The CHIANTI version is read
        from {chianti_path}/VERSION and pn.atomicData is refreshed.
        Existing Atom objects keep their already-loaded data; build new ones to
        use the new path.

        Parameters:
            chianti_path: the CHIANTI database directory, or None to unset

        """
        if chianti_path is None:
            os.environ.pop('XUVTOP', None)
            self.INSTALLED['Chianti'] = False
            self.Chianti_path = None
            self.Chianti_version = None
            self.Chianti_version_main = None
        else:
            chianti_path = os.path.abspath(os.path.expanduser(chianti_path))
            if not os.path.isdir(chianti_path):
                self.log_.error('Directory {0} does not exist'.format(chianti_path),
                                calling=self.calling)
            version, version_main = self._read_chianti_version(chianti_path)
            if version is None:
                self.log_.error('Cannot read {0}/VERSION; not a valid CHIANTI directory'.format(chianti_path),
                                calling=self.calling)
            os.environ['XUVTOP'] = chianti_path
            self.INSTALLED['Chianti'] = True
            self.Chianti_path = chianti_path
            self.Chianti_version = version
            self.Chianti_version_main = version_main
        self._refresh_chianti()

    def _refresh_chianti(self):
        # lazy imports: Config.py must not import pn_chianti/manage_atomic_data at
        # module top (circular); also tolerate being called before pyneb finished importing
        try:
            import pyneb as pn
            from . import pn_chianti
            pn_chianti._refresh_chianti_tools()
            pn.atomicData._initChianti()
        except (ImportError, AttributeError):
            pass

    def get_chianti_path(self):
        return self.Chianti_path

    def set_stout_path(self, stout_path=None):
        """
        Set (or unset, with None) the Stout database directory at runtime,
        as an alternative to the STOUT_DIR environment variable.
        Side effect: sets/removes the STOUT_DIR environment variable for this
        process. pn.atomicData is refreshed.
        Existing Atom objects keep their already-loaded data; build new ones to
        use the new path.

        Parameters:
            stout_path: the Stout database directory, or None to unset

        """
        if stout_path is None:
            os.environ.pop('STOUT_DIR', None)
            self.INSTALLED['Stout'] = False
            self.Stout_dir = None
        else:
            stout_path = os.path.abspath(os.path.expanduser(stout_path))
            if not os.path.isdir(stout_path):
                self.log_.error('Directory {0} does not exist'.format(stout_path),
                                calling=self.calling)
            os.environ['STOUT_DIR'] = stout_path
            self.INSTALLED['Stout'] = True
            self.Stout_dir = stout_path
        try:
            import pyneb as pn
            pn.atomicData._initStout()
        except (ImportError, AttributeError):
            pass

    def get_stout_path(self):
        return self.Stout_dir

    def _get_PypicPath(self):
        pypic_path = self.__pypic_path
        if pypic_path is None:
            pypic_path = './pypics/'
            if os.path.exists(pypic_path):
                if not os.access(pypic_path, os.W_OK):
                    self.log_.warn('Directory {0} not writable'.format(pypic_path),
                                      calling=self.calling)
                    pypic_path = None
            else:
                try:
                    os.makedirs(pypic_path)
                    self.log_.message('Directory {0} created'.format(pypic_path),
                                      calling=self.calling)
                except:
                    pypic_path = None
            if pypic_path is None:
                pypic_path = '/tmp/pypics/'
                if os.path.exists(pypic_path):
                    if not os.access(pypic_path, os.W_OK):
                        self.log_.warn('Directory {0} not writable'.format(pypic_path),
                                          calling=self.calling)                                   
                        pypic_path = None 
                else:
                    try:
                        os.makedirs(pypic_path)
                        self.log_.message('Directory {0} created'.format(pypic_path),
                                          calling=self.calling)
                    except:
                        pypic_path = None
        else:
            if os.path.exists(pypic_path):
                if not os.access(pypic_path, os.W_OK):
                    self.log_.warn('Directory {0} not writable'.format(pypic_path),
                                      calling=self.calling)
                    pypic_path = None
            else:
                try:
                    os.makedirs(pypic_path)
                    self.log_.message('Directory {0} created'.format(pypic_path),
                                      calling=self.calling)
                except:
                    self.log_.warn('Unable to create directory {0}'.format(pypic_path),
                                      calling=self.calling)
                    
                    pypic_path = None
            
        self.log_.message('pypic_path set to {0}'.format(pypic_path),
                                          calling=self.calling)
        self.__pypic_path = pypic_path
        
        return self.__pypic_path
    
    def _set_PypicPath(self, value):            
        self.__pypic_path = value
        test = self.pypic_path
        
    pypic_path = property(_get_PypicPath, _set_PypicPath, None, None)
            
    def use_multiprocs(self):
        self._use_mp = True
    
    def unuse_multiprocs(self):
        self._use_mp = False
        
            