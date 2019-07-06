import os
from .misc import execution_path
from .logging import my_logging

class _Config(object):
    """
    This is the place where to put stuff that any module may need, a kind of COMMON.
    An instantiation is done in the main __init__ file using the "config" name.

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
            try:
                xuvtop = os.environ['XUVTOP']
                vFileName = os.path.join(xuvtop, 'VERSION')
                vFile = open(vFileName)
                versionStr = vFile.readline()
                vFile.close()
                self.Chianti_version = versionStr.strip()
                self.Chianti_version_main = self.Chianti_version.split('.')[0]
            except:
                self.Chianti_version = None
                self.Chianti_version_main = None
        else:
            self.INSTALLED['Chianti'] = False
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
        
        self.DataFiles = {}
            
        self.unuse_multiprocs()
        
        self.kappa = None
        self.set_noExtrapol(False)
                    
        self.__pypic_path = None
        
        self.vactoair_low_wl = 2000. # UV in vacuum
        self.vactoair_high_wl = 1e30 # no upper limit, IR in air!!!
                
    def set_noExtrapol(self, value):
        self._noExtrapol = bool(value)
        
    def get_noExtrapol(self):
        return self._noExtrapol
        
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
        
        
            