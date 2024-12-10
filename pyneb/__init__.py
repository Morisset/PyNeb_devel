"""
PyNeb - python package for the analysis of emission lines
"""
"""
PyNeb - python package for the analysis of emission lines
"""
import sys

#from .version import __version__
from importlib.metadata import version
__version__ = version("PyNeb")

from .utils.Config import _Config
config = _Config()
log_ = config.log_

from .utils.manage_atomic_data import _ManageAtomicData
atomicData = _ManageAtomicData()
atomicData.defaultDict = 'PYNEB_23_01'
atomicData.resetDataFileDict()
from .core.pynebcore import Atom, RecAtom, getAtomDict, getRecEmissivity
from .core.pynebcore import EmissionLine, Observation, parseLineLabel, isValid
from .core.icf import ICF 
from .core.diags import Diagnostics, diags_dict
from .core.emisGrid import EmisGrid, getEmisGridDict
from .core.continuum import Continuum
from .plot.plotAtomicData import DataPlot
from .utils.physics import CST, print_IPs
from .utils.saverestore import save, restore
from .utils.init import LINE_LABEL_LIST, BLEND_LIST, label2levelDict
from .utils.manage_atomic_data import getLevelsNIST
from .utils.misc import ROOT_DIR
from .extinction.red_corr import RedCorr

__all__ = ['Atom', 'RecAtom', 'getAtomDict', 'getRecEmissivity', 'EmissionLine', 
           'Observation', 'parseLineLabel', 'isValid',
           'ICF', 'Diagnostics', 'diags_dict',
           'EmisGrid', 'getEmisGridDict',
           'Continuum',
           'DataPlot',
           'CST', 'print_IPs',
           'save', 'restore',
           'LINE_LABEL_LIST', 'BLEND_LIST', 'label2levelDict',
           'getLevelsNIST',
           'ROOT_DIR',
           'RedCorr']

__pyversion__ = sys.version_info[0]

log_.message('Starting PyNeb version %s' % __version__, calling='pyneb')

if sys.version_info[:2] < (2, 6):
    log_.warn('Python version >= 2.6 needed, seems you have {0}'.format(sys.version_info), calling='PyNeb')

log_.message('PyNeb ready.', calling='PyNeb')
log_.timer('Starting PyNeb version %s' % __version__, quiet=True, calling='PyNeb')

__info__ = """PyNeb version: {}.
PyNeb is cited using: Luridiana, V.; Morisset, C.; Shaw, R. A. 2015, A&A, 573, 42
bibcode: 2015A&A...573A..42L
Default data dictionary: {}. This set of data can be changed by the user. 
Do not forget to cite the papers corresponding to the atomic data you are using.
PyNeb website: https://pypi.python.org/pypi/PyNeb
PyNeb discussion and hotline group email: pyneb@googlegroups.com
""".format(__version__, atomicData.defaultDict)
