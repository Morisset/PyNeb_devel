"""
Test reading atomic data in CHIANTI format, using the small CHIANTI tree
bundled in tests/CHIANTI (O III only, CHIANTI version 10.0.1).

pyneb configures its CHIANTI support at import time from the XUVTOP
environment variable, so this test cannot rely on it having been set:
it points XUVTOP at the bundled tree, redoes the import-time
initializations, and restores everything afterwards. It therefore works
both on a machine with a real CHIANTI installation and on CI without one.
"""
import os

with open('.env', 'w') as f:
    f.write('XUVTOP={}\n'.format(os.path.join(os.path.dirname(__file__), 'CHIANTI')))

from dotenv import load_dotenv

load_dotenv()

import pyneb as pn
import pyneb.utils.pn_chianti as pn_chianti

CHIANTI_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CHIANTI')


class TestChianti:

    def test_O3(self):
        old_xuvtop = os.environ.get('XUVTOP')
        old_installed = pn.config.INSTALLED['Chianti']
        old_version = pn.config.Chianti_version
        old_version_main = pn.config.Chianti_version_main

        os.environ['XUVTOP'] = CHIANTI_DIR
        pn.config.INSTALLED['Chianti'] = True
        with open(os.path.join(CHIANTI_DIR, 'VERSION')) as f:
            pn.config.Chianti_version = f.readline().strip()
        pn.config.Chianti_version_main = pn.config.Chianti_version.split('.')[0]
        # if pyneb was imported without XUVTOP set, the version-dependent
        # _chianti_tools import in pn_chianti was skipped: bind it now
        if not hasattr(pn_chianti, '_chianti_tools'):
            from pyneb.utils import _chianti_tools_9 as _chianti_tools
            pn_chianti._chianti_tools = _chianti_tools

        try:
            pn.atomicData._initChianti()
            pn.atomicData.setDataFile('o_iii_atom.chianti')
            pn.atomicData.setDataFile('o_iii_coll.chianti')
            O3_chianti = pn.Atom('O', 3, NLevels=19)
            assert O3_chianti.atom == 'O3'
            assert O3_chianti.elem == 'O'
            assert O3_chianti.spec == 3
            assert O3_chianti.name == 'oxygen'
            assert O3_chianti.Z == 8
            assert O3_chianti.is_valid is True
            assert O3_chianti.NLevels == 19
            assert (os.path.normpath(O3_chianti.atomFitsPath) ==
                    os.path.normpath(os.path.join(CHIANTI_DIR, 'o', 'o_3')))
        finally:
            pn.atomicData.resetDataFileDict()
            if old_xuvtop is None:
                del os.environ['XUVTOP']
            else:
                os.environ['XUVTOP'] = old_xuvtop
            pn.config.INSTALLED['Chianti'] = old_installed
            pn.config.Chianti_version = old_version
            pn.config.Chianti_version_main = old_version_main
            pn.atomicData._initChianti()
