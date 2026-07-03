"""
Test reading atomic data in Stout format, using the small Stout tree bundled
in tests/STOUT (O III only, truncated to 5 levels from Cloudy c23.01), and
the runtime pn.config.set_stout_path() / set_chianti_path() API.
"""
import os

import numpy as np
import pyneb as pn

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
STOUT_DIR = os.path.join(TESTS_DIR, 'STOUT')
CHIANTI_DIR = os.path.join(TESTS_DIR, 'CHIANTI')


class TestStout:

    def test_O3_stout(self):
        old_path = pn.config.get_stout_path()
        pn.config.set_stout_path(STOUT_DIR)
        try:
            assert pn.config.INSTALLED['Stout'] is True
            assert 'o_3' in pn.atomicData.StoutIONS['atom']
            assert 'o_3' in pn.atomicData.StoutIONS['coll']
            pn.atomicData.setDataFile('o_iii_atom.stout')
            pn.atomicData.setDataFile('o_iii_coll.stout')
            O3 = pn.Atom('O', 3)
            assert O3.is_valid is True
            # the fixture has transitions up to level 5, but the maximum
            # *lower* level is 4: NLevels must come from the upper levels
            assert O3.collNLevels == 5
            assert O3.NLevels == 5
            # tabulated collision strength round-trip: Omega(2->1) at the
            # third tabulated temperature (9e3 K) is 5.20e-01 in the file
            assert np.isclose(O3.getOmegaArray(2, 1)[2], 5.20e-01)
            # A(5->4) is 8.54e-01 in the file
            assert np.isclose(O3.getA(5, 4), 8.54e-01)
            assert O3.getEmissivity(1e4, 1e2, wave=5007) > 0.
        finally:
            pn.atomicData.resetDataFileDict()
            pn.config.set_stout_path(old_path)
        assert pn.config.INSTALLED['Stout'] == (old_path is not None)

    def test_unset_reverts_cleanly(self):
        old_stout = pn.config.get_stout_path()
        old_chianti = pn.config.get_chianti_path()
        try:
            pn.config.set_stout_path(STOUT_DIR)
            pn.config.set_stout_path(None)
            assert pn.config.INSTALLED['Stout'] is False
            assert pn.config.get_stout_path() is None
            assert 'STOUT_DIR' not in os.environ
            assert pn.atomicData.StoutIONS == {'atom': [], 'coll': []}
            assert pn.atomicData.Stout_path is None

            pn.config.set_chianti_path(CHIANTI_DIR)
            pn.config.set_chianti_path(None)
            assert pn.config.INSTALLED['Chianti'] is False
            assert pn.config.get_chianti_path() is None
            assert pn.config.Chianti_version is None
            assert 'XUVTOP' not in os.environ
            assert pn.atomicData.ChiantiIONS == {'atom': [], 'coll': []}
        finally:
            pn.config.set_stout_path(old_stout)
            pn.config.set_chianti_path(old_chianti)
