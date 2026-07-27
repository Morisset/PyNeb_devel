"""
Test reading atomic data in Stout format, using the small Stout tree bundled
in tests/STOUT (O III and Ne III, truncated to 5 levels), and the runtime
pn.config.set_stout_path() / set_chianti_path() API.
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

    def test_Ne3_stout_temperature_blocks(self):
        """
        A Stout coll file can hold several TEMP blocks, each with its own
        temperature grid and its own subset of the transitions. The bundled
        ne_3 fixture has two: levels 1-3 on a 27-point grid, levels 1-5 on a
        10-point one.
        """
        old_path = pn.config.get_stout_path()
        pn.config.set_stout_path(STOUT_DIR)
        try:
            pn.atomicData.setDataFile('ne_iii_atom.stout')
            pn.atomicData.setDataFile('ne_iii_coll.stout')
            Ne3 = pn.Atom('Ne', 3)
            assert Ne3.is_valid is True
            assert Ne3.CollData.N_temps == [27, 10]
            # NLevels is the highest upper level over *all* the blocks: the
            # first one alone only reaches level 3
            assert Ne3.collNLevels == 5
            assert Ne3.NLevels == 5
            # each transition is tabulated on the grid of the block it
            # belongs to, 1->3 in the first one and 4->5 in the second
            assert len(Ne3.getTemArray(lev_i=3, lev_j=1)) == 27
            assert len(Ne3.getTemArray(lev_i=5, lev_j=4)) == 10
            assert np.isclose(Ne3.getOmegaArray(3, 1)[0], 1.57e-01)
            assert np.isclose(Ne3.getOmegaArray(5, 4)[0], 2.09e-01)
            assert Ne3.getEmissivity(1e4, 1e2, wave=3869) > 0.
        finally:
            pn.atomicData.resetDataFileDict()
            pn.config.set_stout_path(old_path)

    def test_Ne3_stout_temperature_blocks_NLevels(self):
        """
        Asking for fewer levels than the file holds must truncate the atom,
        not raise: with several TEMP blocks the per-block level lists are
        ragged, so the number of levels cannot be taken over them at once.
        """
        old_path = pn.config.get_stout_path()
        pn.config.set_stout_path(STOUT_DIR)
        try:
            pn.atomicData.setDataFile('ne_iii_atom.stout')
            pn.atomicData.setDataFile('ne_iii_coll.stout')
            Ne3 = pn.Atom('Ne', 3, NLevels=4)
            assert Ne3.NLevels == 4
            assert Ne3.collNLevels == 4
            assert Ne3.getPopulations(1e4, 1e2).size == 4
            assert Ne3.getEmissivity(1e4, 1e2, 4, 1) > 0.
            # more levels than available: capped at what the file holds
            assert pn.Atom('Ne', 3, NLevels=8).NLevels == 5
        finally:
            pn.atomicData.resetDataFileDict()
            pn.config.set_stout_path(old_path)

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
