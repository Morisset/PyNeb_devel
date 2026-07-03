"""
Test reading atomic data in CHIANTI format, using the small CHIANTI tree
bundled in tests/CHIANTI (O III only, CHIANTI version 10.0.1).

The database location is set at runtime with pn.config.set_chianti_path(),
so the test works both on a machine with a real CHIANTI installation
(XUVTOP set) and on CI without one; the previous state is restored at the end.
"""
import os

import pyneb as pn

CHIANTI_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CHIANTI')


class TestChianti:

    def test_O3(self):
        old_path = pn.config.get_chianti_path()
        pn.config.set_chianti_path(CHIANTI_DIR)
        try:
            assert pn.config.INSTALLED['Chianti'] is True
            assert pn.config.Chianti_version == '10.0.1'
            assert 'o_3' in pn.atomicData.ChiantiIONS['atom']
            assert 'o_3' in pn.atomicData.ChiantiIONS['coll']
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
            pn.config.set_chianti_path(old_path)
        assert pn.config.INSTALLED['Chianti'] == (old_path is not None)
