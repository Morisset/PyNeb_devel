import numpy as np
import pytest
import pyneb as pn
import os
#%%

class TestChianti:
    def test_O3(self):
        os.environ['XUVTOP'] = './tests/CHIANTI'
        #print(os.getcwd())
        pn.atomicData.setDataFile('o_iii_atom.chianti')
        pn.atomicData.setDataFile('o_iii_coll.chianti')
        O3_chianti = pn.Atom('O', 3, NLevels=19)
        pn.atomicData.resetDataFileDict()        
        assert O3_chianti.atom == 'O3'
        assert O3_chianti.elem == 'O'
        assert O3_chianti.spec == 3
        assert O3_chianti.name == 'oxygen'
        assert O3_chianti.Z == 8
        assert O3_chianti.is_valid is True
        assert O3_chianti.NLevels == 19
        assert O3_chianti.atomFitsPath == './tests/CHIANTI/o/o_3/'