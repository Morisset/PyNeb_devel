#%%

import unittest 
import pyneb as pn
import numpy as np
import os
#%%

class TestAtom(unittest.TestCase):
    def test_O3(self):
        os.environ['XUVTOP'] = './tests/CHIANTI'
        print(os.getcwd())
        pn.atomicData.setDataFile('o_iii_atom.chianti')
        pn.atomicData.setDataFile('o_iii_coll.chianti')

        O3 = pn.Atom('O', 3, NLevels=19)