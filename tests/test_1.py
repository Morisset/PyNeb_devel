import unittest 
import pyneb as pn
import numpy as np

class TestAtom(unittest.TestCase):
    def test_Atom_O3(self):
        O3 = pn.Atom('O',3)
        with self.subTest():
            self.assertEqual(O3.atom, 'O3')
        with self.subTest():
            self.assertEqual(O3.elem, 'O')
        with self.subTest():
            self.assertEqual(O3.spec, 3)
        with self.subTest():
            self.assertEqual(O3.name, 'oxygen')
        with self.subTest():
            self.assertEqual(O3.Z, 8)
        with self.subTest():
            self.assertIsInstance(O3.NLevels, np.int64)
        with self.subTest():
            self.assertEqual(O3.energy_eV[0,0], np.inf)

    def test_Atom_N2(self):
        N2 = pn.Atom('N',2)
        with self.subTest():
            self.assertEqual(N2.atom, 'N2')
        with self.subTest():
            self.assertEqual(N2.elem, 'N')
        with self.subTest():
            self.assertEqual(N2.spec, 2)
        with self.subTest():
            self.assertEqual(N2.name, 'nitrogen')
        with self.subTest():
            self.assertEqual(N2.Z, 7)
        with self.subTest():
            self.assertIsInstance(N2.NLevels, np.int64)
        with self.subTest():
            self.assertEqual(N2.energy_eV[0,0], np.inf)

class TestRecAtom(unittest.TestCase):

    def test_RecAtom_atom(self):
        H1 = pn.RecAtom('H',1)
        self.assertEqual(H1.atom, 'H1')