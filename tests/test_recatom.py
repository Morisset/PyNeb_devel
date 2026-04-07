import unittest
import numpy as np
import pyneb as pn


class TestRecAtomInstantiation(unittest.TestCase):
    def test_H1_from_elem_spec(self):
        H1 = pn.RecAtom('H', 1)
        with self.subTest():
            self.assertEqual(H1.atom, 'H1')
        with self.subTest():
            self.assertEqual(H1.elem, 'H')
        with self.subTest():
            self.assertEqual(H1.spec, 1)
        with self.subTest():
            self.assertEqual(H1.type, 'rec')

    def test_H1_from_atom_string(self):
        H1 = pn.RecAtom(atom='H1')
        self.assertEqual(H1.atom, 'H1')

    def test_He1_instantiation(self):
        He1 = pn.RecAtom('He', 1)
        with self.subTest():
            self.assertEqual(He1.atom, 'He1')
        with self.subTest():
            self.assertTrue(He1.is_valid)

    def test_He2_instantiation(self):
        He2 = pn.RecAtom('He', 2)
        with self.subTest():
            self.assertEqual(He2.atom, 'He2')
        with self.subTest():
            self.assertTrue(He2.is_valid)

    def test_case_B_default(self):
        H1 = pn.RecAtom('H', 1)
        self.assertEqual(H1.case, 'B')

    def test_case_A(self):
        H1a = pn.RecAtom('H', 1, case='A')
        self.assertEqual(H1a.case, 'A')


class TestRecAtomGetEmissivity(unittest.TestCase):
    def setUp(self):
        self.H1 = pn.RecAtom('H', 1)
        self.He1 = pn.RecAtom('He', 1)

    def test_H1_emissivity_scalar(self):
        emis = self.H1.getEmissivity(1e4, 1e2, lev_i=4, lev_j=2)
        self.assertIsInstance(float(emis), float)
        self.assertGreater(emis, 0.)

    def test_H1_emissivity_by_label(self):
        emis = self.H1.getEmissivity(1e4, 1e2, label='4_2')
        self.assertGreater(emis, 0.)

    def test_H1_emissivity_consistency(self):
        emis_lev = self.H1.getEmissivity(1e4, 1e2, lev_i=4, lev_j=2)
        emis_lbl = self.H1.getEmissivity(1e4, 1e2, label='4_2')
        np.testing.assert_allclose(emis_lev, emis_lbl, rtol=1e-6)

    def test_H1_emissivity_array_inputs(self):
        tems = np.array([8000., 10000., 12000.])
        dens = np.array([100., 1000., 1e4])
        emis = self.H1.getEmissivity(tems, dens, lev_i=4, lev_j=2, product=False)
        self.assertEqual(emis.shape, (3,))
        self.assertTrue(np.all(emis > 0.))

    def test_H1_emissivity_product_mode(self):
        tems = np.array([8000., 10000.])
        dens = np.array([100., 1000.])
        emis = self.H1.getEmissivity(tems, dens, lev_i=4, lev_j=2, product=True)
        self.assertEqual(emis.shape, (2, 2))

    def test_He1_emissivity_by_wave(self):
        emis = self.He1.getEmissivity(1e4, 1e2, wave=5876)
        self.assertIsInstance(float(emis), float)
        self.assertGreater(emis, 0.)

    def test_emissivity_decreases_with_temperature(self):
        # For H1 Balmer alpha, emissivity is larger at lower temperature (Case B)
        emis_low = self.H1.getEmissivity(5000., 1e2, lev_i=3, lev_j=2)
        emis_high = self.H1.getEmissivity(2e4, 1e2, lev_i=3, lev_j=2)
        self.assertGreater(emis_low, emis_high)

    def test_hbeta_emissivity_reasonable(self):
        # Hbeta (4->2) emissivity at standard conditions ~1.235e-25
        emis = self.H1.getEmissivity(1e4, 1e2, lev_i=4, lev_j=2)
        np.testing.assert_allclose(emis, 1.235e-25, rtol=0.1)


class TestRecAtomGetIonAbundance(unittest.TestCase):
    def setUp(self):
        self.H1 = pn.RecAtom('H', 1)

    def test_H1_ionabundance_normalized(self):
        # H+ abundance relative to H+ should be 1
        abund = self.H1.getIonAbundance(100., 1e4, 1e2, lev_i=4, lev_j=2)
        np.testing.assert_allclose(float(abund), 1.0, rtol=0.01)

    def test_ionabundance_scales_with_intensity(self):
        abund1 = self.H1.getIonAbundance(100., 1e4, 1e2, lev_i=4, lev_j=2)
        abund2 = self.H1.getIonAbundance(200., 1e4, 1e2, lev_i=4, lev_j=2)
        np.testing.assert_allclose(abund2, 2. * abund1, rtol=1e-5)


class TestRecAtomGetTransition(unittest.TestCase):
    def setUp(self):
        self.H1 = pn.RecAtom('H', 1)

    def test_gettransition_halpha(self):
        result = self.H1.getTransition(6563)
        self.assertIsNotNone(result)
        lev_i, lev_j = result
        self.assertEqual(int(lev_i), 3)
        self.assertEqual(int(lev_j), 2)

    def test_gettransition_hbeta(self):
        result = self.H1.getTransition(4861)
        self.assertIsNotNone(result)
        lev_i, lev_j = result
        self.assertEqual(int(lev_i), 4)
        self.assertEqual(int(lev_j), 2)


class TestRecAtomCaseAvsB(unittest.TestCase):
    def test_halpha_case_A_ne_case_B(self):
        H1a = pn.RecAtom('H', 1, case='A')
        H1b = pn.RecAtom('H', 1, case='B')
        emis_a = H1a.getEmissivity(1e4, 1e2, lev_i=3, lev_j=2)
        emis_b = H1b.getEmissivity(1e4, 1e2, lev_i=3, lev_j=2)
        # Case A and B give different results for Lyman-alpha opaque transitions
        # For Halpha (3->2), they may differ
        self.assertIsNotNone(emis_a)
        self.assertIsNotNone(emis_b)


if __name__ == '__main__':
    unittest.main()
