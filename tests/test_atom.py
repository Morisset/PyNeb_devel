import unittest
import numpy as np
import pyneb as pn


class TestAtomInstantiation(unittest.TestCase):
    def test_atom_from_elem_spec(self):
        O3 = pn.Atom('O', 3)
        with self.subTest():
            self.assertEqual(O3.atom, 'O3')
        with self.subTest():
            self.assertEqual(O3.elem, 'O')
        with self.subTest():
            self.assertEqual(O3.spec, 3)

    def test_atom_from_atom_string(self):
        N2 = pn.Atom(atom='N2')
        with self.subTest():
            self.assertEqual(N2.atom, 'N2')
        with self.subTest():
            self.assertEqual(N2.elem, 'N')
        with self.subTest():
            self.assertEqual(N2.spec, 2)

    def test_atom_is_valid(self):
        O3 = pn.Atom('O', 3)
        self.assertTrue(O3.is_valid)

    def test_atom_nlevels_positive(self):
        O3 = pn.Atom('O', 3)
        self.assertGreater(O3.NLevels, 0)

    def test_repr(self):
        O3 = pn.Atom('O', 3)
        rep = repr(O3)
        self.assertIn('O3', rep)


class TestAtomGetA(unittest.TestCase):
    def setUp(self):
        self.O3 = pn.Atom('O', 3)

    def test_getA_specific_transition(self):
        a = self.O3.getA(4, 3)
        self.assertIsInstance(a, float)
        self.assertGreater(a, 0.)

    def test_getA_by_wave(self):
        a = self.O3.getA(wave=5007)
        self.assertIsInstance(a, float)
        self.assertGreater(a, 0.)

    def test_getA_all_transitions(self):
        a_all = self.O3.getA()
        # getA() returns a (NLevels+1) x (NLevels+1) array (1-indexed storage)
        self.assertEqual(a_all.ndim, 2)
        self.assertGreater(a_all.shape[0], self.O3.NLevels - 1)


class TestAtomGetStatWeight(unittest.TestCase):
    def setUp(self):
        self.O3 = pn.Atom('O', 3)

    def test_statweight_ground_state(self):
        w = self.O3.getStatWeight(1)
        self.assertIsInstance(float(w), float)
        self.assertGreater(w, 0.)

    def test_statweight_all_levels(self):
        for lev in range(1, self.O3.NLevels + 1):
            with self.subTest(level=lev):
                w = self.O3.getStatWeight(lev)
                self.assertGreater(w, 0.)


class TestAtomGetEnergy(unittest.TestCase):
    def setUp(self):
        self.O3 = pn.Atom('O', 3)

    def test_energy_ground_state_zero(self):
        e = self.O3.getEnergy(1, unit='1/Ang')
        self.assertEqual(e, 0.)

    def test_energy_eV_increases(self):
        energies = [self.O3.getEnergy(i, unit='eV') for i in range(1, self.O3.NLevels + 1)]
        for i in range(len(energies) - 1):
            with self.subTest(level=i + 1):
                self.assertLess(energies[i], energies[i + 1])

    def test_energy_cm1_units(self):
        e = self.O3.getEnergy(2, unit='cm-1')
        self.assertGreater(e, 0.)


class TestAtomGetOmega(unittest.TestCase):
    def setUp(self):
        self.O3 = pn.Atom('O', 3)

    def test_getOmega_scalar(self):
        omega = self.O3.getOmega(1e4, 4, 3)
        self.assertIsInstance(float(omega), float)
        self.assertGreater(omega, 0.)

    def test_getOmega_array(self):
        tems = np.array([8000., 10000., 12000.])
        omega = self.O3.getOmega(tems, 4, 3)
        self.assertEqual(omega.shape, (3,))
        self.assertTrue(np.all(omega > 0.))


class TestAtomGetCollRates(unittest.TestCase):
    def setUp(self):
        self.O3 = pn.Atom('O', 3)

    def test_getcollrates_shape(self):
        rates = self.O3.getCollRates(1e4)
        self.assertEqual(rates.shape, (self.O3.NLevels, self.O3.NLevels))

    def test_getcollrates_dtype(self):
        rates = self.O3.getCollRates(1e4)
        self.assertEqual(rates.dtype, np.float64)

    def test_getcollrates_nonnegative(self):
        rates = self.O3.getCollRates(1e4)
        self.assertTrue(np.all(rates >= 0.))


class TestAtomGetTransition(unittest.TestCase):
    def setUp(self):
        self.O3 = pn.Atom('O', 3)
        self.S2 = pn.Atom('S', 2)

    def test_gettransition_5007(self):
        lev_i, lev_j = self.O3.getTransition(5007)
        self.assertEqual(int(lev_i), 4)
        self.assertEqual(int(lev_j), 3)

    def test_gettransition_4363(self):
        lev_i, lev_j = self.O3.getTransition(4363)
        self.assertEqual(int(lev_i), 5)
        self.assertEqual(int(lev_j), 4)

    def test_gettransition_density_diagnostic(self):
        lev_i, lev_j = self.S2.getTransition(6716)
        self.assertIsNotNone(lev_i)
        self.assertGreater(int(lev_i), int(lev_j))


class TestAtomGetPopulations(unittest.TestCase):
    def setUp(self):
        self.O3 = pn.Atom('O', 3)

    def test_populations_sum_to_one(self):
        pop = self.O3.getPopulations(1e4, 100.)
        np.testing.assert_allclose(pop.sum(), 1.0, rtol=1e-5)

    def test_populations_nonnegative(self):
        pop = self.O3.getPopulations(1e4, 100.)
        self.assertTrue(np.all(pop >= 0.))

    def test_populations_shape(self):
        pop = self.O3.getPopulations(1e4, 100.)
        self.assertEqual(pop.shape, (self.O3.NLevels,))

    def test_populations_product_mode(self):
        tems = np.array([8000., 10000., 12000.])
        dens = np.array([100., 1000.])
        pop = self.O3.getPopulations(tems, dens, product=True)
        # shape: (NLevels, n_tem, n_den)
        self.assertEqual(pop.shape, (self.O3.NLevels, 3, 2))

    def test_populations_nonproduct_mode(self):
        tems = np.array([8000., 10000., 12000.])
        dens = np.array([100., 500., 1000.])
        pop = self.O3.getPopulations(tems, dens, product=False)
        self.assertEqual(pop.shape, (self.O3.NLevels, 3))


class TestAtomGetCritDensity(unittest.TestCase):
    def setUp(self):
        self.O3 = pn.Atom('O', 3)

    def test_critdensity_shape(self):
        crit = self.O3.getCritDensity(1e4)
        self.assertEqual(crit.shape, (self.O3.NLevels,))

    def test_critdensity_positive_for_excited(self):
        crit = self.O3.getCritDensity(1e4)
        # Excited levels should have positive critical density
        self.assertTrue(np.any(crit > 0.))


class TestAtomGetEmissivity(unittest.TestCase):
    def setUp(self):
        self.O3 = pn.Atom('O', 3)

    def test_emissivity_scalar(self):
        emis = self.O3.getEmissivity(1e4, 100., wave=5007)
        self.assertIsInstance(float(emis), float)
        self.assertGreater(emis, 0.)

    def test_emissivity_by_levels(self):
        emis = self.O3.getEmissivity(1e4, 100., lev_i=4, lev_j=3)
        self.assertGreater(emis, 0.)

    def test_emissivity_consistency(self):
        emis_wave = self.O3.getEmissivity(1e4, 100., wave=5007)
        emis_lev = self.O3.getEmissivity(1e4, 100., lev_i=4, lev_j=3)
        np.testing.assert_allclose(emis_wave, emis_lev, rtol=1e-6)

    def test_emissivity_array_temps(self):
        tems = np.array([8000., 10000., 12000.])
        dens = np.array([100., 500., 1000.])
        emis = self.O3.getEmissivity(tems, dens, wave=5007, product=False)
        self.assertEqual(emis.shape, (3,))
        self.assertTrue(np.all(emis > 0.))

    def test_emissivity_product_mode(self):
        tems = np.array([8000., 10000.])
        dens = np.array([100., 1000.])
        emis = self.O3.getEmissivity(tems, dens, wave=5007, product=True)
        self.assertEqual(emis.shape, (2, 2))

    def test_emissivity_all_transitions(self):
        emis_all = self.O3.getEmissivity(1e4, 100.)
        self.assertEqual(emis_all.shape,
                         (self.O3.NLevels, self.O3.NLevels))


class TestAtomGetTemDen(unittest.TestCase):
    def setUp(self):
        self.O3 = pn.Atom('O', 3)
        self.S2 = pn.Atom('S', 2)

    def test_temperature_from_ratio(self):
        tem = self.O3.getTemDen(150., den=100., wave1=5007, wave2=4363)
        self.assertGreater(tem, 5000.)
        self.assertLess(tem, 30000.)

    def test_density_from_ratio(self):
        den = self.S2.getTemDen(1.0, tem=1e4, wave1=6731, wave2=6716)
        self.assertGreater(den, 10.)
        self.assertLess(den, 1e7)

    def test_temden_array_input(self):
        tems = self.O3.getTemDen([150., 200.], den=[100., 200.], wave1=5007, wave2=4363)
        self.assertEqual(tems.shape, (2,))
        self.assertTrue(np.all(tems > 5000.))

    def test_temden_roundtrip(self):
        true_tem = 10000.
        emis_ratio = (self.O3.getEmissivity(true_tem, 100., wave=4363) /
                      self.O3.getEmissivity(true_tem, 100., wave=5007))
        recovered_tem = self.O3.getTemDen(emis_ratio, den=100., wave1=4363, wave2=5007)
        np.testing.assert_allclose(recovered_tem, true_tem, rtol=1e-2)


class TestAtomGetIonAbundance(unittest.TestCase):
    def setUp(self):
        self.O3 = pn.Atom('O', 3)

    def test_ionabundance_returns_float(self):
        abund = self.O3.getIonAbundance(100., 1.5e4, 100., wave=5007)
        self.assertIsInstance(float(abund), float)
        self.assertGreater(abund, 0.)

    def test_ionabundance_by_levels(self):
        abund = self.O3.getIonAbundance(100., 1.5e4, 100., lev_i=4, lev_j=3)
        self.assertGreater(abund, 0.)

    def test_ionabundance_to_eval(self):
        abund = self.O3.getIonAbundance(130., 1.5e4, 100., to_eval='I(4,3) + I(4,2)')
        self.assertGreater(abund, 0.)


class TestAtomDensityRatios(unittest.TestCase):
    def setUp(self):
        self.O3 = pn.Atom('O', 3)

    def test_lowdens_ratio_finite(self):
        low = self.O3.getLowDensRatio(wave1=5007, wave2=4363)
        self.assertTrue(np.isfinite(low))
        self.assertGreater(low, 0.)

    def test_highdens_ratio_finite(self):
        high = self.O3.getHighDensRatio(wave1=5007, wave2=4363)
        self.assertTrue(np.isfinite(high))
        self.assertGreater(high, 0.)

    def test_lowdens_greater_than_highdens(self):
        # For temperature-sensitive O3 ratio, low density limit > high density limit
        low = self.O3.getLowDensRatio(wave1=5007, wave2=4363)
        high = self.O3.getHighDensRatio(wave1=5007, wave2=4363)
        self.assertGreater(low, high)


if __name__ == '__main__':
    unittest.main()
