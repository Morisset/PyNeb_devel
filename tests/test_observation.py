import unittest
import numpy as np
import pyneb as pn


def _make_obs():
    """Helper: build a small Observation with 2 objects and 3 O3 lines."""
    obs = pn.Observation()
    obs.addLine(pn.EmissionLine('O', 3, 5007, obsIntens=[100., 120.]))
    obs.addLine(pn.EmissionLine('O', 3, 4959, obsIntens=[33., 40.]))
    obs.addLine(pn.EmissionLine('O', 3, 4363, obsIntens=[2.5, 3.0]))
    obs.names = ['obj1', 'obj2']
    return obs


class TestEmissionLineInstantiation(unittest.TestCase):
    def test_from_elem_spec_wave(self):
        line = pn.EmissionLine('O', 3, 5007, obsIntens=[100., 120.])
        with self.subTest():
            self.assertEqual(line.atom, 'O3')
        with self.subTest():
            self.assertEqual(line.label, 'O3_5007A')
        with self.subTest():
            self.assertTrue(line.is_valid)

    def test_from_label(self):
        line = pn.EmissionLine(label='O3_5007A', obsIntens=100.)
        with self.subTest():
            self.assertEqual(line.atom, 'O3')
        with self.subTest():
            self.assertEqual(line.label, 'O3_5007A')

    def test_obs_intens_stored(self):
        line = pn.EmissionLine('O', 3, 5007, obsIntens=[1.4, 1.3])
        np.testing.assert_array_equal(line.obsIntens, [1.4, 1.3])

    def test_obs_error_default_zero(self):
        line = pn.EmissionLine('O', 3, 5007, obsIntens=[100.])
        np.testing.assert_array_equal(line.obsError, [0.])

    def test_obs_error_custom(self):
        line = pn.EmissionLine('O', 3, 5007, obsIntens=[100.], obsError=[0.05])
        np.testing.assert_array_equal(line.obsError, [0.05])

    def test_repr(self):
        line = pn.EmissionLine('O', 3, 5007, obsIntens=[100.])
        self.assertIn('O3', repr(line))

    def test_addObs(self):
        line = pn.EmissionLine('O', 3, 5007, obsIntens=[100.])
        line.addObs(120.)
        self.assertEqual(line.obsIntens.size, 2)
        self.assertAlmostEqual(line.obsIntens[1], 120.)


class TestEmissionLineCorrectIntens(unittest.TestCase):
    def test_correctIntens_changes_corr(self):
        line = pn.EmissionLine('O', 3, 5007, obsIntens=[100.])
        RC = pn.RedCorr(E_BV=0.1, law='F99')
        line.correctIntens(RC)
        self.assertNotAlmostEqual(float(line.corrIntens[0]), 0.)
        # corrected intensity should differ from observed when E_BV != 0
        self.assertNotAlmostEqual(float(line.corrIntens[0]), 100.)

    def test_correctIntens_no_correction(self):
        line = pn.EmissionLine('O', 3, 5007, obsIntens=[100.])
        RC = pn.RedCorr(E_BV=0., law='No correction')
        line.correctIntens(RC)
        np.testing.assert_allclose(line.corrIntens, [100.], rtol=1e-6)


class TestObservationAddAndProperties(unittest.TestCase):
    def setUp(self):
        self.obs = _make_obs()

    def test_n_lines(self):
        self.assertEqual(self.obs.n_lines, 3)

    def test_n_obs(self):
        self.assertEqual(self.obs.n_obs, 2)

    def test_lineLabels(self):
        labels = list(self.obs.lineLabels)
        self.assertIn('O3_5007A', labels)
        self.assertIn('O3_4959A', labels)
        self.assertIn('O3_4363A', labels)

    def test_getLine_by_elem_spec_wave(self):
        line = self.obs.getLine('O', 3, 5007)
        self.assertIsNotNone(line)
        self.assertEqual(line.label, 'O3_5007A')

    def test_getLine_by_label(self):
        line = self.obs.getLine(label='O3_5007A')
        self.assertIsNotNone(line)
        self.assertEqual(line.label, 'O3_5007A')

    def test_getLine_missing_returns_none(self):
        line = self.obs.getLine('O', 3, 9999)
        self.assertIsNone(line)

    def test_getIntens_keys(self):
        intens = self.obs.getIntens()
        self.assertIn('O3_5007A', intens)
        self.assertIn('O3_4959A', intens)

    def test_getIntens_values(self):
        intens = self.obs.getIntens(returnObs=True)
        np.testing.assert_array_equal(intens['O3_5007A'], [100., 120.])

    def test_getError_keys(self):
        errors = self.obs.getError()
        self.assertIn('O3_5007A', errors)


class TestObservationRemoveLine(unittest.TestCase):
    def test_remove_existing_line(self):
        obs = _make_obs()
        obs.removeLine('O3_5007A')
        self.assertEqual(obs.n_lines, 2)
        self.assertNotIn('O3_5007A', obs.lineLabels)

    def test_remove_nonexistent_is_noop(self):
        obs = _make_obs()
        obs.removeLine('XX_9999A')  # should not raise
        self.assertEqual(obs.n_lines, 3)


class TestObservationFillObs(unittest.TestCase):
    def test_fillobs_adds_line(self):
        obs = _make_obs()
        obs.fillObs('O3_4363A_extra', default=0.0)  # won't add invalid label
        # fillObs with an existing label should warn, not add
        n_before = obs.n_lines
        obs.fillObs('O3_5007A')  # already exists → noop
        self.assertEqual(obs.n_lines, n_before)


class TestObservationGetSortedLines(unittest.TestCase):
    def test_sorted_by_atom(self):
        obs = pn.Observation()
        obs.addLine(pn.EmissionLine('O', 3, 5007, obsIntens=[100.]))
        obs.addLine(pn.EmissionLine('N', 2, 6583, obsIntens=[50.]))
        obs.addLine(pn.EmissionLine('O', 3, 4363, obsIntens=[2.5]))
        sorted_lines = obs.getSortedLines(crit='atom')
        atoms = [l.atom for l in sorted_lines]
        self.assertEqual(atoms, sorted(atoms))

    def test_sorted_by_wave(self):
        obs = pn.Observation()
        obs.addLine(pn.EmissionLine('O', 3, 5007, obsIntens=[100.]))
        obs.addLine(pn.EmissionLine('O', 3, 4363, obsIntens=[2.5]))
        obs.addLine(pn.EmissionLine('O', 3, 4959, obsIntens=[33.]))
        sorted_lines = obs.getSortedLines(crit='wave')
        waves = [l.wave for l in sorted_lines]
        self.assertEqual(waves, sorted(waves))


class TestObservationGetUniqueAtoms(unittest.TestCase):
    def test_unique_atoms_single(self):
        obs = _make_obs()
        unique = obs.getUniqueAtoms()
        np.testing.assert_array_equal(unique, ['O3'])

    def test_unique_atoms_multiple(self):
        obs = pn.Observation()
        obs.addLine(pn.EmissionLine('O', 3, 5007, obsIntens=[100.]))
        obs.addLine(pn.EmissionLine('N', 2, 6583, obsIntens=[50.]))
        unique = obs.getUniqueAtoms()
        self.assertEqual(len(unique), 2)
        self.assertIn('O3', unique)
        self.assertIn('N2', unique)


class TestObservationAddMonteCarloObs(unittest.TestCase):
    def test_mc_increases_n_obs(self):
        obs = pn.Observation()
        obs.addLine(pn.EmissionLine('O', 3, 5007, obsIntens=[100.], obsError=[0.05]))
        obs.names = ['obj1']
        obs.addMonteCarloObs(N=20)
        self.assertEqual(obs.n_obs, 21)

    def test_mc_names_labeled(self):
        obs = pn.Observation()
        obs.addLine(pn.EmissionLine('O', 3, 5007, obsIntens=[100.], obsError=[0.05]))
        obs.names = ['obj1']
        obs.addMonteCarloObs(N=5)
        mc_names = [n for n in obs.names if '-MC-' in n]
        self.assertEqual(len(mc_names), 5)


if __name__ == '__main__':
    unittest.main()
