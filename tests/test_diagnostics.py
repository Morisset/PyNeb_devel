import unittest
import numpy as np
import pyneb as pn


class TestDiagnosticsInstantiation(unittest.TestCase):
    def test_empty_diagnostics(self):
        diags = pn.Diagnostics()
        self.assertEqual(len(list(diags.getDiagLabels())), 0)

    def test_addAll_in_constructor(self):
        diags = pn.Diagnostics(addAll=True)
        labels = list(diags.getDiagLabels())
        self.assertGreater(len(labels), 0)

    def test_addAll_count_matches_getAllDiags(self):
        diags = pn.Diagnostics(addAll=True)
        all_labels = list(diags.getAllDiagLabels())
        my_labels = list(diags.getDiagLabels())
        self.assertEqual(len(my_labels), len(all_labels))


class TestDiagnosticsAddAndGet(unittest.TestCase):
    def setUp(self):
        self.diags = pn.Diagnostics()
        self.diags.addDiag('[OIII] 4363/5007')
        self.diags.addDiag('[SII] 6731/6716')

    def test_getDiagLabels_contains_added(self):
        labels = list(self.diags.getDiagLabels())
        self.assertIn('[OIII] 4363/5007', labels)
        self.assertIn('[SII] 6731/6716', labels)

    def test_getDiagFromLabel_returns_tuple(self):
        result = self.diags.getDiagFromLabel('[OIII] 4363/5007')
        self.assertIsNotNone(result)
        self.assertIsInstance(result, tuple)
        # tuple: (atom, eval_expression, error_expression)
        self.assertEqual(result[0], 'O3')

    def test_getDiagFromLabel_missing_returns_none(self):
        result = self.diags.getDiagFromLabel('[NONEXISTENT] 9999/8888')
        self.assertIsNone(result)

    def test_getAllDiags_returns_non_empty(self):
        all_diags = self.diags.getAllDiags()
        self.assertGreater(len(all_diags), 100)

    def test_getAllDiagLabels_returns_keys(self):
        all_labels = list(self.diags.getAllDiagLabels())
        self.assertIn('[OIII] 4363/5007', all_labels)


class TestDiagnosticsAddCustom(unittest.TestCase):
    def test_addDiag_custom_tuple(self):
        diags = pn.Diagnostics()
        diags.addDiag('my_diag', ('O3', 'L(5007)/L(4363)', 'RMS([E(5007), E(4363)])'))
        result = diags.getDiagFromLabel('my_diag')
        self.assertIsNotNone(result)
        self.assertEqual(result[0], 'O3')

    def test_addDiag_by_atom(self):
        diags = pn.Diagnostics()
        diags.addDiag(atom='O3')
        labels = list(diags.getDiagLabels())
        # All O3 diagnostics should be added
        for label in labels:
            self.assertEqual(diags.getDiagFromLabel(label)[0], 'O3')

    def test_addDiag_list(self):
        diags = pn.Diagnostics()
        diags.addDiag(['[OIII] 4363/5007', '[SII] 6731/6716'])
        self.assertIn('[OIII] 4363/5007', list(diags.getDiagLabels()))
        self.assertIn('[SII] 6731/6716', list(diags.getDiagLabels()))


class TestDiagnosticsDelDiag(unittest.TestCase):
    def test_delDiag_removes_entry(self):
        diags = pn.Diagnostics()
        diags.addDiag('[OIII] 4363/5007')
        diags.delDiag('[OIII] 4363/5007')
        self.assertNotIn('[OIII] 4363/5007', list(diags.getDiagLabels()))

    def test_delDiag_nonexistent_does_not_raise(self):
        diags = pn.Diagnostics()
        diags.delDiag('[NONEXISTENT] 9999/8888')  # should not raise


class TestDiagnosticsSetAtoms(unittest.TestCase):
    def test_setAtoms_injects_object(self):
        diags = pn.Diagnostics()
        diags.addDiag('[OIII] 4363/5007')
        O3 = pn.Atom('O', 3)
        diags.setAtoms({'O3': O3})
        self.assertIn('O3', diags.atomDict)
        self.assertIs(diags.atomDict['O3'], O3)


class TestDiagnosticsAddClabel(unittest.TestCase):
    def test_addClabel_extends_tuple(self):
        diags = pn.Diagnostics()
        diags.addDiag('[OIII] 4363/5007')
        diags.addClabel('[OIII] 4363/5007', 'My custom label')
        result = diags.getDiagFromLabel('[OIII] 4363/5007')
        self.assertEqual(len(result), 4)
        self.assertEqual(result[3], 'My custom label')


class TestDiagnosticsGetCrossTemDen(unittest.TestCase):
    def setUp(self):
        self.diags = pn.Diagnostics()

    def test_basic_cross_temden(self):
        tem, den = self.diags.getCrossTemDen(
            '[OIII] 4363/5007', '[SII] 6731/6716',
            value_tem=0.005, value_den=1.0
        )
        with self.subTest('temperature in physical range'):
            self.assertGreater(tem, 5000.)
            self.assertLess(tem, 30000.)
        with self.subTest('density in physical range'):
            self.assertGreater(den, 10.)
            self.assertLess(den, 1e7)

    def test_cross_temden_with_known_values(self):
        # These values give roughly T~9250 K, N~690 cm^-3
        tem, den = self.diags.getCrossTemDen(
            '[OIII] 4363/5007', '[SII] 6731/6716',
            value_tem=0.005, value_den=1.0
        )
        np.testing.assert_allclose(tem, 9250., rtol=0.05)
        np.testing.assert_allclose(den, 690., rtol=0.1)


class TestDiagnosticsEvalDiag(unittest.TestCase):
    def test_eval_diag_returns_correct_ratio(self):
        diags = pn.Diagnostics()
        diags.addDiag('[OIII] 4363/5007')
        obs = pn.Observation()
        obs.addLine(pn.EmissionLine('O', 3, 5007, obsIntens=[100.], corrected=True))
        obs.addLine(pn.EmissionLine('O', 3, 4363, obsIntens=[0.5], corrected=True))
        obs.names = ['obj1']
        val = diags.eval_diag('[OIII] 4363/5007', obs)
        # 4363/5007 = 0.5/100 = 0.005
        np.testing.assert_allclose(val, [0.005], rtol=1e-6)


if __name__ == '__main__':
    unittest.main()
