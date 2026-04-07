import unittest
import numpy as np
import pyneb as pn


class TestICFInstantiation(unittest.TestCase):
    def test_instantiation(self):
        icf = pn.ICF()
        self.assertIsNotNone(icf)

    def test_all_icfs_not_empty(self):
        icf = pn.ICF()
        self.assertGreater(len(icf.all_icfs), 0)


class TestICFGetAvailableICFs(unittest.TestCase):
    def setUp(self):
        self.icf = pn.ICF()

    def test_returns_dict(self):
        result = self.icf.getAvailableICFs()
        self.assertIsInstance(result, dict)

    def test_contains_oxygen(self):
        result = self.icf.getAvailableICFs()
        self.assertIn('O', result)

    def test_filter_by_element(self):
        result = self.icf.getAvailableICFs('S')
        self.assertIn('S', result)
        self.assertEqual(list(result.keys()), ['S'])

    def test_direct_icf_present(self):
        result = self.icf.getAvailableICFs('O')
        all_icfs = result['O']
        self.assertIn('direct_O.23', all_icfs)


class TestICFAddAndDel(unittest.TestCase):
    def setUp(self):
        self.icf = pn.ICF()

    def test_addICF_adds_entry(self):
        self.icf.addICF(
            label='TEST_ADD_1',
            elem='Ne',
            atom='abun["Ne2"] + abun["Ne3"]',
            icf='(abun["O2"] + abun["O3"]) / abun["O2"]',
            type_='PNe',
            ref='Test reference',
            url='http://example.com'
        )
        self.assertIn('TEST_ADD_1', self.icf.all_icfs)

    def test_addICF_duplicate_warns(self):
        self.icf.addICF('TEST_DUP_1', 'O', 'abun["O2"]', '1')
        # Adding again should not raise, just warn
        self.icf.addICF('TEST_DUP_1', 'O', 'abun["O2"]', '1')
        # Still only one entry
        count = sum(1 for k in self.icf.all_icfs if k == 'TEST_DUP_1')
        self.assertEqual(count, 1)

    def test_delICF_removes_entry(self):
        self.icf.addICF('TEST_DEL_1', 'Ne', 'abun["Ne2"]', '1')
        self.icf.delICF('TEST_DEL_1')
        self.assertNotIn('TEST_DEL_1', self.icf.all_icfs)


class TestICFGetExpression(unittest.TestCase):
    def setUp(self):
        self.icf = pn.ICF()

    def test_getExpression_direct_returns_string(self):
        expr = self.icf.getExpression('direct_O.23')
        self.assertIsInstance(expr, str)

    def test_getExpression_contains_element(self):
        expr = self.icf.getExpression('direct_O.23')
        self.assertIn('O', expr)

    def test_getExpression_custom(self):
        self.icf.addICF('MYTEST_1', 'Ne',
                        'abun["Ne2"] + abun["Ne3"]',
                        '(abun["O2"] + abun["O3"]) / abun["O2"]')
        expr = self.icf.getExpression('MYTEST_1')
        self.assertIn('Ne', expr)
        self.icf.delICF('MYTEST_1')


class TestICFGetReference(unittest.TestCase):
    def setUp(self):
        self.icf = pn.ICF()

    def test_getReference_exists_for_known_icf(self):
        # 'direct_O.23' belongs to the 'direct' paper entry
        ref = self.icf.getReference('direct_O.23')
        self.assertIsNotNone(ref)

    def test_getURL_exists_for_known_icf(self):
        url = self.icf.getURL('direct_O.23')
        # May be None or a string, just should not raise
        self.assertTrue(url is None or isinstance(url, str))

    def test_getReference_added_icf(self):
        self.icf.addICF('REFTEST_1', 'O', 'abun["O2"]', '1',
                        ref='My Paper 2025', url='http://ads.example.com')
        ref = self.icf.getReference('REFTEST_1')
        self.assertEqual(ref, 'My Paper 2025')
        self.icf.delICF('REFTEST_1')

    def test_getURL_added_icf(self):
        self.icf.addICF('URLTEST_1', 'O', 'abun["O2"]', '1',
                        ref='My Paper', url='http://test.url')
        url = self.icf.getURL('URLTEST_1')
        self.assertEqual(url, 'http://test.url')
        self.icf.delICF('URLTEST_1')


class TestICFGetElemAbundance(unittest.TestCase):
    def setUp(self):
        self.icf = pn.ICF()

    def test_direct_O23_simple(self):
        atom_abun = {'O2': 0.001, 'O3': 0.002}
        result = self.icf.getElemAbundance(atom_abun, icf_list=['direct_O.23'])
        self.assertIn('direct_O.23', result)
        np.testing.assert_allclose(result['direct_O.23'], 0.003, rtol=1e-6)

    def test_result_is_sum_when_icf_is_one(self):
        atom_abun = {'O2': 1e-4, 'O3': 2e-4}
        result = self.icf.getElemAbundance(atom_abun, icf_list=['direct_O.23'])
        np.testing.assert_allclose(result['direct_O.23'], 3e-4, rtol=1e-6)

    def test_missing_ion_gives_nan(self):
        # When ions are not provided, result should be NaN
        atom_abun = {}
        result = self.icf.getElemAbundance(atom_abun, icf_list=['direct_O.23'])
        self.assertTrue(np.isnan(result['direct_O.23']))

    def test_elem_abun_stored_as_attribute(self):
        atom_abun = {'O2': 0.001, 'O3': 0.002}
        self.icf.getElemAbundance(atom_abun, icf_list=['direct_O.23'])
        self.assertIn('direct_O.23', self.icf.elem_abun)


class TestICFGetType(unittest.TestCase):
    def test_getType_returns_valid(self):
        icf = pn.ICF()
        for label in list(icf.all_icfs.keys())[:5]:
            with self.subTest(label=label):
                t = icf.getType(label)
                self.assertIn(t, ['HII', 'PNe', 'All'])


if __name__ == '__main__':
    unittest.main()
