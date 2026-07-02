"""
Tests for pn.ICF (ionization correction factors).

Each test builds its own ICF instance (the icf dictionary is per-instance,
so addICF/delICF cannot leak between tests). icf_list is always passed
explicitly: evaluating all ~200 ICFs is slow and floods the log.
"""
import numpy as np
import pytest

import pyneb as pn


ATOM_ABUN = {'O2': 1e-5, 'O3': 2e-5, 'N2': 1e-6, 'Ne3': 1.2e-5}


class TestAvailableICFs:

    def test_all_icfs(self):
        icf = pn.ICF()
        avail = icf.getAvailableICFs()
        assert 'O' in avail
        assert 'direct_O.23' in avail['O']

    def test_filtered_by_type(self):
        icf = pn.ICF()
        avail = icf.getAvailableICFs('S', type_=['HII'])
        for elem in avail:
            for label in avail[elem]:
                assert icf.getType(label) == 'HII'


class TestICFMetadata:

    def test_getExpression(self):
        icf = pn.ICF()
        assert icf.getExpression('direct_O.23') == 'O = O2 + O3'

    def test_getType(self):
        icf = pn.ICF()
        assert icf.getType('TPP77_14') == 'PNe'

    def test_getComment(self):
        icf = pn.ICF()
        assert icf.getComment('direct_O.23') == 'just summing visible ions'

    def test_getReference(self):
        icf = pn.ICF()
        assert 'Torres-Peimbert' in icf.getReference('TPP77_14')

    def test_getURL(self):
        icf = pn.ICF()
        url = icf.getURL('TPP77_14')
        assert isinstance(url, str)
        assert url.startswith('http')


class TestGetElemAbundance:

    def test_direct_sum(self):
        icf = pn.ICF()
        res = icf.getElemAbundance(ATOM_ABUN, icf_list=['direct_O.23'])
        assert res['direct_O.23'] == pytest.approx(3e-5)

    def test_icf_factor(self):
        icf = pn.ICF()
        res = icf.getElemAbundance(ATOM_ABUN, icf_list=['TPP77_14', 'TPP77_15'])
        # N = N2 * (O2 + O3) / O2 -> 1e-6 * 3
        assert res['TPP77_14'] == pytest.approx(3e-6)
        # Ne = Ne3 * (O2 + O3) / O3 -> 1.2e-5 * 1.5
        assert res['TPP77_15'] == pytest.approx(1.8e-5)
        # an ICF must be >= 1 (it corrects for unseen ionization stages)
        assert icf.icf_value['TPP77_14'] == pytest.approx(3.0)
        assert icf.icf_value['TPP77_14'] >= 1.
        assert icf.icf_value['TPP77_15'] >= 1.

    def test_array_input(self):
        icf = pn.ICF()
        abun = {key: np.array([value, 2 * value]) for key, value in ATOM_ABUN.items()}
        res = icf.getElemAbundance(abun, icf_list=['direct_O.23', 'TPP77_14'])
        assert res['direct_O.23'] == pytest.approx([3e-5, 6e-5])
        assert res['TPP77_14'] == pytest.approx([3e-6, 6e-6])

    def test_absent_ion_gives_nan(self):
        icf = pn.ICF()
        res = icf.getElemAbundance({'O2': 1e-5}, icf_list=['direct_O.23'])
        assert np.isnan(res['direct_O.23'])


class TestAddDelICF:

    def test_addICF_delICF(self):
        icf = pn.ICF()
        icf.addICF(label='TEST_1',
                   elem='Ne',
                   atom='abun["Ne2"] + abun["Ne3"]',
                   icf='(abun["O2"] + abun["O3"]) / abun["O2"]',
                   type_='PNe',
                   comment='test icf')
        assert 'TEST_1' in icf.all_icfs
        abun = {'O2': 1e-5, 'O3': 2e-5, 'Ne2': 1e-5, 'Ne3': 1.2e-5}
        res = icf.getElemAbundance(abun, icf_list=['TEST_1'])
        assert res['TEST_1'] == pytest.approx(2.2e-5 * 3.)
        icf.delICF('TEST_1')
        assert 'TEST_1' not in icf.all_icfs
