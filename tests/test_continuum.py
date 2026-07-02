"""
Tests for pn.Continuum: nebular continuum emission (bound-free, free-free,
two-photon) and the Balmer-jump temperature diagnostics.
"""
import numpy as np
import pytest

import pyneb as pn


@pytest.fixture(scope="module")
def C():
    return pn.Continuum()


class TestGetContinuum:

    def test_unnormalized_scalar(self, C):
        cont = C.get_continuum(tem=1e4, den=1e2, wl=np.array([3600., 3700.]), HI_label=None)
        assert cont.shape == (2,)
        assert np.all(cont > 0)
        # order of magnitude of the H recombination continuum (erg/s.cm3/A)
        assert 1e-28 < cont[0] < 1e-25

    def test_balmer_jump_sign(self, C):
        # the continuum just below the Balmer jump exceeds the one just above
        cont = C.get_continuum(tem=1e4, den=1e2, wl=np.array([3643., 3861.]), HI_label=None)
        assert cont[0] > cont[1]

    def test_normalized_tem_array(self, C):
        tem = np.array([8000., 10000., 12000.])
        cont = C.get_continuum(tem=tem, den=1e2, He1_H=0.1, He2_H=0.01, HI_label='11_2')
        # default wl has 5 values
        assert cont.shape == (5, 3)
        assert np.all(cont > 0)


class TestBalmerJump:

    def test_BJ_HI_positive_decreasing_with_T(self, C):
        bj_cool = C.BJ_HI(tem=5000., den=1e2, He1_H=0.1, He2_H=0.01)
        bj_hot = C.BJ_HI(tem=20000., den=1e2, He1_H=0.1, He2_H=0.01)
        assert bj_cool > 0
        assert bj_hot > 0
        # the Balmer jump decreases with temperature: this is the basis of the T_BJ diagnostic
        assert bj_cool > bj_hot

    def test_T_BJ_roundtrip(self, C):
        bj = C.BJ_HI(tem=1e4, den=1e2, He1_H=0.1, He2_H=0.01)
        tem = C.T_BJ(bj, den=1e2, He1_H=0.1, He2_H=0.01)
        assert tem == pytest.approx(1e4, rel=1e-3)

    def test_T_BJ_Liu(self, C):
        expected = 368 * (1 + 0.259 * 0.1 + 3.409 * 0.01) * 100. ** (-3. / 2)
        assert C.T_BJ_Liu(BJ_H11=100., He1_H=0.1, He2_H=0.01) == pytest.approx(expected, rel=1e-10)


class TestComponents:

    def test_two_photon(self, C):
        wl = np.array([1000., 2000., 4000.])
        tp = C.two_photon(1e4, 1e2, wl)
        # no two-photon emission below Ly-alpha
        assert tp[0] == 0.
        assert tp[1] > 0
        assert tp[2] > 0

    def test_two_photon_density_suppression(self, C):
        wl = np.array([2000., 4000.])
        tp_low = C.two_photon(1e4, 1e2, wl)
        tp_high = C.two_photon(1e4, 1e6, wl)
        # collisional de-excitation of 2s suppresses the two-photon continuum
        assert np.all(tp_high < tp_low)

    def test_gff_order_unity(self, C):
        g = C.gff(1., 1e4, np.array([4000., 6000.]))
        assert np.all(g > 0.5)
        assert np.all(g < 2.0)

    def test_FreeFree_increases_with_He(self, C):
        wl = np.array([4000.])
        ff_H = C.FreeFree(1e4, wl)
        ff_He = C.FreeFree(1e4, wl, He1_H=0.1, He2_H=0.05)
        assert ff_He[0] > ff_H[0]

    def test_make_cont_Ercolano_invalid(self, C):
        assert np.all(np.isnan(C.make_cont_Ercolano(1e4, 'XX', np.array([4000.]))))
        assert np.all(np.isnan(C.make_cont_Ercolano(50., 'H', np.array([4000.]))))
