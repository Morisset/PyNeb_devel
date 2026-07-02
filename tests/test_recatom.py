"""
Comprehensive tests for pyneb.RecAtom (recombination-line ions).

Covers: constructor, getWave, getEmissivity (numerical values,
        Balmer decrement, array inputs), density insensitivity,
        and the module-level getRecEmissivity helper.
"""
import numpy as np
import pytest
import pyneb as pn


# ---------------------------------------------------------------------------
# A — Constructor
# ---------------------------------------------------------------------------

class TestRecAtomConstructor:
    def test_H1_basic_attributes(self, H1):
        assert H1.atom == 'H1'
        assert H1.elem == 'H'
        assert H1.spec == 1
        assert H1.name == 'hydrogen'
        assert H1.is_valid is True
        assert H1.type == 'rec'

    def test_He2_basic_attributes(self, He2):
        assert He2.atom == 'He2'
        assert He2.elem == 'He'
        assert He2.spec == 2

    def test_atom_keyword_H1(self):
        a = pn.RecAtom(atom='H1')
        assert a.atom == 'H1'
        assert a.elem == 'H'
        assert a.spec == 1

    def test_atom_keyword_He2(self):
        a = pn.RecAtom(atom='He2')
        assert a.atom == 'He2'


# ---------------------------------------------------------------------------
# B — Wavelength lookup
# ---------------------------------------------------------------------------

class TestGetWave:
    def test_H1_halpha(self, H1):
        assert H1.getWave(3, 2) == pytest.approx(6562.8, rel=1e-3)

    def test_H1_hbeta(self, H1):
        assert H1.getWave(4, 2) == pytest.approx(4861.3, rel=1e-3)

    def test_He2_4686(self, He2):
        assert He2.getWave(4, 3) == pytest.approx(4685.5, rel=1e-3)


# ---------------------------------------------------------------------------
# C — Emissivity
# ---------------------------------------------------------------------------

class TestRecAtomGetEmissivity:
    def test_H1_hbeta_value(self, H1):
        em = H1.getEmissivity(1e4, 1e2, lev_i=4, lev_j=2)
        assert em == pytest.approx(1.235e-25, rel=5e-3)

    def test_H1_halpha_value(self, H1):
        em = H1.getEmissivity(1e4, 1e2, lev_i=3, lev_j=2)
        assert em == pytest.approx(3.536e-25, rel=5e-3)

    def test_balmer_decrement_canonical(self, H1):
        Ha = H1.getEmissivity(1e4, 1e2, lev_i=3, lev_j=2)
        Hb = H1.getEmissivity(1e4, 1e2, lev_i=4, lev_j=2)
        assert Ha / Hb == pytest.approx(2.863, rel=1e-3)

    def test_hgamma_hbeta_ratio(self, H1):
        Hg = H1.getEmissivity(1e4, 1e2, lev_i=5, lev_j=2)
        Hb = H1.getEmissivity(1e4, 1e2, lev_i=4, lev_j=2)
        assert Hg / Hb == pytest.approx(0.468, rel=0.01)

    def test_balmer_decrement_decreases_with_temperature(self, H1):
        tems = np.array([5e3, 1e4, 1.5e4, 2e4])
        Ha = H1.getEmissivity(tems, 1e2, lev_i=3, lev_j=2)
        Hb = H1.getEmissivity(tems, 1e2, lev_i=4, lev_j=2)
        ratios = Ha / Hb
        assert np.all(np.diff(ratios) < 0)

    def test_array_tem_shape(self, H1):
        tems = np.array([8e3, 1e4, 1.2e4])
        em = H1.getEmissivity(tems, 1e2, lev_i=4, lev_j=2)
        assert em.shape == (3,)
        assert np.all(em > 0)

    def test_hbeta_decreases_with_temperature(self, H1):
        tems = np.array([5e3, 1e4, 2e4])
        em = H1.getEmissivity(tems, 1e2, lev_i=4, lev_j=2)
        assert np.all(np.diff(em) < 0)

    def test_He2_4686_value(self, He2):
        em = He2.getEmissivity(1e4, 1e2, lev_i=4, lev_j=3)
        assert em == pytest.approx(1.516e-24, rel=0.01)

    def test_H1_density_insensitive(self, H1):
        dens = np.array([1e2, 1e3, 1e4])
        em = np.array([H1.getEmissivity(1e4, d, lev_i=3, lev_j=2) for d in dens])
        variation = (em.max() - em.min()) / em.mean()
        assert variation < 0.01  # < 1% variation


# ---------------------------------------------------------------------------
# D — Module-level getRecEmissivity helper
# ---------------------------------------------------------------------------

class TestGetRecEmissivity:
    def test_hbeta_value(self):
        em = pn.getRecEmissivity(1e4, 1e2, 4, 2, atom='H1')
        assert em == pytest.approx(1.235e-25, rel=5e-3)

    def test_matches_RecAtom(self, H1):
        em_func = pn.getRecEmissivity(1e4, 1e2, 4, 2, atom='H1')
        em_atom = H1.getEmissivity(1e4, 1e2, lev_i=4, lev_j=2)
        assert em_func == pytest.approx(em_atom, rel=1e-6)
