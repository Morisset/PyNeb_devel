"""
Comprehensive tests for pyneb.RedCorr (reddening correction).

Covers: constructor defaults, E_BV/cHbeta roundtrip, getCorr,
        getCorrHb, setCorr, available laws.
"""
import numpy as np
import pytest
import pyneb as pn


# ---------------------------------------------------------------------------
# A — Constructor
# ---------------------------------------------------------------------------

class TestRedCorrConstructor:
    def test_defaults(self):
        rc = pn.RedCorr()
        assert rc.E_BV == pytest.approx(0.0)
        assert rc.R_V == pytest.approx(3.1)
        assert rc.law == 'No correction'

    def test_constructor_E_BV(self):
        rc = pn.RedCorr(E_BV=0.5, law='CCM89')
        assert rc.E_BV == pytest.approx(0.5)
        assert rc.cHbeta > 0.

    def test_constructor_cHbeta(self):
        rc = pn.RedCorr(cHbeta=0.5, law='CCM89')
        assert rc.cHbeta == pytest.approx(0.5)
        assert rc.E_BV > 0.

    def test_E_BV_cHbeta_roundtrip(self):
        rc = pn.RedCorr(law='CCM89', E_BV=0.5)
        chb = rc.cHbeta
        rc2 = pn.RedCorr(law='CCM89', cHbeta=chb)
        assert rc2.E_BV == pytest.approx(0.5, rel=1e-4)


# ---------------------------------------------------------------------------
# B — getCorr
# ---------------------------------------------------------------------------

class TestGetCorr:
    def test_no_correction_law_returns_unity(self):
        rc = pn.RedCorr(law='No correction', E_BV=1.0)
        assert rc.getCorr(5007.) == pytest.approx(1.0)
        assert rc.getCorr(4861.) == pytest.approx(1.0)
        assert rc.getCorr(6563.) == pytest.approx(1.0)

    def test_zero_E_BV_returns_unity_CCM89(self):
        rc = pn.RedCorr(law='CCM89', E_BV=0.0)
        assert rc.getCorr(5007.) == pytest.approx(1.0, abs=1e-10)

    def test_zero_E_BV_returns_unity_F99(self):
        rc = pn.RedCorr(law='F99', E_BV=0.0)
        assert rc.getCorr(5007.) == pytest.approx(1.0, abs=1e-10)

    def test_positive_E_BV_gives_correction_above_one(self):
        rc = pn.RedCorr(law='CCM89', E_BV=0.5)
        assert rc.getCorr(4861.) > 1.0

    def test_shorter_wavelength_needs_more_correction(self):
        rc = pn.RedCorr(law='CCM89', E_BV=0.5)
        c_uv = rc.getCorr(4000.)
        c_opt = rc.getCorr(5007.)
        c_red = rc.getCorr(6563.)
        assert c_uv > c_opt > c_red

    def test_array_wavelength_input(self):
        rc = pn.RedCorr(law='CCM89', E_BV=0.5)
        waves = np.array([4000., 5007., 6563.])
        corr = rc.getCorr(waves)
        assert corr.shape == (3,)
        assert np.all(corr > 0)


# ---------------------------------------------------------------------------
# C — getCorrHb
# ---------------------------------------------------------------------------

class TestGetCorrHb:
    def test_at_Hbeta_is_approximately_one(self):
        rc = pn.RedCorr(law='CCM89', E_BV=1.0)
        assert rc.getCorrHb(4861.) == pytest.approx(1.0, rel=1e-3)

    def test_at_Halpha_less_than_one(self):
        rc = pn.RedCorr(law='CCM89', E_BV=0.5)
        assert rc.getCorrHb(6563.) < 1.0

    def test_at_OII_greater_than_one(self):
        rc = pn.RedCorr(law='CCM89', E_BV=0.5)
        assert rc.getCorrHb(3727.) > 1.0

    def test_monotonic_with_wavelength(self):
        rc = pn.RedCorr(law='CCM89', E_BV=0.5)
        c_uv = rc.getCorrHb(3727.)
        c_hb = rc.getCorrHb(4861.)
        c_ha = rc.getCorrHb(6563.)
        assert c_uv > c_hb > c_ha


# ---------------------------------------------------------------------------
# D — setCorr
# ---------------------------------------------------------------------------

class TestSetCorr:
    def test_theoretical_ratio_gives_zero_extinction(self):
        rc = pn.RedCorr(law='CCM89')
        rc.setCorr(2.85 / 2.85, 6563., 4861.)
        assert rc.E_BV == pytest.approx(0.0, abs=1e-5)

    def test_reddened_ratio_gives_positive_extinction(self):
        rc = pn.RedCorr(law='CCM89')
        rc.setCorr(6.5 / 2.85, 6563., 4861.)
        assert rc.E_BV > 0.

    def test_setCorr_roundtrip(self):
        # setCorr takes obs_over_theo = (obs_ratio) / (theoretical_ratio)
        # obs_Ha/obs_Hb = 2.85 * getCorr(Hb) / getCorr(Ha)
        # obs_over_theo = getCorr(Hb) / getCorr(Ha)
        rc = pn.RedCorr(law='CCM89', E_BV=0.5)
        ha_corr = rc.getCorr(6563.)
        hb_corr = rc.getCorr(4861.)
        obs_over_theo = hb_corr / ha_corr  # = obs_Ha/obs_Hb / 2.85

        rc2 = pn.RedCorr(law='CCM89')
        rc2.setCorr(obs_over_theo, 6563., 4861.)
        assert rc2.E_BV == pytest.approx(0.5, rel=0.01)


# ---------------------------------------------------------------------------
# E — Available laws
# ---------------------------------------------------------------------------

class TestAvailableLaws:
    def test_printLaws_no_exception(self):
        rc = pn.RedCorr()
        rc.printLaws()

    def test_known_laws_available(self):
        rc = pn.RedCorr()
        laws = rc.getLaws()
        for expected in ['CCM89', 'F99', 'Cal00', 'No correction', 'S79 H83 CCM89']:
            assert expected in laws

    def test_at_least_12_laws(self):
        rc = pn.RedCorr()
        assert len(rc.getLaws()) >= 12

    def test_law_switching_gives_different_results(self):
        rc = pn.RedCorr(E_BV=0.5)
        rc.law = 'CCM89'
        c1 = rc.getCorr(5007.)
        rc.law = 'F99'
        c2 = rc.getCorr(5007.)
        assert c1 != pytest.approx(c2, rel=1e-6)
