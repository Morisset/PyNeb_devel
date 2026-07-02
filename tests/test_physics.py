"""
Tests for pyneb.utils.physics: air/vacuum wavelength conversion, ground-state
configurations, physical constants, and the Hbeta emissivity fit.
"""
import numpy as np
import pytest

import pyneb as pn
from pyneb.core.pynebcore import getHbEmissivity
from pyneb.utils.physics import airtovac, vactoair, gsFromAtom, CST


class TestAirVac:

    def test_airtovac_reference_value(self):
        # Reference value from the IDL astrolib AIRTOVAC documentation
        assert airtovac(np.array([6056.125]))[0] == pytest.approx(6057.8019, abs=0.01)

    def test_roundtrip(self):
        w = np.array([3000., 5000., 7000., 9000.])
        assert vactoair(airtovac(w)) == pytest.approx(w, rel=1e-5)

    def test_vacuum_larger_than_air(self):
        w = np.array([4000., 6000., 8000.])
        assert np.all(airtovac(w) > w)

    def test_below_conversion_window_unchanged(self):
        low_wl = pn.config.vactoair_low_wl
        w = np.array([low_wl / 2.])
        assert airtovac(w)[0] == pytest.approx(w[0])
        assert vactoair(w)[0] == pytest.approx(w[0])


class TestGsFromAtom:

    def test_ground_states(self):
        assert gsFromAtom('O3') == 'p2'
        assert gsFromAtom('O2') == 'p3'
        assert gsFromAtom('N2') == 'p2'
        assert gsFromAtom('S2') == 'p3'


class TestCST:

    def test_hbeta_wavelength(self):
        assert CST.HBETA == pytest.approx(4861.33, abs=0.01)

    def test_rydberg(self):
        assert CST.RYD_EV == pytest.approx(13.6057, rel=1e-4)
        assert CST.RYD_ANG == pytest.approx(911.267, rel=1e-4)

    def test_boltzmann(self):
        assert CST.BOLTZMANN_eVK == pytest.approx(8.6173e-5, rel=1e-3)


class TestGetHbEmissivity:

    def test_value_at_1e4(self):
        # Aller (1984) fit: 1.387e-25 / t**0.983 / 10**(-0.0424/t) with t = 1
        expected = 1.387e-25 / 10 ** (0.0424)
        assert getHbEmissivity(tem=1e4) == pytest.approx(expected, rel=1e-3)

    def test_against_recatom(self, H1):
        # The Aller fit should agree with the H I recombination data within a few %
        fit = getHbEmissivity(tem=1e4)
        data = H1.getEmissivity(1e4, 1e2, label='4_2')
        assert fit == pytest.approx(data, rel=0.05)
