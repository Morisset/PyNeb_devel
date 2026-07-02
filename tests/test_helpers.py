"""
Tests for the module-level helper functions of pyneb: line-label parsing,
getAtomDict, and the save/restore pickle utilities.
"""
import numpy as np
import pytest

import pyneb as pn
from pyneb.core.pynebcore import getLineLabel


class TestGetLineLabel:

    def test_simple_line(self):
        atom, wave_label, label = getLineLabel('O', 3, 5007)
        assert atom == 'O3'
        assert wave_label == '5007A'
        assert label == 'O3_5007A'

    def test_recombination_line(self):
        atom, wave_label, label = getLineLabel('H', 1, 4861, blend=False)
        assert label in ('H1_4861A', 'H1r_4861A')
        assert wave_label == '4861A'

    def test_blend(self):
        atom, wave_label, label = getLineLabel('O', 2, 3727, blend=True)
        assert label.endswith('+')
        assert atom == 'O2'


class TestParseLineLabel:

    def test_simple_line(self):
        elem, spec, atom, wave, wave_label, blend = pn.parseLineLabel('O3_5007A')
        assert elem == 'O'
        assert int(spec) == 3
        assert atom == 'O3'
        assert wave == pytest.approx(5007.)
        assert wave_label == '5007A'
        assert blend is False

    def test_blend(self):
        elem, spec, atom, wave, wave_label, blend = pn.parseLineLabel('O2_7319A+')
        assert atom == 'O2'
        assert wave == pytest.approx(7319.)
        assert blend is True

    def test_micron_wavelength(self):
        elem, spec, atom, wave, wave_label, blend = pn.parseLineLabel('S3_18.7m')
        assert atom == 'S3'
        # 18.7 microns expressed in Angstroms
        assert wave == pytest.approx(187000.)
        assert wave_label == '18.7m'

    def test_recombination_atom(self):
        elem, spec, atom, wave, wave_label, blend = pn.parseLineLabel('H1r_4861A')
        assert elem == 'H'
        assert atom == 'H1r'
        assert wave == pytest.approx(4861.)


class TestIsValid:

    def test_valid_line(self):
        assert pn.isValid('O3_5007A') is True

    def test_invalid_line(self):
        assert pn.isValid('O3_9999A') is False

    def test_valid_blend(self):
        assert pn.isValid('O2_3727A+') is True


class TestGetAtomDict:

    def test_atom_list(self):
        d = pn.getAtomDict(atom_list=['O3', 'N2'])
        # recombination counterparts (e.g. 'O3r') are added when data exist
        assert {'O3', 'N2'} <= set(d.keys())
        assert set(d.keys()) <= {'O3', 'N2', 'O3r', 'N2r'}
        assert isinstance(d['O3'], pn.Atom)
        assert d['O3'].atom == 'O3'
        assert d['N2'].elem == 'N'
        assert d['N2'].spec == 2
        if 'O3r' in d:
            assert isinstance(d['O3r'], pn.RecAtom)


class TestSaveRestore:

    def test_kwargs_roundtrip(self, tmp_path):
        file_ = str(tmp_path / 'data.pypic')
        arr = np.arange(3.)
        pn.save(file_, a=1, arr=arr)
        d = pn.restore(file_)
        assert d['a'] == 1
        assert np.array_equal(d['arr'], arr)

    def test_positional_args(self, tmp_path):
        file_ = str(tmp_path / 'data.pypic')
        my_var = 42
        pn.save(file_, 'my_var')
        d = pn.restore(file_)
        assert d['my_var'] == 42

    def test_positional_unknown_name(self, tmp_path):
        file_ = str(tmp_path / 'data.pypic')
        with pytest.raises(NameError):
            pn.save(file_, 'no_such_variable_here')
