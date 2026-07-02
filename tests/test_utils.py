"""
Tests for pyneb utility functions.

Covers: int_to_roman, roman_to_int (from pyneb/utils/misc.py),
        parseAtom, and ROOT_DIR.
"""
import os
import numpy as np
import pytest

from pyneb.utils.misc import int_to_roman, roman_to_int, parseAtom, ROOT_DIR


# ---------------------------------------------------------------------------
# A — int_to_roman
# ---------------------------------------------------------------------------

class TestIntToRoman:
    def test_basic_values(self):
        assert int_to_roman(1) == 'I'
        assert int_to_roman(2) == 'II'
        assert int_to_roman(3) == 'III'
        assert int_to_roman(4) == 'IV'
        assert int_to_roman(5) == 'V'
        assert int_to_roman(9) == 'IX'
        assert int_to_roman(10) == 'X'
        assert int_to_roman(40) == 'XL'
        assert int_to_roman(50) == 'L'
        assert int_to_roman(90) == 'XC'
        assert int_to_roman(100) == 'C'
        assert int_to_roman(400) == 'CD'
        assert int_to_roman(500) == 'D'
        assert int_to_roman(900) == 'CM'
        assert int_to_roman(1000) == 'M'
        assert int_to_roman(1999) == 'MCMXCIX'
        assert int_to_roman(2000) == 'MM'
        assert int_to_roman(3999) == 'MMMCMXCIX'

    def test_out_of_range_returns_None(self):
        assert int_to_roman(0) is None
        assert int_to_roman(-1) is None
        assert int_to_roman(4000) is None

    def test_non_integer_returns_None(self):
        assert int_to_roman(1.5) is None
        assert int_to_roman('III') is None

    def test_roundtrip(self):
        for n in range(1, 4000):
            assert roman_to_int(int_to_roman(n)) == n


# ---------------------------------------------------------------------------
# B — roman_to_int
# ---------------------------------------------------------------------------

class TestRomanToInt:
    def test_basic_values(self):
        assert roman_to_int('I') == 1
        assert roman_to_int('II') == 2
        assert roman_to_int('III') == 3
        assert roman_to_int('IV') == 4
        assert roman_to_int('V') == 5
        assert roman_to_int('IX') == 9
        assert roman_to_int('X') == 10
        assert roman_to_int('XL') == 40
        assert roman_to_int('L') == 50
        assert roman_to_int('C') == 100
        assert roman_to_int('M') == 1000
        assert roman_to_int('MCMXCIX') == 1999

    def test_case_insensitive(self):
        assert roman_to_int('iii') == 3
        assert roman_to_int('iv') == 4
        assert roman_to_int('xlii') == 42

    def test_invalid_returns_None(self):
        assert roman_to_int('IIII') is None   # invalid repetition
        assert roman_to_int('A') is None       # invalid character
        assert roman_to_int(1) is None         # non-string


# ---------------------------------------------------------------------------
# C — parseAtom
# ---------------------------------------------------------------------------

class TestParseAtom:
    def test_basic_single_char_elements(self):
        assert parseAtom('O3') == ('O', '3')
        assert parseAtom('N2') == ('N', '2')
        assert parseAtom('S2') == ('S', '2')
        assert parseAtom('H1') == ('H', '1')

    def test_multi_char_elements(self):
        assert parseAtom('He1') == ('He', '1')
        assert parseAtom('He2') == ('He', '2')
        assert parseAtom('Ar4') == ('Ar', '4')
        assert parseAtom('Fe3') == ('Fe', '3')

    def test_case_normalization(self):
        elem, spec = parseAtom('o3')
        assert elem == 'O'
        assert spec == '3'

        elem2, spec2 = parseAtom('n2')
        assert elem2 == 'N'

    def test_returns_strings(self):
        elem, spec = parseAtom('O3')
        assert isinstance(elem, str)
        assert isinstance(spec, str)


# ---------------------------------------------------------------------------
# D — ROOT_DIR
# ---------------------------------------------------------------------------

class TestRootDir:
    def test_is_existing_directory(self):
        assert os.path.isdir(ROOT_DIR)

    def test_points_to_pyneb_package(self):
        assert os.path.basename(ROOT_DIR) == 'pyneb'

    def test_contains_init_py(self):
        assert os.path.isfile(os.path.join(ROOT_DIR, '__init__.py'))
