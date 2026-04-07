import unittest
import numpy as np
import pyneb as pn
from pyneb.utils.misc import (int_to_roman, roman_to_int, parseAtom,
                               strExtract, quiet_divide)
from pyneb.utils.physics import vactoair, airtovac


class TestIntToRoman(unittest.TestCase):
    _cases = [
        (1, 'I'),
        (4, 'IV'),
        (9, 'IX'),
        (14, 'XIV'),
        (40, 'XL'),
        (90, 'XC'),
        (400, 'CD'),
        (900, 'CM'),
        (1999, 'MCMXCIX'),
        (2000, 'MM'),
        (3999, 'MMMCMXCIX'),
    ]

    def test_known_conversions(self):
        for n, expected in self._cases:
            with self.subTest(n=n):
                self.assertEqual(int_to_roman(n), expected)


class TestRomanToInt(unittest.TestCase):
    def test_known_conversions(self):
        cases = [('I', 1), ('IV', 4), ('IX', 9), ('XIV', 14),
                 ('XL', 40), ('XC', 90), ('CD', 400), ('CM', 900),
                 ('MCMXCIX', 1999), ('MM', 2000)]
        for roman, expected in cases:
            with self.subTest(roman=roman):
                self.assertEqual(roman_to_int(roman), expected)

    def test_roundtrip_range(self):
        for n in range(1, 100):
            with self.subTest(n=n):
                self.assertEqual(roman_to_int(int_to_roman(n)), n)

    def test_non_string_returns_none(self):
        self.assertIsNone(roman_to_int(4))

    def test_invalid_roman_returns_none(self):
        self.assertIsNone(roman_to_int('IIII'))

    def test_lowercase_accepted(self):
        # roman_to_int upper-cases input
        self.assertEqual(roman_to_int('iv'), 4)


class TestParseAtom(unittest.TestCase):
    def test_parse_O3(self):
        elem, spec = parseAtom('O3')
        self.assertEqual(elem, 'O')
        self.assertEqual(spec, '3')

    def test_parse_N2(self):
        elem, spec = parseAtom('N2')
        self.assertEqual(elem, 'N')
        self.assertEqual(spec, '2')

    def test_parse_He1(self):
        elem, spec = parseAtom('He1')
        self.assertEqual(elem, 'He')
        self.assertEqual(spec, '1')

    def test_parse_S2(self):
        elem, spec = parseAtom('S2')
        self.assertEqual(elem, 'S')
        self.assertEqual(spec, '2')

    def test_parse_lowercase_normalised(self):
        elem, spec = parseAtom('o3')
        self.assertEqual(elem, 'O')

    def test_parse_preserves_spec(self):
        _, spec = parseAtom('O3')
        self.assertEqual(spec, '3')


class TestStrExtract(unittest.TestCase):
    def test_extract_between_strings(self):
        result = strExtract('test123', 'e', '1')
        self.assertEqual(result, 'st')

    def test_extract_after_string(self):
        result = strExtract('test123', 'st')
        self.assertEqual(result, '123')

    def test_extract_with_int_start(self):
        result = strExtract('test123', 2)
        self.assertEqual(result, 'st123')

    def test_extract_with_int_stop(self):
        result = strExtract('test123', 0, 4)
        self.assertEqual(result, 'test')

    def test_extract_full_string_no_params(self):
        result = strExtract('test123')
        self.assertEqual(result, 'test123')


class TestQuietDivide(unittest.TestCase):
    def test_normal_division(self):
        result = quiet_divide(np.array([6., 9.]), np.array([2., 3.]))
        np.testing.assert_allclose(result, [3., 3.])

    def test_division_by_zero_gives_inf(self):
        result = quiet_divide(np.array([1., 2.]), np.array([1., 0.]))
        self.assertTrue(np.isinf(result[1]))

    def test_zero_divided_by_zero_gives_nan(self):
        result = quiet_divide(np.array([0.]), np.array([0.]))
        self.assertTrue(np.isnan(result[0]))

    def test_no_exception_raised(self):
        try:
            quiet_divide(np.array([1., 0.]), np.array([0., 0.]))
        except ZeroDivisionError:
            self.fail('quiet_divide raised ZeroDivisionError')


class TestVactoair(unittest.TestCase):
    def test_vactoair_optical_range(self):
        wave_vac = 5007.
        wave_air = vactoair(wave_vac)
        self.assertLess(wave_air, wave_vac)

    def test_vactoair_below_2000_unchanged(self):
        # Below pn.config.vactoair_low_wl, no correction applied
        wave = 1500.
        result = vactoair(wave)
        self.assertAlmostEqual(float(result), wave, places=3)

    def test_vactoair_array(self):
        waves = np.array([4000., 5007., 6563.])
        result = vactoair(waves)
        self.assertEqual(result.shape, waves.shape)

    def test_vactoair_value(self):
        np.testing.assert_allclose(vactoair(5007.), 5005.60, rtol=1e-3)


class TestAirtovac(unittest.TestCase):
    def test_airtovac_optical_range(self):
        wave_air = 5007.
        wave_vac = airtovac(wave_air)
        self.assertGreater(wave_vac, wave_air)

    def test_airtovac_below_2000_unchanged(self):
        wave = 1500.
        result = airtovac(wave)
        self.assertAlmostEqual(float(result), wave, places=3)

    def test_airtovac_array(self):
        waves = np.array([4000., 5007., 6563.])
        result = airtovac(waves)
        self.assertEqual(result.shape, waves.shape)

    def test_roundtrip_vactoair_airtovac(self):
        waves = np.array([4000., 5007., 6563.])
        np.testing.assert_allclose(airtovac(vactoair(waves)), waves, rtol=1e-4)


if __name__ == '__main__':
    unittest.main()
