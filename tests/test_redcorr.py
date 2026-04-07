import unittest
import numpy as np
import pyneb as pn


class TestRedCorrInstantiation(unittest.TestCase):
    def test_default_law(self):
        RC = pn.RedCorr()
        self.assertEqual(RC.law, 'No correction')

    def test_custom_law_F99(self):
        RC = pn.RedCorr(E_BV=0.3, law='F99')
        self.assertEqual(RC.law, 'F99')

    def test_custom_law_CCM89(self):
        RC = pn.RedCorr(E_BV=0.3, law='CCM89')
        self.assertEqual(RC.law, 'CCM89')

    def test_ebv_set(self):
        RC = pn.RedCorr(E_BV=0.5)
        self.assertAlmostEqual(RC.E_BV, 0.5)

    def test_rv_default(self):
        RC = pn.RedCorr()
        self.assertAlmostEqual(RC.R_V, 3.1)

    def test_chbeta_init(self):
        RC = pn.RedCorr(cHbeta=0.3, law='F99')
        self.assertAlmostEqual(RC.cHbeta, 0.3)


class TestRedCorrGetLaws(unittest.TestCase):
    def test_getlaws_returns_keys(self):
        RC = pn.RedCorr()
        laws = list(RC.getLaws())
        self.assertIn('No correction', laws)
        self.assertIn('F99', laws)
        self.assertIn('CCM89', laws)

    def test_getlaws_count(self):
        RC = pn.RedCorr()
        laws = list(RC.getLaws())
        self.assertGreater(len(laws), 5)


class TestRedCorrGetCorr(unittest.TestCase):
    def test_no_correction_returns_one(self):
        RC = pn.RedCorr(E_BV=0.5, law='No correction')
        corr = RC.getCorr(5007.)
        self.assertAlmostEqual(float(corr), 1.0, places=6)

    def test_F99_positive_ebv_gt_one(self):
        RC = pn.RedCorr(E_BV=0.1, law='F99')
        corr = RC.getCorr(5007.)
        self.assertGreater(corr, 1.0)

    def test_correction_increases_for_bluer_wave(self):
        # Extinction is stronger at shorter wavelengths
        RC = pn.RedCorr(E_BV=0.3, law='F99')
        corr_blue = RC.getCorr(4000.)
        corr_red = RC.getCorr(7000.)
        self.assertGreater(corr_blue, corr_red)

    def test_array_input(self):
        RC = pn.RedCorr(E_BV=0.1, law='CCM89')
        waves = np.array([4000., 5007., 6563.])
        corr = RC.getCorr(waves)
        self.assertEqual(corr.shape, (3,))
        self.assertTrue(np.all(corr > 0.))

    def test_zero_ebv_returns_one(self):
        RC = pn.RedCorr(E_BV=0., law='F99')
        corr = RC.getCorr(5007.)
        self.assertAlmostEqual(float(corr), 1.0, places=6)


class TestRedCorrGetCorrHb(unittest.TestCase):
    def test_getCorrHb_at_hbeta_is_one(self):
        for law in ['CCM89', 'F99', 'S79 H83 CCM89']:
            with self.subTest(law=law):
                RC = pn.RedCorr(E_BV=0.5, law=law)
                corr = RC.getCorrHb(4861.)
                np.testing.assert_allclose(corr, 1.0, atol=1e-2)

    def test_getCorrHb_bluer_gt_one(self):
        RC = pn.RedCorr(E_BV=0.3, law='F99')
        corr = RC.getCorrHb(4000.)
        self.assertGreater(corr, 1.0)

    def test_getCorrHb_redder_lt_one(self):
        RC = pn.RedCorr(E_BV=0.3, law='F99')
        corr = RC.getCorrHb(7000.)
        self.assertLess(corr, 1.0)


class TestRedCorrGetErrCorr(unittest.TestCase):
    def test_errcorr_positive(self):
        RC = pn.RedCorr(E_BV=0.3, law='F99')
        err = RC.getErrCorr(5007., 0.01)
        self.assertGreater(err, 0.)

    def test_errcorr_scales_with_ebv_error(self):
        RC = pn.RedCorr(E_BV=0.3, law='F99')
        err1 = RC.getErrCorr(5007., 0.01)
        err2 = RC.getErrCorr(5007., 0.02)
        np.testing.assert_allclose(err2, 2. * err1, rtol=1e-3)


class TestRedCorrSetCorr(unittest.TestCase):
    def test_setCorr_from_Balmer_decrement(self):
        RC = pn.RedCorr(law='F99')
        # Theoretical Ha/Hb = 2.85
        RC.setCorr(obs_over_theo=2.85 / 2.85, wave1=6563., wave2=4861.)
        # No reddening: E(B-V) should be close to 0
        np.testing.assert_allclose(RC.E_BV, 0., atol=0.05)

    def test_setCorr_reddened(self):
        RC = pn.RedCorr(law='F99')
        # Observed Ha/Hb > theoretical: positive reddening
        RC.setCorr(obs_over_theo=3.5 / 2.85, wave1=6563., wave2=4861.)
        self.assertGreater(RC.E_BV, 0.)


class TestRedCorrConversions(unittest.TestCase):
    def test_chbeta_from_ebv_roundtrip(self):
        RC = pn.RedCorr(law='F99')
        ebv = 0.4
        chb = RC.cHbetaFromEbv(ebv)
        recovered = RC.EbvFromCHbeta(chb)
        np.testing.assert_allclose(recovered, ebv, rtol=1e-6)

    def test_ebv_from_chbeta_roundtrip(self):
        RC = pn.RedCorr(law='CCM89')
        chb = 0.3
        ebv = RC.EbvFromCHbeta(chb)
        recovered = RC.cHbetaFromEbv(ebv)
        np.testing.assert_allclose(recovered, chb, rtol=1e-6)

    def test_AV_from_EBV(self):
        RC = pn.RedCorr(E_BV=1.0)
        # AV = R_V * E_BV = 3.1 * 1.0
        np.testing.assert_allclose(RC.AV, 3.1, rtol=1e-6)

    def test_EBV_from_AV(self):
        RC = pn.RedCorr()
        RC.AV = 3.1
        np.testing.assert_allclose(RC.E_BV, 1.0, rtol=1e-6)


class TestRedCorrPropertySetters(unittest.TestCase):
    def test_set_law_changes_correction(self):
        RC = pn.RedCorr(E_BV=0.5)
        RC.law = 'CCM89'
        corr_ccm = RC.getCorr(5007.)
        RC.law = 'F99'
        corr_f99 = RC.getCorr(5007.)
        # Different laws give different corrections
        self.assertNotAlmostEqual(corr_ccm, corr_f99, places=3)

    def test_set_ebv_changes_correction(self):
        RC = pn.RedCorr(E_BV=0., law='F99')
        corr0 = RC.getCorr(5007.)
        RC.E_BV = 0.5
        corr1 = RC.getCorr(5007.)
        self.assertGreater(corr1, corr0)


if __name__ == '__main__':
    unittest.main()
