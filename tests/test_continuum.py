import unittest
import numpy as np
import pyneb as pn


WL = np.array([3500., 4000., 4500., 5000., 6000.])


class TestContinuumGff(unittest.TestCase):
    def setUp(self):
        self.cont = pn.Continuum()

    def test_gff_positive(self):
        g = self.cont.gff(1., 1e4, WL)
        self.assertTrue(np.all(g > 0.))

    def test_gff_array_shape(self):
        g = self.cont.gff(1., 1e4, WL)
        self.assertEqual(g.shape, WL.shape)

    def test_gff_nonzero(self):
        g = self.cont.gff(1., 1e4, WL)
        # Gaunt factor should be around 1 for optical wavelengths
        self.assertTrue(np.all(g > 0.5))
        self.assertTrue(np.all(g < 10.))

    def test_gff_scalar_wavelength(self):
        g = self.cont.gff(1., 1e4, np.array([5000.]))
        self.assertGreater(float(g), 0.)


class TestContinuumTwoPhoton(unittest.TestCase):
    def setUp(self):
        self.cont = pn.Continuum()

    def test_two_photon_nonnegative(self):
        tp = self.cont.two_photon(1e4, 100., WL)
        self.assertTrue(np.all(tp >= 0.))

    def test_two_photon_array_shape(self):
        tp = self.cont.two_photon(1e4, 100., WL)
        self.assertEqual(tp.shape, WL.shape)

    def test_two_photon_zero_below_lyman_alpha(self):
        # Two-photon emission is cut off for λ < Ly-α (1215.7 Å),
        # where y = 1215.7/λ > 1. The code sets A[mask]=0 for y > 1.
        wl_below = np.array([500., 900., 1200.])
        tp = self.cont.two_photon(1e4, 100., wl_below)
        np.testing.assert_array_equal(tp, 0.)

    def test_two_photon_decreases_with_density(self):
        # At very high density, two-photon is collisionally suppressed
        tp_low = self.cont.two_photon(1e4, 1., WL)
        tp_high = self.cont.two_photon(1e4, 1e8, WL)
        self.assertTrue(np.all(tp_low >= tp_high))


class TestContinuumFreeFree(unittest.TestCase):
    def setUp(self):
        self.cont = pn.Continuum()

    def test_freefree_positive(self):
        ff = self.cont.FreeFree(1e4, WL)
        self.assertTrue(np.all(ff > 0.))

    def test_freefree_array_shape(self):
        ff = self.cont.FreeFree(1e4, WL)
        self.assertEqual(ff.shape, WL.shape)

    def test_freefree_larger_with_helium(self):
        ff_no_he = self.cont.FreeFree(1e4, WL)
        ff_with_he = self.cont.FreeFree(1e4, WL, He1_H=0.1, He2_H=0.01)
        self.assertTrue(np.all(ff_with_he >= ff_no_he))


class TestContinuumGetContinuum(unittest.TestCase):
    def setUp(self):
        self.cont = pn.Continuum()

    def test_get_continuum_returns_array(self):
        gc = self.cont.get_continuum(1e4, 100., wl=WL)
        self.assertEqual(gc.shape, WL.shape)

    def test_get_continuum_positive(self):
        gc = self.cont.get_continuum(1e4, 100., wl=WL)
        self.assertTrue(np.all(gc >= 0.))

    def test_get_continuum_increases_with_helium(self):
        gc_no_he = self.cont.get_continuum(1e4, 100., wl=WL)
        gc_with_he = self.cont.get_continuum(1e4, 100., He1_H=0.1,
                                              He2_H=0.01, wl=WL)
        self.assertTrue(np.all(gc_with_he >= gc_no_he))

    def test_get_continuum_default_wavelengths(self):
        # Default wl covers the Balmer continuum region
        gc = self.cont.get_continuum(1e4, 100.)
        self.assertGreater(gc.size, 0)
        self.assertTrue(np.all(gc >= 0.))


if __name__ == '__main__':
    unittest.main()
