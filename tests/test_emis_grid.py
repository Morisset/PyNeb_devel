import unittest
import numpy as np
import pyneb as pn


class TestEmisGridInstantiation(unittest.TestCase):
    def test_basic_creation(self):
        grid = pn.EmisGrid('O', 3, n_tem=10, n_den=10,
                           tem_min=5000., tem_max=20000.)
        with self.subTest():
            self.assertEqual(grid.elem, 'O')
        with self.subTest():
            self.assertEqual(grid.spec, 3)

    def test_tem_array_shape(self):
        n = 10
        grid = pn.EmisGrid('O', 3, n_tem=n, n_den=10)
        self.assertEqual(grid.tem.shape, (n,))

    def test_den_array_shape(self):
        n = 8
        grid = pn.EmisGrid('O', 3, n_tem=10, n_den=n)
        self.assertEqual(grid.den.shape, (n,))

    def test_tem_range(self):
        grid = pn.EmisGrid('O', 3, n_tem=10, n_den=10,
                           tem_min=5000., tem_max=20000.)
        np.testing.assert_allclose(grid.tem.min(), 5000., rtol=1e-3)
        np.testing.assert_allclose(grid.tem.max(), 20000., rtol=1e-3)

    def test_den_range(self):
        grid = pn.EmisGrid('O', 3, n_tem=10, n_den=10,
                           den_min=10., den_max=1e6)
        np.testing.assert_allclose(grid.den.min(), 10., rtol=1e-3)
        np.testing.assert_allclose(grid.den.max(), 1e6, rtol=1e-3)

    def test_emis_grid_shape(self):
        nt, nd = 10, 8
        O3 = pn.Atom('O', 3)
        grid = pn.EmisGrid('O', 3, n_tem=nt, n_den=nd)
        nl = O3.NLevels
        self.assertEqual(grid.emis_grid.shape, (nl, nl, nt, nd))

    def test_creation_from_atom_object(self):
        O3 = pn.Atom('O', 3)
        grid = pn.EmisGrid(n_tem=10, n_den=10, atomObj=O3)
        self.assertEqual(grid.elem, 'O')
        self.assertEqual(grid.spec, 3)


class TestEmisGridGetGrid(unittest.TestCase):
    def setUp(self):
        self.nt = 10
        self.nd = 8
        self.grid = pn.EmisGrid('O', 3, n_tem=self.nt, n_den=self.nd)

    def test_getGrid_by_levels(self):
        g = self.grid.getGrid(lev_i=4, lev_j=3)
        self.assertEqual(g.shape, (self.nt, self.nd))

    def test_getGrid_by_wave(self):
        g = self.grid.getGrid(wave=5007)
        self.assertEqual(g.shape, (self.nt, self.nd))

    def test_getGrid_levels_equals_wave(self):
        g_lev = self.grid.getGrid(lev_i=4, lev_j=3)
        g_wave = self.grid.getGrid(wave=5007)
        np.testing.assert_array_equal(g_lev, g_wave)

    def test_getGrid_nonnegative(self):
        g = self.grid.getGrid(lev_i=4, lev_j=3)
        self.assertTrue(np.all(g >= 0.))

    def test_getGrid_to_eval_ratio(self):
        g = self.grid.getGrid(to_eval='L(5007)/L(4959)')
        self.assertEqual(g.shape, (self.nt, self.nd))
        # Ratio should be positive and finite where both lines are defined
        self.assertTrue(np.all(np.isfinite(g)))
        self.assertTrue(np.all(g > 0.))

    def test_getGrid_to_eval_sum(self):
        g_sum = self.grid.getGrid(to_eval='L(5007) + L(4959)')
        g_5007 = self.grid.getGrid(wave=5007)
        g_4959 = self.grid.getGrid(wave=4959)
        np.testing.assert_allclose(g_sum, g_5007 + g_4959, rtol=1e-10)


if __name__ == '__main__':
    unittest.main()
