"""
Tests for pn.EmisGrid and pn.getEmisGridDict.

Uses the small session-scoped O3_grid fixture (10x8 points); tests must not
mutate it. Files are written to tmp_path only.
"""
import os

import numpy as np
import pytest

import pyneb as pn

GRID_KWARGS = dict(n_tem=10, n_den=8, tem_min=8000., tem_max=12000.,
                   den_min=1e2, den_max=1e4)


class TestEmisGridAttributes:

    def test_atom(self, O3_grid):
        assert O3_grid.elem == 'O'
        assert O3_grid.spec == 3
        assert O3_grid.atom.atom == 'O3'

    def test_grid_shapes(self, O3_grid):
        n_lev = O3_grid.atom.NLevels
        assert O3_grid.emis_grid.shape == (n_lev, n_lev, 10, 8)
        assert O3_grid.tem2D.shape == (10, 8)
        assert O3_grid.den2D.shape == (10, 8)

    def test_tem_den_ranges(self, O3_grid):
        assert O3_grid.tem.size == 10
        assert O3_grid.den.size == 8
        assert O3_grid.tem.min() == pytest.approx(8000.)
        assert O3_grid.tem.max() == pytest.approx(12000.)
        assert O3_grid.den.min() == pytest.approx(1e2)
        assert O3_grid.den.max() == pytest.approx(1e4)


class TestGetGrid:

    def test_by_wave(self, O3_grid):
        g = O3_grid.getGrid(wave=5007)
        assert g.shape == (10, 8)
        assert np.all(np.isfinite(g))
        assert np.all(g > 0)
        i, j = O3_grid.atom.getTransition(5007)
        assert np.array_equal(g, O3_grid.emis_grid[i - 1, j - 1])

    def test_by_levels(self, O3_grid):
        i, j = O3_grid.atom.getTransition(5007)
        assert np.array_equal(O3_grid.getGrid(lev_i=i, lev_j=j),
                              O3_grid.getGrid(wave=5007))

    def test_ratio_5007_4959_constant(self, O3_grid):
        # pure ratio of transition probabilities: independent of (T, n)
        r = O3_grid.getGrid(to_eval='L(5007)/L(4959)')
        assert np.allclose(r, r.mean(), rtol=1e-4)
        assert r.mean() == pytest.approx(2.98, rel=0.02)

    def test_4363_5007_increases_with_T(self, O3_grid):
        r = O3_grid.getGrid(to_eval='L(4363)/L(5007)')
        # axis 0 is temperature
        assert np.all(np.diff(r[:, 0]) > 0)


class TestSaveRestore:

    def test_roundtrip(self, O3_grid, tmp_path):
        file_ = str(tmp_path / 'O3.pypic')
        O3_grid.save(file_)
        g2 = pn.EmisGrid('O', 3, restore_file=file_, **GRID_KWARGS)
        assert g2.compute_new_grid is False
        assert np.allclose(g2.emis_grid, O3_grid.emis_grid)
        assert np.allclose(g2.tem, O3_grid.tem)

    def test_mismatch_recomputes(self, O3_grid, tmp_path):
        file_ = str(tmp_path / 'O3.pypic')
        O3_grid.save(file_)
        kwargs = dict(GRID_KWARGS, n_tem=5)
        g2 = pn.EmisGrid('O', 3, restore_file=file_, **kwargs)
        assert g2.compute_new_grid is True
        assert g2.tem.size == 5


class TestGetEmisGridDict:

    def test_compute_and_restore(self, tmp_path):
        old_pypic_path = pn.config.pypic_path
        try:
            d = pn.getEmisGridDict(atom_list=['O3'], pypic_path=str(tmp_path), **GRID_KWARGS)
            assert 'O3' in d
            assert isinstance(d['O3'], pn.EmisGrid)
            pypics = [f for f in os.listdir(str(tmp_path)) if f.endswith('.pypic')]
            assert len(pypics) == 1
            # second call restores from the saved file instead of recomputing
            d2 = pn.getEmisGridDict(atom_list=['O3'], pypic_path=str(tmp_path), **GRID_KWARGS)
            assert d2['O3'].compute_new_grid is False
        finally:
            pn.config.pypic_path = old_pypic_path


class TestEmisGridPlots:

    @pytest.fixture(autouse=True)
    def close_figures(self):
        import matplotlib.pyplot as plt
        yield
        plt.close('all')

    def test_plotImage(self, O3_grid):
        import matplotlib.pyplot as plt
        O3_grid.plotImage(to_eval='L(4363)/L(5007)')
        assert len(plt.gcf().axes) > 0

    def test_plotContours(self, O3_grid):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        O3_grid.plotContours('L(4363)/L(5007)', ax=ax)
        assert len(ax.collections) > 0

    def test_plotLineRatio(self, O3_grid):
        import matplotlib.pyplot as plt
        O3_grid.plotLineRatio('L(4363)/L(5007)', par='den')
        assert len(plt.gcf().axes) > 0
