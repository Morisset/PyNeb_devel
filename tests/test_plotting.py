"""
Smoke tests for the plotting facilities, run on the non-interactive Agg
backend (set in conftest.py). Beyond "does not crash", each test checks
that something was actually drawn.
"""
import matplotlib.pyplot as plt
import numpy as np
import pytest

import pyneb as pn


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close('all')


class TestAtomPlots:

    def test_plotEmiss(self, O3):
        fig, ax = plt.subplots()
        O3.plotEmiss(tem_min=8000, tem_max=12000, ax=ax)
        assert len(ax.lines) > 0

    def test_plotGrotrian(self, O3):
        O3.plotGrotrian(tem=1e4, den=1e2)
        fig = plt.gcf()
        assert len(fig.axes) > 0
        ax = fig.axes[0]
        assert len(ax.lines) + len(ax.collections) + len(ax.texts) > 0


class TestRedCorrPlot:

    def test_plot(self):
        RC = pn.RedCorr(law='CCM89', E_BV=1.)
        fig, ax = plt.subplots()
        RC.plot(w_inf=3000., w_sup=9000., ax=ax)
        assert len(ax.lines) > 0


class TestDataPlot:

    def test_plotA(self):
        dp = pn.DataPlot('O', 3)
        dp.plotA()
        fig = plt.gcf()
        assert len(fig.axes) > 0


class TestDiagnosticsPlot:

    def test_plot(self, O3_grid):
        # build a corrected observation whose ratio is guaranteed to cross the grid
        ratio = O3_grid.getGrid(to_eval='L(4363)/L(5007)')[5, 4]
        obs = pn.Observation(corrected=True)
        obs.addLine(pn.EmissionLine(label='O3_5007A', obsIntens=[100.],
                                    obsError=[0.05], corrected=True))
        obs.addLine(pn.EmissionLine(label='O3_4363A', obsIntens=[100. * ratio],
                                    obsError=[0.05], corrected=True))
        obs.names = ['obj1']
        diags = pn.Diagnostics()
        diags.addDiag('[OIII] 4363/5007')
        fig, ax = plt.subplots()
        diags.plot({'O3': O3_grid}, obs, ax=ax)
        assert len(ax.collections) + len(ax.lines) > 0
