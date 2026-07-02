"""
Session-scoped fixtures for PyNeb test suite.

Using session scope so each Atom/RecAtom is built only once per test run,
avoiding repeated ~1-second data-loading costs.
"""
import os

import matplotlib
matplotlib.use('Agg')

import pytest
import pyneb as pn


@pytest.fixture(scope="session")
def O3():
    return pn.Atom('O', 3)


@pytest.fixture(scope="session")
def N2():
    return pn.Atom('N', 2)


@pytest.fixture(scope="session")
def S2():
    return pn.Atom('S', 2)


@pytest.fixture(scope="session")
def O2():
    return pn.Atom('O', 2)


@pytest.fixture(scope="session")
def H1():
    return pn.RecAtom('H', 1)


@pytest.fixture(scope="session")
def He2():
    return pn.RecAtom('He', 2)


@pytest.fixture(scope="session")
def O3_grid():
    """
    Small O III emissivity grid, shared read-only by emisGrid and plotting
    tests. The default 100x100 grid is far too slow for a test suite.
    """
    return pn.EmisGrid('O', 3, n_tem=10, n_den=8,
                       tem_min=8000., tem_max=12000.,
                       den_min=1e2, den_max=1e4)


@pytest.fixture(scope="session")
def smc24_path():
    """
    Path to the smc24.dat observation file (lines_in_rows format).
    A copy lives in tests/ because sample_scripts is not part of the
    installed package, against which CI runs the tests.
    """
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), 'smc24.dat')
