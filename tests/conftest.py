"""
Session-scoped fixtures for PyNeb test suite.

Using session scope so each Atom/RecAtom is built only once per test run,
avoiding repeated ~1-second data-loading costs.
"""
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
