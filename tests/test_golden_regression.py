"""
Numerical regression tests against golden reference values.

The goldens (tests/golden/golden_refs.npz) were produced by
tests/generate_golden.py on a known-good revision. Every test recomputes
the same quantity with the current code and compares to the stored value.

If these tests fail because the default atomic data files changed on
purpose, regenerate the goldens:

    python tests/generate_golden.py
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from generate_golden import (GOLDEN_PATH, LOOSE_PREFIXES, ZERO_PATTERN_PREFIXES,
                             compute_all, get_datafiles, pin_datafiles)

REGEN_MSG = ('Golden reference mismatch: the atomic data files in use differ from '
             'the ones recorded in golden_refs.npz. If this is intentional, '
             'regenerate the goldens with: python tests/generate_golden.py')

if os.path.exists(GOLDEN_PATH):
    _golden = np.load(GOLDEN_PATH)
    GOLDEN_KEYS = sorted(k for k in _golden.files if not k.startswith('__'))
else:
    _golden = None
    GOLDEN_KEYS = []


@pytest.fixture(scope='module')
def computed():
    if _golden is None:
        pytest.fail('tests/golden/golden_refs.npz is missing. '
                    'Generate it with: python tests/generate_golden.py')
    recorded = {k[len('__datafiles_'):]: str(_golden[k])
                for k in _golden.files if k.startswith('__datafiles_')}
    pin_datafiles(recorded)
    if get_datafiles() != recorded:
        pytest.fail(REGEN_MSG)
    return compute_all()


def test_golden_file_exists():
    assert os.path.exists(GOLDEN_PATH), (
        'tests/golden/golden_refs.npz is missing. '
        'Generate it with: python tests/generate_golden.py')


def test_key_sets_match(computed):
    assert sorted(computed) == GOLDEN_KEYS


@pytest.mark.parametrize('key', GOLDEN_KEYS)
def test_golden(key, computed):
    new = computed[key]
    ref = _golden[key]
    assert new.shape == ref.shape, f'{key}: shape {new.shape} != {ref.shape}'
    rtol = 1e-5 if any(p in key for p in LOOSE_PREFIXES) else 1e-10
    np.testing.assert_allclose(new, ref, rtol=rtol, atol=0., err_msg=key)
    if any(p in key for p in ZERO_PATTERN_PREFIXES):
        np.testing.assert_array_equal(new == 0., ref == 0.,
                                      err_msg=f'{key}: structural zero pattern changed')
