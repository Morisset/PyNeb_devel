"""
Comprehensive tests for pyneb.Diagnostics.

Covers: constructor, addDiag, delDiag, getAllDiags, getDiagFromLabel,
        getUniqueAtoms, addAll, and getCrossTemDen (roundtrip + sanity).

Each test uses a function-scoped fixture so addDiag mutations don't leak
between tests.
"""
import numpy as np
import pytest
import pyneb as pn


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def diags():
    """Fresh Diagnostics object for each test (function scope)."""
    return pn.Diagnostics()


# ---------------------------------------------------------------------------
# A — Constructor
# ---------------------------------------------------------------------------

class TestDiagnosticsConstructor:
    def test_empty_after_construction(self, diags):
        assert list(diags.getDiagLabels()) == []

    def test_atomDict_empty_after_construction(self, diags):
        assert diags.atomDict == {}

    def test_getAllDiags_returns_large_dict(self, diags):
        all_diags = diags.getAllDiags()
        assert len(all_diags) >= 100

    def test_known_labels_present_in_getAllDiags(self, diags):
        all_diags = diags.getAllDiags()
        assert '[OIII] 4363/5007' in all_diags
        assert '[NII] 5755/6548' in all_diags
        assert '[SII] 6731/6716' in all_diags


# ---------------------------------------------------------------------------
# B — addDiag / delDiag
# ---------------------------------------------------------------------------

class TestAddDelDiag:
    def test_addDiag_single_label(self, diags):
        diags.addDiag('[NII] 5755/6548')
        assert '[NII] 5755/6548' in list(diags.getDiagLabels())
        assert len(list(diags.getDiagLabels())) == 1

    def test_addDiag_list_of_labels(self, diags):
        diags.addDiag(['[NII] 5755/6548', '[SII] 6731/6716'])
        labels = list(diags.getDiagLabels())
        assert '[NII] 5755/6548' in labels
        assert '[SII] 6731/6716' in labels
        assert len(labels) == 2

    def test_addDiag_by_atom(self, diags):
        diags.addDiag(atom='N2')
        labels = list(diags.getDiagLabels())
        assert len(labels) > 0
        # All returned labels should belong to N2
        for lbl in labels:
            info = diags.getDiagFromLabel(lbl)
            assert info[0] == 'N2'

    def test_addAll(self, diags):
        diags.addAll()
        assert len(list(diags.getDiagLabels())) >= 100

    def test_delDiag(self, diags):
        diags.addDiag('[NII] 5755/6548')
        diags.delDiag('[NII] 5755/6548')
        assert list(diags.getDiagLabels()) == []

    def test_getDiagFromLabel_structure(self, diags):
        diags.addDiag('[NII] 5755/6548')
        info = diags.getDiagFromLabel('[NII] 5755/6548')
        assert info is not None
        assert info[0] == 'N2'
        assert 'L(5755)' in info[1]
        assert len(info) == 3

    def test_getDiagFromLabel_missing_returns_None(self, diags):
        assert diags.getDiagFromLabel('nonexistent_label_xyz') is None


# ---------------------------------------------------------------------------
# C — getUniqueAtoms
# ---------------------------------------------------------------------------

class TestGetUniqueAtoms:
    def test_unique_atoms_N2_S2(self, diags):
        diags.addDiag('[NII] 5755/6548')
        diags.addDiag('[SII] 6731/6716')
        atoms = diags.getUniqueAtoms()
        assert set(atoms) == {'N2', 'S2'}

    def test_unique_atoms_single(self, diags):
        diags.addDiag('[OIII] 4363/5007')
        atoms = diags.getUniqueAtoms()
        assert 'O3' in atoms


# ---------------------------------------------------------------------------
# D — getCrossTemDen
# ---------------------------------------------------------------------------

class TestGetCrossTemDen:
    def test_O3_S2_roundtrip(self):
        """Compute synthetic ratios at known T,n and recover them."""
        O3 = pn.Atom('O', 3)
        S2 = pn.Atom('S', 2)
        true_tem = 1e4
        true_den = 500.

        r_O3 = (O3.getEmissivity(true_tem, true_den, wave=4363) /
                O3.getEmissivity(true_tem, true_den, wave=5007))
        r_S2 = (S2.getEmissivity(true_tem, true_den, wave=6731) /
                S2.getEmissivity(true_tem, true_den, wave=6716))

        diags = pn.Diagnostics()
        diags.addDiag(['[OIII] 4363/5007', '[SII] 6731/6716'])
        tem, den = diags.getCrossTemDen('[OIII] 4363/5007', '[SII] 6731/6716',
                                        r_O3, r_S2)
        assert float(tem) == pytest.approx(true_tem, rel=3e-3)
        assert float(den) == pytest.approx(true_den, rel=5e-3)

    def test_N2_S2_physical_range(self):
        """Representative observed ratios should yield physically sensible T, n."""
        diags = pn.Diagnostics()
        diags.addDiag(['[NII] 5755/6548', '[SII] 6731/6716'])
        tem, den = diags.getCrossTemDen('[NII] 5755/6548', '[SII] 6731/6716',
                                        0.05, 1.1)
        assert 5000 < float(tem) < 30000
        assert 10 < float(den) < 1e6

    def test_atomDict_populated_after_getCrossTemDen(self):
        diags = pn.Diagnostics()
        diags.addDiag(['[OIII] 4363/5007', '[SII] 6731/6716'])
        diags.getCrossTemDen('[OIII] 4363/5007', '[SII] 6731/6716', 0.006, 0.92)
        assert 'O3' in diags.atomDict
        assert 'S2' in diags.atomDict
