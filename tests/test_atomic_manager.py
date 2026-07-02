"""
Tests for the pn.atomicData manager (_ManageAtomicData singleton).

The manager holds global state: the only mutating test restores the
original data file in a finally block.
"""
import pytest

import pyneb as pn


class TestReadOnly:

    def test_getDataFile_with_type(self):
        file_ = pn.atomicData.getDataFile('O3', 'atom')
        assert 'o_iii' in file_
        assert file_.endswith('.dat') or file_.endswith('.fits') or file_.endswith('.chianti')

    def test_getDataFile_all_types(self):
        files = pn.atomicData.getDataFile('O3')
        # (atom, coll, rec, trc)
        assert any(f is not None and 'o_iii_atom' in f for f in files)
        assert any(f is not None and 'o_iii_coll' in f for f in files)

    def test_getAllAvailableFiles(self):
        files = pn.atomicData.getAllAvailableFiles('O3')
        assert len(files) > 0
        assert all('o_iii' in f for f in files)
        # currently used files are marked with '*'
        current_atom = pn.atomicData.getDataFile('O3', 'atom')
        assert '* ' + current_atom in files

    def test_getAllAtoms(self):
        atoms = pn.atomicData.getAllAtoms()
        for atom in ('O3', 'N2', 'S2'):
            assert atom in atoms

    def test_getAllAtoms_rec(self):
        rec_atoms = pn.atomicData.getAllAtoms(coll=False, rec=True)
        assert 'H1' in rec_atoms

    def test_getPredefinedDataFileDict(self):
        all_dicts = pn.atomicData.getAllPredefinedDict()
        assert 'PYNEB_23_01' in all_dicts
        predef = pn.atomicData.getPredefinedDataFileDict('PYNEB_23_01')
        assert 'atom' in predef['O3']
        assert 'coll' in predef['O3']


class TestSetDataFile:

    def test_setDataFile_roundtrip(self):
        orig = pn.atomicData.getDataFile('O3', 'coll')
        alts = [f.strip('* ') for f in
                pn.atomicData.getAllAvailableFiles('O3', data_type='coll')]
        # avoid chianti files, which require the XUVTOP environment variable
        alt = next(f for f in alts if f != orig and f.endswith('.dat'))
        try:
            pn.atomicData.setDataFile(alt)
            assert pn.atomicData.getDataFile('O3', 'coll') == alt
            O3_alt = pn.Atom('O', 3)
            assert O3_alt.collFitsFile == alt
        finally:
            pn.atomicData.setDataFile(orig)
        assert pn.atomicData.getDataFile('O3', 'coll') == orig
