import pytest 
import pyneb as pn

def test_Atom_atom():
    O3 = pn.Atom('O',3)
    assert O3.atom == 'O3'

def test_Atom_elem():
    O3 = pn.Atom('O',3)
    assert O3.elem == 'O'

def test_RecAtom_atom():
    H1 = pn.RecAtom('H',1)
    assert H1.atom == 'H1'

