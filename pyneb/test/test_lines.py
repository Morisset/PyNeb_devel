import pyneb as pn
import numpy as np

def test_line_lists():
    """
    Check the validity of all the lines in diags_dict by comparing the wavelengths
    """
    def L(wave):
        if wave > 10000:
            return 1
        if pn.isValid('{0}_{1}A'.format(atom,wave)):
            return 1
        else:
            print('Wrong: {0}_{1}A'.format(atom,wave))
            return 1
        
    def I(i, j):
        return 1
    
    def B(wave):
        if pn.isValid('{0}_{1}'.format(atom,wave)):
            return 1
        else:
            print('Wrong: {0}_{1}'.format(atom,wave))
            return 1
    
    for diag in pn.diags_dict:
        atom = pn.diags_dict[diag][0]
        #print atom
        atom_obj = pn.Atom(atom = atom)
        eval(pn.diags_dict[diag][1])


def test_diags():
    """

    """
    pn.atomicData.setDataFile('o_iii_coll_AK99.dat')
    tem = 1e4 
    den = 1e2
    L = lambda wave: atom.getEmissivity(tem, den, wave = wave)
    I = lambda i, j: atom.getEmissivity(tem, den, i, j)
    def B(wave):
        if pn.isValid('{0}_{1}'.format(atom.atom,wave)):
            return 1
        else:
            print('!!! Wrong: {0}_{1}'.format(atom.atom,wave))
            return 1
        
    for diag in pn.diags_dict:
        print diag
        atom = pn.Atom(atom = pn.diags_dict[diag][0])
        to_eval = pn.diags_dict[diag][1]
        try:
            lines_rat = eval(to_eval)
        except:
            print('{} failed'.format(to_eval))
            
def test_labels():
    """
    This tool scans all the lines from pn.LINE_LABEL_LIST, extract the wavelength and compare with the 
    wavelength given by the corresponding atom for the corresponding transition. Is the difference is greater 
    than 1AA for optical or 0.1mu for IR, it prints the information
    """
    for atom_str in pn.LINE_LABEL_LIST:
        try:
            if atom_str[-1] == 'r':
                atom = pn.RecAtom(atom=atom_str)
            else:
                atom = pn.Atom(atom=atom_str)
        except:
            atom = None
        if atom is not None:
            for line in pn.LINE_LABEL_LIST[atom_str]:
                if line[-1] == 'm':
                    wavelength = float(line[:-1])*1e4
                else:
                    wavelength = float(line[:-1])
                i, j = atom.getTransition(wavelength)
                wave = atom.wave_Ang[i-1, j-1]
                if wave < 1e4:
                    dif = np.abs(wavelength-wave)
                else:
                    dif = np.abs(wavelength-wave)/1e3
                #print atom_str
                if dif > 1.0:
                    print('{0} {1:.2f} =? {2:.2f}  {3:.2f}'.format(atom_str, wavelength, wave, dif))
                    

def test_atom_data():
    """
    This tool try to instantiate all the atoms using all the atomic data.
    """
    atoms = pn.atomicData.getAllAtoms()
    i = 0
    for atom in atoms:
        for atom_file in pn.atomicData.getAllAvailableFiles(atom, 'atom'):
            for coll_file in pn.atomicData.getAllAvailableFiles(atom, 'coll'):
                print atom, atom_file, coll_file
                pn.atomicData.setDataFile(atom_file)
                pn.atomicData.setDataFile(coll_file)
                atm = pn.Atom(atom=atom)
                print atm
                i += 1
    print('Just did {} atoms for the {} ions'.format(i, len(atoms)))
    # Just did 263 atoms for the 55 ions

def compare_line_list():
    OLD_LINE_LABEL_LIST = pn.utils.init.OLD_LINE_LABEL_LIST
    LINE_LABEL_LIST = pn.LINE_LABEL_LIST
    
    for atom in np.sort(LINE_LABEL_LIST.keys()):
        if atom in OLD_LINE_LABEL_LIST:
            print atom + ' NEW',
            for line in LINE_LABEL_LIST[atom]:
                if line not in OLD_LINE_LABEL_LIST[atom]:
                    print line,
            print ''
            print atom + ' OLD',
            for line in OLD_LINE_LABEL_LIST[atom]:
                if line not in LINE_LABEL_LIST[atom]:
                    print line,
            print ''
        else:
            print atom + ' NEW ATOM'
            print LINE_LABEL_LIST[atom]
        print '==================='
        
        
                    