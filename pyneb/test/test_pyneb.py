import pyneb as pn
import numpy as np
import matplotlib.pyplot as plt

"""
Usage:
pytest
pytest test_lines.py::test_atom_data


"""


verbose = False

def rel_error(calculated, expected):
    return (abs((expected - calculated) * 1.0) / (expected * 1.0))

def test_plot_data(save_pdf = False, verbose=verbose):
    """
    This function try to make the dataplot (As and Omegas) for any of the available atoms
    """

    plt.ion()
    atom_list = pn.atomicData.getAllAtoms()
    atom_list.remove('3He2')
    for atom in atom_list:
        if verbose:
            print('Doing {0}'.format(atom))
        dp = pn.DataPlot(atom=atom, OmegaInterp='Linear')
        
        dp.atom_n_max = np.min((dp.atom_n_max, 7))
        if len(dp.atom_data) > 1:
            if verbose:
                print('Plotting {0}'.format(atom))
            dp.plotAllA(save = save_pdf)
            plt.close()
          
        dp.coll_n_max = np.min((dp.coll_n_max, 7))
        if len(dp.coll_data) > 1:
            if verbose:
                print('Plotting {0}'.format(atom))
            dp.plotOmega(save = save_pdf)
            plt.close()
        
def test_diag_list():
    
    diags = ['[ArIII] 5192/7136', 
             '[ArIII] 9.0m/21.8m',
             '[ArIV] 4740/4711',
             '[CIII] 1909/1907',
             '[ClIII] 5538/5518',
             '[NII] 121m/20.5m',
             '[NII] 5755/6584',
             '[NeIII] 15.6m/36.0m',
             '[NeIII] 3343/3930+',
             '[NeV] 14.3m/24.2m',
             '[NeV] 2973/3370+',
             '[OI] 5577/6300+',
             '[OI] 63m/147m',
             '[OII] 3726/3729',
             '[OII] 3727+/7325+',
             '[OIII] 1664+/5007',
             '[OIII] 5007/88m', 
             '[OIII] 4363/5007',
             '[OIV] 1400+/25.9m',
             '[OIV] 1401/1405',
             '[SII] 4072+/6720+',
             '[SII] 4069/4076', 
             '[SII] 6731/6716',
             '[SIII] 6312/9069',
             '[SIII] 18.7m/33.5m'
             ] 
    # Build the 25 subplots
    f, axes = plt.subplots(5, 5, figsize=(20,20))
    # loop on the diagnostics
    for i, d in enumerate(diags):
        # extract from the diagnostic dictionnary 
        # the expression to be evaluated
        to_eval = pn.diags_dict[d][1]
        atom = pn.diags_dict[d][0]
        # split e.g. obtain 'O', 3 from 'O3'
        elem, spec = pn.utils.misc.parseAtom(atom)
        # instantiate the EmisGrid object for the given ion 
        EG = pn.EmisGrid(elem, spec)
        # select the axis in which to plot
        ax = axes.ravel()[i]
        # make the plot
        EG.plotContours(to_eval = to_eval, ax=ax)
   
    plt.close()

def test_line_labels():
    """
    This tool scans all the lines from pn.LINE_LABEL_LIST, extract the wavelength and compare with the 
    wavelength given by the corresponding atom for the corresponding transition. Is the difference is greater 
    than 1AA for optical or 0.1mu for IR, it prints the information
    """
    for atom_str in pn.LINE_LABEL_LIST:
        try:
            if atom_str[0] == '3':
                atom = None
            elif atom_str[-1] == 'r':
                atom = pn.RecAtom(atom=atom_str)
            else:
                atom = pn.Atom(atom=atom_str)
        except:
            atom = None
        if atom is not None:
            if atom.NLevels != 0:
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
                    if dif > 1.0:
                        print(atom, line)
                        assert dif < 1.0

def test_atom_instanciation():
    """
    This tool try to instantiate all the atoms using all the atomic data.
    """
    atoms = pn.atomicData.getAllAtoms()
    atoms.remove('3He2')
    for atom in atoms:
        for atom_file in pn.atomicData.getAllAvailableFiles(atom, 'atom'):
            for coll_file in pn.atomicData.getAllAvailableFiles(atom, 'coll'):
                #print atom, atom_file, coll_file
                try:
                    pn.atomicData.setDataFile(atom_file)
                    pn.atomicData.setDataFile(coll_file)
                    atm = pn.Atom(atom=atom)
                    assert atm.atom == atom
                except:
                    assert None == atom
        for rec_file in pn.atomicData.getAllAvailableFiles(atom, 'rec'):
                try:
                    pn.atomicData.setDataFile(rec_file)
                    atm = pn.Atom(atom=atom)
                    assert atm.atom == atom
                except:
                    assert None == atom
            
def test_atom_props():
    
    O3 = pn.Atom("O", 3)
    
    assert O3.elem == "O"
    assert O3.spec == 3
    assert O3.atom == "O3"
    assert O3.getTransition(5007) == (4, 3)
    
def test_atom_emissivity_values(error=5e-2):
    
    pn.atomicData.setDataFile('o_iii_atom_SZ00-WFD96.dat')
    pn.atomicData.setDataFile('o_iii_coll_SSB14.dat')
    atom = pn.Atom("O", 3)

    result = atom.getEmissivity(20000, 1e2)
    emi = [9.55e-22, 4.00e-28, 7.77e-22, 1.42e-24, 3.98e-21, 1.19e-20, 9.49e-23, 
        2.78e-25, 4.01e-22]
    a = []
    for i in result:
        for j in i:
            if (j > 0) :
                a.append(j)

    for aux in range(len(a) - 1):
        result = (a[aux] / emi[aux])
        assert (rel_error(a[aux], emi[aux]) < error)
    
def test_atom_emissivity_shape():
    pn.atomicData.setDataFile('o_iii_atom_SZ00-WFD96.dat')
    pn.atomicData.setDataFile('o_iii_coll_AK99.dat')
    atom = pn.Atom("O", 3)
    result = atom.getEmissivity(20000, 1e2)
    assert (result.shape == (6, 6))
    result = atom.getEmissivity(20000, [1e2, 1e3])
    assert (result.shape == (6, 6, 2))
    result = atom.getEmissivity([20000, 10000], 1e2)
    assert (result.shape == (6, 6, 2))
    result = atom.getEmissivity([20000, 10000], [1e2, 1e3])        
    assert (result.shape == (6, 6, 2, 2))
    result = atom.getEmissivity([20000, 10000], [1e2, 1e3], product = True)        
    assert (result.shape == (6, 6, 2, 2))
    result = atom.getEmissivity([20000, 10000], [1e2, 1e3], product = False)        
    assert (result is None)
    
    result = atom.getEmissivity(20000, 1e2, wave = 5007)
    assert (result.shape == ())
    result = atom.getEmissivity(20000, [1e2, 1e3], wave = 5007)
    assert (result.shape == (2,))
    result = atom.getEmissivity([20000, 10000], 1e2, wave = 5007)
    assert (result.shape == (2,)) # the bug I'm looking for
    result = atom.getEmissivity([20000, 10000], [1e2, 1e3], wave = 5007)        
    assert (result.shape == (2, 2))
    result = atom.getEmissivity([20000, 10000], [1e2, 1e3], product = True, wave = 5007)        
    assert (result.shape == (2, 2))
    
def test_HbEmissivity(max_error=5e-4):
    H1 = pn.RecAtom('H', 1)
    tem = np.array([.5, 1, 1.5, 2., 3.]) * 1e4
    den = 1e2
    Hb_expected = np.array([2.2e-25, 1.235e-25, 8.6e-26, 6.58e-26, 4.440e-26])
    Hb_computed = H1.getEmissivity(tem, den, label='4_2')
    for i in range(len(tem)):
        assert(rel_error(H1.getEmissivity(tem[i], den, label='4_2'), Hb_expected[i]) < max_error)
        assert(rel_error((Hb_computed[i]), Hb_expected[i]) < max_error)
    
def test_ion_abundance(max_error=5e-3):
    pn.atomicData.setDataFile('o_iii_atom_SZ00-WFD96.dat')
    pn.atomicData.setDataFile('o_iii_coll_AK99.dat')
    atom = pn.Atom("O", 3)
    assert(rel_error(atom.getIonAbundance(10000, 10000, 100, 4, 3), 3.62e-3) < max_error)
    assert(rel_error(atom.getIonAbundance(10000, 20000, 100, 4, 3), 5.662e-4) < max_error)
    assert(rel_error(atom.getIonAbundance(10000, 10000, 100, wave=5007), 3.62e-3) < max_error)
    assert(rel_error(atom.getIonAbundance(10000, 20000, 100, wave=5007), 5.662e-4) < max_error)
    pn.atomicData.setDataFile('s_ii_atom_PKW09.dat')
    pn.atomicData.setDataFile('s_ii_coll_TZ10.dat')
    atom = pn.Atom("S", 2)
    assert(rel_error(atom.getIonAbundance(200, 20000, 100, 3, 1), 2.295e-6) < max_error)
    pn.atomicData.setDataFile('o_ii_atom_Z82-WFD96.dat')
    pn.atomicData.setDataFile('o_ii_coll_P06-T07.dat')
    atom = pn.Atom("O", 2)
    assert(rel_error(atom.getIonAbundance(200, 20000, 1000, 3, 1), 1.548e-5) < max_error)

def test_getA(max_error=5e-3):
    pn.atomicData.setDataFile('o_iii_atom_SZ00-WFD96.dat')
    pn.atomicData.setDataFile('o_iii_coll_AK99.dat')
    atom = pn.Atom("O", 3)
    assert(rel_error(atom.getA(1, 1) + 1, 1.0) < max_error)
    assert(rel_error(atom.getA(2, 1), 2.62e-05) < max_error)
    assert(rel_error(atom.getA(3, 1), 3.17e-11) < max_error)
    assert(rel_error(atom.getA(4, 1), 2.41e-06) < max_error)
    assert(rel_error(atom.getA(3, 2), 9.76e-05) < max_error)
    assert(rel_error(atom.getA(4, 2), 0.0068) < max_error)
    assert(rel_error(atom.getA(5, 2), 0.215) < max_error)

def test_getOmega(max_error=5e-3):
    pn.atomicData.setDataFile('o_iii_atom_SZ00-WFD96.dat')
    pn.atomicData.setDataFile('o_iii_coll_AK99.dat')
    atom = pn.Atom("O", 3)
    assert(rel_error(atom.getOmega(10000, 2, 1), 0.522) < max_error)
    assert(rel_error(atom.getOmega(10000, 3, 1), 0.2573) < max_error)
    assert(rel_error(atom.getOmega(10000, 4, 1), 0.2434) < max_error)
    assert(rel_error(atom.getOmega(10000, 5, 1), 0.0321) < max_error)
    assert(rel_error(atom.getOmega(10000, 3, 2), 1.232) < max_error)
    assert(rel_error(atom.getOmega(10000, 4, 2), 0.73) < max_error)
    assert(rel_error(atom.getOmega(10000, 5, 2), 0.0962) < max_error)
    assert(rel_error(atom.getOmega(10000, wave=5007), 1.217) < max_error)
    
def test_getTemDen(max_error=5e-3):
    pn.atomicData.setDataFile('o_iii_atom_SZ00-WFD96.dat')
    pn.atomicData.setDataFile('o_iii_coll_AK99.dat')
    atom = pn.Atom('O', 3)
    assert(rel_error(atom.getTemDen(100., den=10000., wave1=5007, wave2=4363), 11308.) < max_error)
    assert(rel_error(atom.getTemDen(1., tem=10000., wave1=88e4, wave2=52e4), 183) < max_error)
    pn.atomicData.setDataFile('o_ii_atom_Z82-WFD96.dat')
    pn.atomicData.setDataFile('o_ii_coll_P06-T07.dat')
    atom = pn.Atom('O', 2)
    assert(rel_error(atom.getTemDen(1.5, tem=10000., wave1=3726, wave2=3729), 1496) < max_error)

                   