"""
this file contains the definitions of the functions used to answer the questions of
exercise  2.1 of  PLASMA DIAGNOSTICS

to execute these functions, one first needs to
import a few libraries as listed at the beginning of the file run_ex.py
as well as this file.
the functions can then be called  as explained in the file run_ex.py

the functions defined here can be modified. 
In that case, is is necessary, before using the modified version, 
to do "reload(ex1_1)" from the terminal
sometimes this does not work well. 
it is then recommended to quit python (ctrl D) and enter again (ipython --pylab)
"""
import numpy as np
import matplotlib.pyplot as plt
import pyneb as pn
from pyneb.utils.misc import parseAtom

#pn.atomicData.setDataFile('cl_iii_atom_M83-KS86.fits')

def p1(ion):
        # split ion into elem and spec, e.g 'O3' into 'O' and 3
        elem, spec = parseAtom(ion)
        # instanciate the corresponding Atom object
        atom = pn.Atom(elem, spec)
        # print information including transition probabilities
        #atom.printIonic(printA = True)
        # prepare a new figure
        plt.figure()
        # plot energy levels
        atom.plotGrotrian()


def p2(diag):
    # get the ion, and diagnostic description from the dictionary:
    ion, diag_eval, err = pn.diags_dict[diag]
    # split ion into elem and spec, e.g 'O3' into 'O' and 3
    elem, spec = parseAtom(ion)
    # prepare a new figure
    plt.figure()
    # create a grid of emissivities 
    #NB: one can use a pypic file containing all the emissivities, if already made
    # in that case set restore_file to the name of the pypic file.
    grid = pn.EmisGrid(elem, spec, restore_file=None, OmegaInterp='Linear')
    # plot the contours
    grid.plotContours(to_eval=diag_eval, low_level=None, high_level=None, n_levels=20, 
                      linestyles='-', clabels=True, log_levels=True, 
                      title='{0} {1}'.format(ion, diag_eval))
    # save the plot into pdf files                  
    plt.savefig('{0}_{1}.pdf'.format(ion, diag_eval.replace('/', '_')))


# the following is to plot all the possible diagnostic ratiios available in pyneb    
def plot_all(save=False):
    pn.log_.level=1
    AA = pn.getAtomDict(OmegaInterp='Linear')
    for diag in pn.diags_dict:
        atom, diag_eval, err = pn.diags_dict[diag]
        if atom in AA:
            plt.figure()
            grid = pn.EmisGrid(atomObj=AA[atom])
            grid.plotContours(to_eval=diag_eval)
            if save:
                plt.savefig('{0}_{1}.pdf'.format(atom, diag_eval.replace('/', '_')))
        
        
