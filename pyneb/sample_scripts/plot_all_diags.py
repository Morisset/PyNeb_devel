# Plot the contour diagram of the selected diagnostics
# If SAVE, store them in files

import matplotlib.pyplot as plt
import pyneb as pn
from pyneb.core.diags import diags_dict

def plot_all(save=False):
    pn.log_.level=1
    AA = pn.getAtomDict(OmegaInterp='Linear')
    # Loop over all the diags stored in pn.core.diags.diags_dict
    for diag in diags_dict:
        atom, diag_eval, err = diags_dict[diag]
        # Skip Fe III as they are so many
        if (atom in AA) and (atom != 'Fe3'):
            print atom
            plt.figure()
            grid = pn.EmisGrid(atomObj=AA[atom])
            grid.plotContours(to_eval=diag_eval)
            if save:
                plt.savefig('{0}_{1}.pdf'.format(atom, diag_eval.replace('/', '_')))

if __name__ == '__main__':

    plot_all(False)
    plt.show()