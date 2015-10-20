# Analysis of a simple two-component model, meant to illustrate the bias arising from assuming
# that the region is homogeneous in density
# First, an emission region made up of two different subregions is modelled,
# each with a different mass and density. The resulting overall emissivity is computed
# Second, the region is analyzed as if it were a homogeneous region

import pyneb as pn
import matplotlib.pyplot as plt
from pyneb.utils.misc import parseAtom


def plot_2comp(tem1=1e4, tem2=1e4, dens1=3e2, dens2=5e5, mass1=1, mass2=5e-4):
        
    # List of diagnostics used to analyze the region
    diags = pn.Diagnostics()
    
    for diag in pn.diags_dict:
        if diag[0:7] != '[FeIII]':
            diags.addDiag(diag)
    diags.addClabel('[SIII] 6312/9069', '[SIII]A')
    diags.addClabel('[OIII] 4363/5007', '[OIII]A')
    
    # Define all the ions that are involved in the diagnostics
    all_atoms = diags.atomDict
    
    pn.log_.message('Atoms built')
    
    obs = pn.Observation(corrected = True)
    for atom in all_atoms:
        # Computes all the intensities of all the lines of all the ions considered
        for wavelength in all_atoms[atom].lineList:
            elem, spec = parseAtom(atom)
            intens1 = all_atoms[atom].getEmissivity(tem1, dens1, wave = wavelength) * dens1 * mass1
            intens2 = all_atoms[atom].getEmissivity(tem2, dens2, wave = wavelength) * dens2 * mass2
            obs.addLine(pn.EmissionLine(elem, spec, wavelength,
                                         obsIntens=[intens1, intens2, intens1+intens2], 
                                         obsError=[0.0, 0.0, 0.0]))
    
    pn.log_.message('Virtual observations computed')
    
    emisgrids = pn.getEmisGridDict(atomDict = all_atoms)
    
    pn.log_.message('EmisGrids available')
    
    # Produce a diagnostic plot for each of the two regions and another one for the (misanalyzed) overall region
    
    plt.subplot(2,2,1)
    diags.plot(emisgrids, obs, i_obs=0)

    plt.subplot(2,2,2)
    diags.plot(emisgrids, obs, i_obs=1)
    
    plt.subplot(2,1,2)
    pn.log_.level=3
    diags.plot(emisgrids, obs, i_obs=2)

if __name__ == '__main__':

    plot_2comp(tem1=1e4, tem2=1e4, dens1=3e2, dens2=5e5, mass1=1, mass2=5e-4)
    plt.show()
