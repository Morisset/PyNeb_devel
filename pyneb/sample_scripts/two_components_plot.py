# Analysis plot of a simple two-component model, meant to illustrate the bias arising 
# from assuming hat the region is homogeneous in density
# First, an emission region made up of two different subregions is modelled,
# each with a different mass and density. The resulting overall emissivity is computed
# Second, the region is analyzed as if it were a homogeneous region

import pyneb as pn
import matplotlib.pyplot as plt
from pyneb.utils.misc import parseAtom


def plot_2comp(tem1=1e4, tem2=1e4, dens1=3e2, dens2=5e5, mass1=1, mass2=5e-4):
        
    # List of diagnostics used to analyze the region
    diags = pn.Diagnostics()
    
    diags.addDiag([
                    '[NI] 5198/5200',
                    '[NII] 5755/6548',
                    '[OII] 3726/3729',
                    '[OII] 3727+/7325+',
                    '[OIII] 4363/5007', 
                    '[ArIII] 5192/7136',
                    '[ArIII] 5192/7300+',
                    '[ArIV] 4740/4711',
                    '[ArIV] 7230+/4720+',
                    '[SII] 6731/6716', 
                    '[SII] 4072+/6720+',
                    '[SIII] 6312/9069',
                    '[ClIII] 5538/5518'
                    ])
    """    
    for diag in pn.diags_dict:
        if diag[0:7] != '[FeIII]':
            diags.addDiag(diag)
            print('Adding', diag)
    diags.addClabel('[SIII] 6312/9069', '[SIII]A')
    diags.addClabel('[OIII] 4363/5007', '[OIII]A')
    """    
    # Define all the ions that are involved in the diagnostics
    adict = diags.atomDict
    pn.log_.message('Atoms built')
    
    obs = pn.Observation(corrected = True)
    for atom in adict:
        # Computes all the intensities of all the lines of all the ions considered
        for line in pn.LINE_LABEL_LIST[atom]:
            if line[-1] == 'm':
                wavelength = float(line[:-1])*1e4
            else:
                wavelength = float(line[:-1])
            elem, spec = parseAtom(atom)
            intens1 = adict[atom].getEmissivity(tem1, dens1, wave = wavelength) * dens1 * mass1
            intens2 = adict[atom].getEmissivity(tem2, dens2, wave = wavelength) * dens2 * mass2
            obs.addLine(pn.EmissionLine(elem, spec, wavelength,
                                         obsIntens=[intens1, intens2, intens1+intens2], 
                                         obsError=[0.0, 0.0, 0.0]))
    
    pn.log_.message('Virtual observations computed')
    emisgrids = pn.getEmisGridDict(atomDict=adict)
    
    pn.log_.message('EmisGrids available')
    
    # Produce a diagnostic plot for each of the two regions and another one for the 
    # (misanalyzed) overall region
    f, axes = plt.subplots(2,2)

    diags.plot(emisgrids, obs, i_obs=0, ax=axes[0,0])
    diags.plot(emisgrids, obs, i_obs=1, ax=axes[0,1])
    diags.plot(emisgrids, obs, i_obs=2, ax=axes[1,0])

if __name__ == '__main__':

    plot_2comp(tem1=1e4, tem2=1e4, dens1=3e2, dens2=5e5, mass1=1, mass2=5e-4)
    plt.show()
