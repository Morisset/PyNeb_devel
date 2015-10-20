# plot_Fe3 example script
# Plots diagnostic diagrams from several Fe III lines, using two different sets of atomic data

import numpy as np
import pyneb as pn
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

pn.config.set_noExtrapol(True)
pn.log_.level = 3
# Read the observed intensities, reddening corrected, from a file
obs = pn.Observation('Fe3.dat', fileFormat='lines_in_rows', err_default = 0.01, corrected = True)

# Build an atom with default values from Quinet 1996, two kinds of interpolation
pn.atomicData.setDataFile('fe_iii_atom_Q96_J00.dat')
pn.atomicData.setDataFile('fe_iii_coll_Z96.dat')
Fe3_Z = pn.Atom('Fe', 3)

Fe3EM_Q_C = {'Fe3': pn.EmisGrid(n_tem = 100, n_den = 100, tem_min=5000., tem_max=15000., 
                 den_min=100., den_max=1.e5, atomObj= Fe3_Z)}
Fe3EM_Q_L = {'Fe3': pn.EmisGrid(n_tem = 100, n_den = 100, tem_min=5000., tem_max=15000., 
                 den_min=100., den_max=1.e5, OmegaInterp='Linear', atomObj= Fe3_Z)}

# Change to alternate set of data from Bautista, Ballance and Quinet 2010
pn.atomicData.setDataFile('fe_iii_atom_BBQ10.dat')
pn.atomicData.setDataFile('fe_iii_coll_BBQ10.dat')

Fe3_B = pn.Atom('Fe', 3)
# Build an atom with alternate set of data, two kinds of interpolation
Fe3EM_B_C = {'Fe3': pn.EmisGrid(n_tem = 100, n_den = 100, tem_min=5000., tem_max=15000., 
                 den_min=100., den_max=1.e5, atomObj= Fe3_B)}
Fe3EM_B_L = {'Fe3': pn.EmisGrid(n_tem = 100, n_den = 100, tem_min=5000., tem_max=15000., 
                 den_min=100., den_max=1.e5, OmegaInterp='Linear', atomObj= Fe3_B)}

pn.atomicData.setDataFile('fe_iii_atom_BBQ10.fits')
pn.atomicData.setDataFile('fe_iii_coll_BBQ10.fits')
Fe3_F = pn.Atom('Fe', 3)
# Build an atom with alternate set of data, two kinds of interpolation
Fe3EM_F_C = {'Fe3': pn.EmisGrid(n_tem = 100, n_den = 100, tem_min=5000., tem_max=15000., 
                 den_min=100., den_max=1.e5, atomObj= Fe3_F)}
Fe3EM_F_L = {'Fe3': pn.EmisGrid(n_tem = 100, n_den = 100, tem_min=5000., tem_max=15000., 
                 den_min=100., den_max=1.e5, OmegaInterp='Linear', atomObj= Fe3_F)}

pn.atomicData.setDataFile('fe_iii_atom_Q96_J00.dat')
pn.atomicData.setDataFile('fe_iii_coll_BB14.dat')
Fe3_Ba = pn.Atom('Fe', 3)
# Build an atom with alternate set of data, two kinds of interpolation
Fe3EM_Ba_C = {'Fe3': pn.EmisGrid(n_tem = 100, n_den = 100, tem_min=5000., tem_max=15000., 
                 den_min=100., den_max=1.e5, atomObj= Fe3_Ba)}
Fe3EM_Ba_L = {'Fe3': pn.EmisGrid(n_tem = 100, n_den = 100, tem_min=5000., tem_max=15000., 
                 den_min=100., den_max=1.e5, OmegaInterp='Linear', atomObj= Fe3_Ba)}

diags = pn.Diagnostics()
diags.addDiagsFromObs(obs)
for i, label in enumerate(diags.getDiagLabels()):
    diags.addClabel(label, str(i))


plt.figure()
diags.plot(Fe3EM_Q_L, obs)
plt.title('Zhang, Quinet 1996')

plt.figure()
diags.plot(Fe3EM_B_L, obs)
plt.title('Bautista et al 2010')
"""
plt.figure()
diags.plot(Fe3EM_F_L, obs)
plt.title('Bautista et al 2010, fits')
"""
plt.figure()
diags.plot(Fe3EM_Ba_L, obs)
plt.title('Badnell and Ballance 2014')

plt.show()

"""
# Uncomment this if you want the numerical solution of two specific diagnostics

tem, den = diags.getCrossTemDen('[FeIII] 4659/4009', '[FeIII] 4703/4735', obs=obs)
print tem, den
"""


def plot_As():
    f, axes = plt.subplots(2, 2)
    
    im1 = axes[0, 0].imshow(np.log10(Fe3_Ba.getA()), interpolation='None', vmin=-7, vmax=1.5)
    axes[0, 0].set_title('Quinet, P. 1996')
    
    im2 = axes[0, 1].imshow(np.log10(Fe3_B.getA()), interpolation='None', vmin=-7, vmax=1.5)
    axes[0, 1].set_title('Bautista, Ballance and Quinet 2010')
    divider = make_axes_locatable(axes[0,1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb2 = f.colorbar(im2, cax=cax)
    cb2.set_label('log(A)')
    
    im3 = axes[1, 0].imshow(np.log10(Fe3_B.getA() / Fe3_Ba.getA()), interpolation='None', vmin=-1, vmax=1)
    axes[1, 0].set_title('Q96/BBQ10')
    divider = make_axes_locatable(axes[1,0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb3 = f.colorbar(im3, cax=cax)
    
def plot_Omegas():
    f, axes = plt.subplots(1, 3, figsize=(20, 7))
    
    im1 = axes[0].imshow(np.log10(Fe3_Ba.getOmega(1e4)), interpolation='None', vmin=-4, vmax=1.)
    axes[0].set_title('Badnell and Ballance 2014')

    im2 = axes[1].imshow(np.log10(Fe3_B.getOmega(1e4)), interpolation='None', vmin=-4, vmax=1.)
    axes[1].set_title('Bautista, Ballance and Quinet 2010')

    im3 = axes[2].imshow(np.log10(Fe3_Z.getOmega(1e4)), interpolation='None', vmin=-4, vmax=1.)
    axes[2].set_title('Zhang, H. 1996')
    
    divider = make_axes_locatable(axes[2])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = f.colorbar(im3, cax=cax)
    cb.set_label('log(Omega)')
    
# Reset default dataset
#pn.atomicData.setDataFileDict(pn.atomicData.defaultDict)