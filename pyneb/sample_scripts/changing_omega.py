# Sample script
# Produces contour plots of OIII and SII line ratios
# Two versions of the SII plot are produced with two different sets of collision strengths
# A colormap and a plot of the ratio of the two SII determinations is also shown

import numpy as np
import matplotlib.pyplot as plt
import pyneb as pn

# Compute emission grids
EM_O3 = pn.EmisGrid('O', 3, n_den=100, n_tem=100)
EM_S2 = pn.EmisGrid('S', 2, n_den=100, n_tem=100)
pn.atomicData.setDataFile('s_ii_coll_RBS96.fits')
EM_S2_GR = pn.EmisGrid('S', 2, n_den=100, n_tem=100)
pn.atomicData.resetDataFileDict()

# Print populations a la ionic
EM_O3.atom.printIonic()

# Contour plots for OIII and the two versions of the SII ratios
EM_O3.plotContours('L(4363)/L(4959)')
plt.show()
EM_S2.plotContours('L(6731)/L(6716)', linestyles='--', low_level= -0.3, high_level=0.15)
plt.show()
EM_S2_GR.plotContours('L(6731)/L(6716)', linestyles=':', low_level= -0.3, high_level=0.15)
plt.show()

# Image of the ratio between the two SII computations
plt.imshow(EM_S2.getGrid(to_eval='L(6731) / L(6716)') / 
           EM_S2_GR.getGrid(to_eval='L(6731) / L(6716)'))
plt.colorbar()
plt.show()

# Plot of the ratio between the two SII computations
plt.semilogx(EM_S2.den, (EM_S2.getGrid(to_eval='L(6731) / L(6716)') / EM_S2_GR.getGrid(to_eval='L(6731) / L(6716)'))[50, :])
plt.xlabel('N$_e$')
plt.ylabel('N$_e,1$/N$_e,2$')
plt.show()

EM_S2.save('em_S2.pypics')
EM_S2_bis = pn.EmisGrid(restore_file='em_S2.pypics')
