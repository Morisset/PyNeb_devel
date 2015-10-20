# EmisGrid example script
# Shows how to produce, save, restore and plot emission maps

import pyneb as pn
import matplotlib.pyplot as plt

plt.figure()
O3map = pn.EmisGrid('O', 3, n_tem=30, n_den=30) # Instantiate an O III atom and computes the emissivity of all its lines in a 30x30 grid
print "Plot the 4363/5007 ratio map with a 30x30 Ne-Te grid"
O3map.plotImage(wave1=4363, wave2=5007)  # Plot

plt.figure()
O3map = pn.EmisGrid('O', 3, n_tem=100, n_den=100) # Increase the grid resolution
print "Increase the grid resolution to 100x100"
O3map.plotImage(wave1=4363, wave2=5007)  # Plot

plt.figure()
O3map = pn.EmisGrid('O', 3, n_tem=100, n_den=100, tem_max=1.5e4, den_max=1e5) # Change the axes limits
print "Change the axes limits"
O3map.plotImage(wave1=4363, wave2=5007)  # Plot
O3map.save('O3_1k.pypic') # Save for later use

plt.figure()
O3map = pn.EmisGrid(restore_file='O3_1k.pypic') # Recover the map from a previous computation
print "Plot with a different color map"
O3map.plotImage(wave1=4363, wave2=5007, cmap=plt.cm.bone) # Plot with a different color map

plt.figure()
print "Plot the 4363/88m line ratio"
O3map.plotImage(wave1=4363, wave2=884000) # Plot a different line ratio 
