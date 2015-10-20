"""
Sample PyNeb script
Plots the [O III] 4363/5007 ratio as a function of Te for several Ne values
"""

# Import relevant packages
import numpy as np
import matplotlib.pyplot as plt
import pyneb as pn

# Set high verbosity level to keep track of atom creation
pn.log_.level = 3

# Create a collection of atoms - a bit overkill if we just need O III
adict = pn.getAtomDict()

# Lower verbosity level
pn.log_.level = 2

# Function to compute line ratio
def line_ratio(atom, wave1, wave2, tem, den):
    emis1 = adict[atom].getEmissivity(tem, den, wave = wave1)
    emis2 = adict[atom].getEmissivity(tem, den, wave = wave2)
    return emis1 / emis2

# Define array of Te 
tem = np.arange(5000, 18000, 30)

# Plot
plt.figure(1)
for den in [1e2, 1e3, 1e4, 1e5]:
    plt.semilogy(tem, line_ratio('O3', 4363, 5007, tem, den), label = 'Ne={0:.0e}'.format(den))
plt.xlabel('T$_e$ [K]')
plt.ylabel(r'[OIII] 4363/5007 $\AA$')
plt.legend(loc=2)

plt.show()


