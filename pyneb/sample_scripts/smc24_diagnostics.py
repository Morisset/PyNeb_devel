# Diagnostic plot

# Imports
import pyneb as pn
import matplotlib.pyplot as plt

### General settings
# Setting verbosity level. Enter pn.my_logging? for details
pn.log_.level = 3

# Set to True if the emission maps have already been generated
restore = True

# Define where emission maps are to be stored (restore = False) or read from (restore = True)
pypic_path = './pypics/'

# Adopt an extinction law
extinction_law = 'CCM89'

# Define the data file
obs_data = 'smc24.dat'

# Define the plot title
title = 'SMC 24'

### Read and deredden observational data
# define an Observation object and assign it to name 'obs'
obs = pn.Observation()

# fill obs with data read from file obs_data, with lines varying across rows and a default percent error on line intensities
obs.readData(obs_data, fileFormat='lines_in_rows', err_default=0.05)

# deredden data with Cardelli's law
obs.extinction.law = extinction_law
obs.correctData()

### Include the diagnostics of interest
# instantiate the Diagnostics class
diags = pn.Diagnostics()
# include in diags the relevant line ratios
diags.addDiag([
              '[NII] 5755/6584', 
              '[OII] 3726/3729', 
              '[OIII] 4363/5007', 
              '[SII] 6731/6716', 
              '[SII] 4072+/6720+',
              '[SIII] 6312/18.7m', 
              '[NeIII] 3930+/15.6m', 
              ])
diags.addClabel('[SII] 6731/6716', '[SII]a')
diags.addClabel('[SII] 4072+/6720+', '[SII]b')

# Create the emission maps to be compared to the observation data (some overkill here)
emisgrids = pn.getEmisGridDict(atom_list=diags.getUniqueAtoms(), den_max=1e6)

### Plot
# Create the contour plot as the intersection of tem-den emission maps with dereddened line ratios
diags.plot(emisgrids, obs)

# Place the title
plt.title(title)

# Display the plot
plt.show()
