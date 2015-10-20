# Sample extinction in PyNeb
# shows how to display available extinction laws, select one or define a new one,
# and do some simple dereddening calculations
# Further examples can be found in other sample scripts

import numpy as np
import matplotlib.pyplot as plt
import pyneb as pn

# Convert wavelength to x
def x(wave):
    return 10000. / wave

# Define an extinction law (to be used below)
def my_X(wave, par=0):
    x = 10000. / wave
    Rv = 3.1
    X_lin = x/2. # linear part of the extinction law
    X_bump = 0.5*x**2. -6*x + 20. # bump part of the extinction law
    return Rv*np.where(x<5., X_lin, X_bump)

# Define a reddening correction object
RC = pn.RedCorr()

# List the available laws
RC.printLaws()

# Plot the available laws
RC.plot(laws='all')
plt.show()

# Choose the one we intend to use 
RC.law = 'CCM 89'
# or define a new one
RC.UserFunction = my_X
RC.law = 'user'

# Plot the selected law as a function of x
# Define an array in lambda to do the plot
wave= np.logspace(2.5, 5, 100)
# Plot commands
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_ylim([0, 15])
ax.set_xlim([0, 10])
ax.plot(x(wave), my_X(wave), label='%s' % (RC.law))
plt.xlabel('1/$\lambda$ ($\mu^{-1}$)')
plt.ylabel('A$_V$/E(B-V)')
plt.legend(loc='upper left')
plt.show()

# Correct observed line ratios
wave1 = 5007
I_obs1 = 4.0
wave2 = 4686
I_obs2 = 0.10

# Correct based on the given law and the observed Ha/Hb ratio
RC = pn.RedCorr(law='CCM 89')
I_obs_HaHb = 3.5 
I_theo_HaHb = 2.86 
RC.setCorr(I_obs_HaHb / I_theo_HaHb, 6563., 4861.)
print 'Correct based on the given law and the observed Ha/Hb ratio:'
print str(wave1) + ': I_obs =', I_obs1, ' I_dered =', I_obs1 * RC.getCorrHb(wave1)
print str(wave2) + ': I_obs =', I_obs2, ' I_dered =', I_obs2 * RC.getCorrHb(wave2)

# Correct based on the given law and c(Hb)
RC = pn.RedCorr(law='CCM 89', cHbeta=0.3)
print '\nCorrect based on the given law and the observed Ha/Hb ratio:'
print str(wave1) + ': I_obs =', I_obs1, ' I_dered =', I_obs1 * RC.getCorrHb(wave1)
print str(wave2) + ': I_obs =', I_obs2, ' I_dered =', I_obs2 * RC.getCorrHb(wave2)

