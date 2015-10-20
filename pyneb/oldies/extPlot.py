# Simple plot of extinction functions

import numpy as np
import matplotlib.pyplot as plt

import Extinction as ext

e = ext.ExtClass()

flam = e.extFnDict['GAL_CCM']
plt.plot(flam[:,0],flam[:,1], 'k-')

flam = e.extFnDict['GAL_GCC']
plt.plot(flam[:,0],flam[:,1], 'k-')

flam = e.extFnDict['GAL_SM']
plt.plot(flam[:,0],flam[:,1], 'k-')

flam = e.extFnDict['GAL_JBK']
plt.plot(flam[:,0],flam[:,1], 'k-')

flam = e.extFnDict['LMC_HOWARTH']
plt.plot(flam[:,0],flam[:,1], 'r-')

flam = e.extFnDict['LMC_GCMLW']
plt.plot(flam[:,0],flam[:,1], 'r-')

flam = e.extFnDict['SMC_GCMLW']
plt.plot(flam[:,0],flam[:,1], 'b-')

flam = e.extFnDict['SMC_PREVOT']
plt.plot(flam[:,0],flam[:,1], 'b-')

plt.xlabel(r'Wave number ($\mu$m$^{-1}$)')
plt.ylabel(r'F$_{\lambda}$')
plt.title('Extinction functions')
plt.show()
