"""
this file contains the definitions of the functions used to answer the questions of
exercise  1.1 of  LINE EMISSIVITIES

to execute these functions, one first needs to
import a few libraries as listed at the beginning of the file run_ex.py
as well as this file.
the functions can then be called  as explained in the file run_ex.py

the functions defined here can be modified. 
In that case, is is necessary, before using the modified version, 
to do "reload(ex1_1)" from the terminal
sometimes this does not work well. 
it is then recommended to quit python (ctrl D) and enter again (ipython --pylab)
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm # imports colormaps (which are not automatically imported with matplotlib)
import pyneb as pn


def p1():

    O3 = pn.Atom('O',3)  # define the Atom object for O++ and  fills it with data
    
    O3.printIonic() # print for each level the possible transitions and their wavelengths 
    
    O3.printIonic(printA=True) # print, in addition the transition probabilities


def p23(den=1e2, style='-', legend=True, ): # den=1e2, style='-', legend=True are the default options
	# other possible line styles for calling the function are '--','-.',':', and '.'
    O3 = pn.Atom('O',3) 	# define the Atom object for O++ and  fills it with data
    # O3.plotEmiss? shows all the default parameters 
    O3.plotEmiss(den=den, style=style, legend=legend) # function to plot the emissivities  
    # O3.plotEmiss? shows all the default parameters 
    plt.ylim((-29, -19))	# change the y limits on the current plot
    
def p3( ):
    O3 = pn.Atom('O',3) 	# define the Atom object for O++ and  fills it with data
    O3.printIonic(tem=10000., printCrit=True) # prints the critical densities for each transition
    
def p4():
    O3 = pn.Atom('O',3)
    print(O3)				# Print the element and ionisation, with the names of the atomic data files
    print(O3.atomFile)	# print the name of the atomic data file "atom" (containing Energies and transitions proba)
    print(O3.collFile)	# print the name of the atomic data file "coll" (containg the Omegas)
    print(O3.atomPath)	# print the directory where the atomFitsFile is read from
    print(O3.collPath)	# print the directory where the collFitsFile is rwad from
    O3.printSources()		# print the bibliographical sources from the atom and coll fits files.


def p5(den=1e2, style='-', legend=True, coeff=1.): # in parenthesis are the default options

    O3 = pn.Atom('O',3)
    if coeff != 1.:
        O3._A *= coeff 		# multiplies all the A values by coeff 
    O3.plotEmiss(den=den, style=style, legend=legend) # defines a plot of line emissivities with updated A 
    plt.ylim((-29, -19))
    
def p6(den=1e2, style='-', legend=True, coeff=1.):

    O3 = pn.Atom('O',3)
    if coeff != 1.:
        Ntrans = O3.collNLevels * (O3.collNLevels - 1) / 2
        for i, rec in enumerate(O3._CollData):
            for j in np.arange(Ntrans-1)+1:
                rec[j] = rec[j] * coeff  # multiplies all the omega values by coeff 
            O3._CollData[i] = rec
        O3.initOmegas()	# update some values depending on Omegas
    O3.plotEmiss(den=den, style=style, legend=legend) # defines a plot of line emissivities with updated omegas 
    plt.ylim((-29, -19))

def p7a(den=1e2, style='-', legend=True):
    
    atom = pn.Atom('O',3)
    tem_min = 1e2
    tem_max = 1.e7
    ionic_abund=1.
    atom.plotEmiss(tem_min=tem_min, tem_max=tem_max, den=den, ionic_abund=ionic_abund, style=style, legend=legend)
    plt.ylim((-25, -19))

def p7b():     
    dataplot = pn.DataPlot('O', 3)	#	defines the ion for which atomic data data will be plotted
    dataplot.plotOmega()      	# plots the  omega values in the range where they are calculated by atomic physics
 
def p7c(den=1e2, style='-', legend=True, OmegaInterp='Cheb'):
    
    # OmegaInterp is the method for interpolation and extrapolation of omegas.
    # it can take the values 'Cheb' (default) or 'Linear'
    atom = pn.Atom('O',3, OmegaInterp=OmegaInterp)
    tem_min = 1e2
    tem_max = 1.e7
    ionic_abund=1.
    atom.plotEmiss(tem_min=tem_min, tem_max=tem_max, den=den, ionic_abund=ionic_abund, style=style, 
                   legend=legend, temLog=True)
    plt.ylim((-25, -19))

def p8(elem = 'O', den=1e2, style='-', legend=True, OmegaInterp='Cheb'):
    
    # OmegaInterp is the method for interpolation and extrapolation of omegas.
    # it can take the values 'Cheb' (default) or 'Linear'
    atom = pn.Atom(elem,3, OmegaInterp=OmegaInterp)
    tem_min = 1e2
    tem_max = 1.e7
    ionic_abund=1.
    atom.plotEmiss(tem_min=tem_min, tem_max=tem_max, den=den, ionic_abund=ionic_abund, style=style, 
                   legend=legend, temLog=True)
    plt.ylim((-25, -19))