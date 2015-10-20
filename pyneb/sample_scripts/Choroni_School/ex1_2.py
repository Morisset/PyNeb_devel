"""
this file contains the definitions of the functions used to answer the questions of
exercise  1.2 of  LINE EMISSIVITIES

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
from matplotlib import cm
import pyneb as pn


# define some tools common to all the questions

tem = np.logspace(np.log10(10), np.log10(30000), 1000)

#------------------------------------------------------------------------------------------------------------------
"""
question 1
"""
def alpha_B(Te): # Recomb. coefficient, case B

    T4 = Te/1e4
    return 2.6e-13/T4

def Gphot(Te, Tstar): # Gains for a H nebula at Te, photoionized by a star at Tstar

    return 3./2. * pn.CST.BOLTZMANN * Tstar * alpha_B(Te)

def plot_Gain(Tstar=50000, linestyle='-'):
    plt.plot(tem, np.log10(Gphot(tem, Tstar)), label='Gain', linewidth=3, linestyle=linestyle, color='blue')

#------------------------------------------------------------------------------------------------------------------
"""
question 2
"""
def Lff(Te):	# H free-free losses
    
    Z = 1
    gff = 1.3
    #formula for loss due to free-free radiation
    #Lff =  (32. * np.pi * pn.CST.ECHARGE**6 * Z**2 / (3.**(1.5) * pn.CST.HPLANCK * pn.CST.EMASS * pn.CST.CLIGHT**3) *
    #        (2 * np.pi * pn.CST.BOLTZMANN * Te / pn.CST.EMASS)**0.5 * gff)
    
    Lff = 1.42e-27 * Z**2 * Te**0.5 * gff
    return Lff

def beta_B(Te):	#Energy averaged recomb. coeff.
    """
    We now take = alpha, but will need to be updated ***
     # TODO: find expression for Beta_B
    """
    return alpha_B(Te)

def Lrec(Te):	# H losses due to recombinationfree-free losses
    return  pn.CST.BOLTZMANN * Te * beta_B(Te)

def LtotalH(Te):	# total H losses (collisional excitation of Lya not considered)
    return Lrec(Te) + Lff(Te)

def plot_LossH():	# plot H energylosses
    plt.plot(tem, np.log10(Lff(tem)), label='H Free-Free loss', color='black', linestyle='--')
    plt.plot(tem, np.log10(Lrec(tem)), label='H recomb. loss', color='black')
    plt.plot(tem, np.log10(LtotalH(tem)), label='TOTAL H-loss', linewidth=3, color='black')
    
#------------------------------------------------------------------------------------------------------------------
"""
question 3
"""

def plot_LossO(atom, OoH=4e-4, den=1e2):	#energylosses due to O3 for an abundance O/H of 4e-4
    # O3 has been defined at the beginning as O3 = pn.Atom('O', 3, OmegaInterp='Linear')
    print(atom)
    atom.plotEmiss(den=den, ionic_abund=OoH, plot_total = True, legend=True, total_color='green')

#------------------------------------------------------------------------------------------------------------------
"""
question 4
"""   

def plot_Loss(atom, OoH=4e-4, den=1e2, onlyTotal=False, linestyle='-'):

    total_Hloss = LtotalH(tem)
    
    total_loss = LtotalH(tem) # !!!Don't use total_loss = total_Hloss!!!
    for wave in atom.lineList:
        total_loss += OoH * atom.getEmissivity(tem, den, wave=wave) # adds the losses for all the lines of O3 
        
    plt.plot(tem, np.log10(total_loss), color='red', label = 'TOTAL loss', linewidth=3, linestyle=linestyle)
    if not onlyTotal:
        plt.plot(tem, np.log10(total_Hloss), label='TOTAL H-loss', linewidth=1, color='black')
        atom.plotEmiss(den=den, ionic_abund=OoH, plot_total=True, plot_only_total=True, 
                   legend=True, total_color='green', total_label='TOTAL O-loss')

