'''
Created on Jun 8, 2016

@author: morisset
'''
import numpy as np
import pyneb as pn

o3 = pn.Atom('O',3)
tem = np.linspace(5000, 15000, 50)
den = np.logspace(1, 5, 50)
em = o3.getEmissivity(tem, den, wave=5007, product=True)
print em
