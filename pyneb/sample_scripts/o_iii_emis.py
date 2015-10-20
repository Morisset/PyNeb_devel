import pyneb as pn
import numpy as np
import matplotlib.pyplot as plt

o3=pn.Atom('O', 3)
tem=np.arange(100)*300+300
den = 1000
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_ylim([1.e-30, 5e-20])
lineList=o3.lineList
lineList=[1666, 4363, 4959, 5007, 518000, 880000]
for line in lineList:
    y=o3.getEmissivity(tem, den, wave=line)
    plt.semilogy(tem, y,  label="{:.0f}".format(line))

plt.xlabel('T$_e$ [K]')
plt.ylabel("j(T) [erg cm$^{-3}$ s${-1}$]")
plt.legend(loc='lower right')
plt.title('[O III] emissivities @ N$_e$={:.0f}'.format(den))
plt.show()
