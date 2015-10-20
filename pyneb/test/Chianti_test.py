import pyneb as pn
import matplotlib.pyplot as plt

for atom in pn.atomicData.getAllAtoms():
    print('Doing {0}'.format(atom))
    dp = pn.DataPlot(atom=atom, OmegaInterp='Linear')
    
    dp.atom_n_max = np.min((dp.atom_n_max, 7))
    if len(dp.atom_data) > 1:
        print('Plotting {0}'.format(atom))
        dp.plotAllA(save = True)
        plt.close()
      
    dp.coll_n_max = np.min((dp.coll_n_max, 7))
    if len(dp.coll_data) > 1:
        print('Plotting {0}'.format(atom))
        dp.plotOmega(save = True)
        plt.close()
