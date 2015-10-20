import pyneb as pn
import numpy as np
import matplotlib.pyplot as plt

def p1(obs, do_plot=True):    
    
    # Change the atomic data:
    pn.atomicData.setDataFile('s_iii_coll_TG99.fits')
    pn.atomicData.setDataFile('s_iii_atom_MZ82b-HSC95-LL93.fits')
    
    #Instanciate the Diagnostic object
    diags = pn.Diagnostics()
    # We will use the following diagnostics
    diags.addDiag([
                '[NII] 5755/6584',
                #'[OII] 3727+/7325+',
                '[OIII] 4363/5007', 
                '[SII] 6731/6716',
                '[SIII] 6312/9069',
#                '[SII] 4072+/6720+',
                ])
    
    # we can choose to do the diagnostic plots of the whole set of observations
    if do_plot:
        emisgrids = pn.getEmisGridDict(atomDict=diags.atomDict, restore=False, save=False,
                                       n_tem=150, n_den=150, tem_min=5000., tem_max=20000., 
                                       den_min=10., den_max=1e5)
    
        plt.figure(figsize=(30, 30))
        for i, obs_name in enumerate(obs.names):
            plt.subplot(6, 5, i+1)
            diags.plot(emisgrids, obs, i_obs=i)
            plt.title(obs_name)
        plt.savefig('diags.pdf')
        plt.close()
    
    # Determination of Te and Ne by intersection of 2 diagnostics:
    pn.log_.level = 3
    temp_O3, dens_S2a = diags.getCrossTemDen('[OIII] 4363/5007', '[SII] 6731/6716', obs=obs)
    temp_S3, dens_S2b = diags.getCrossTemDen('[SIII] 6312/9069', '[SII] 6731/6716', obs=obs)
    temp_N2, dens_S2c = diags.getCrossTemDen('[NII] 5755/6584', '[SII] 6731/6716', obs=obs)
#    temp_S2, dens_S2d = diags.getCrossTemDen('[SII] 4072+/6720+', '[SII] 6731/6716', obs=obs)
    
    pn.log_.level = 1
    #In the cases the density is not defined, we choose to set it to 10 and compute new Te
    tt = np.isnan(dens_S2a)
    dens_S2a[tt] = 10.
    O3 = pn.Atom('O', 3)
    temp_O3[tt] = O3.getTemDen((obs.getLine(label='O3_4363A').corrIntens/
                                obs.getLine(label='O3_5007A').corrIntens)[tt], 
                               den=10, wave1=4363, wave2=5007)
    tt = np.isnan(dens_S2b)
    dens_S2b[tt] = 10.
    S3 = pn.Atom('S', 3)
    temp_S3[tt] = S3.getTemDen((obs.getLine(label='S3_6312A').corrIntens/
                                obs.getLine(label='S3_9069A').corrIntens)[tt], 
                               den=10, wave1=6312, wave2=9069)
    tt = np.isnan(dens_S2c)
    dens_S2c[tt] = 10.
    N2 = pn.Atom('N', 2)
    temp_N2[tt] = N2.getTemDen((obs.getLine(label='N2_5755A').corrIntens/
                                obs.getLine(label='N2_6583A').corrIntens)[tt], 
                               den=10, wave1=5755, wave2=6583)
    
    # Here we adopt an average value of the density
    mean_dens = (dens_S2a + dens_S2b + dens_S2c) / 3
    
    
    O2 = pn.Atom('O', 2)
    temp_O2 = O2.getTemDen((obs.getLine(label='O2_3727A+').corrIntens/obs.getLine(label='O2_7325A+').corrIntens),
                           den = mean_dens, to_eval = '(L(3726)+L(3729))/(I(4,2)+I(5,2)+I(4,3)+I(5,3))')
    S2 = pn.Atom('S', 2)
    temp_S2 = S2.getTemDen((obs.getLine(label='S2_4072A+').corrIntens/obs.getLine(label='S2_6716A').corrIntens),
                           den = mean_dens, to_eval = '(L(4076)+L(4069))/L(6716)')

    return mean_dens, temp_O2, temp_S2, temp_N2, temp_O3, temp_S3

def p2(temp_S3, temp_O3):
    tS3 = temp_S3.copy()
    tO3 = temp_O3.copy()
    # we define 3 Te for low, mid and high ionization regions
    Tmid = tS3
    tt = np.isnan(tS3)
    Tmid[tt] = tO3[tt] * 0.83 + 1700
    
    Thigh = tO3
    tt = np.isnan(tO3)
    Thigh[tt] = (tS3[tt] - 1700) / 0.83
    
    #Tlow = (temp_S2+temp_O2) / 2.    
    Tlow = Thigh*0.7 + 3000

    return Tlow, Tmid, Thigh
    
def p2_P(temp_N2, temp_O3):
    tN2 = temp_N2.copy()
    tO3 = temp_O3.copy()
    # we define 3 Te for low, mid and high ionization regions
    Thigh = tO3
    Tlow = 3140 + 0.672 * tO3
    tt = np.isnan(tO3)
    Tlow[tt] = tN2[tt]
    Thigh[tt] = (tN2[tt] - 3140 ) / 0.672
    return Tlow, Thigh
    
def p3(obs, Tlow, Tmid, Thigh, mean_dens, verbose=True):
    #we associate the Te with each ion
    Te_dic = {'N2' : Tlow,
              'O2' : Tlow,
              'S2' : Tlow,
              'S3' : Tmid,
              'Ar3' : Tmid,
              'O3'  : Thigh,
              'Ne3' : Thigh}
    
    all_atoms = pn.getAtomDict(atom_list=obs.getUniqueAtoms())
    
    # This is the dictionnary which will contain the ionic abundances
    ab_dic = {}
    
    # we  use the following lines to determine the ionic abundances
    ab_labels = ['N2_6584A', 'O2_3727A+', 'O3_5007A', 'S2_6716A', 'S3_9069A', 'Ar3_7136A', 'Ne3_3869A']
    for line in obs.getSortedLines():
        if line.label in ab_labels:
            ab = all_atoms[line.atom].getIonAbundance(line.corrIntens, Te_dic[line.atom], mean_dens, 
                                                      to_eval=line.to_eval, Hbeta=100)
            if verbose:
                print '{0:9s}'.format(line.label) + ' '.join(['{0:>8.2f}'.format(t) for t in 12 + np.log10(ab)])
            ab_dic[line.atom] = ab
    return ab_dic

if __name__ == '__main__':
    plt.figure()
    #Read the observations. Notice that in this file the errors are absolutes
    obs = pn.Observation('NGC300.dat', corrected=True, errIsRelative=False)
    # the galactocentric distance has been also read and the values are in "obsIntens"
    Rgal = obs.getLine(label='DIST').obsIntens
    # Compute the densities and temperatures
    mean_dens, temp_O2, temp_S2, temp_N2, temp_O3, temp_S3 = p1(obs)
    # Instantiate a table to write the results in anm ascii file
    # compute Te for the low, mid and high ionisation regions
    
    Tlow, Tmid, Thigh = p2(temp_S3, temp_O3)
    # Compute the ionic abundances
    ab_dic = p3(obs, Tlow, Tmid, Thigh, mean_dens)

    OsH = ab_dic['O2'] + ab_dic['O3']
    NsH = ab_dic['N2'] * OsH / ab_dic['O2']
    NesH = ab_dic['Ne3'] * OsH / ab_dic['O3']
    
    plt.scatter(Rgal, np.log10(OsH) + 12, label = 'O/H', color='red')
    plt.scatter(Rgal, np.log10(NsH) + 12, label = 'N/H', color='blue')
    plt.scatter(Rgal, np.log10(NesH) + 12, label = 'Ne/H', color='green')
    
    Tlow2, Thigh2 = p2_P(temp_N2, temp_O3)
    Tmid = Tlow2

    ab_dic2 = p3(obs, Tlow2, Tmid, Thigh2, mean_dens)
    
    OsH2 = ab_dic2['O2'] + ab_dic2['O3']
    NsH2 = ab_dic2['N2'] * OsH2 / ab_dic2['O2']
    NesH2 = ab_dic2['Ne3'] * OsH2 / ab_dic2['O3']
    
    plt.scatter(Rgal, np.log10(OsH2) + 12, label = 'O/H', color='red', marker='o', alpha=.3)
    plt.scatter(Rgal, np.log10(NsH2) + 12, label = 'N/H', color='blue', marker='o', alpha=.3)
    plt.scatter(Rgal, np.log10(NesH2) + 12, label = 'Ne/H', color='green', marker='o', alpha=.3)
        
    plt.legend()
    plt.xlabel('RGal')
    plt.ylabel('12+log(X/H)')
    plt.show()
    
