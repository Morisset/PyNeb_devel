import numpy as np
import pyneb as pn
from pyneb.utils.physics import IP
from pyneb.utils.misc import parseAtom

def getIonAb(obs_file, Ne_O2, Ne_Ar4, T_N2, T_O3, printIonAb = True):
    ### General settings
    # define an Observation object and assign it to name 'obs'
    obs = pn.Observation()
    # fill obs with data read from file obs_data, with lines varying across rows and a default percent error on line intensities
    obs.readData(obs_file, fileFormat='lines_in_rows', err_default=0.05, corrected=True)
    
    if 'H1_4.1m' in obs.lineLabels:
        temp = 14000
        dens = 10**3.5
        # Compute theoretical H5-4/Hbeta ratio from Hummer and Storey
        # Instatiate the H1 recombination atom
        H1 = pn.RecAtom('H', 1)
        IR2Opt_theo = H1.getEmissivity(temp, dens, label='5_4') / H1.getEmissivity(temp, dens, label='4_2')
        IR2Opt_obs = obs.getLine(label='H1_4.1m').corrIntens / obs.getLine(label='H1_4861A').corrIntens
        for line in obs.lines: 
            if line.label[-1] == 'm':
                line.corrIntens *= IR2Opt_theo/IR2Opt_obs
                
    # instanciation of all the needed Atom objects
    all_atoms = pn.getAtomDict(atom_list=obs.getUniqueAtoms())
    
    # define a dictionnary for the ionic abundances
    ab_dic = {}
    # we  use the following lines to determine the ionic abundances
    ab_labels = ['C3_1907A', 'C3_1909A',
                 'N2_6583A', 'N3_57.4m', 'N4_1487A', 
                 'O2_3726A', 'O2_3729A', 'O3_5007A', 'O4_25.9m',
                 'S2_6731A', 'S3_9069A', 'S4_10.5m',
                 'Ar3_7136A', 'Ar4_4740A', 'Ar5_7006A', 'Ar6_4.5m', 
                 'Ne2_12.8m', 'Ne3_3869A', 'Ne5_3426A', 'Ne6_7.6m', 
                 'Cl3_5518A', 'Cl3_5538A', 'Cl4_5323A', 
                 'Mg4_4.5m', 'Mg5_5.6m']
    # loop on the observed lines to determine the corresponding ionic abundances
    for line in obs.getSortedLines():
        # this is one way to define temp and dens in each zone
        # must be adapted to each case
        if (line.atom in all_atoms) and (line.label in ab_labels):
            IP_cut = 30. #
            if IP[line.atom] > IP_cut:
                temp = T_O3
                dens = Ne_Ar4
                IP_used = 'H'
            else:
                temp = T_N2
                dens = Ne_O2
                IP_used = 'L'                 
            ab = all_atoms[line.atom].getIonAbundance(line.corrIntens, temp, dens, to_eval=line.to_eval, Hbeta=100)
            if printIonAb:
                print('{0:13s} {1} '.format(line.label, IP_used) + ' '.join(['{0:>8.4f}'.format(t) for t in ab * 1e6]))
            if line.atom not in ab_dic:
                ab_dic[line.atom] = []
            ab_dic[line.atom].append(ab)
    
    He1 = pn.RecAtom('He', 1)
    He2 = pn.RecAtom('He', 2)
    ab_dic['He2']= He1.getIonAbundance(obs.getLine(label='He1_5876A').corrIntens, 
                                 T_N2, Ne_O2, wave=5876.0)
    ab_dic['He3'] = He2.getIonAbundance(obs.getLine(label='He2_4686A').corrIntens, 
                                 T_O3, Ne_Ar4, lev_i= 4, lev_j= 3)
    for atom in ab_dic:
        mean = np.mean(np.asarray(ab_dic[atom]))
        ab_dic[atom] = mean
    
    for ion in np.sort(ab_dic.keys()):
        if printIonAb:
            print('{0} = {1:.3}'.format(ion, np.log10(ab_dic[ion])+12))

    return ab_dic

def getElemAb(ab_dic, printAb = True):
    # Instantiation of the ICF object
    icf = pn.ICF()
    # Computing the elemental abundances from all ICF rules.
    icf_labels = []
    #for label in icf_labels:
    #    print('{0} is a rule for {1}'.format(label, icf.all_icfs[label]['elem']))
    
    # The following computes the elemental abundances from the rules.
    elem_abun = icf.getElemAbundance(ab_dic, icf_labels)
    if printAb:
        for icf_ref in np.sort(elem_abun.keys()):
            if np.log10(elem_abun[icf_ref]) > -10:
                print('{0} {1}={3:.3} using {2}'.format(icf_ref, 
                                              icf.all_icfs[icf_ref]['elem'],
                                              icf.getExpression(icf_ref), 
                                              np.log10(elem_abun[icf_ref])+12))

    return elem_abun
    
def p5():

    # this will plot all the available data in pyneb for the omegas
    ions = ['C2','C3', 'N2', 'N3','O2', 'O3', 'O4','Ne2','Ne3', 'Ne5', 'S2', 'S3', 'S4','Cl3', 'Ar2','Ar3','Ar4','Ar5',]

    for ion in ions:
        # split ion into elem and spec, e.g 'O3' into 'O' and 3
        elem, spec = parseAtom(ion)
        # instanciate the corresponding Atom object
        atom = pn.Atom(elem, spec)
        # print information including transition probabilities
        #atom.printIonic(printA = True)
        dp = pn.DataPlot(elem, spec)
        dp.plotOmega(save=True)

        # prepare a new figure
        #plt.figure()
        # plot energy levels
        #dataplot.plotOmega()
        #plt.savefig('{0}_{1}_omegas.pdf'.format(elem, spec.replace('/', '_')))
        
     
