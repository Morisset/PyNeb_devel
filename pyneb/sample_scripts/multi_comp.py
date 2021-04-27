import pyneb as pn
import numpy as np

"""
This is a toy-program to test the effect of a 2-densities region on the nebular analysis of the observations.
The 2 regions have different temperature (tem1 and tem2), densities (dens1 and dens2), masses (mass1 and mass2, 
whatever unit, the important being the ratio). Intensities for various lines from various ions (included Hbeta) are computed: 
emis1 from region 1, emis2 from region2 and emis 3 from the sum of the 2 regions.
Then Te and Ne are computed from these pseudo-observations, using various diagnostic line ratios.
Finally, ionic fractions are obtained and compared to the original ionic fractions used to compute the intensities.
"""

pn.log_.level = 3

# Parameters of the toy model:
tem1 = 1.e4
tem2 = 1.e4
dens1 = 1e3
dens2 = 5.e4
mass1 = 1
mass2 = 1e-2
Sp = -6 #log S+/H+
Np = -5 #log N+/H+
Op = -4 #log O+/H+
Opp = -4 #log O++/H+

# The atoms and lines we will need:
S2 = pn.Atom('S', 2)
N2 = pn.Atom('N', 2)
O2 = pn.Atom('O', 2)
O3 = pn.Atom('O', 3)
S2_lambda = (4069, 4076, 6716, 6731)
N2_lambda = (6583, 5755)
O2_lambda = (3726, 3729)
O3_lambda = (5007, 4363)

# computing line intensities:
emis1 = {}
emis2 = {}
emis3 = {}
for line in S2_lambda:
    emis1['SII_'+str(line)] = S2.getEmissivity(tem1, dens1, wave = line) * 10.**Sp * dens1 * mass1
    emis2['SII_'+str(line)] = S2.getEmissivity(tem2, dens2, wave = line) * 10.**Sp * dens2 * mass2
    emis3['SII_'+str(line)] = emis1['SII_'+str(line)] + emis2['SII_'+str(line)]
for line in N2_lambda:
    emis1['NII_'+str(line)] = N2.getEmissivity(tem1, dens1, wave = line) * 10.**Np * dens1 * mass1
    emis2['NII_'+str(line)] = N2.getEmissivity(tem2, dens2, wave = line) * 10.**Np * dens2 * mass2
    emis3['NII_'+str(line)] = emis1['NII_'+str(line)] + emis2['NII_'+str(line)]
for line in O3_lambda:
    emis1['OIII_'+str(line)] = O3.getEmissivity(tem1, dens1, wave = line) * 10.**Opp * dens1 * mass1
    emis2['OIII_'+str(line)] = O3.getEmissivity(tem2, dens2, wave = line) * 10.**Opp * dens2 * mass2
    emis3['OIII_'+str(line)] = emis1['OIII_'+str(line)] + emis2['OIII_'+str(line)]
for line in O2_lambda:
    emis1['OII_'+str(line)] = O2.getEmissivity(tem1, dens1, wave = line) * 10.**Op * dens1 * mass1
    emis2['OII_'+str(line)] = O2.getEmissivity(tem2, dens2, wave = line) * 10.**Op * dens2 * mass2
    emis3['OII_'+str(line)] = emis1['OII_'+str(line)] + emis2['OII_'+str(line)]
        
emis1['Hbeta'] = pn.getHbEmissivity(tem1) * dens1 * mass1
emis2['Hbeta'] = pn.getHbEmissivity(tem2) * dens2 * mass2
emis3['Hbeta'] = emis1['Hbeta'] + emis2['Hbeta'] 

# determining Te and Ne for different regions and diagnostics:
# Region 1 only
temp_1, dens_S2_1 = pn.getCrossTemDen('N', 2, 'L(6583) / L(5755)', 'S', 2, 'L(6716)/L(6731)', 
                                      emis1['NII_6583']/emis1['NII_5755'], 
                                      emis1['SII_6716']/emis1['SII_6731'], 
                                      guess_tem=10000, tol_tem = 1., tol_den = 1., max_iter = 5)
# Region 2 only
temp_2, dens_S2_2 = pn.getCrossTemDen('N', 2, 'L(6583) / L(5755)','S', 2, 'L(6716)/L(6731)', 
                                      emis2['NII_6583']/emis2['NII_5755'], 
                                      emis2['SII_6716']/emis2['SII_6731'], 
                                      guess_tem=10000, tol_tem = 1., tol_den = 1., max_iter = 5)
# sum of region 1 and 2, using Te(NII) and Ne(SII)
temp_3, dens_S2_3 = pn.getCrossTemDen('N', 2, 'L(6583) / L(5755)', 'S', 2, 'L(6716)/L(6731)', 
                                      emis3['NII_6583']/emis3['NII_5755'], 
                                      emis3['SII_6716']/emis3['SII_6731'], 
                                      guess_tem=10000, tol_tem = 1., tol_den = 1., max_iter = 5)
# sum of region 1 and 2, using Te(NII) and Ne(SII 4 lines)
temp_4, dens_S2_4 = pn.getCrossTemDen('N', 2, 'L(6583) / L(5755)',
                                      'S', 2, '(L(4069) + L(4076)) / (L(6716)+L(6731))', 
                                      emis3['NII_6583']/emis3['NII_5755'],
                                      (emis3['SII_4069']+emis3['SII_4076']) / (emis3['SII_6716']+emis3['SII_6731']), 
                                      guess_tem=10000, tol_tem = 1., tol_den = 1., max_iter = 5)
# sum of region 1 and 2, using Te(OIII) and Ne(SII)
temp_5, dens_S2_5 = pn.getCrossTemDen('O', 3, 'L(5007) / L(4363)', 'S', 2, 'L(6716)/L(6731)', 
                                      emis3['OIII_5007']/emis3['OIII_4363'], 
                                      emis3['SII_6716']/emis3['SII_6731'], 
                                      guess_tem=10000, tol_tem = 1., tol_den = 1., max_iter = 5)

def ab_ion5(atom, line, wave):
    """
    this return ionic abundances in 5 cases: region 1, 2 and the sum of the 2 regions (using different Te and Ne)
    """
    ab_1 = np.log10(atom.getIonAbundance(emis1[line]/emis1['Hbeta']*100, temp_1, dens_S2_1, wave= wave))
    ab_2 = np.log10(atom.getIonAbundance(emis2[line]/emis2['Hbeta']*100, temp_2, dens_S2_2, wave= wave))
    ab_3 = np.log10(atom.getIonAbundance(emis3[line]/emis3['Hbeta']*100, temp_3, dens_S2_3, wave= wave))
    ab_4 = np.log10(atom.getIonAbundance(emis3[line]/emis3['Hbeta']*100, temp_4, dens_S2_4, wave= wave))
    ab_5 = np.log10(atom.getIonAbundance(emis3[line]/emis3['Hbeta']*100, temp_5, dens_S2_5, wave= wave))
    return ab_1, ab_2, ab_3, ab_4, ab_5

print('   Dens1 = %.0f, Dens2 = %.0f,  Dens3 = %.0f,  Dens4 = %.0f,  Dens5 = %.0f' % \
    (dens_S2_1, dens_S2_2, dens_S2_3, dens_S2_4, dens_S2_5))
print('   Temp1 = %.0f, Temp2 = %.0f, Temp3 = %.0f, Temp4 = %.0f, Temp5 = %.0f' % \
    (temp_1, temp_2, temp_3, temp_4, temp_5))

for line in S2_lambda:
    ab_Sp1, ab_Sp2, ab_Sp3, ab_Sp4, ab_Sp5 = ab_ion5(S2, 'SII_'+str(line), line)        
    print(' S+%i : %.2f          %.2f         %.2f            %.2f           %.2f' % \
        (line,ab_Sp1-Sp, ab_Sp2-Sp, ab_Sp3-Sp, ab_Sp4-Sp, ab_Sp5-Sp))

ab_Np1, ab_Np2, ab_Np3, ab_Np4, ab_Np5 = ab_ion5(N2, 'NII_6583', 6583)  
print(' N+%i : %.2f          %.2f         %.2f            %.2f           %.2f' % \
        (6583,ab_Np1-Np, ab_Np2-Np, ab_Np3-Np, ab_Np4-Np, ab_Np5-Np))

ab_Op1, ab_Op2, ab_Op3, ab_Op4, ab_Op5 = ab_ion5(O2, 'OII_3726', 3726)  
print(' O+%i : %.2f          %.2f         %.2f            %.2f           %.2f' % \
        (3726,ab_Op1-Op, ab_Op2-Op, ab_Op3-Op, ab_Op4-Op, ab_Op5-Op))
ab_Op1, ab_Op2, ab_Op3, ab_Op4, ab_Op5 = ab_ion5(O2, 'OII_3729', 3729)  
print(' O+%i : %.2f          %.2f         %.2f            %.2f           %.2f' % \
        (3729,ab_Op1-Op, ab_Op2-Op, ab_Op3-Op, ab_Op4-Op, ab_Op5-Op))

ab_Opp1, ab_Opp2, ab_Opp3, ab_Opp4, ab_Opp5 = ab_ion5(O3, 'OIII_5007', 5007)  
print('O++%i : %.2f          %.2f         %.2f            %.2f           %.2f' % \
        (5007,ab_Opp1-Opp, ab_Opp2-Opp, ab_Opp3-Opp, ab_Opp4-Opp, ab_Opp5-Opp))

