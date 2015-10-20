# Analysis of a simple two-component model, meant to illustrate the bias arising from assuming
# that the region is homogeneous in density
# First, an emission region made up of two different subregions is modelled,
# each with a different mass and density. The resulting overall emissivity is computed
# Second, the region is analyzed as if it were a homogeneous region

import numpy as np
import pyneb as pn

# Create atoms
O2 = pn.Atom('O', 2)
O3 = pn.Atom('O', 3)
N2 = pn.Atom('N', 2)
S2 = pn.Atom('S', 2)
S3 = pn.Atom('S', 3)

# Set abundances
O_H = 1e-4 # O/H
N_H = 1e-4
S_H = 1e-4

O2_1 = 0.3 * O_H # O+/O * O/H
O3_1 = 0.7 * O_H
N2_1 = 0.3 * N_H
S2_1 = 0.3 * S_H
S3_1 = 0.7 * S_H

O2_2 = 0.9 * O_H
O3_2 = 0.1 * O_H
N2_2 = 0.9 * N_H
S2_2 = 0.9 * S_H
S3_2 = 0.1 * S_H

# Set properties of two regions
tem1 = 1.e4
tem2 = 1.2e4
den1 = 1.e2
den2 = 5.e5
vol1 = 5.e8
vol2 = 1.

# These factors are necessary to compute luminosities from specific emissivities
fac1 = vol1 * den1**2
fac2 = vol2 * den2**2

# Compute emissivi
I1_3729 = O2.getEmissivity(tem1, den1, wave = 3729) * O2_1 * fac1
I2_3729 = O2.getEmissivity(tem2, den2, wave = 3729) * O2_2 * fac2
I_3729 = I1_3729 + I2_3729
I1_3726 = O2.getEmissivity(tem1, den1, wave = 3726) * O2_1 * fac1
I2_3726 = O2.getEmissivity(tem2, den2, wave = 3726) * O2_2 * fac2
I_3726 = I1_3726 + I2_3726
I1_7325 = (O2.getEmissivity(tem1, den1, wave = 7318) + O2.getEmissivity(tem1, den1, wave = 7319) +
           O2.getEmissivity(tem1, den1, wave = 7329) + O2.getEmissivity(tem1, den1, wave = 7330)) * O2_1 * fac1
I2_7325 = (O2.getEmissivity(tem2, den2, wave = 7318) + O2.getEmissivity(tem2, den2, wave = 7319) +
           O2.getEmissivity(tem2, den2, wave = 7329) + O2.getEmissivity(tem2, den2, wave = 7330)) * O2_2 * fac2
I_7325 = I1_7325 + I2_7325

I1_5007 =  O3.getEmissivity(tem1, den1, wave = 5007) * O3_1 * fac1
I2_5007 =  O3.getEmissivity(tem2, den2, wave = 5007) * O3_2 * fac2
I_5007 = I1_5007 + I2_5007
I1_4363 =  O3.getEmissivity(tem1, den1, wave = 4363) * O3_1 * fac1
I2_4363 =  O3.getEmissivity(tem2, den2, wave = 4363) * O3_2 * fac2
I_4363 = I1_4363 + I2_4363

I1_6584 = N2.getEmissivity(tem1, den1, wave = 6584) * N2_1 * fac1
I2_6584 = N2.getEmissivity(tem2, den2, wave = 6584) * N2_2 * fac2
I_6584 = I1_6584 + I2_6584
I1_5755 = N2.getEmissivity(tem1, den1, wave = 5755) * N2_1 * fac1
I2_5755 = N2.getEmissivity(tem2, den2, wave = 5755) * N2_2 * fac2
I_5755 = I1_5755 + I2_5755

I1_6725 = (S2.getEmissivity(tem1, den1, wave = 6717) + S2.getEmissivity(tem1, den1, wave = 6731)) * S2_1 * fac1
I2_6725 = (S2.getEmissivity(tem2, den2, wave = 6717) + S2.getEmissivity(tem2, den2, wave = 6731)) * S2_2 * fac2
I_6725 = I1_6725 + I2_6725
I1_4072 = S2.getEmissivity(tem1, den1, wave = 4072) * S2_1 * fac1
I2_4072 = S2.getEmissivity(tem2, den2, wave = 4072) * S2_2 * fac2
I_4072 = I1_4072 + I2_4072

I1_9069 = S3.getEmissivity(tem1, den1, wave = 9069) * S3_1 * fac1
I2_9069 = S3.getEmissivity(tem2, den2, wave = 9069) * S3_2 * fac2
I_9069 = I1_9069 + I2_9069
I1_6312 = S3.getEmissivity(tem1, den1, wave = 6312) * S3_1 * fac1
I2_6312 = S3.getEmissivity(tem2, den2, wave = 6312) * S3_2 * fac2
I_6312 = I1_6312 + I2_6312

I1_Beta = pn.getRecEmissivity(tem1, den1, 4, 2, 'H1') * fac1
I2_Beta = pn.getRecEmissivity(tem2, den2, 4, 2, 'H1') * fac2
I_Beta = I1_Beta + I2_Beta

# A partir de aqui suponemos que tenemos observaciones:
den_diag = 'L(3729)/L(3726)'
tem_diag_N2 = 'L(5755) / L(6584)'
tem_diag_O2 = '(L(3726)+L(3729)) / (L(7318) + L(7319) + L(7329) + L(7330))'
tem_diag_O3 = 'L(4363) / L(5007)'
tem_diag_S2 = 'L(4072) / (L(6717) + L(6731))'
tem_diag_S3 = 'L(6312) / L(9069)'

Ab_Op_diag = 'L(3729)+L(3726)'
Ab_Opp_diag = 'L(5007)'

computed_den_1 = O2.getTemDen(I1_3729 / I1_3726, tem = 1e4, to_eval = den_diag)
computed_den_2 = O2.getTemDen(I2_3729 / I2_3726, tem = 1e4, to_eval = den_diag)
computed_den = O2.getTemDen(I_3729 / I_3726, tem = 1e4, to_eval = den_diag)

computed_tem_O2_1 = O2.getTemDen((I1_3729 + I1_3726) / I1_7325, den = computed_den_1, to_eval = tem_diag_O2)
computed_tem_O2_2 = O2.getTemDen((I2_3729 + I2_3726) / I2_7325, den = computed_den_2, to_eval = tem_diag_O2)
computed_tem_O2 = O2.getTemDen((I_3729 + I_3726) / I_7325, den = computed_den, to_eval = tem_diag_O2)

computed_tem_O3_1 = O3.getTemDen(I1_4363 / I1_5007, den = computed_den_1, to_eval = tem_diag_O3)
computed_tem_O3_2 = O3.getTemDen(I2_4363 / I2_5007, den = computed_den_2, to_eval = tem_diag_O3)
computed_tem_O3 = O3.getTemDen(I_4363 / I_5007, den = computed_den, to_eval = tem_diag_O3)

computed_tem_N2_1 = N2.getTemDen(I1_5755 / I1_6584, den = computed_den_1, to_eval = tem_diag_N2)
computed_tem_N2_2 = N2.getTemDen(I2_5755 / I2_6584, den = computed_den_2, to_eval = tem_diag_N2)
computed_tem_N2 = N2.getTemDen(I_5755 / I_6584, den = computed_den, to_eval = tem_diag_N2)

computed_tem_S2_1 = S2.getTemDen(I1_4072 / I1_6725, den = computed_den_1, to_eval = tem_diag_S2)
computed_tem_S2_2 = S2.getTemDen(I2_4072 / I2_6725, den = computed_den_2, to_eval = tem_diag_S2)
computed_tem_S2 = S2.getTemDen(I_4072 / I_6725, den = computed_den, to_eval = tem_diag_S2)

computed_tem_S3_1 = S3.getTemDen(I1_6312 / I1_9069, den = computed_den_1, to_eval = tem_diag_S3)
computed_tem_S3_2 = S3.getTemDen(I2_6312 / I2_9069, den = computed_den_2, to_eval = tem_diag_S3)
computed_tem_S3 = S3.getTemDen(I_6312 / I_9069, den = computed_den, to_eval = tem_diag_S3)

computed_O2_1 = np.log10(O2.getIonAbundance((I1_3729 + I1_3726) / I1_Beta * 100, computed_tem_N2_1, computed_den_1, to_eval = Ab_Op_diag))
computed_O2_2 = np.log10(O2.getIonAbundance((I2_3729 + I2_3726) / I2_Beta * 100, computed_tem_N2_2, computed_den_2, to_eval = Ab_Op_diag))
computed_O2 = np.log10(O2.getIonAbundance((I_3729 + I_3726) / I_Beta * 100, computed_tem_N2, computed_den, to_eval = Ab_Op_diag))

computed_O3_1 = np.log10(O3.getIonAbundance(I1_5007 / I1_Beta * 100, computed_tem_O3_1, computed_den_1, to_eval = Ab_Opp_diag))
computed_O3_2 = np.log10(O3.getIonAbundance(I2_5007 / I2_Beta * 100, computed_tem_O3_2, computed_den_2, to_eval = Ab_Opp_diag))
computed_O3 = np.log10(O3.getIonAbundance(I_5007 / I_Beta * 100, computed_tem_O3, computed_den, to_eval = Ab_Opp_diag))

print('\n          Input    Computed')
print('Ne1      {0:6.2e}   {1:6.2e}'.format(den1, computed_den_1))
print('Ne2      {0:6.2e}   {1:6.2e}'.format(den2, computed_den_2))
print('Ne                  {0:6.2e}'.format(computed_den))

print('T1(NII)  {0:6.2e}   {1:6.2e}'.format(tem1, computed_tem_N2_1))
print('T2(NII)  {0:6.2e}   {1:6.2e}'.format(tem2, computed_tem_N2_2))
print('T(NII)              {0:6.2e}'.format(computed_tem_N2))

print('T1(OII)  {0:6.2e}   {1:6.2e}'.format(tem1, computed_tem_O2_1))
print('T2(OII)  {0:6.2e}   {1:6.2e}'.format(tem2, computed_tem_O2_2))
print('T(OII)              {0:6.2e}'.format(computed_tem_O2))

print('T1(OIII) {0:6.2e}   {1:6.2e}'.format(tem1, computed_tem_O3_1))
print('T2(OIII) {0:6.2e}   {1:6.2e}'.format(tem2, computed_tem_O3_2))
print('T(OIII)             {0:6.2e}'.format(computed_tem_O3))

print('T1(SII)  {0:6.2e}   {1:6.2e}'.format(tem1, computed_tem_S2_1))
print('T2(SII)  {0:6.2e}   {1:6.2e}'.format(tem2, computed_tem_S2_2))
print('T(SII)              {0:6.2e}'.format(computed_tem_S2))

print('T1(SIII) {0:6.2e}   {1:6.2e}'.format(tem1, computed_tem_S3_1))
print('T2(SIII) {0:6.2e}   {1:6.2e}'.format(tem2, computed_tem_S3_2))
print('T(SIII)             {0:6.2e}'.format(computed_tem_S3))

print('O+1     {0:6.2f}     {1:6.2f}'.format(np.log10(O2_1), computed_O2_1))
print('O+2     {0:6.2f}     {1:6.2f}'.format(np.log10(O2_2), computed_O2_2))
print('O+                 {0:6.2f}'.format(computed_O2))

print('O++1    {0:6.2f}     {1:6.2f}'.format(np.log10(O3_1), computed_O3_1))
print('O++2    {0:6.2f}     {1:6.2f}'.format(np.log10(O3_2), computed_O3_2))
print('O++                {0:6.2f}'.format(computed_O3))

print('O/H1    {0:6.2f}     {1:6.2f}'.format(np.log10(O2_1 + O3_1), np.log10(10**computed_O2_1 + 10**computed_O3_1)))
print('O/H2    {0:6.2f}     {1:6.2f}'.format(np.log10(O2_2 + O3_2), np.log10(10**computed_O2_1 + 10**computed_O3_1)))
print('O/H                {0:6.2f}'.format(np.log10(10**computed_O2 + 10**computed_O3)))

print("""
We consider 2 regions with different properties. Both have the same O/H = {:.0e}.
The first one has a density of {:.0f} cm-3, an electron temperature of {:.0f} K, an O+/H of {:.2f}, an O++/H of {:.2f} and occupies almost {:.0f}% of the volume.
The second one has a density of {:.0f} cm-3, an electron temperature of {:.0f} K, an O+/H of {:.2f}, an O++/H of {:.2f} and occupies {:.1e}% of the volume.
If we use the [OII] and [OIII] lines to determine the density and temperatures of the whole region (summing the emissions coming from the 2 regions), we
obtain a density of {:.0f} cm-3 and an electron temperature of {:.0f} and {:.0f} K for N+ and O++ respectively.
The determination of O+/H and O++/H are then {:.2f} and {:.2f}, leading to O/H={:.2f} instead of O/H={:.2f}.
""".format(O_H, den1, tem1, np.log10(O2_1), np.log10(O3_1), vol1/(vol1+vol2)*100, den2, tem2, np.log10(O2_2), np.log10(O3_2),
           vol2/(vol1+vol2)*100, computed_den, computed_tem_N2, computed_tem_O3, computed_O2,
           computed_O3, np.log10(10**computed_O2 + 10**computed_O3), -4 ))

