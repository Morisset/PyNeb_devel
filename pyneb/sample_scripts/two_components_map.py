# Simple script which computes the temperature inferred from the line intensities emitted 
# by two regions with different densities and volumes and interpreted as if it were a 
# unique region

# Imports
import numpy as np
import pyneb as pn
import matplotlib.pyplot as plt

# Creation of the necessary atoms
O2 = pn.Atom('O', 2)
O3 = pn.Atom('O', 3)
S2 = pn.Atom('S', 2)
S3 = pn.Atom('S', 3)

# Function to compute temperatures inferred from the sum of the emissions from two different regions
def TOandTS(tem1=1.e4, tem2=1e4, den1=1e2, den2=5e5, vol1=5e8, denDiag='S2', verbose=False):
	# Set total abundances of O, N, and S
    OsH = 1e-4 # O/H
    NsH = 1e-4
    SsH = 1e-4
    # Set ionic abundances of O, N, and S
    Op1 = 0.3 * OsH # O+/O * O/H
    Opp1 = 0.7 * OsH
    Np1 = 0.3 * NsH
    Sp1 = 0.3 * SsH
    Spp1 = 0.7 * SsH
    
    Op2 = 0.9 * OsH
    Opp2 = 0.1 * OsH
    Np2 = 0.9 * NsH
    Sp2 = 0.9 * SsH
    Spp2 = 0.1 * SsH
	
	# Volume of one region (the other is a parameter)
    vol2 = 1.
    
    # Multiplying factors used to compute emissivities
    fac1 = vol1 * den1**2
    fac2 = vol2 * den2**2
    
	# Intensities emitted by each region individually and by the two of them together
    I1_5007 =  O3.getEmissivity(tem1, den1, wave = 5007) * Opp1 * fac1
    I2_5007 =  O3.getEmissivity(tem2, den2, wave = 5007) * Opp2 * fac2
    I_5007 = I1_5007 + I2_5007
    I1_4363 =  O3.getEmissivity(tem1, den1, wave = 4363) * Opp1 * fac1
    I2_4363 =  O3.getEmissivity(tem2, den2, wave = 4363) * Opp2 * fac2
    I_4363 = I1_4363 + I2_4363
            
    I1_9069 = S3.getEmissivity(tem1, den1, wave = 9069) * Spp1 * fac1
    I2_9069 = S3.getEmissivity(tem2, den2, wave = 9069) * Spp2 * fac2
    I_9069 = I1_9069 + I2_9069
    I1_6312 = S3.getEmissivity(tem1, den1, wave = 6312) * Spp1 * fac1
    I2_6312 = S3.getEmissivity(tem2, den2, wave = 6312) * Spp2 * fac2
    I_6312 = I1_6312 + I2_6312
        
    # The following bit simulates an observational analysis
	# The result depends on the density diagnostic assumed, so in each case 
	# it is stored in different variables
    if denDiag == 'O2':
        I1_3729 = O2.getEmissivity(tem1, den1, wave = 3729) * Op1 * fac1
        I2_3729 = O2.getEmissivity(tem2, den2, wave = 3729) * Op2 * fac2
        I_3729 = I1_3729 + I2_3729
        I1_3726 = O2.getEmissivity(tem1, den1, wave = 3726) * Op1 * fac1
        I2_3726 = O2.getEmissivity(tem2, den2, wave = 3726) * Op2 * fac2
        I_3726 = I1_3726 + I2_3726
        Ne_diag = 'L(3729)/L(3726)'    
        d_Ne = O2.getTemDen(I_3729 / I_3726, tem = 1e4, to_eval = Ne_diag)
    elif denDiag == 'S2':
        I1_6717 = (S2.getEmissivity(tem1, den1, wave = 6717)) * Sp1 * fac1
        I2_6717 = (S2.getEmissivity(tem2, den2, wave = 6717)) * Sp2 * fac2
        I_6717 = I1_6717 + I2_6717
        I1_6731 = (S2.getEmissivity(tem1, den1, wave = 6731)) * Sp1 * fac1
        I2_6731 = (S2.getEmissivity(tem2, den2, wave = 6731)) * Sp2 * fac2
        I_6731 = I1_6731 + I2_6731
        Ne_diag = 'L(6717)/L(6731)'
        d_Ne = S2.getTemDen(I_6717 / I_6731, tem = 1e4, to_eval = Ne_diag)
    else:
        pn.log_.error('unknown denDiag: {0}'.format(denDiag), calling = 'two_components_map')
        return None
    
    # Temperature derived from observed intensities in OIII lines
    T_Opp_diag = 'L(4363) / L(5007)'
    d_T_Opp = O3.getTemDen(I_4363 / I_5007, den = d_Ne, to_eval = T_Opp_diag)
        
    # Temperature derived from observed intensities in SIII lines
    T_Spp_diag = 'L(6312) / L(9069)'
    d_T_Spp = S3.getTemDen(I_6312 / I_9069, den = d_Ne, to_eval = T_Spp_diag)
    
    if verbose:
        print('         Input    Determinado')
        print('Ne                  {0:6.2e}'.format(d_Ne))
        print('T(OIII)             {0:6.2e}'.format(d_T_Opp))        
        print('T(SIII)             {0:6.2e}'.format(d_T_Spp))

    return d_T_Opp, d_T_Spp

# Defines a grid in log(den), log(vol). 
# The T diagnostics will be computed for each grid point. 
n_den_points = 40
n_vol_points = 30
den2grid = np.logspace(3, 7, n_den_points)
vol1grid = np.logspace(3, 9, n_vol_points)

T_Opp = np.zeros ((n_vol_points, n_den_points)) 
T_Spp = np.zeros ((n_vol_points, n_den_points))


# Compute T(O III) and T(S III) for all the points of the grid. The density and volume
# can be changed
for i, vol1 in enumerate(vol1grid):
    T_Opp[i], T_Spp[i] = TOandTS(den1=1e2, den2=den2grid, vol1=vol1, denDiag='O2')
    print 'Computing row n.', i, 'of', n_vol_points

plt.figure(1)
plt.subplot(2, 2, 1)
plt.imshow(T_Opp, interpolation = 'none', origin = 'lower', vmin = 0, vmax = 20000,
           extent=[np.log10(min(den2grid)), np.log10(max(den2grid)), np.log10(min(vol1grid)), np.log10(max(vol1grid))], aspect = 'auto')
plt.xlabel('log(den2)')
plt.ylabel('log(vol1)')
plt.title('T(OIII)')
plt.colorbar()

plt.subplot(2, 2, 2)
plt.imshow(T_Spp, interpolation = 'none', origin = 'lower', vmin = 0, vmax = 20000,
           extent=[np.log10(min(den2grid)), np.log10(max(den2grid)), np.log10(min(vol1grid)), np.log10(max(vol1grid))], aspect = 'auto')
plt.xlabel('log(den2)')
plt.ylabel('log(vol1)')
plt.title('T(S III)')
plt.colorbar()

plt.subplot(2, 2, 3)
plt.imshow(T_Opp/T_Spp, interpolation = 'none', origin = 'lower',vmin = 0.9, vmax = 1.1,
           extent=[np.log10(min(den2grid)), np.log10(max(den2grid)), np.log10(min(vol1grid)), np.log10(max(vol1grid))], aspect = 'auto')
plt.xlabel('log(den2)')
plt.ylabel('log(vol1)')
plt.title('T(OIII) / T(SIII)')
plt.colorbar()

plt.show()
