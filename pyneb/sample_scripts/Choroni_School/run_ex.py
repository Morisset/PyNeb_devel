import matplotlib.pyplot as plt
import numpy as np
import pyCloudy as pc
import pyneb as pn

#-------------------------------------- ex1_1 questions 1 to 8 ----------------------------------------------
"""
import ex1_1

ex1_1.p1()  # answer to 1.1.1

# answers to 1.1.2 and 1.2.3
plt.figure(1) 	# open window 1
plt.show() 		# to see the windows containing the figures if not working from an xterm
ex1_1.p23()		# calls p23 with default options
ex1_1.p23(den=1e5, legend=False, style='--') 	# calls p23 with specified options and superimposes the new emissivity plot
plt.savefig('OIII_emissivities.pdf') # saves the figure in a pdf file
ex1_1.p3()

ex1_1.p4()

plt.figure(2) 	# open window 2
ex1_1.p5()		# plots the nominal line emissivities for den=1e2
ex1_1.p5(legend=False, style='--', coeff=5)	# plots the line emissivities after changing the A values by a factor of 5

plt.figure(3) 	# open window 2
ex1_1.p5(den=1e5)		# plots the nominal line emissivities for den=1e5
ex1_1.p5(den=1e5, legend=False, style='--', coeff=5)	# plots the line emissivities after changing the A values by a factor of 5

plt.figure(4) 	# open window 4
ex1_1.p6()
ex1_1.p6(legend=False, style='--', coeff=5) # plots the line emissivities after changing the omegas by a factor of 5

plt.figure(5) 	# open window 5
ex1_1.p6(den=1e5)		# plots the nominal line emissivities for den=1e5)
ex1_1.p6(den=1e5, legend=False, style='--', coeff=5) # plots the line emissivities after changing the omegas by a factor of 5

plt.figure(6) 	# open window 6
ex1_1.p7a()

ex1_1.p7b()

plt.figure(8) 	# open window 8
ex1_1.p7c()
ex1_1.p7c(OmegaInterp='Linear', legend=False, style='--')

plt.figure(9)     # open window 8
ex1_1.p8(elem='S')
ex1_1.p8(elem='S', OmegaInterp='Linear', legend=False, style='--')

pn.atomicData.setDataFile('o_iii_coll_Pal12-AK99.fits')
ex1_1.p8(elem='O')
ex1_1.p8(elem='O', OmegaInterp='Linear', legend=False, style='--')
#-------------------------------------- ex1_2 figura 1 ----------------------------------------------
import ex1_2
atom = pn.Atom('O', 3, OmegaInterp='Linear')

plt.figure(1, figsize=(10,8))
ex1_2.plot_Gain()
plt.xlabel('T [K]')
plt.ylabel(r'log energy rate/(n$_e$ n$_p$) [erg cm$^3$ s$^{-1}$]') # 
plt.legend(loc=4)

ex1_2.plot_LossH()
plt.legend(loc=4)

plt.figure(3, figsize=(10,8))
ex1_2.plot_Gain()
plt.xlabel('T [K]')
plt.ylabel(r'log energy rate/(n$_e$ n$_p$) [erg cm$^3$ s$^{-1}$]') 

ex1_2.plot_LossO(atom)
plt.legend(loc=4)
plt.ylim((-40, -20))
plt.savefig('GainLoss.pdf')

plt.figure(4, figsize=(10,8))
ex1_2.plot_Gain()
ex1_2.plot_Loss(atom)
plt.ylim((-25, -22))

# this will help for question 5:
Teq = 9145
plt.plot([Teq, Teq], [-25, -22], linestyle=':') # plots a vertical dotted line

#-------------------------------------- ex1_2 figura 2 ----------------------------------------------
#*************** a faire tourner pour 1.2.5
plt.figure(5, figsize=(12, 10))
ex1_2.plot_Gain()
ex1_2.plot_Loss(atom)
plt.legend(loc=4)
plt.ylim((-25, -21))
plt.plot((Teq, Teq), (-25, -21), color='black')
plt.plot((Teq - 2000, Teq - 2000), (-25, -21), color='black')
plt.plot((Teq + 2000, Teq + 2000), (-25, -21), color='black')
ex1_2.plot_Loss(atom, OoH=2e-4, onlyTotal=True, linestyle='--')
ex1_2.plot_Loss(atom, OoH=8e-4, onlyTotal=True, linestyle=':')
ex1_2.plot_Gain(30000, linestyle='--')
ex1_2.plot_Gain(80000, linestyle=':')
plt.savefig('EquTermO3.pdf')

atom = pn.Atom('O', 2, OmegaInterp='Linear')
Teq = 11000
plt.figure(7, figsize=(12, 10))
ex1_2.plot_Gain()
ex1_2.plot_Loss(atom)
plt.legend(loc=4)
plt.ylim((-25, -21))
plt.plot((Teq, Teq), (-25, -21), color='black')
plt.plot((Teq - 2000, Teq - 2000), (-25, -21), color='black')
plt.plot((Teq + 2000, Teq + 2000), (-25, -21), color='black')
ex1_2.plot_Loss(atom, OoH=2e-4, onlyTotal=True, linestyle='--')
ex1_2.plot_Loss(atom, OoH=1.1e-3, onlyTotal=True, linestyle=':')
ex1_2.plot_Gain(25000, linestyle='--')
ex1_2.plot_Gain(90000, linestyle=':')
plt.savefig('EquTermO2.pdf')

#-------------------------------------- ex2_1 ----------------------------------------------
import ex2_1

"""
# question 1 #
"""
ions = ['C3', 'N1','N2', 'O1','O2', 'O3', 'Ne3', 'S2', 'S3','Cl3', 'Ar4']

for ion in ions:
    ex2_1.p1(ion)
    plt.savefig('{0}_energies.pdf'.format(ion))

"""
# question 2
"""
# print the list of available plasma diagnostics  
for diag in np.sort(pn.diags_dict.keys()):
    print('{0} -> {1}'.format(diag, pn.diags_dict[diag]))

# select the ones you want:

diags = ['[CIII] 1909/1907', '[NII] 5755/6584', '[OII] 3726/3729', '[OII] 3727+/7325+', '[OIII] 4363/5007',
         '[OIII] 51m/88m', '[OIII] 5007/88m', '[OIII] 1666/5007', '[NeIII] 15.6m/36.0m', '[NeIII] 3869/15.6m',
         '[SII] 6731/6716', '[SII] 4072+/6720+', '[SIII] 18.7m/33.6m', '[SIII] 9069/18.7m', '[SIII] 6312/9069',
         '[ClIII] 5538/5518', '[ArIV] 4740/4711']
for diag in diags:
    ex2_1.p2(diag)


#-------------------------------------- ex2_2 ----------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
import pyCloudy as pc
import pyneb as pn
import ex2_2

ex2_2.p1()
ex2_2.p2()
ex2_2.p3()
# changing the atomic data
pn.atomicData.setDataFileDict('IRAF_09')
ex2_2.p3(fignum=4, pypic_path='/tmp/pypics_IRAF/')
pn.atomicData.resetDataFileDict()

#-------------------------------------- ex3_1 ----------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
import pyCloudy as pc
import pyneb as pn
import ex3_1

"""
# Question 1
"""
Ne_O2 = 10**3.4
Ne_Ar4 = 10**3.8
T_N2 = 13400.
T_O3 = 14500.

"""
# Question 2
"""
ion_ab_dic = ex3_1.getIonAb('IC2165.dat', Ne_O2, Ne_Ar4, T_N2, T_O3, printIonAb = True)

"""
# Question 3
"""
elem_abun = ex3_1.getElemAb(ion_ab_dic, printAb = True)

"""
# Question 4
"""
ion_ab_dic_UVIR = ex3_1.getIonAb('IC2165_IR.dat', Ne_O2, Ne_Ar4, T_N2, T_O3, printIonAb = True)
elem_abun_UVIR = ex3_1.getElemAb(ion_ab_dic_UVIR, printAb = True)
"""
# Question 5
"""
pn.atomicData.setDataFileDict('IRAF_09')
ion_ab_dic_IRAF = ex3_1.getIonAb('IC2165_IR.dat', Ne_O2, Ne_Ar4, T_N2, T_O3, printIonAb = True)
elem_abun_IRAF = ex3_1.getElemAb(ion_ab_dic_IRAF, printAb = True)
pn.atomicData.resetDataFileDict()
for icf in elem_abun_UVIR:
    print('{0}: {1:.3f} {2:.3f}'.format(icf, np.log10(elem_abun_UVIR[icf])+12, np.log10(elem_abun_IRAF[icf])+12))
#-------------------------------------- ex3_2 ----------------------------------------------
import ex3_2
import atpy
import asciitable

#Read the observations. Notice that in this file the errors are absolutes
obs = pn.Observation('NGC300.dat', corrected=True, errIsRelative=False)
# the galactocentric distance has been also read and the values are in "obsIntens"
Rgal = obs.getLine(label='DIST').obsIntens
# Compute the densities and temperatures
mean_dens, temp_O2, temp_S2, temp_N2, temp_O3, temp_S3 = ex3_2.p1(obs)

# Instantiate a table to write the results in anm ascii file
T = atpy.Table()
# Adding the columns with their names into T
T.add_column('NAME', obs.names)
T.add_column('T_O2', temp_O2)
T.add_column('T_S2', temp_S2)
T.add_column('T_N2', temp_N2)
T.add_column('T_O3', temp_O3)
# writing to a file
T.write('temperatures.ascii', type='ascii',
        formats={'NAME':'%5s', 'T_O2': '%7.1f', 'T_N2': '%7.1f', 'T_S2':'%7.1f', 'T_O3':'%7.1f'},
        overwrite=True, Writer=asciitable.FixedWidth, delimiter=None)

T.write('temperatures.tex', type='latex',
        formats={'NAME':'%5s', 'T_O2': '%7.1f', 'T_N2': '%7.1f', 'T_S2':'%7.1f', 'T_O3':'%7.1f'},
        overwrite=True)

# compute Te for the low, mid and high ionisation regions
Tlow, Tmid, Thigh = ex3_2.p2(temp_S3, temp_O3)
# Compute the ionic abundances
ab_dic = ex3_2.p3(obs, Tlow, Tmid, Thigh, mean_dens)

OoH = ab_dic['O2'] + ab_dic['O3']
NoH = ab_dic['N2'] * OoH / ab_dic['O2']
NeoH = ab_dic['Ne3'] * OoH / ab_dic['O3']

Tlow2, Thigh2 = ex3_2.p2_P(temp_N2, temp_O3)
ab_dic2 = ex3_2.p3(obs, Tlow2, Tmid, Thigh2, mean_dens)

OoH2 = ab_dic2['O2'] + ab_dic2['O3']
NoH2 = ab_dic2['N2'] * OoH2 / ab_dic2['O2']
NeoH2 = ab_dic2['Ne3'] * OoH2 / ab_dic2['O3']

plt.scatter(Rgal, np.log10(OoH) + 12, label = 'O/H', color='red')
plt.scatter(Rgal, np.log10(NoH) + 12, label = 'N/H', color='blue')
plt.scatter(Rgal, np.log10(NeoH) + 12, label = 'Ne/H', color='green')
plt.scatter(Rgal, np.log10(OoH2) + 12, label = 'O/H', color='red', alpha=.3)
plt.scatter(Rgal, np.log10(NoH2) + 12, label = 'N/H', color='blue', alpha=.3)
plt.scatter(Rgal, np.log10(NeoH2) + 12, label = 'Ne/H', color='green', alpha=.3)

plt.legend()
plt.xlabel('RGal')
plt.ylabel('12+log(X/H)')

# Print all the atomic data used:
for ion in np.sort(pn.config.DataFiles.keys()):
    print('{0} -> {1}'.format(ion, pn.config.DataFiles[ion]))

#-------------------------------------- ex4 ------------------------------------------------
"""
# STRONG LINE METHOD
"""

import ex4_1

OoH_2 = ex4_1.get_OoH_2()


OoH_m = ex4_1.get_OoH_2('-')
OoH_p = ex4_1.get_OoH_2('+')

plt.scatter(np.log10(OoH)+12, OoH_2)
plt.xlim((8, 9.2))
plt.ylim((8, 9.2))

OoH_3 = ex4_1.get_OoH_3()
plt.scatter(np.log10(OoH)+12, OoH_3, marker = '+')

OoH_4 = ex4_1.get_OoH_4(Tlow, Thigh, mean_dens)
plt.scatter(np.log10(OoH)+12, OoH_4, marker = 'd')
#-------------------------------------- ex5 ------------------------------------------------
"""
# CLOUDY FIRST STEPS
"""
#-------------------------------------- ex6_1 ----------------------------------------------
"""
# CLOUDY WITHOUT PYCLOUDY
"""
import ex6_1
import pyCloudy as pc

"""
# Question 6.1.1
"""
models_dir = 'CloudyModels' #Create this directory if necessary
pc.config.cloudy_exe = '/usr/local/Cloudy/c10.00/cloudy.exe' # point to the location of Cloudy.exe

# Print the Makefile in th emodels_dir directory (only needed one time)
pc.print_make_file(models_dir)
# Print a file containing the line for intensities
ex6_1.print_line_file(models_dir)

# prepare models with different Q(H)
for qH in [43, 44, 45, 46, 47, 48]:
    ex6_1.make_mod(models_dir, name='M_61.1_{0}'.format(qH), Teff=5e4, qH=qH, dens=3, r_in=17)
# run all the models in models_dirs starting with M_61.1
pc.run_cloudy(dir_=models_dir, n_proc=3, use_make=True,  model_name='M_61.1')
M43 = pc.CloudyModel('{0}/M_61.1_43'.format(models_dir), read_lin = True)
# load all the models in the models list
models = pc.load_models('{0}/M_61.1_'.format(models_dir), read_lin = True)

# make some plots
ex6_1.plot_Hb(models)

for m in models:
    # A posteriori cut of the model
    # Line intensities are not changed. Integral of Emissivities are changed.
    m.r_out_cut = 1.5e17
# make the same plots
ex6_1.plot_Hb(models)

# prepare models with different Teff
for i, Teff in enumerate(np.array([40, 50, 80, 130, 200]) * 1e3):
    ex6_1.make_mod(models_dir, name='M_61.1b_{0}'.format(i), Teff=Teff, qH=47, dens=3, r_in=17)
pc.run_cloudy(dir_=models_dir, n_proc=3, use_make=True,  model_name='M_61.1b_')
models = pc.load_models('{0}/M_61.1b_'.format(models_dir), read_lin = True)
ex6_1.plot_Hb_Teff(models)

# Another way is to ask every student to compute a model and to derive Hbeta/Q0 and we see the results.

"""
# Question 6.1.2
"""
# make some models changing input parameters
ex6_1.make_mod(models_dir, name='M_61.2_1', Teff=5e4, qH=45, dens=3, r_in=17)
ex6_1.make_mod(models_dir, name='M_61.2_2', Teff=8e4, qH=45, dens=3, r_in=17)
ex6_1.make_mod(models_dir, name='M_61.2_3', Teff=5e4, qH=47, dens=3, r_in=17)
ex6_1.make_mod(models_dir, name='M_61.2_4', Teff=5e4, qH=45, dens=5, r_in=17)
ex6_1.make_mod(models_dir, name='M_61.2_5', Teff=5e4, qH=45, dens=3, r_in=18)
# run the models
pc.run_cloudy(dir_=models_dir, n_proc=3, use_make=True,  model_name='M_61.2')
# load the models
models = pc.load_models('{0}/M_61.2_'.format(models_dir), read_lin = True)

ex6_1.plot_Hb(models)

# make some plots
plt.figure(figsize=(9,9))
ex6_1.plot_ion(models, elem='O', ion=2)
plt.figure(figsize=(9,9))
ex6_1.plot_ion(models, elem='O', ion=1)

"""
# Question 6.1.3 
"""
plt.figure(figsize=(9,9))
ex6_1.plot_radial2(models[0], 'red')
ex6_1.plot_radial2(models[1], 'green')
ex6_1.plot_radial2(models[2], 'blue')
ex6_1.plot_radial2(models[3], 'cyan')
ex6_1.plot_radial2(models[4], 'magenta')
ex6_1.plot_radial2(models[5], 'yellow')

"""
# Question 6.1.4
"""
qH = 48.23
r_in = 16
ex6_1.make_mod(models_dir, name='M_61.4_1', Teff=5e4, qH=qH, dens=2, r_in=r_in)
ex6_1.make_mod(models_dir, name='M_61.4_2', Teff=5e4, qH=qH, dens=2, r_in=r_in, metals = 0.01)
ex6_1.make_mod(models_dir, name='M_61.4_3', Teff=5e4, qH=qH, dens=2, r_in=r_in, metals = 0.1)
# Take care of the stopping criteria!!
ex6_1.make_mod(models_dir, name='M_61.4_4', Teff=5e4, qH=qH, dens=2, r_in=r_in, metals = 3.0)
pc.run_cloudy(dir_=models_dir, n_proc=3, use_make=True,  model_name='M_61.4')
models = pc.load_models('{0}/M_61.4_'.format(models_dir), read_lin = True)
plt.figure(figsize=(9,9))
ex6_1.plot_metals(models)

"""
# Question 6.1.5
"""
ex6_1.make_mod(models_dir, name='M_61.5_1', Teff=5e4, qH=qH, dens=2, r_in=r_in, nograins=False)
ex6_1.make_mod(models_dir, name='M_61.5_2', Teff=5e4, qH=qH, dens=2, r_in=r_in, metals = 0.01, nograins=False)
ex6_1.make_mod(models_dir, name='M_61.5_3', Teff=5e4, qH=qH, dens=2, r_in=r_in, metals = 0.1, nograins=False)
# Take care of the stopping criteria!!
ex6_1.make_mod(models_dir, name='M_61.5_4', Teff=5e4, qH=qH, dens=2, r_in=r_in, metals = 3.0, nograins=False)
pc.run_cloudy(dir_=models_dir, n_proc=3, use_make=True,  model_name='M_61.5')
modelsG = pc.load_models('{0}/M_61.5_'.format(models_dir), read_lin = True)
plt.figure(figsize=(9,9))
ex6_1.plot_metals(modelsG)


#-------------------------------------- ex6_2 ----------------------------------------------
import ex6_2

"""
# Pregunta 1
"""
print(ex6_2.QH0(U_mean = 1e-2, Ne=1e2, ff=1, Te = 1e4))

"""
# Pregunta 2 : Run Starburst99...
# look at the quanta file:
#  .10100E+07   52.652  -0.360   51.974  -0.850   48.593  -3.980   42.542
#  .50100E+07   51.841  -1.029   50.830  -1.809   48.973  -3.339   42.345
"""
print('QHe1/QH0 at 1Myr = {0:.5}, at 5 Myr = {1:.5}'.format(10**(51.974 - 52.652), 10**(50.830 - 51.841)))

"""
# Pregunta 3
"""

models_dir = 'CloudyModels'
ex6_2.make_model('M_ISB_008_1', models_dir, SED = 'table star "ISB_008.mod"', qH=49.22, SED_params='1000000')
ex6_2.plot_model('M_ISB_008_1', models_dir)

"""
# Pregunta 4
"""
for T in np.linspace(30000., 70000, 9):
    ex6_2.make_model('T_BB_{0:.0f}'.format(T/1e3), models_dir, SED = 'BB', qH=49.22, SED_params=T, n_zones=2, iterate=0)
ex6_2.search_T('T_BB', models_dir)

"""
# Pregunta 5
"""
ex6_2.make_model('M_BB_515', models_dir, SED = 'BB', qH=49.22, SED_params=51500)
ex6_2.plot_model('M_BB_515', models_dir, style='--')

"""
# Pregunta 6
"""
for T in np.linspace(30000., 46000, 9):
    ex6_2.make_model('T_WM_{0:.0f}'.format(T/1e3), models_dir, SED = 'table star "wmbasic.mod"', 
                   qH=49.22, SED_params=[T, 4, -0.3], n_zones=2, iterate=0)
ex6_2.search_T('T_WM', models_dir, SED = 'WM')

"""
# Pregunta 7
"""
ex6_2.make_model('M_WM_420', models_dir, SED = 'table star "wmbasic.mod"', qH=49.22, SED_params=[42000, 4, -0.3])
ex6_2.plot_model('M_WM_420', models_dir, style=':')

"""
# Pregunta 8
"""

ex6_2.print_Xi('M_', models_dir)

"""
# Pregunta 9
"""
ex6_2.plot_SED('M_', models_dir)

#-------------------------------------- ex6_3 ----------------------------------------------

import ex6_3

# Define where the models will be build (location of the Cloudy input and output files)
models_dir = '/Users/christophemorisset/DATAS/Choroni'
# Define where cloudy.exe is
pc.config.cloudy_exe = '/usr/local/Cloudy/c10.00/cloudy.exe' # point to the location of Cloudy.exe
# Write the Makefile with the correct cloudy.exe location
pc.print_make_file(models_dir)
# Print the file telling Cloudy which lines we want to output
ex6_3.print_line_list(models_dir)

"""
# Pregunta 2
"""
# this will run 56 models...
# Do not forget to put the ISB_008.mod files in the models_dir directory.
ex6_3.run_grid(models_dir, n_proc=3)

"""
# Pregunta 3
"""
# this read some SDSS data
ex6_3.plot_obs()

"""
# Pregunta 4
"""
# reading all the models
Ms = pc.load_models('/Users/christophemorisset/DATAS/Choroni/G', 
                    read_lin = True, read_emis = False, read_cont = False, 
                    list_elem = [], read_phy = False, read_rad=False)
ex6_3.plot_grid(Ms)

#-------------------------------------- ex6_4 ----------------------------------------------
import pyCloudy as pc
import ex6_4

pc.log_.level = 2
models_dir2 = '/Users/christophemorisset/DATAS/Choroni'
# define the model name and properties
model_name = 'M64_A'
i = 1
r_in = 15. 
dens = 4. 
Teff = 100000 
Q0 = 47.
distance = 1.0
ab_dict = {'He':-1.0, 'C':-3.10, 'N':-4.30, 'O':-3.0, 'Ne':-4., 'Mg':-4.95,
           'Si':-4.90, 'S':-5.35, 'Cl':-7.00, 'Ar':-6.2, 'Fe':-7.40}

pc.log_.level = 3
# create the object that generates the input files
Min = ex6_4.In(models_dir2, '{0}_{1}'.format(model_name, i), r_in, dens, Teff, Q0,
               ab_dict, distance)
Min.print_model()

# run the models
pc.run_cloudy(dir_=models_dir2, n_proc=3, use_make=True,  model_name=model_name)

# read the models
pc.log_.level = 2
Mouts = ex6_4.Outs(models_dir2, model_name)
# output the parameters and line intensities, with the observations in 1rst column
Mouts.print_res()

#-------------------------------------- ex6_5 ----------------------------------------------
import ex6_5
import pyCloudy as pc
models_dir = 'CloudyModels'

model_name = 'IC418_BB'
Min = ex6_5.In('{0}/{1}'.format(models_dir, model_name))
# Here you can change things in the model before writing the input file
# Example:
# run a Blackbody model
Min.Teff = 39500.
Min.SED = 'BB'
Min.print_model()

# change the model to a CMFGEN stellar atomsphere model
# do not forget to create the mod103.mod using:
# echo 'compile stars "mod103.ascii"' | cloudy.exe

model_name = 'IC418_CMF'
Min.model_name = '{0}/{1}'.format(models_dir, model_name)
Min.SED = 'STAR'
Min.print_model()

# you can comment the following to avoid re-run the model
pc.run_cloudy(dir_=models_dir, n_proc=2, use_make=True,  model_name='IC418_')

# reading the outputs, 
# reading the observations,
# doing the cross_calibration of the UV and IR relative to the optical, using the model results

# The 2 previous models are read:
MoutBB = ex6_5.Out('{0}/IC418_BB'.format(models_dir))
MoutCMF = ex6_5.Out('{0}/IC418_CMF'.format(models_dir))

# pretty print the results:
MoutBB.print_res()
MoutCMF.print_res()

# using a 3D models to take the slit into account
MoutCMF.set_3D(use=True)
MoutCMF.print_res()

# Draw a 3 colors image:
plt.figure()
MoutCMF.plot_3col()

# Draw an image and the slit:
MoutCMF.plot_im('H__1__4861A')

# PLot the continuum: top: The stellar continuum in the ionizing region
# bottom: the nebular+stellar in the UV+Opt+IR region
MoutCMF.plot_cont()
#-------------------------------------- ex7_3 ----------------------------------------------
import ex6_5

models_dir = 'CloudyModels'
Mout = ex6_5.Out('{0}/IC418_CMF'.format(models_dir))
plt.figure()

plt.subplot(2,2,1)
Mout.define_profiles(size_spectrum = 41,vel_max = 40, vel_params=[20.,0.,0])
Mout.plot_profile()

plt.subplot(2,2,2)
Mout.define_profiles(size_spectrum = 41,vel_max = 40, vel_params=[0.,20.,0])
Mout.plot_profile()

plt.subplot(2,2,3)
Mout.define_profiles(size_spectrum = 41,vel_max = 40, vel_params=[0.,0.,20])
Mout.plot_profile()

plt.subplot(2,2,4)
Mout.plot_profile(75, 75)

#-------------------------------------- ex7_3 ----------------------------------------------
import ex7_3
ex7_3.plot_2comp(tem1=1e4, tem2=1e4, dens1=3e2, dens2=5e5, mass1=1, mass2=5e-4)
        

"""