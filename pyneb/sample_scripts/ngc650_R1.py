# Sample script of complete analysis of an object, including physical conditions and ionic abundances
# Read line intensity from an ascii file
# Compute desired diagnostics
# Compute abundances 

import matplotlib.pyplot as plt
import pyneb as pn
import os


### General settings
# Setting verbosity level. Enter pn.log_? for details
pn.log_.level = 2
# Set to True if the emission maps have already been generated
restore = True 
# Adopt IRAF atomic data set
pn.atomicData.setDataFileDict('IRAF_09')
# Define where emission maps are to be stored (restore = False) or read from (restore = True)
pypic_path = './pypics/'
# Adopt an extinction law
extinction_law = 'CCM 89'
# Instantiate the Diagnostics class
diags = pn.Diagnostics()
# Include in diags the relevant line ratios
diags.addDiag([
				'[NI] 5198/5200',
				'[NII] 5755/6548',
				'[OIII] 4363/5007',
				'[ArIII] 5192/7136',
				'[ArIII] 5192/7300+',
				'[ArIV] 4740/4711',
				'[ArIV] 7230+/4720+',
				'[SII] 6731/6716',
				'[NeIII] 15.6m/36.0m',
				'[NeV] 14.3m/24.2m'
				])
# Create the emission maps to be compared to the observation data (some overkill here)
emisgrids = pn.getEmisGridDict(atom_list=diags.getUniqueAtoms(), den_max=1e6, pypic_path=pypic_path)


# Object to be analyzed (could be a list with more than one element)
ngc_list = ['ngc650_R1']

# Loops over items in list (only one in this example)
for ngc in ngc_list:
	
	# Define the data file
	obs_data = ngc +'.dat'
	# Define title
	title = ngc.upper()
	 
	### Read and deredden observational data

	# Define an Observation object and assign it to name 'obs'
	obs = pn.Observation()
	# Fill obs with data read from file obs_data, with lines varying across rows and 
	# a default percent error on line intensities
	obs.readData(obs_data, fileFormat='lines_in_rows', err_default=0.05)
	# Deredden data with Cardelli's law
	obs.extinction.law = extinction_law
	# Apply correction
	obs.correctData()

	#Create file path for data if it does not already exist
	filepath = './'+str(ngc)+'/'
	if not os.path.exists(filepath):
		os.mkdir(filepath)
		
	### Plot
	# Create the contour plot as the intersection of tem-den emission maps with dereddened line ratios
	diags.plot(emisgrids, obs)
	# Place the title
	plt.title(title)
	# Display the plot
	plt.show()
	# To save, click on the icon
	plt.close()

	# Define all atoms to make calculations (etAtomDict is overkill, but worth 
	# given the large number of ions required)
	all_atoms = pn.getAtomDict()
	
	# Simultaneously determine T(OIII) and N(SII)
	try:
		tem_O3, den_S2 = diags.getCrossTemDen('[OIII] 4363/5007', '[SII] 6731/6716', obs=obs)
	except:
		tem_O3 = 'NA'
		den_S2 = 'NA'
		pass
		
	# Simultaneously determine T(NII) and N(SII)
	try:
		tem_N2, den_tmp = diags.getCrossTemDen('[NII] 5755/6548', '[SII] 6731/6716', obs=obs)
	except:
		tem_N2 = 'NA'
		pass

	# Simultaneously determine T(OIII) and N(ArIV)		
	try:
		tem, den_Ar4 = diags.getCrossTemDen('[OIII] 4363/5007', '[ArIV] 4740/4711', obs=obs)
	except:
		den_Ar4 = 'NA'
		pass
		
	# Simultaneously determine T(OIII) and N(SIII)
	try:	
		tem, den_S3 = diags.getCrossTemDen('[OIII] 4363/5007', '[SIII] 18.7m/33.6m', obs=obs)
	except:
		den_S3='NA'
		pass

	# Simultaneously determine T(OIII) and N(NeIII)
	try:
		tem, den_Ne3 = diags.getCrossTemDen('[OIII] 4363/5007', '[NeIII] 15.6m/36.0m', obs=obs) 
	except:
		den_Ne3='NA'
		pass
		
	# Intensities
	try:
		for line in obs.lines:
			if line.label == 'He2_4686A':
				i4686 = line.corrIntens
			if line.label == 'He1_5876A':
				i5876 = line.corrIntens
			if line.label == 'He1_6678A':
				i6678 = line.corrIntens
			if line.label == 'H1_4861A':
				i4861 = line.corrIntens
			if line.label == 'O3_5007A':
				i5007 = line.corrIntens
			if line.label == 'O3_4363A':
				i4363 = line.corrIntens
	except (RuntimeError, TypeError, NameError):
		pass
				
	# Alternate way of computing T(OIII)
	if tem_O3 == 'NA':
		tem_O3 = all_atoms['O3'].getTemDen(i5007/i4363, den=100., wave1=5007, wave2=4363)

	# Printout of physical conditions
	print('tem_O3: ', tem_O3)
	print('tem_N2: ', tem_N2)
	print('den_S2: ', den_S2)
	print('den_Ar4: ', den_Ar4)	
	print('den_S3:', den_S3)
	print('den_Ne3: ', den_Ne3)
	print('i4686: ', i4686)
	print('i5876: ', i5876)
	print('i6678: ', i6678)
	print('i4861: ', i4861)

	# Calculation and printout of abundances
	try:
		for line in obs.lines:
			if line.atom in all_atoms:
				ab = all_atoms[line.atom].getIonAbundance(line.corrIntens, tem_O3, den_S2, to_eval=line.to_eval)
				ab2 = all_atoms[line.atom].getIonAbundance(line.corrIntens, tem_N2, den_S2, to_eval=line.to_eval)
				ab3 = all_atoms[line.atom].getIonAbundance(line.corrIntens, tem_O3, den_Ar4, to_eval=line.to_eval)
				ab4 = all_atoms[line.atom].getIonAbundance(line.corrIntens, tem_O3, den_Ne3, to_eval=line.to_eval)
				print('{0:9s}'.format(line.label) + '  '.join(['{0:>20.10e}'.format(t) for t in (ab)]))
				print('{0:9s}'.format(line.label) + '  '.join(['{0:>20.10e}'.format(t) for t in  (ab2)]))
				print('{0:9s}'.format(line.label) + '  '.join(['{0:>20.10e}'.format(t) for t in  (ab3)]))
				print('{0:9s}'.format(line.label) + '  '.join(['{0:>20.10e}'.format(t) for t in  (ab4)]))
			else:
				pn.log_.warn('line from %s not used because ion not found' % line.atom, calling='ngc605_R1.py')
		pn.log_.timer('Ending ngc605_R1.py', calling='ngc605_R1.py')
	except (RuntimeError, TypeError, NameError):
		pass

	# Write result in a file
	fout = open(filepath + ngc +'.dat', "w")
	fout.write('tem_03 ' + str(tem_O3)+"\n")
	fout.write('tem_N2 ' + str(tem_N2)+"\n")
	fout.write('den_S2 ' + str(den_S2)+"\n")
	fout.write('den_Ar4 ' + str(den_Ar4)+"\n")
	fout.close()
