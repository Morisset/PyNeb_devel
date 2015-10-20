#######################################################################
# GETTING STARTED
# To get started, move into directory where PyNeb resides and enter python

# import code and modules
import pyneb as pn
import numpy as np

#######################################################################
# DEFINING ATOMS

# define an OII atom
O2 = pn.Atom("O", 2)

# alternate syntax to define an atom (spec is a string)
N2 = pn.Atom("N", "2")

# check atom definition
O2.elem
O2.spec
O2.atom
O2.name

# inspect the header of atomic data file
# print 'Atomic data header:'
print O2.AtomData.AtomHeader 

# inspect the header of collisional data file
# print 'Collisional header:'
print O2.CollData.CollHeader 

# explore the atom: builtin data
O2.gs # ground-state configuration
# array of stat weights
O2.getStatWeight()
# stat weight of a given level
lev_i = 2
O2.getStatWeight(lev_i)


# explore the atom: adopted atomic data
pn.atomicData.getPredefinedDataFileDict() # we suggest using the tab for this command...
pn.atomicData.getDirForFile('o_ii_atom_WFD96.fits') # wanna know where the file lies?
O2.printSources() # print bibliographic references of data used to build O2

O2.NLevels # number of levels in the selected data
O2.getEnergy(2) # energy of first excited level (ground = 1) in Angstrom^-1
O2.getA(2,1) # transition probability of 2->1
# set temperature and density
tem = 15000.
den = 1000.
O2.getPopulations(tem, den) # compute populations
O2.getCritDensity(tem, level=2) # critical density of level 2 at tem
O2.getOmega(tem, 2, 1) # effective collision strength of transition 2->1 at T=10000K
O2.getOmegaArray(2, 1) # array of effective collision strengths for 2->1 as a function of T
O2.getTemArray() # print array of temperatures of tabulated Omegas
O2.getCollRates(tem) # print collisional Rates at T=tem


# This bit calls the script DataPlot.py to plot atomic data. 
dataplot = pn.DataPlot('O', 3)
dataplot.plotA() # transition probabilities plot 
dataplot.plotRelA() # relative transition probabilities plot
dataplot.plotOmega() # collision strength plot    
        

# customize atomic data 
# First step: check which directories are searched for atomic data files
pn.atomicData.getAllDataFilePaths()
# Add your selected directory to the list
pn.atomicData.addDataFilePath('./')
# Check if it's been added
pn.atomicData.getAllDataFilePaths()
# Remove it if you gave the wrong dir
pn.atomicData.removeDataFilePath('./')
# Set 'o_iii_TEST.fits' to be the OIII atom file
pn.atomicData.setDataFile('o_iii_TEST.fits', 'O3', 'atom')
# Check that it is read from the right place
pn.atomicData.getDirForFile('o_ii_atom_TEST.fits')
# Summarize all data used
pn.atomicData.printDirForAllFiles()

# define an atom with the new data
O2test = pn.Atom("O", 2)

# define all atoms at once and put them in a dictionary
# (all of them defined with the latest dataset selected)
atoms = pn.getAtomDict() #this generates a lot of warnings as not all element-spectrum combination exist
# see what atoms have been built
atoms
# build only some atoms                       
atoms = pn.getAtomDict(atom_list=['O1', 'O2', 'O3', 'N2', 'N3'])

# explore some specific atom in the atoms collection
atoms['N2'].elem

# if you want to be able to access them directly rather than through a dictionary:
for key in atoms.keys():
    vars()[key]=atoms[key]

# for example
O2.elem

# list all atom features
dir(O2)

#######################################################################
# MAKING CALCULATIONS

# set temperature and density
tem = 15000.
den = 1000.

# compute populations
O2.getPopulations(tem, den)

# compute emissivity of transition (lev_i, lev_j)
O2.getEmissivity(tem, den, 3, 2)

# also works if tem is an array
tem = np.array([10000, 20000, 30000])
O2.getOmega(tem, 2, 1)
O2.getCollRates(tem)
O2.getPopulations(tem, den)

# tem and den can be arrays as well as single numbers
tem = np.array([10000, 12000, 13000]) 
den = np.array([100, 200, 300])
O2.getPopulations(tem, den) # returns the n_tem x n_den x n_levels array of populations
O2.getPopulations(tem, den, product=False) #  element-by-element multiplication of tem and den (no scalar product: returns [pop(tem_1, den_1), pop(tem_2, den_2), ... pop(tem_n, den_n)] 

# find transition corresponding to given wavelength
N2.printTransition(6584)

# temperature determination from an intensity ratio
N2.getTemDen(0.01, den=1000., wave1=5755, wave2=6584)

# same as above, by specifying the levels involved
N2.getTemDen(0.01, den=1000., lev_i1=5, lev_j1=4, lev_i2=4, lev_j2=3)

# same as above, by specifying the levels involved
N2.getTemDen(0.01, den=1000., to_eval = 'L(5755) / L(6584)')

# no formal difference between temperature and density diagnostics, so beware of what you do
N2.getTemDen(0.01, tem=8782., wave1=5755, wave2=6584)
N2.getTemDen(0.01, tem=8882., wave1=5755, wave2=6584)

# ionic abundance (intensity, temperature, density, transition)
O2.getIonAbundance(100, 1.5e4, 100., wave=3727)

# printout as in old nebular
O2.printIonic() # only prints transitions and corresponding wavelengths. Useful for a quick glance at the atom.
O2.printIonic(tem=10000, den=100) # also prints line emissivities
O2.printIonic(tem=10000, den=100, printA=True, printPop=True, printCrit=True) # also prints populations and critical densities

# Compute Hb emissivity at T=10000K
pn.getRecEmissivity(10000, 1e2, 4, 2, atom='H1')

# simultaneously compute temperature and density from pairs of line ratios
# First of all, a Diagnostics object must be created and initialized with the relevant diagnostics.
diags = pn.Diagnostics()   # this creates the object
diags.getAllDiags()  # see what Diagnostics exist
tem, den = diags.getCrossTemDen('[NII] 5755/6548', '[SII] 6731/6716', 50, 1.0, 
                                guess_tem=10000, tol_tem = 1., tol_den = 1., max_iter = 5)

#TO BE CONTINUED FROM HERE
print tem, den

#######################################################################
# HANDLING OBSERVATIONS

# explore the line label list to write an observation record
pn.LINE_LABEL_LIST
# if only Ar IV is needed
pn.LINE_LABEL_LIST['Ar4']

# list the available extinction laws
pn.RedCorr().printLaws()
