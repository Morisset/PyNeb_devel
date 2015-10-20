import numpy as np
import pyCloudy as pc
import matplotlib.pyplot as plt

def print_line_list(models_dir):
    line_list = open('{0}/line_list.dat'.format(models_dir), 'w')
    txt = """
H  1  4861A 
H  1  6563A 
N  2  6584A 
O  3  5007A 
"""
    line_list.write(txt)
    line_list.close()

def make_model(name, logU, logZ, models_dir='./'):
    """
    
    Z : log of Z/Zsol
    """
    pc.log_.level=3

    abund_AGSS09 = {'He' : 10.93, 'C' : 8.43, 'N' : 7.83, 'O' : 8.69, 'Ne' : 7.93, 'Mg' : 7.6,
             'S' : 7.12, 'Ar' : 6.40, 'Fe' : 7.5, 'Cl' : 5.5, 'Si' : 7.51}
    for elem in abund_AGSS09:
        abund_AGSS09[elem] -= 12
        if elem != 'He':
            abund_AGSS09[elem] += logZ

    options = ('no molecules',
               'no level2 lines',
               'no fine opacities',
               'atom h-like levels small',
               'atom he-like levels small',
               'COSMIC RAY BACKGROUND',
               'element limit off -8',
               )
    
    c_input = pc.CloudyInput('{0}/{1}'.format(models_dir, name))
    c_input.set_star(SED = 'table star "ISB_008.mod"', SED_params = 1000000, 
                     lumi_unit = 'ionization parameter', lumi_value=logU)
    # Defining the density. You may also use set_dlaw(parameters) if you have a density law defined in dense_fabden.cpp.
    c_input.set_cste_density(2)
    c_input.set_abund(ab_dict = abund_AGSS09, nograins = True)
    c_input.set_other(options)
    c_input.set_line_file('line_list.dat')
    c_input.set_iterate() # (0) for no iteration, () for one iteration, (N) for N iterations.
    c_input.set_sphere(False) # () or (True) : sphere, or (False): open geometry.
    c_input.set_distance(dist=10., unit='Mpc', linear=True) # unit can be 'kpc', 'Mpc', 'parsecs', 'cm'. If linear=False, the distance is in log.
    c_input.print_input()
    
def run_grid(models_dir, n_proc):
    pc.print_make_file(models_dir)

    # To have the full list of possible extensions:
    # print pc.config.SAVE_LIST
    # Here we set to nothing the list of cloudy's output files.
    # Only the .lin files will be output
    pc.config.SAVE_LIST = []
    pc.config.SAVE_LIST_ELEMS = []

    # Metallicity table
    Zs = np.array([0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.5])
    # Ionization parameter table
    logUs = [-2, -2.3, -2.7, -3, -3.3, -3.7, -4]

    # loop on the tables to write the input files
    for Z in Zs:
        for logU in logUs:
            make_model(name='G_{0:.1f}_{1:.1f}'.format(logU, Z), logU=logU,
                       logZ=np.log10(Z), models_dir=models_dir)
            
    # Run all the models
    pc.run_cloudy(dir_ = models_dir, n_proc = n_proc)

def plot_obs():
    obs = np.genfromtxt('BPT4Graz_f4.dat', names=True)
    plt.plot(obs['xN2Ha'], obs['yO3Hb'], 'y,', linestyle='None')    
        
def plot_grid(Ms):
    # This plots the results of the grid.

    # A small function to extract a line intensity from all the models 
    extract_line = lambda label:np.array([M.get_line(label) for M in Ms])
    
    Ha = extract_line('H__1__6563A')
    Hb = extract_line('H__1__4861A')
    O3 = extract_line('O__3__5007A')
    N2 = extract_line('N__2__6584A')

    # Recover the list of unique values of the input parameters from the name of the model:
    Z = np.array([float(M.model_name_s.split('_')[2]) for M in Ms])
    Z_u = np.sort(np.unique(Z))
    logU = np.array([float(M.model_name_s.split('_')[1]) for M in Ms])
    logU_u = np.sort(np.unique(logU))
    
    for logU1 in logU_u:
        s = np.where(logU == logU1)[0] # find the indices where logU is the current logU1
        indx = s[Z[s].argsort()] #find and sort the indices of the models with the 
        plt.plot(np.log10(N2/Ha)[indx], np.log10(O3/Hb)[indx], linestyle='-', label = 'logU={0:.2f}'.format(logU1))
    for Z1 in Z_u:
        s = np.where(Z == Z1)[0] # 
        indx = s[logU[s].argsort()]
        plt.plot(np.log10(N2/Ha)[indx], np.log10(O3/Hb)[indx], linestyle='--', label = 'Zneb={0:.2f}'.format(Z1))
    plt.legend()
    
    plt.xlim((-1.5, 1))
    plt.ylim((-1.5, 1.5))
