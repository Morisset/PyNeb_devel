import numpy as np
import pyCloudy as pc
import matplotlib.pyplot as plt
from pyneb.utils.physics import IP

# TODO: Add comments
"""
Pregunta 1
"""
def alpha_B(Te):
    """
    Recomb. coefficient, case B
    """
    T4 = Te/1e4
    return 2.6e-13/T4


def U_mean_def(QH0, Ne, Rstr):
    """
    \int_0^{Rstr}{U.dV} / \int_0^{Rstr}{dV}
    Return the mean over Stromgren volume of U
    """
    return 3.* QH0 / (4. * np.pi * pc.CST.CLIGHT * Ne * Rstr**2)


def QH0_def(Rstr, Ne, ff, Te = 1e4):
    """
    Volume of 
    """
    return 4. / 3. * np.pi * Rstr**3 * Ne**2 * ff * alpha_B(Te)

def Rstr(QH0, Ne, ff, Te = 1e4):
    
    return (3. * QH0 / (4. * np.pi * ff * alpha_B(Te) * Ne**2))**(1./3.)

def U_mean(QH0, Ne, ff, Te = 1e4):
    
    return (Ne * QH0 * ff**2 * 3 / (4. * np.pi) * alpha_B(Te)**2)**(1./3.) / pc.CST.CLIGHT

def QH0(U_mean, Ne, ff, Te = 1e4):
    
    return U_mean**3 * pc.CST.CLIGHT**3 * 4. * np.pi / 3. / (Ne * ff**2 * alpha_B(Te)**2)  

# -------------------------------------- 
"""
Pregunta 3
"""
def make_model(name, models_dir='./', SED='BB', qH=None, SED_params=None, n_zones = None, iterate=1):
    pc.log_.level=3

    abund_AGSS09 = {'He' : 10.93, 'C' : 8.43, 'N' : 7.83, 'O' : 8.69, 'Ne' : 7.93, 'Mg' : 7.6,
             'S' : 7.12, 'Ar' : 6.40, 'Fe' : 7.5, 'Cl' : 5.5, 'Si' : 7.51}
    for elem in abund_AGSS09:
        abund_AGSS09[elem] -= 12
        if elem != 'He':
            abund_AGSS09[elem] -= 0.3

    options = ('no molecules',
               'no level2 lines',
               'no fine opacities',
               'atom h-like levels small',
               'atom he-like levels small',
               'COSMIC RAY BACKGROUND',
               'element limit off -8',
               )
    
    c_input = pc.CloudyInput('{0}/{1}'.format(models_dir, name))
    if SED == 'BB':
        c_input.set_BB(Teff = SED_params, lumi_unit = 'q(H)', lumi_value = qH)
    else:
        c_input.set_star(SED = SED, SED_params = SED_params, lumi_unit = 'q(H)', lumi_value=qH)
    # Defining the density. You may also use set_dlaw(parameters) if you have a density law defined in dense_fabden.cpp.
    c_input.set_cste_density(2, ff = 1.)
    # Defining the inner radius. A second parameter would be the outer radius (matter-bounded nebula).
    c_input.set_radius(r_in = np.log10(pc.CST.PC/10))
    c_input.set_abund(ab_dict = abund_AGSS09, nograins = True)
    c_input.set_other(options)
    c_input.set_iterate(iterate) # (0) for no iteration, () for one iteration, (N) for N iterations.
    c_input.set_sphere() # () or (True) : sphere, or (False): open geometry.
    c_input.set_distance(dist=1., unit='kpc', linear=True) # unit can be 'kpc', 'Mpc', 'parsecs', 'cm'. If linear=False, the distance is in log.
    if n_zones is not None:
        c_input.set_stop('zones {0}'.format(n_zones))
    c_input.print_input()
    
    c_input.run_cloudy()
    
def plot_model(name, models_dir = './', style='-', fig_num = 1):
    pc.log_.level=3

    M = pc.CloudyModel('{0}/{1}'.format(models_dir, name), read_emis = False)
    X = M.radius/1e19
    colors = ['r', 'g', 'b', 'y', 'm', 'c']
    plt.figure(fig_num)
    
    plt.subplot(3, 3, 1)
    plt.plot(X, M.get_ionic('H', 0), label='H0', linestyle=style, c= colors[0])
    plt.plot(X, M.get_ionic('H', 1), label='H+', linestyle=style, c= colors[1])
    plt.plot(X, M.get_ionic('He', 0), label='He0', linestyle=style, c= colors[2])
    plt.plot(X, M.get_ionic('He', 1), label='He+', linestyle=style, c= colors[3])
    plt.plot(X, M.get_ionic('He', 2), label='He++', linestyle=style, c= colors[4])
    if style== '-':
        plt.legend()
    plt.title(name)
    
    for i_plot, elem in enumerate(['N', 'O', 'Ne', 'S', 'Ar']):
        plt.subplot(3, 3, i_plot + 2)
        for i in np.arange(4):
            plt.plot(X, M.get_ionic(elem, i), linestyle=style, c=colors[i])
        plt.text(np.max(X)/2, 0.9, elem)
        if i_plot == 0:
            plt.title(M.date_model)
    
    plt.subplot(3, 3, 7)
    plt.plot(X, M.ne, label=r'N$_e$', linestyle=style, c='blue')
    plt.plot(X, M.nH, label='N$_H$', linestyle=style, c='red')
    if style== '-':
        plt.legend(loc=3)
    plt.xlabel(r'R [10$^{19}$cm]')
    
    plt.subplot(3, 3, 8)
    plt.plot(X, M.te, label=r'T$_e$', linestyle=style, c='blue')
    if style== '-':
        plt.legend(loc=3)
    
    plt.subplot(3, 3, 9)
    plt.plot(X, M.log_U, label='log U', c='blue')
    if style== '-':
        plt.legend()
     
def search_T(name, models_dir = './', SED = 'BB'):
    
    Ms = pc.load_models('{0}/{1}'.format(models_dir, name), read_emis = False)
    if SED == 'BB':
        T = np.array([float(pc.sextract(M.out['Blackbody'], 'dy ', '*')) for M in Ms])
    elif SED == 'WM':
        T = np.array([float(pc.sextract(M.out['table star'], 'mod" ', '4.0')) for M in Ms])
    QH0 = np.array([M.Q0 for M in Ms])
    QHe0 = np.array([M.Q[1::].sum() for M in Ms])
    
    plt.plot(T/1e3, QHe0/QH0)
    plt.xlabel('T [kK]')
    plt.ylabel('QHe0/QH0')
        
def print_Xi(name, models_dir = './'):
    Ms = pc.load_models('{0}/{1}'.format(models_dir, name), read_emis = False)
    names = [M.model_name_s for M in Ms]
    print(names)
    print('H0/H:   {0:.2e} {1:.2e} {2:.2e}'.format(Ms[0].get_ab_ion_vol('H', 0), 
                                                   Ms[1].get_ab_ion_vol('H', 0), 
                                                   Ms[2].get_ab_ion_vol('H', 0)))
    
    print('H1/H:   {0:.2e} {1:.2e} {2:.2e}'.format(Ms[0].get_ab_ion_vol('H', 1), 
                                                   Ms[1].get_ab_ion_vol('H', 1), 
                                                   Ms[2].get_ab_ion_vol('H', 1)))
    
    print('He0/H:  {0:.2e} {1:.2e} {2:.2e}'.format(Ms[0].get_ab_ion_vol('He', 0), 
                                                   Ms[1].get_ab_ion_vol('He', 0), 
                                                   Ms[2].get_ab_ion_vol('He', 0)))
    
    print('He1/H:  {0:.2e} {1:.2e} {2:.2e}'.format(Ms[0].get_ab_ion_vol('He', 1), 
                                                   Ms[1].get_ab_ion_vol('He', 1), 
                                                   Ms[2].get_ab_ion_vol('He', 1)))
    
    print('He2/H:  {0:.2e} {1:.2e} {2:.2e}'.format(Ms[0].get_ab_ion_vol('He', 2), 
                                                   Ms[1].get_ab_ion_vol('He', 2), 
                                                   Ms[2].get_ab_ion_vol('He', 2)))
    for elem in ['N', 'O', 'Ne', 'S', 'Ar']:
        for i in np.arange(4):
            print('{0:2s}{1}/H:  {2:.2e} {3:.2e} {4:.2e}'.format(elem, i, Ms[0].get_ab_ion_vol(elem, i), 
                                                             Ms[1].get_ab_ion_vol(elem, i), 
                                                             Ms[2].get_ab_ion_vol(elem, i)))
        
def plot_SED(name, models_dir = './', unit='Jy'):
    Ms = pc.load_models('{0}/{1}'.format(models_dir, name), read_emis = False)
    plt.figure()
    plt.subplot(2, 1, 1)
    for M in Ms:
        plt.plot(M.get_cont_x(unit = 'eV'), np.log10(M.get_cont_y(unit = 'esHz')), label=M.model_name_s)
    plt.xlim((10., 60))
    plt.ylim((18, 24))
    plt.ylabel('log [erg.s-1.Hz-1]')
    plt.legend(loc=3)
    
    plt.subplot(2, 1, 2)
    for M in Ms:
        plt.plot(M.get_cont_x(unit = 'eV'), np.log10(M.get_cont_y(unit = 'Q')), label=M.model_name_s)
    plt.xlim((10., 60))
    plt.ylim((42., 50))
    plt.xlabel('E [eV]')
    plt.ylabel('QH0(E)')
    # TODO: avoid overlap
    for ip in IP:
        plt.plot([IP[ip], IP[ip]], [49, 50])
        plt.text(IP[ip], 48, ip)
    
