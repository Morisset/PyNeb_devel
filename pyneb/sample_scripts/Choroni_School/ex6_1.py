import numpy as np
import matplotlib.pyplot as plt
import pyCloudy as pc
from pyCloudy.utils.misc import sextract

def make_mod(models_dir, name, Teff, qH, dens, r_in, ff=1.0, metals=None, nograins=True):
    
    pc.log_.level=3

    options = ('no molecules',
               'no level2 lines',
               'no fine opacities',
               'atom h-like levels small',
               'atom he-like levels small',
               'COSMIC RAY BACKGROUND',
               'element limit off -8',
               )
    
    c_input = pc.CloudyInput('{0}/{1}'.format(models_dir, name))
    c_input.set_BB(Teff = Teff, lumi_unit = 'q(H)', lumi_value = qH)
    # Defining the density. You may also use set_dlaw(parameters) if you have a density law defined in dense_fabden.cpp.
    c_input.set_cste_density(dens, ff = ff)
    # Defining the inner radius. A second parameter would be the outer radius (matter-bounded nebula).
    c_input.set_radius(r_in = r_in)
    if nograins:
        # If nograins is set, then the metals option only acts on the metallicity
        c_input.set_abund(predef = 'ism', nograins = nograins, metals=metals)
    else:
        # If nograins is False (ism grains will thus be used), the dust abundance follows the metallicity
        c_input.set_abund(predef = 'ism', nograins = nograins, metalsgrains=metals)
    c_input.set_other(options)
    c_input.set_iterate() # (0) for no iteration, () for one iteration, (N) for N iterations.
    c_input.set_sphere() # () or (True) : sphere, or (False): open geometry.
    c_input.set_distance(dist=1., unit='kpc', linear=True) # unit can be 'kpc', 'Mpc', 'parsecs', 'cm'. If linear=False, the distance is in log.
    c_input.set_line_file(line_file = 'lines.dat', absolute=True)
    c_input.set_emis_tab(emis_tab_str = ['H  1  4861A'])
    c_input.set_stop(['temperature off', 'pfrac 0.02'])
    c_input.print_input()

def print_line_file(models_dir, name='lines.dat'):
    
    line_file = open('{0}/{1}'.format(models_dir, name), 'w')
    line_file.write('H  1  4861A \n')
    line_file.write('TOTL  3727A \n')
    line_file.write('O  3  5007A \n')
    line_file.write('O  3 88.33m \n')
    line_file.close()
    
def plot_Hb(models):
    
    Q0 = []
    for M in models:
        Q0.append(M.Q0)
    Q0 = np.asarray(Q0)
    
    # Equivalent compact way:
    Q0 = np.asarray([M.Q0 for M in models])
    
    Hbeta = np.asarray([M.get_line('H__1__4861A') for M in models])
    Hbeta2 = np.asarray([np.log10(M.get_emis_vol('H__1__4861A')) for M in models])
    r_out_cut = np.asarray([M.r_out_cut for M in models])
    plt.figure()
    
    plt.subplot(2,2,1)
    plt.scatter(np.log10(Q0), Hbeta)
    plt.xlabel('log(Q0)')
    plt.ylabel(r'log(H$\beta$)')

    plt.subplot(2,2,2)
    plt.scatter(np.log10(Q0), Hbeta2)
    plt.xlabel('log(Q0)')
    plt.ylabel(r'log(H$\beta$)')

    plt.subplot(2,2,3)
    plt.scatter(np.log10(Q0), np.log10(r_out_cut))
    plt.xlabel('log(Q0)')
    plt.ylabel(r'log(Rout cut)')
    
def plot_Hb_Teff(models):
    
    Teff = np.asarray([float(sextract(M.out['Blackbody'], 'body', '*')) for M in models])
    
    Hbeta = np.asarray([M.get_line('H__1__4861A') for M in models])
    Hbeta2 = np.asarray([np.log10(M.get_emis_vol('H__1__4861A')) for M in models])
    r_out_cut = np.asarray([M.r_out_cut for M in models])
    plt.figure()
    
    plt.subplot(2,2,1)
    plt.scatter(Teff, Hbeta)
    plt.xlabel('Teff')
    plt.ylabel(r'log(H$\beta$)')

    plt.subplot(2,2,2)
    plt.scatter(Teff, Hbeta2)
    plt.xlabel('Teff')
    plt.ylabel(r'log(H$\beta$)')

    plt.subplot(2,2,3)
    plt.scatter(Teff, np.log10(r_out_cut))
    plt.xlabel('Teff')
    plt.ylabel(r'log(Rout cut)')
    
def plot_ion(models, elem='O', ion=2):
    colors = ['red', 'green', 'blue', 'yellow', 'magenta', 'cyan']
    markers=['o', '+', '<', '>', 'd', '.']
    for i, M in enumerate(models):
        DR =(M.r_out-M.r_in)/M.r_in 
        plt.scatter(M.log_U, M.get_ionic(elem, ion), color = colors[i], marker=markers[i],
                    label='Model {0} DR={1:5.2f}'.format(i, np.log10(DR)))
    plt.xlabel('log_U')
    plt.ylabel('{0}{1}+/H+'.format(elem, ion))
    plt.legend(loc=2)

def plot_radial(M, color):
    
    plt.subplot(2,2,1)
    plt.plot(M.radius, M.te, color = color)
    plt.title('Te')

    plt.subplot(2,2,2)
    plt.plot(M.radius, M.get_ionic('O', 1), color = color)
    plt.plot(M.radius, M.get_ionic('O', 2), color = color, linestyle=':')
    plt.title('O+ and O++ (dot)')
    
    plt.subplot(2,2,3)
    plt.plot(M.radius, M.get_ionic('H', 1), color = color)
    plt.title('H+')
    
    plt.subplot(2,2,4)
    plt.plot(M.radius, M.get_ionic('N', 1), color = color)
    plt.plot(M.radius, M.get_ionic('N', 2), color = color, linestyle=':')
    plt.plot(M.radius, M.get_ionic('N', 3), color = color, linestyle='--')
    plt.title('N+, N++ (dot), and N+++ (dashed)')

def plot_radial2(M, color):
    norm_radius = (M.depth/M.depth[-1])
    plt.subplot(2,2,1)
    plt.plot(norm_radius, M.te, color = color)
    plt.title('Te')

    plt.subplot(2,2,2)
    plt.plot(norm_radius, M.get_ionic('O', 1), color = color)
    plt.plot(norm_radius, M.get_ionic('O', 2), color = color, linestyle=':')
    plt.title('O+ and O++ (dot)')
    
    plt.subplot(2,2,3)
    plt.plot(norm_radius, M.get_ionic('H', 1), color = color)
    plt.title('H+')
    
    plt.subplot(2,2,4)
    plt.plot(norm_radius, M.get_ionic('N', 1), color = color)
    plt.plot(norm_radius, M.get_ionic('N', 2), color = color, linestyle=':')
    plt.plot(norm_radius, M.get_ionic('N', 3), color = color, linestyle='--')
    plt.title('N+, N++ (dot), and N+++ (dashed)')

def plot_metals(models):
    
    Z = np.asarray([M.abund['O'] for M in models])
    r_out = np.asarray([np.log10(M.r_out) for M in models])
    Op = np.asarray([M.get_ab_ion_vol_ne('O', 1) for M in models])
    Opp = np.asarray([M.get_ab_ion_vol_ne('O', 2) for M in models])
    TOp = np.asarray([M.get_T0_ion_vol_ne('O', 1) for M in models])
    TOpp = np.asarray([M.get_T0_ion_vol_ne('O', 2) for M in models])
    O2 = np.asarray([M.get_line('TOTL__3727A') - M.get_line('H__1__4861A') for M in models])
    O3 = np.asarray([M.get_line('O__3__5007A') - M.get_line('H__1__4861A') for M in models])
    O3ir = np.asarray([M.get_line('O__3_8833M') - M.get_line('H__1__4861A') for M in models])
    
    plt.subplot(2,2,1)
    plt.scatter(Z, r_out, label='Rout')
    plt.legend()
    plt.ylim((18.4, 18.85))

    plt.subplot(2,2,2)
    plt.scatter(Z, Op, label='O+')
    plt.scatter(Z, Opp, label='O++', marker='+', color='red')
    plt.ylim(0.3, 0.65)
    plt.legend(loc=2)
    
    plt.subplot(2,2,3)
    plt.scatter(Z, TOp, label='T(O+)')
    plt.scatter(Z, TOpp, label='T(O++)', marker='+', color='red')
    plt.ylim((00, 18000))
    plt.legend()
    
    plt.subplot(2,2,4)
    plt.scatter(Z, O2, label='[OII]/Hb')
    plt.scatter(Z, O3, label='[OIII]/Hb', marker='+', color='red')
    plt.scatter(Z, O3ir, label='[OIII]IR/Hb', marker='d', color='magenta')
    plt.ylim(-2.5, 0.5)
    plt.legend(loc=2)
    
