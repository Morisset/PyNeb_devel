from math import sqrt
import numpy as np

class CST(object):
    BOLTZMANN = 1.3806488e-16 # erg/K - NIST 2010
    CLIGHT = 2.99792458e10 # cm/s - NIST 2010
    HPLANCK = 6.62606957e-27 # erg s - NIST 2010
    EMASS = 9.10938291e-28 # g - NIST 2010
    ECHARGE = 1.602176565e-19 # Electron charge in Coulomb - NIST 2010
    PI = 3.141592653589793238462643
    BOLTZMANN_ANGK = (BOLTZMANN) / (HPLANCK * CLIGHT * 1.e8) # Boltzmann constant in (Ang * Kelvin) ** -1
    RYD = 109737.31568539 # Rydberg constant in cm^-1 - NIST 2010 
    RYD_EV = HPLANCK * CLIGHT * RYD * 1.e-7 / ECHARGE # infinite mass Rydberg in eV
    RYD_ANG = 1.e8 / RYD # infinite mass Rydberg in A
    RYD2ERG = HPLANCK * CLIGHT * RYD
    KCOLLRATE = sqrt(2 * PI / BOLTZMANN) * (HPLANCK / (2 * PI))**2 / EMASS**1.5 # constant of collisional rate equation
    BOLTZMANN_eVK = 8.617343e-5 # Boltzmann constant in eV/K
    
    HBETA = 4861.3316598713955

Z = {}
Z["H"] = 1
Z["He"] = 2
Z["Li"] = 3
Z["Be"] = 4
Z["B"] = 5
Z["C"] = 6
Z["N"] = 7
Z["O"] = 8
Z["F"] = 9
Z["Ne"] = 10
Z["Na"] = 11
Z["Mg"] = 12
Z["Al"] = 13
Z["Si"] = 14
Z["P"] = 15
Z["S"] = 16
Z["Cl"] = 17
Z["Ar"] = 18
Z["K"] = 19
Z["Ca"] = 20
Z["Sc"] = 21
Z["Ti"] = 22
Z["V"] = 23
Z["Cr"] = 24
Z["Mn"] = 25
Z["Fe"] = 26
Z["Co"] = 27
Z["Ni"] = 28
Z["Cu"] = 29
Z["Zn"] = 30
Z["Ga"] = 31
Z["Ge"] = 32
Z["As"] = 33
Z["Se"] = 34
Z["Br"] = 35
Z["Kr"] = 36
Z["Rb"] = 37
Z["Sr"] = 38
Z["Y"] = 39
Z["Zr"] = 40
Z["Nb"] = 41
Z["Mo"] = 42
Z["Tc"] = 43
Z["Ru"] = 44
Z["Rh"] = 45
Z["Pd"] = 46
Z["Ag"] = 47
Z["Cd"] = 48
Z["In"] = 49
Z["Sn"] = 50
Z["Sb"] = 51
Z["Te"] = 52
Z["I"] = 53
Z["Xe"] = 54
Z["Cs"] = 55
Z["Ba"] = 56
Z["La"] = 57
Z["Ce"] = 58
Z["Pr"] = 59
Z["Nd"] = 60
Z["Pm"] = 61
Z["Sm"] = 62
Z["Eu"] = 63
Z["Gd"] = 64
Z["Tb"] = 65
Z["Dy"] = 66
Z["Ho"] = 67
Z["Er"] = 68
Z["Tm"] = 69
Z["Yb"] = 70
Z["Lu"] = 71
Z["Hf"] = 72
Z["Ta"] = 73
Z["W"] = 74
Z["Re"] = 75
Z["Os"] = 76
Z["Ir"] = 77
Z["Pt"] = 78
Z["Au"] = 79
Z["Hg"] = 80
Z["Tl"] = 81
Z["Pb"] = 82
Z["Bi"] = 83
Z["Po"] = 84
Z["At"] = 85
Z["Rn"] = 86
Z["Fr"] = 87
Z["Ra"] = 88
Z["Ac"] = 89
Z["Th"] = 90
Z["Pa"] = 91
Z["U"] = 92
Z["Np"] = 93
Z["Pu"] = 94
Z["Am"] = 95
Z["Cm"] = 96
Z["Bk"] = 97
Z["Cf"] = 98
Z["Es"] = 99
Z["Fm"] = 100
Z["Md"] = 101
Z["No"] = 102
Z["Lr"] = 103
Z["Rf"] = 104
Z["Db"] = 105
Z["Sg"] = 106
Z["Bh"] = 107
Z["Hs"] = 108
Z["Mt"] = 109


IP = {}
IP['H2'] = 13.598
IP['He1'] = 0.0
IP['He2'] = 24.6
IP['He3'] = 54.4
IP['H1r'] = 13.598
IP['He1r'] = 24.6
IP['He2r'] = 54.4
IP['C1'] = 0.0
IP['C2'] = 11.26
IP['C3'] = 24.38
IP['C4'] = 47.89
IP['C5'] = 64.49
IP['N1'] = 0.0
IP['N2'] = 14.53
IP['N3'] = 29.60
IP['N4'] = 47.45
IP['N5'] = 77.47
IP['O1'] = 0.0
IP['O2'] = 13.62
IP['O3'] = 35.12
IP['O4'] = 54.94
IP['S1'] = 0.0
IP['S2'] = 10.36
IP['S3'] = 23.34
IP['S4'] = 34.79
IP['S5'] = 47.22
IP['Ne1'] = 0.0
IP['Ne2'] = 21.56
IP['Ne3'] = 40.96
IP['Ne4'] = 63.45
IP['Ne5'] = 97.12
IP['Ne6'] = 126.21
IP['Cl1'] = 0.0
IP['Cl2'] = 12.97
IP['Cl3'] = 23.81
IP['Cl4'] = 39.61
IP['Ar1'] = 0.0
IP['Ar2'] = 15.76
IP['Ar3'] = 27.63
IP['Ar4'] = 40.74
IP['Ar5'] = 59.81
IP['Ba1'] = 0.00
IP['Ba2'] = 5.21
IP['Ba3'] = 10.003826
IP['Ba4'] = 35.84
IP['Xe1'] = 0.00
IP['Xe3'] = 20.98
IP['Xe4'] = 31.05
IP['Xe5'] = 42.20
IP['Xe6'] = 54.1
IP['Mg1'] = 0.00
IP['Mg2'] = 7.65
IP['Mg3'] = 15.03
IP['Mg4'] = 80.14
IP['Mg5'] = 109.26
IP['Fe1'] = 0.0
IP['Fe2'] = 7.90 
IP['Fe3'] = 16.188
IP['Fe4'] = 30.652
IP['Fe5'] = 54.8
IP['Fe6'] = 75.0
IP['Fe7'] = 98.99
IP['Ni1'] = 0.
IP['Ni2'] = 7.639877
IP['Ni3'] = 18.168837
IP['Ni4'] = 35.187
IP['Kr1'] = 0.0
IP['Kr2'] = 13.9996
IP['Kr3'] = 24.35984
IP['Kr4'] = 35.67
IP['Kr5'] = 50.85
IP['Na1'] = 0.0
IP['Na2'] = 5.1390767
IP['Na3'] = 47.28636
IP['Na4'] = 71.6200
IP['Na5'] = 98.936
IP['Rb1'] = 0.0
IP['Rb2'] = 4.177128
IP['Rb3'] = 27.28954
IP['Rb4'] = 39.2472
IP['Rb5'] = 52.20
IP['K1'] = 0.0
IP['K2'] = 4.34066354
IP['K3'] = 31.62500
IP['K4'] = 45.8031
IP['K5'] = 60.917
IP['K6'] = 82.66

sym2name = {
            'H': 'hydrogen',
            'He': 'helium',
            '3He': 'helium',
            'B': 'boron',
            'C': 'carbon',
            'N': 'nitrogen',
            'O': 'oxygen',
            'F': 'fluorine',
            'Ne': 'neon',
            'Na': 'sodium',
            'Mg': 'magnesium',
            'Al': 'aluminium',
            'Si': 'silicon',
            'P': 'phosphorus',
            'S': 'sulfur',
            'Cl': 'chlorine',
            'Ar': 'argon',
            'K': 'potassium',
            'Ca': 'calcium',
            'Cr': 'chromium',
            'Fe': 'iron',
            'Co': 'cobalt',
            'Ni': 'nickel',
            'Cu': 'copper',
            'Zn': 'zinc',
            'Se': 'selenium',
            'Br': 'bromine',
            'Kr': 'krypton',
            'Ag': 'silver',
            'Au': 'gold',
            'Ba': 'barium',
            'Xe': 'xenon',
            'Rb': 'rubidium'
            }

gsDict = {
            'p1': ['C2', 'N3', 'O4', 'Si2', 'S4', 'Se4', 'Xe6'],
            'p2': ['N2', 'O3', 'Ne5', 'S3', 'Cl4', 'Ar5', 'Se3', 'Kr5', 'Xe5', 'Rb6'],
            'p3': ['N1', 'O2', 'Ne4', 'Na5', 'S2', 'Cl3', 'Ar4', 'K5', 'Se2', 'Kr4', 'Xe4', 'Rb5', 'Br3'],
            'p4': ['O1', 'Ne3', 'Na4', 'Mg5', 'Cl2', 'Ar3', 'K4', 'Ca5', 'Kr3', 'Xe3', 'Rb4'],
            'p5': ['Ne2', 'Na3', 'Si6', 'Cl1', 'Ar2', 'K3', 'Ba4'],
            's1': ['C4', 'N5', 'Ba2'],
            's2': ['C3', 'N4', 'O5', 'Si3'],
            'd2': ['Fe7'],
            'd3': ['Fe6'],
            'd4': ['Fe5'],
            'd5': ['Fe4'],
            'd6': ['Fe3'],
            'd8': ['Ni3']
            }

def gsFromAtom(atom):
    result = 'unknown'
    for gs in gsDict:
        if atom in gsDict[gs]:
            result = gs
    return result
    
gsLevelDict = {
            'p1': ['$^2$P$_{1/2}$', '$^2$P$_{3/2}$', '$^4$P$_{1/2}$', '$^4$P$_{3/2}$', '$^4$P$_{5/2}$', '$^2$D$_{3/2}$', '$^2$D$_{5/2}$', '$^2$S$_{1/2}$'],
            'p2': ['$^3$P$_0$', '$^3$P$_1$', '$^3$P$_2$', '$^1$D$_2$', '$^1$S$_0$', '$^5$S$_2$'],
            'p3': ['$^4$S$_{3/2}$', '$^2$D$_{3/2}$', '$^2$D$_{5/2}$', '$^2$P$_{1/2}$', '$^2$P$_{3/2}$', '$^4$P$_{5/2}$', '$^4$P$_{3/2}$', '$^4$P$_{1/2}$'],
            'p4': ['$^3$P$_0$', '$^3$P$_1$', '$^3$P$_2$', '$^1$D$_2$', '$^1$S$_0$'],
            'p5': ['$^2$P$_{3/2}$', '$^2$P$_{1/2}$'],
            's1': ['$^2$S$_{1/2}$', '$^2$P$_{1/2}$', '$^2$P$_{3/2}$'],
            's2': ['$^1$S$_0$', '$^3$P$_0$', '$^3$P$_1$', '$^3$P$_2$', '$^1$P$_1$'],
            'd2': ['$^3$F$_{2}$', '$^3$F$_{3}$', '$^3$F$_{4}$','$^1$D$_{2}$', '$^3$P$_{0}$', '$^3$P$_{1}$', '$^3$P$_{2}$', '$^1$G$_{4}$', '$^1$S$_{0}$'], 
            'd3': ['$^4$F$_{3/2}$', '$^4$F$_{5/2}$', '$^4$F$_{7/2}$','$^4$F$_{9/2}$', '$^4$P$_{1/2}$', '$^4$P$_{3/2}$', '$^4$P$_{5/2}$', '$^2$G$_{7/2}$', '$^2$G$_{9/2}$', '$^2$P$_{3/2}$', '$^2$P$_{1/2}$', '$^2$D2$_{5/2}$', '$^2$D2$_{3/2}$', '$^2$H$_{9/2}$', '$^2$H$_{11/2}$', '$^2$F$_{7/2}$', '$^2$F$_{5/2}$', '$^2$D1$_{5/2}$', '$^2$D1$_{3/2}$'], 
            'd6': ['$^5$D$_4$', '$^5$D$_3$', '$^5$D$_2$', '$^5$D$_1$', '$^5$D$_0$', '$^3$P2$_2$', '$^3$H$_6$'], 
            'd8': ['$3$F$_4$', '$^3$F$_3$', '$^3$F$_2$', '$^1$D$_2$', '$^3$P$_2$', '$^3$P2$_1$', '$^3$P$_0$']
            }


_predefinedDataFileDict = {'IRAF_09_orig':
                           {'Al2': {'atom': 'al_ii_atom.fits', 'coll': 'al_ii_coll.fits'},
                            'Ar2': {'atom': 'ar_ii_atom.fits', 'coll': 'ar_ii_coll.fits'},
                            'Ar3': {'atom': 'ar_iii_atom.fits', 'coll': 'ar_iii_coll.fits'},
                            'Ar4': {'atom': 'ar_iv_atom.fits', 'coll': 'ar_iv_coll.fits'},
                            'Ar5': {'atom': 'ar_v_atom.fits', 'coll': 'ar_v_coll.fits'},
                            'C1': {'atom': 'c_i_atom.fits', 'coll': 'c_i_coll.fits'},
                            'C2': {'atom': 'c_ii_atom.fits', 'coll': 'c_ii_coll.fits'},
                            'C3': {'atom': 'c_iii_atom.fits', 'coll': 'c_iii_coll.fits'},
                            'Ca5': {'atom': 'ca_v_atom.fits', 'coll': 'ca_v_coll.fits'},
                            'Cl3': {'atom': 'cl_iii_atom.fits', 'coll': 'cl_iii_coll.fits'},
                            'Cl4': {'atom': 'cl_iv_atom.fits', 'coll': 'cl_iv_coll.fits'},
                            'K4': {'atom': 'k_iv_atom.fits', 'coll': 'k_iv_coll.fits'},
                            'K5': {'atom': 'k_v_atom.fits', 'coll': 'k_v_coll.fits'},
                            'Mg5': {'atom': 'mg_v_atom.fits', 'coll': 'mg_v_coll.fits'},
                            'Mg7': {'atom': 'mg_vii_atom.fits', 'coll': 'mg_vii_coll.fits'},
                            'N1': {'atom': 'n_i_atom.fits', 'coll': 'n_i_coll.fits'},
                            'N2': {'atom': 'n_ii_atom.fits', 'coll': 'n_ii_coll.fits'},
                            'N3': {'atom': 'n_iii_atom.fits', 'coll': 'n_iii_coll.fits'},
                            'N4': {'atom': 'n_iv_atom.fits', 'coll': 'n_iv_coll.fits'},
                            'Na4': {'atom': 'na_iv_atom.fits', 'coll': 'na_iv_coll.fits'},
                            'Na6': {'atom': 'na_vi_atom.fits', 'coll': 'na_vi_coll.fits'},
                            'Ne2': {'atom': 'ne_ii_atom.fits', 'coll': 'ne_ii_coll.fits'},
                            'Ne3': {'atom': 'ne_iii_atom.fits', 'coll': 'ne_iii_coll.fits'},
                            'Ne4': {'atom': 'ne_iv_atom.fits', 'coll': 'ne_iv_coll.fits'},
                            'Ne5': {'atom': 'ne_v_atom.fits', 'coll': 'ne_v_coll.fits'},
                            'Ne6': {'atom': 'ne_vi_atom.fits', 'coll': 'ne_vi_coll.fits'},
                            'O1': {'atom': 'o_i_atom.fits', 'coll': 'o_i_coll.fits'},
                            'O2': {'atom': 'o_ii_atom.fits', 'coll': 'o_ii_coll.fits'},
                            'O3': {'atom': 'o_iii_atom.fits', 'coll': 'o_iii_coll.fits'},
                            'O4': {'atom': 'o_iv_atom.fits', 'coll': 'o_iv_coll.fits'},
                            'O5': {'atom': 'o_v_atom.fits', 'coll': 'o_v_coll.fits'},
                            'S2': {'atom': 's_ii_atom.fits', 'coll': 's_ii_coll.fits'},
                            'S3': {'atom': 's_iii_atom.fits', 'coll': 's_iii_coll.fits'},
                            'S4': {'atom': 's_iv_atom.fits', 'coll': 's_iv_coll.fits'},
                            'Si2': {'atom': 'si_ii_atom.fits', 'coll': 'si_ii_coll.fits'},
                            'Si3': {'atom': 'si_iii_atom.fits', 'coll': 'si_iii_coll.fits'},
                            },
                         'IRAF_09':
                           {'H1': {'rec': 'h_i_rec_SH95.fits'},
                            'He1': {'rec': 'he_i_rec_Pal12-Pal13.fits'},
                            'He2': {'rec': 'he_ii_rec_SH95.fits'},
                            'Al2': {'atom': 'al_ii_atom_JSP86-HK87-VVF96-KS86.fits', 'coll': 'al_ii_coll_KHAF92-TBK85-TBK84.fits'},
                            'Ar2': {'atom': 'ar_ii_atom_Bal06.fits', 'coll': 'ar_ii_coll_PB95.fits'},
                            'Ar3': {'atom': 'ar_iii_atom_B60-M83-KS86.fits', 'coll': 'ar_iii_coll_GMZ95.fits'},
                            'Ar4': {'atom': 'ar_iv_atom_B60-MZ82a-KS86.fits', 'coll': 'ar_iv_coll_ZBL87.fits'},
                            'Ar5': {'atom': 'ar_v_atom_B60-LL93-MZ82-KS86.fits', 'coll': 'ar_v_coll_GMZ95.fits'},
                            'C1': {'atom': 'c_i_atom_B60-NS84-FFS85.fits', 'coll': 'c_i_coll_JBK87-PA76.fits'},
                            'C2': {'atom': 'c_ii_atom_MG93-PP95-NS81-GMZ98.fits', 'coll': 'c_ii_coll_BP92.fits'},
                            'C3': {'atom': 'c_iii_atom_BK93-G83-NS78-WFD96.fits', 'coll': 'c_iii_coll_Bal85.fits'},
                            'Ca5': {'atom': 'ca_v_atom_B60-M83-KS86.fits', 'coll': 'ca_v_coll_GMZ95.fits'},
                            'Cl3': {'atom': 'cl_iii_atom_M83-KS86.fits', 'coll': 'cl_iii_coll_BZ89.fits'},
                            'Cl4': {'atom': 'cl_iv_atom_B60-H85-KS86-MZ82-EM84.fits', 'coll': 'cl_iv_coll_GMZ95.fits'},
                            'K4': {'atom': 'k_iv_atom_B60-M83-KS86.fits', 'coll': 'k_iv_coll_GMZ95.fits'},
                            'K5': {'atom': 'k_v_atom_B60-M83-KS86.fits', 'coll': 'k_v_coll_BZL88.fits'},
                            'Mg5': {'atom': 'mg_v_atom_Bal06-B60-GMZ97.fits', 'coll': 'mg_v_coll_BZ94.fits'},
                            'Mg7': {'atom': 'mg_vii_atom_Bal06-B60-BD95-GMZ97.fits', 'coll': 'mg_vii_coll_LB94-U.fits'},
                            'N1': {'atom': 'n_i_atom_B60-KS86-F75-WFD96.fits', 'coll': 'n_i_coll_PA76-DMR76.fits'},
                            'N2': {'atom': 'n_ii_atom_MG93-B60-PP95-GMZ97-WFD96.fits', 'coll': 'n_ii_coll_LB94.fits'},
                            'N3': {'atom': 'n_iii_atom_MG93-BFFJ95-BP92-GMZ98.fits', 'coll': 'n_iii_coll_BP92.fits'},
                            'N4': {'atom': 'n_iv_atom_NS79-WFD96.fits', 'coll': 'n_iv_coll_RBHB94.fits'},
                            'Na4': {'atom': 'na_iv_atom_Bal06-B60-GMZ97.fits', 'coll': 'na_iv_coll_BZ94.fits'},
                            'Na6': {'atom': 'na_vi_atom_Bal06-B60-GMZ97.fits', 'coll': 'na_vi_coll_LB94.fits'},
                            'Ne2': {'atom': 'ne_ii_atom_Bal06.fits', 'coll': 'ne_ii_coll_GMB01.fits'},
                            'Ne3': {'atom': 'ne_iii_atom_Bal06-B60-GMZ97.fits', 'coll': 'ne_iii_coll_BZ94.fits'},
                            'Ne4': {'atom': 'ne_iv_atom_E84-BBZ89-BK88.fits', 'coll': 'ne_iv_coll_G81.fits'},
                            'Ne5': {'atom': 'ne_v_atom_Bal06-B60-M83-GMZ97-U-BD93.fits', 'coll': 'ne_v_coll_LB94.fits'},
                            'Ne6': {'atom': 'ne_vi_atom_Bal06-KS86-MVGK95.fits', 'coll': 'ne_vi_coll_ZGP94.fits'},
                            'O1': {'atom': 'o_i_atom_M83b-WFD96.fits', 'coll': 'o_i_coll_BK95.fits'},
                            'O2': {'atom': 'o_ii_atom_B60-KS86-F75-WFD96.fits', 'coll': 'o_ii_coll_P76-McLB93-v1.fits'},
                            'O3': {'atom': 'o_iii_atom_MG93-B60-M85-GMZ97-WFD96.fits', 'coll': 'o_iii_coll_LB94.fits'},
                            'O4': {'atom': 'o_iv_atom_U-U-GMZ98b.fits', 'coll': 'o_iv_coll_BP92.fits'},
                            'O5': {'atom': 'o_v_atom_BJ68-B80-H80-NS79.fits', 'coll': 'o_v_coll_BBDK85.fits'},
                            'S2': {'atom': 's_ii_atom_B60-VVF96-KHOC93.fits', 'coll': 's_ii_coll_RBS96.fits'},
                            'S3': {'atom': 's_iii_atom_B60-Sal84-LL93-HSC95-MZ82b-KS86.fits', 'coll': 's_iii_coll_GMZ95.fits'},
                            'S4': {'atom': 's_iv_atom_BDF80-JKD86-DHKD82.fits', 'coll': 's_iv_coll_DHKD82.fits'},
                            'Si2': {'atom': 'si_ii_atom_Dal91-BL93-CSB93-N77.fits', 'coll': 'si_ii_coll_DK91.fits'},
                            'Si3': {'atom': 'si_iii_atom_WL95-M83-OKH88-FW90-KS86.fits', 'coll': 'si_iii_coll_DK94.fits'},
                            },
                           'PYNEB_13_01':
                            {'H1': {'rec': 'h_i_rec_SH95.fits'},
                             'He1': {'rec': 'he_i_rec_Pal12-Pal13.fits'},
                             'He2': {'rec': 'he_ii_rec_SH95.fits'},
                             'Al2': {'atom': 'al_ii_atom_JSP86-HK87-VVF96-KS86.fits', 'coll': 'al_ii_coll_KHAF92-TBK85-TBK84.fits'},
                             'Ar2': {'atom': 'ar_ii_atom_Bal06.fits', 'coll': 'ar_ii_coll_PB95.fits'},
                             'Ar3': {'atom': 'ar_iii_atom_B60-M83-KS86.fits', 'coll': 'ar_iii_coll_GMZ95.fits'},
                             'Ar4': {'atom': 'ar_iv_atom_MZ82.fits', 'coll': 'ar_iv_coll_RB97.fits'},
                             'Ar5': {'atom': 'ar_v_atom_B60-LL93-MZ82-KS86.fits', 'coll': 'ar_v_coll_GMZ95.fits'},
                             'Ba2': {'atom': 'ba_ii_atom_C04.fits', 'coll': 'ba_ii_coll_SB98.fits'},
                             'Ba4': {'atom': 'ba_iv_atom_SC10-BHQZ95.fits', 'coll': 'ba_iv_coll_SB98.fits'},
                             'C1': {'atom': 'c_i_atom_B60-NS84-FFS85.fits', 'coll': 'c_i_coll_JBK87-PA76.fits'},
                             'C2': {'atom': 'c_ii_atom_MG93-PP95-NS81-GMZ98.fits', 'coll': 'c_ii_coll_BP92.fits'},
                             'C3': {'atom': 'c_iii_atom_BK93-G83-NS78-WFD96.fits', 'coll': 'c_iii_coll_Bal85.fits'},
                             'C4': {'atom': 'c_iv_atom_M83-WFD96.fits', 'coll': 'c_iv_coll_AK04.fits'},
                             'Ca5': {'atom': 'ca_v_atom_B60-M83-KS86.fits', 'coll': 'ca_v_coll_GMZ95.fits'},
                             'Cl2': {'atom': 'cl_ii_atom_RK74-MZ83.fits', 'coll': 'cl_ii_coll_T04.fits'},
                             'Cl3': {'atom': 'cl_iii_atom_M83-KS86.fits', 'coll': 'cl_iii_coll_BZ89.fits'},
                             'Cl4': {'atom': 'cl_iv_atom_B60-H85-KS86-MZ82-EM84.fits', 'coll': 'cl_iv_coll_GMZ95.fits'},
                             'Fe3': {'atom': 'fe_iii_atom_Q96.fits', 'coll': 'fe_iii_coll_Q96.fits'},
                             'K4': {'atom': 'k_iv_atom_B60-M83-KS86.fits', 'coll': 'k_iv_coll_GMZ95.fits'},
                             'K5': {'atom': 'k_v_atom_B60-M83-KS86.fits', 'coll': 'k_v_coll_BZL88.fits'},
                             'Mg5': {'atom': 'mg_v_atom_Bal06-B60-GMZ97.fits', 'coll': 'mg_v_coll_BZ94.fits'},
                             'Mg7': {'atom': 'mg_vii_atom_Bal06-B60-BD95-GMZ97.fits', 'coll': 'mg_vii_coll_LB94-U.fits'},
                             'N1': {'atom': 'n_i_atom_B60-KS86-F75-WFD96.fits', 'coll': 'n_i_coll_PA76-DMR76.fits'},
                             'N2': {'atom': 'n_ii_atom_MG93-B60-PP95-GMZ97-WFD96.fits', 'coll': 'n_ii_coll_T11.fits'},
                             'N3': {'atom': 'n_iii_atom_MG93-BFFJ95-BP92-GMZ98.fits', 'coll': 'n_iii_coll_BP92.fits'},
                             'N4': {'atom': 'n_iv_atom_NS79-WFD96.fits', 'coll': 'n_iv_coll_RBHB94.fits'},
                             'Na4': {'atom': 'na_iv_atom_Bal06-B60-GMZ97.fits', 'coll': 'na_iv_coll_BZ94.fits'},
                             'Na6': {'atom': 'na_vi_atom_Bal06-B60-GMZ97.fits', 'coll': 'na_vi_coll_LB94.fits'},
                             'Ne2': {'atom': 'ne_ii_atom_Bal06.fits', 'coll': 'ne_ii_coll_GMB01.fits'},
                             'Ne3': {'atom': 'ne_iii_atom_Bal06-B60-GMZ97.fits', 'coll': 'ne_iii_coll_McLB00.fits'},
                             'Ne4': {'atom': 'ne_iv_atom_E84-BBZ89-BK88.fits', 'coll': 'ne_iv_coll_G81.fits'},
                             'Ne5': {'atom': 'ne_v_atom_Bal06-B60-M83-GMZ97-U-BD93.fits', 'coll': 'ne_v_coll_LB94.fits'},
                             'Ne6': {'atom': 'ne_vi_atom_Bal06-KS86-MVGK95.fits', 'coll': 'ne_vi_coll_ZGP94.fits'},
                             'O1': {'atom': 'o_i_atom_M83b-WFD96.fits', 'coll': 'o_i_coll_BK95.fits'},
                             'O2': {'atom': 'o_ii_atom_B60-KS86-F75-Z82-WFD96.fits', 'coll': 'o_ii_coll_P06-T07.fits'},
                             'O3': {'atom': 'o_iii_atom_B60-M85-WFD96-SZ00-WFD96.fits', 'coll': 'o_iii_coll_AK99.fits'},
                             'O4': {'atom': 'o_iv_atom_U-U-GMZ98b.fits', 'coll': 'o_iv_coll_BP92.fits'},
                             'O5': {'atom': 'o_v_atom_BJ68-B80-H80-NS79.fits', 'coll': 'o_v_coll_BBDK85.fits'},
                             'S2': {'atom': 's_ii_atom_B60-VVF96-PKW09.fits', 'coll': 's_ii_coll_TZ10.fits'},
                             'S3': {'atom': 's_iii_atom_B60-Sal84-PKW09.fits', 'coll': 's_iii_coll_GMZ95.fits'},
                             'S4': {'atom': 's_iv_atom_BDF80-JKD86-DHKD82.fits', 'coll': 's_iv_coll_DHKD82.fits'},
                             'Si2': {'atom': 'si_ii_atom_Dal91-BL93-CSB93-N77.fits', 'coll': 'si_ii_coll_DK91.fits'},
                             'Si3': {'atom': 'si_iii_atom_WL95-M83-OKH88-FW90-KS86.fits', 'coll': 'si_iii_coll_DK94.fits'},
                             'Xe3': {'atom': 'xe_iii_atom_M93-BHQZ95.fits', 'coll': 'xe_iii_coll_SB98.fits'},
                             'Xe4': {'atom': 'xe_iv_atom_S04-BHQZ95.fits', 'coll': 'xe_iv_coll_SB98.fits'}
                             },
                           'PYNEB_14_01':
                            {'H1': {'rec': 'h_i_rec_SH95.fits'},
                             'He1': {'rec': 'he_i_rec_Pal12-Pal13.fits'},
                             'He2': {'rec': 'he_ii_rec_SH95.fits'},
                             'Al2': {'atom': 'al_ii_atom_JSP86-HK87-VVF96-KS86.fits', 'coll': 'al_ii_coll_KHAF92-TBK85-TBK84.fits'},
                             'Ar2': {'atom': 'ar_ii_atom_Bal06.fits', 'coll': 'ar_ii_coll_PB95.fits'},
                             'Ar3': {'atom': 'ar_iii_atom_B60-M83-KS86.fits', 'coll': 'ar_iii_coll_GMZ95.fits'},
                             'Ar4': {'atom': 'ar_iv_atom_MZ82.fits', 'coll': 'ar_iv_coll_RB97.fits'},
                             'Ar5': {'atom': 'ar_v_atom_B60-LL93-MZ82-KS86.fits', 'coll': 'ar_v_coll_GMZ95.fits'},
                             #'Ba2': {'atom': 'ba_ii_atom_C04.fits', 'coll': 'ba_ii_coll_SB98.fits'},
                             #'Ba4': {'atom': 'ba_iv_atom_SC10-BHQZ95.fits', 'coll': 'ba_iv_coll_SB98.fits'},
                             'C1': {'atom': 'c_i_atom_B60-NS84-FFS85.fits', 'coll': 'c_i_coll_JBK87-PA76.fits'},
                             'C2': {'atom': 'c_ii_atom_MG93-PP95-NS81-GMZ98.fits', 'coll': 'c_ii_coll_BP92.fits'},
                             'C3': {'atom': 'c_iii_atom_BK93-G83-NS78-WFD96.fits', 'coll': 'c_iii_coll_Bal85.fits'},
                             'C4': {'atom': 'c_iv_atom_M83-WFD96.fits', 'coll': 'c_iv_coll_AK04.fits'},
                             'Ca5': {'atom': 'ca_v_atom_B60-M83-KS86.fits', 'coll': 'ca_v_coll_GMZ95.fits'},
                             'Cl2': {'atom': 'cl_ii_atom_RK74-MZ83.fits', 'coll': 'cl_ii_coll_T04.fits'},
                             'Cl3': {'atom': 'cl_iii_atom_M83-KS86.fits', 'coll': 'cl_iii_coll_BZ89.fits'},
                             'Cl4': {'atom': 'cl_iv_atom_B60-H85-KS86-MZ82-EM84.fits', 'coll': 'cl_iv_coll_GMZ95.fits'},
                             #'Fe3': {'atom': 'fe_iii_atom_Z96_Q96_J00.dat', 'coll': 'fe_iii_coll_Z96.dat'},
                             'K4': {'atom': 'k_iv_atom_B60-M83-KS86.fits', 'coll': 'k_iv_coll_GMZ95.fits'},
                             'K5': {'atom': 'k_v_atom_B60-M83-KS86.fits', 'coll': 'k_v_coll_BZL88.fits'},
                             'Mg5': {'atom': 'mg_v_atom_Bal06-B60-GMZ97.fits', 'coll': 'mg_v_coll_BZ94.fits'},
                             'Mg7': {'atom': 'mg_vii_atom_Bal06-B60-BD95-GMZ97.fits', 'coll': 'mg_vii_coll_LB94-U.fits'},
                             'N1': {'atom': 'n_i_atom_B60-KS86-F75-WFD96.fits', 'coll': 'n_i_coll_PA76-DMR76.fits'},
                             'N2': {'atom': 'n_ii_atom_MG93-B60-PP95-GMZ97-WFD96.fits', 'coll': 'n_ii_coll_T11.fits'},
                             'N3': {'atom': 'n_iii_atom_MG93-BFFJ95-BP92-GMZ98.fits', 'coll': 'n_iii_coll_BP92.fits'},
                             'N4': {'atom': 'n_iv_atom_NS79-WFD96.fits', 'coll': 'n_iv_coll_RBHB94.fits'},
                             'Na4': {'atom': 'na_iv_atom_Bal06-B60-GMZ97.fits', 'coll': 'na_iv_coll_BZ94.fits'},
                             'Na6': {'atom': 'na_vi_atom_Bal06-B60-GMZ97.fits', 'coll': 'na_vi_coll_LB94.fits'},
                             'Ne2': {'atom': 'ne_ii_atom_Bal06.fits', 'coll': 'ne_ii_coll_GMB01.fits'},
                             'Ne3': {'atom': 'ne_iii_atom_Bal06-B60-GMZ97.fits', 'coll': 'ne_iii_coll_McLB00.fits'},
                             'Ne4': {'atom': 'ne_iv_atom_E84-BBZ89-BK88.fits', 'coll': 'ne_iv_coll_G81.fits'},
                             'Ne5': {'atom': 'ne_v_atom_Bal06-B60-M83-GMZ97-U-BD93.fits', 'coll': 'ne_v_coll_LB94.fits'},
                             'Ne6': {'atom': 'ne_vi_atom_Bal06-KS86-MVGK95.fits', 'coll': 'ne_vi_coll_ZGP94.fits'},
                             'O1': {'atom': 'o_i_atom_M83b-WFD96.fits', 'coll': 'o_i_coll_BK95.fits'},
                             'O2': {'atom': 'o_ii_atom_B60-KS86-F75-Z82-WFD96.fits', 'coll': 'o_ii_coll_P06-T07.fits'},
                             'O3': {'atom': 'o_iii_atom_B60-M85-WFD96-SZ00-WFD96.fits', 'coll': 'o_iii_coll_AK99.fits'},
                             'O4': {'atom': 'o_iv_atom_U-U-GMZ98b.fits', 'coll': 'o_iv_coll_BP92.fits'},
                             'O5': {'atom': 'o_v_atom_BJ68-B80-H80-NS79.fits', 'coll': 'o_v_coll_BBDK85.fits'},
                             'S2': {'atom': 's_ii_atom_B60-VVF96-PKW09.fits', 'coll': 's_ii_coll_TZ10.fits'},
                             'S3': {'atom': 's_iii_atom_B60-Sal84-PKW09.fits', 'coll': 's_iii_coll_GMZ95.fits'},
                             'S4': {'atom': 's_iv_atom_BDF80-JKD86-DHKD82.fits', 'coll': 's_iv_coll_DHKD82.fits'},
                             'Si2': {'atom': 'si_ii_atom_Dal91-BL93-CSB93-N77.fits', 'coll': 'si_ii_coll_DK91.fits'},
                             'Si3': {'atom': 'si_iii_atom_WL95-M83-OKH88-FW90-KS86.fits', 'coll': 'si_iii_coll_DK94.fits'},
                             'Xe3': {'atom': 'xe_iii_atom_BHQZ95.dat', 'coll': 'xe_iii_coll_SB98.fits'},
                             'Xe4': {'atom': 'xe_iv_atom_BHQZ95.dat', 'coll': 'xe_iv_coll_SB98.fits'},
                             'Xe6': {'atom': 'xe_vi_atom_BHQZ95.dat', 'coll': 'xe_vi_coll_SB98.fits'},
                             'Kr3': {'atom': 'kr_iii_atom_BH86.dat', 'coll':'kr_iii_coll_S97.dat'},
                             'Kr4': {'atom': 'kr_iv_atom_BH86.dat', 'coll': 'kr_iv_coll_S97.dat'},
                             'Kr5': {'atom': 'kr_v_atom_BH86.dat', 'coll': 'kr_v_coll_S97.dat'},
                             'Se3': {'atom': 'se_iii_atom_BH86.dat', 'coll': 'se_iii_coll_S97.dat'},
                             'Se4': {'atom': 'se_iv_atom_B05.dat', 'coll': 'se_iv_coll_B05.dat'},
                             'Br3': {'atom': 'br_iii_atom_BH86.dat', 'coll': 'br_iii_coll_S97.dat'},
                             'Br4': {'atom': 'br_iv_atom_BH86.dat', 'coll': 'br_iv_coll_S97.dat'},
                             'Rb4': {'atom': 'rb_iv_atom_BH86.dat', 'coll': 'rb_iv_coll_S97.dat'},
                             'Rb5': {'atom': 'rb_v_atom_BH86.dat', 'coll': 'rb_v_coll_S97.dat'},
                             'Rb6': {'atom': 'rb_vi_atom_BH86.dat', 'coll': 'rb_vi_coll_S97.dat'}
                             },
                           'PYNEB_14_02':
                            {'H1': {'rec': 'h_i_rec_SH95.fits'},
                             'He1': {'rec': 'he_i_rec_Pal12-Pal13.fits'},
                             'He2': {'rec': 'he_ii_rec_SH95.fits'},
                             'Al2': {'atom': 'al_ii_atom_JSP86-HK87-VVF96-KS86.dat', 'coll': 'al_ii_coll_KHAF92-TBK85-TBK84.dat'},
                             'Ar2': {'atom': 'ar_ii_atom_Bal06.dat', 'coll': 'ar_ii_coll_PB95.dat'},
                             'Ar3': {'atom': 'ar_iii_atom_M83-KS86.dat', 'coll': 'ar_iii_coll_GMZ95.dat'},
                             'Ar4': {'atom': 'ar_iv_atom_MZ82.dat', 'coll': 'ar_iv_coll_RB97.dat'},
                             'Ar5': {'atom': 'ar_v_atom_LL93-MZ82-KS86.dat', 'coll': 'ar_v_coll_GMZ95.dat'},
                             'Ba2': {'atom': 'ba_ii_atom_C04.dat', 'coll': 'ba_ii_coll_SB98.dat'},
                             'Ba4': {'atom': 'ba_iv_atom_BHQZ95.dat', 'coll': 'ba_iv_coll_SB98.dat'},
                             'C1': {'atom': 'c_i_atom_FFS85.dat', 'coll': 'c_i_coll_JBK87-PA76.dat'},
                             'C2': {'atom': 'c_ii_atom_GMZ98.dat', 'coll': 'c_ii_coll_BP92.dat'},
                             'C3': {'atom': 'c_iii_atom_G83-NS78-WFD96.dat', 'coll': 'c_iii_coll_Bal85.dat'},
                             'C4': {'atom': 'c_iv_atom_WFD96.dat', 'coll': 'c_iv_coll_AK04.dat'},
                             'Ca5': {'atom': 'ca_v_atom_M83-KS86.dat', 'coll': 'ca_v_coll_GMZ95.dat'},
                             'Cl2': {'atom': 'cl_ii_atom_MZ83.dat', 'coll': 'cl_ii_coll_T04.dat'},
                             'Cl3': {'atom': 'cl_iii_atom_M83-KS86.dat', 'coll': 'cl_iii_coll_BZ89.dat'},
                             'Cl4': {'atom': 'cl_iv_atom_KS86-MZ82-EM84.dat', 'coll': 'cl_iv_coll_GMZ95.dat'},
                             'Fe3': {'atom': 'fe_iii_atom_Q96_J00.dat', 'coll': 'fe_iii_coll_Z96.dat'},
                             'K4': {'atom': 'k_iv_atom_M83-KS86.dat', 'coll': 'k_iv_coll_GMZ95.dat'},
                             'K5': {'atom': 'k_v_atom_M83-KS86.dat', 'coll': 'k_v_coll_BZL88.dat'},
                             'Mg5': {'atom': 'mg_v_atom_GMZ97.dat', 'coll': 'mg_v_coll_BZ94.dat'},
                             'Mg7': {'atom': 'mg_vii_atom_GMZ97.dat', 'coll': 'mg_vii_coll_LB94-U.dat'},
                             'N1': {'atom': 'n_i_atom_KS86-WFD96.dat', 'coll': 'n_i_coll_PA76-DMR76.dat'},
                             'N2': {'atom': 'n_ii_atom_GMZ97-WFD96.dat', 'coll': 'n_ii_coll_T11.dat'},
                             'N3': {'atom': 'n_iii_atom_GMZ98.dat', 'coll': 'n_iii_coll_BP92.dat'},
                             'N4': {'atom': 'n_iv_atom_WFD96.dat', 'coll': 'n_iv_coll_RBHB94.dat'},
                             'Na4': {'atom': 'na_iv_atom_GMZ97.dat', 'coll': 'na_iv_coll_BZ94.dat'},
                             'Na6': {'atom': 'na_vi_atom_GMZ97.dat', 'coll': 'na_vi_coll_LB94.dat'},
                             'Ne2': {'atom': 'ne_ii_atom_Bal06.dat', 'coll': 'ne_ii_coll_GMB01.dat'},
                             'Ne3': {'atom': 'ne_iii_atom_GMZ97.dat', 'coll': 'ne_iii_coll_McLB00.dat'},
                             'Ne4': {'atom': 'ne_iv_atom_BBZ89-BK88.dat', 'coll': 'ne_iv_coll_G81.dat'},
                             'Ne5': {'atom': 'ne_v_atom_GMZ97-U-BD93.dat', 'coll': 'ne_v_coll_LB94.dat'},
                             'Ne6': {'atom': 'ne_vi_atom_GMZ98.dat', 'coll': 'ne_vi_coll_ZGP94.dat'},
                             'Ni3': {'atom': 'ni_iii_atom_B01.dat', 'coll': 'ni_iii_coll_B01.dat'},
                             'O1': {'atom': 'o_i_atom_WFD96.dat', 'coll': 'o_i_coll_BK95.dat'},
                             'O2': {'atom': 'o_ii_atom_Z82-WFD96.dat', 'coll': 'o_ii_coll_P06-T07.dat'},
                             'O3': {'atom': 'o_iii_atom_SZ00-WFD96.dat', 'coll': 'o_iii_coll_AK99.dat'},
                             'O4': {'atom': 'o_iv_atom_GMZ98.dat', 'coll': 'o_iv_coll_BP92.dat'},
                             'O5': {'atom': 'o_v_atom_H80-NS79.dat', 'coll': 'o_v_coll_BBDK85.dat'},
                             'S2': {'atom': 's_ii_atom_PKW09.dat', 'coll': 's_ii_coll_TZ10.dat'},
                             'S3': {'atom': 's_iii_atom_PKW09.dat', 'coll': 's_iii_coll_GMZ95.dat'},
                             'S4': {'atom': 's_iv_atom_JKD86-DHKD82.dat', 'coll': 's_iv_coll_DHKD82.dat'},
                             'Si2': {'atom': 'si_ii_atom_BL93-CSB93-N77.dat', 'coll': 'si_ii_coll_DK91.dat'},
                             'Si3': {'atom': 'si_iii_atom_M83-OKH88-FW90-KS86.dat', 'coll': 'si_iii_coll_DK94.dat'},
                             'Xe3': {'atom': 'xe_iii_atom_BHQZ95.dat', 'coll': 'xe_iii_coll_SB98.dat'},
                             'Xe4': {'atom': 'xe_iv_atom_BHQZ95.dat', 'coll': 'xe_iv_coll_SB98.dat'},
                             'Xe6': {'atom': 'xe_vi_atom_BHQZ95.dat', 'coll': 'xe_vi_coll_SB98.dat'},
                             'Kr3': {'atom': 'kr_iii_atom_BH86.dat', 'coll':'kr_iii_coll_S97.dat'},
                             'Kr4': {'atom': 'kr_iv_atom_BH86.dat', 'coll': 'kr_iv_coll_S97.dat'},
                             'Kr5': {'atom': 'kr_v_atom_BH86.dat', 'coll': 'kr_v_coll_S97.dat'},
                             'Se3': {'atom': 'se_iii_atom_BH86.dat', 'coll': 'se_iii_coll_S97.dat'},
                             'Se4': {'atom': 'se_iv_atom_B05.dat', 'coll': 'se_iv_coll_B05.dat'},
                             'Br3': {'atom': 'br_iii_atom_BH86.dat', 'coll': 'br_iii_coll_S97.dat'},
                             'Br4': {'atom': 'br_iv_atom_BH86.dat', 'coll': 'br_iv_coll_S97.dat'},
                             'Rb4': {'atom': 'rb_iv_atom_BH86.dat', 'coll': 'rb_iv_coll_S97.dat'},
                             'Rb5': {'atom': 'rb_v_atom_BH86.dat', 'coll': 'rb_v_coll_S97.dat'},
                             'Rb6': {'atom': 'rb_vi_atom_BH86.dat', 'coll': 'rb_vi_coll_S97.dat'}
                             },
                           'PYNEB_14_03':
                            {'H1': {'rec': 'h_i_rec_SH95.fits', 'trc': 'h_i_trc_SH95-caseB.dat'},
                             'He1': {'rec': 'he_i_rec_Pal12-Pal13.fits'},
                             'He2': {'rec': 'he_ii_rec_SH95.fits', 'trc': 'he_ii_trc_SH95-caseB.dat'},
                             '3He2': {'atom': '3he_ii_atom_cloudy.dat', 'coll': '3he_ii_coll_cloudy.dat'},
                             'Al2': {'atom': 'al_ii_atom_JSP86-HK87-VVF96-KS86.dat', 'coll': 'al_ii_coll_KHAF92-TBK85-TBK84.dat'},
                             'Ar2': {'atom': 'ar_ii_atom_Bal06.dat', 'coll': 'ar_ii_coll_PB95.dat'},
                             'Ar3': {'atom': 'ar_iii_atom_M83-KS86.dat', 'coll': 'ar_iii_coll_GMZ95.dat'},
                             #'Ar3': {'atom': 'ar_iii_atom_MB09.dat', 'coll': 'ar_iii_coll_MB09.dat'},
                             'Ar4': {'atom': 'ar_iv_atom_MZ82.dat', 'coll': 'ar_iv_coll_RB97.dat'},
                             'Ar5': {'atom': 'ar_v_atom_LL93-MZ82-KS86.dat', 'coll': 'ar_v_coll_GMZ95.dat'},
                             'Ba2': {'atom': 'ba_ii_atom_C04.dat', 'coll': 'ba_ii_coll_SB98.dat'},
                             'Ba4': {'atom': 'ba_iv_atom_BHQZ95.dat', 'coll': 'ba_iv_coll_SB98.dat'},
                             'C1': {'atom': 'c_i_atom_FFS85.dat', 'coll': 'c_i_coll_JBK87-PA76.dat'},
                             'C2': {'atom': 'c_ii_atom_GMZ98.dat', 'coll': 'c_ii_coll_BP92.dat'},
                             'C3': {'atom': 'c_iii_atom_G83-NS78-WFD96.dat', 'coll': 'c_iii_coll_Bal85.dat'},
                             'C4': {'atom': 'c_iv_atom_WFD96.dat', 'coll': 'c_iv_coll_AK04.dat'},
                             'Ca5': {'atom': 'ca_v_atom_M83-KS86.dat', 'coll': 'ca_v_coll_GMZ95.dat'},
                             'Cl2': {'atom': 'cl_ii_atom_MZ83.dat', 'coll': 'cl_ii_coll_T04.dat'},
                             'Cl3': {'atom': 'cl_iii_atom_M83-KS86.dat', 'coll': 'cl_iii_coll_BZ89.dat'},
                             'Cl4': {'atom': 'cl_iv_atom_KS86-MZ82-EM84.dat', 'coll': 'cl_iv_coll_GMZ95.dat'},
                             'Fe3': {'atom': 'fe_iii_atom_Q96_J00.dat', 'coll': 'fe_iii_coll_Z96.dat'},
                             'Fe4': {'atom': 'fe_iv_atom_FFRR08.dat', 'coll': 'fe_iv_coll_ZP97.dat'},
                             'Fe5': {'atom': 'fe_v_atom_Nal00.dat', 'coll': 'fe_v_coll_BGMcL07.dat'},
                             'Fe6': {'atom': 'fe_vi_atom_CP00.dat', 'coll': 'fe_vi_coll_CP99.dat'},
                             'Fe7': {'atom': 'fe_vii_atom_WB08.dat', 'coll': 'fe_vii_coll_WB08.dat'},
                             'K4': {'atom': 'k_iv_atom_M83-KS86.dat', 'coll': 'k_iv_coll_GMZ95.dat'},
                             'K5': {'atom': 'k_v_atom_M83-KS86.dat', 'coll': 'k_v_coll_BZL88.dat'},
                             'Mg5': {'atom': 'mg_v_atom_GMZ97.dat', 'coll': 'mg_v_coll_BZ94.dat'},
                             'Mg7': {'atom': 'mg_vii_atom_GMZ97.dat', 'coll': 'mg_vii_coll_LB94-U.dat'},
                             'N1': {'atom': 'n_i_atom_KS86-WFD96.dat', 'coll': 'n_i_coll_PA76-DMR76.dat'},
                             'N2': {'atom': 'n_ii_atom_GMZ97-WFD96.dat', 'coll': 'n_ii_coll_T11.dat'},
                             'N3': {'atom': 'n_iii_atom_GMZ98.dat', 'coll': 'n_iii_coll_BP92.dat'},
                             'N4': {'atom': 'n_iv_atom_WFD96.dat', 'coll': 'n_iv_coll_RBHB94.dat'},
                             'Na4': {'atom': 'na_iv_atom_GMZ97.dat', 'coll': 'na_iv_coll_BZ94.dat'},
                             'Na6': {'atom': 'na_vi_atom_GMZ97.dat', 'coll': 'na_vi_coll_LB94.dat'},
                             'Ne2': {'atom': 'ne_ii_atom_Bal06.dat', 'coll': 'ne_ii_coll_GMB01.dat'},
                             'Ne3': {'atom': 'ne_iii_atom_GMZ97.dat', 'coll': 'ne_iii_coll_McLB00.dat'},
                             'Ne4': {'atom': 'ne_iv_atom_BBZ89-BK88.dat', 'coll': 'ne_iv_coll_G81.dat'},
                             'Ne5': {'atom': 'ne_v_atom_GMZ97-U-BD93.dat', 'coll': 'ne_v_coll_LB94.dat'},
                             'Ne6': {'atom': 'ne_vi_atom_GMZ98.dat', 'coll': 'ne_vi_coll_ZGP94.dat'},
                             'Ni3': {'atom': 'ni_iii_atom_B01.dat', 'coll': 'ni_iii_coll_B01.dat'},
                             'O1': {'atom': 'o_i_atom_WFD96.dat', 'coll': 'o_i_coll_BK95.dat'},
                             'O2': {'atom': 'o_ii_atom_Z82-WFD96.dat', 'coll': 'o_ii_coll_P06-T07.dat'},
                             'O3': {'atom': 'o_iii_atom_SZ00-WFD96.dat', 'coll': 'o_iii_coll_AK99.dat'},
                             'O4': {'atom': 'o_iv_atom_GMZ98.dat', 'coll': 'o_iv_coll_BP92.dat'},
                             'O5': {'atom': 'o_v_atom_H80-NS79.dat', 'coll': 'o_v_coll_BBDK85.dat'},
                             'S2': {'atom': 's_ii_atom_PKW09.dat', 'coll': 's_ii_coll_TZ10.dat'},
                             'S3': {'atom': 's_iii_atom_PKW09.dat', 'coll': 's_iii_coll_TG99.dat'},
                             'S4': {'atom': 's_iv_atom_JKD86-DHKD82.dat', 'coll': 's_iv_coll_DHKD82.dat'},
                             'Si2': {'atom': 'si_ii_atom_BL93-CSB93-N77.dat', 'coll': 'si_ii_coll_DK91.dat'},
                             'Si3': {'atom': 'si_iii_atom_M83-OKH88-FW90-KS86.dat', 'coll': 'si_iii_coll_DK94.dat'},
                             'Xe3': {'atom': 'xe_iii_atom_BHQZ95.dat', 'coll': 'xe_iii_coll_SB98.dat'},
                             'Xe4': {'atom': 'xe_iv_atom_BHQZ95.dat', 'coll': 'xe_iv_coll_SB98.dat'},
                             'Xe6': {'atom': 'xe_vi_atom_BHQZ95.dat', 'coll': 'xe_vi_coll_SB98.dat'},
                             'Kr3': {'atom': 'kr_iii_atom_BH86.dat', 'coll':'kr_iii_coll_S97.dat'},
                             'Kr4': {'atom': 'kr_iv_atom_BH86.dat', 'coll': 'kr_iv_coll_S97.dat'},
                             'Kr5': {'atom': 'kr_v_atom_BH86.dat', 'coll': 'kr_v_coll_S97.dat'},
                             'Se3': {'atom': 'se_iii_atom_BH86.dat', 'coll': 'se_iii_coll_S97.dat'},
                             'Se4': {'atom': 'se_iv_atom_B05.dat', 'coll': 'se_iv_coll_B05.dat'},
                             'Br3': {'atom': 'br_iii_atom_BH86.dat', 'coll': 'br_iii_coll_S97.dat'},
                             'Br4': {'atom': 'br_iv_atom_BH86.dat', 'coll': 'br_iv_coll_S97.dat'},
                             'Rb4': {'atom': 'rb_iv_atom_BH86.dat', 'coll': 'rb_iv_coll_S97.dat'},
                             'Rb5': {'atom': 'rb_v_atom_BH86.dat', 'coll': 'rb_v_coll_S97.dat'},
                             'Rb6': {'atom': 'rb_vi_atom_BH86.dat', 'coll': 'rb_vi_coll_S97.dat'}
                             }
                           } 

_predefinedDataFileDict['PYNEB_16_01'] = _predefinedDataFileDict['PYNEB_14_03']
#_predefinedDataFileDict['PYNEB_16_01']['Fe2'] = {'atom': 'fe_ii_atom_NS88G62.dat', 'coll': 'fe_ii_coll_BP96ZP95.dat'}
 
def airtovac(wave):
    """
Convert air wavelengths to vacuum wavelengths
Parameters
----------
wave: float, array
The wavelength in air [Angstrom]
Returns
-------
Wavelength: array
Wavelength in vacuum [Angstrom]
Notes
-----

.. note:: This function was ported from the IDL Astronomy User's Library.
:IDL - Documentation:
NAME:
AIRTOVAC
PURPOSE:
Convert air wavelengths to vacuum wavelengths
EXPLANATION:
Wavelengths are corrected for the index of refraction of air under
standard conditions. Wavelength values below 2000 A will not be
altered. Uses the IAU standard for conversion given in Morton
(1991 Ap.J. Suppl. 77, 119)
CALLING SEQUENCE:
AIRTOVAC, WAVE
INPUT/OUTPUT:
WAVE - Wavelength in Angstroms, scalar or vector
WAVE should be input as air wavelength(s), it will be
returned as vacuum wavelength(s). WAVE is always converted to
double precision upon return.
EXAMPLE:
If the air wavelength is W = 6056.125 (a Krypton line), then
AIRTOVAC, W yields an vacuum wavelength of W = 6057.8019
METHOD:
See Morton (Ap. J. Suppl. 77, 119) for the formula used
REVISION HISTORY
Written W. Landsman November 1991
Converted to IDL V5.0 W. Landsman September 1997
    """
    
    sigma2 = (1.e4/wave)**2. #Convert to wavenumber squared
    
    # Compute conversion factor
    
    fact = 1. + 6.4328e-5 + 2.94981e-2/(146.-sigma2) + 2.5540e-4/(41.-sigma2)
    
    fact = fact*(wave >= 2000.) + 1.*(wave < 2000.)
    
    wave = wave*fact #Convert Wavelength
    
    return wave



def vactoair(wave):
    """
Convert vacuum wavelengths to air wavelengths

Parameters
----------
wave: float, array
The wavelength in vacuum [Angstrom]
Returns
-------
Wavelength: array,
Wavelength in air [Angstrom]
Notes
-----

.. note:: This function was ported from the IDL Astronomy User's Library.
:IDL - Documentation:

NAME:
VACTOAIR
PURPOSE:
Convert vacuum wavelengths to air wavelengths
EXPLANATION:
Corrects for the index of refraction of air under standard conditions.
Wavelength values below 2000 A will not be altered. Accurate to
about 0.005 A

CALLING SEQUENCE:
VACTOAIR, WAVE

INPUT/OUTPUT:
WAVE - Wavelength in Angstroms, scalar or vector
WAVE should be input as vacuum wavelength(s), it will be
returned as air wavelength(s). WAVE is always converted to
double precision

EXAMPLE:
If the vacuum wavelength is W = 2000, then

IDL> VACTOAIR, W

yields an air wavelength of W = 1999.353 Angstroms

METHOD:
An approximation to the 4th power of inverse wavenumber is used
See IUE Image Processing Manual Page 6-15.

REVISION HISTORY
Written, D. Lindler 1982
Documentation W. Landsman Feb. 1989
Converted to IDL V5.0 W. Landsman September 1997
    """

    wave2 = wave**2.
    fact = 1. + 2.735182e-4 + 131.4182/wave2 + 2.76249e8/(wave2**2.)
    fact = fact * ( wave >= 2000. ) + 1.*( wave < 2000. )
    
    # Convert wavelengths
    
    wave = wave/fact
    
    return wave

