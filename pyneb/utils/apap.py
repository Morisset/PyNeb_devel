import re
import gzip
import os
import numpy as np
import pyneb as pn
from misc import int_to_roman 
from physics import sym2name
from manage_atomic_data import getLevelsNIST

term2_dic = {'0':'S', '1':'P', '2':'D', '3':'F', '4':'G', '5':'H', '6':'I'}

def conv_apap(str_in):
    """
    transform "1.3-01" into 1.3e-01 
    """
    a, b = re.split('\+|-',str_in)
    if '-' in str_in:
        return float(a) / 10**float(b)
    else:
        return float(a) * 10**float(b)

def read_apap(apap_file, max_levels=10):
    if max_levels is None:
        max_levels = np.Inf
        
    term1 = []
    term2 = []
    energy = []
    J = []
    if apap_file.split('.')[-1] == 'gz':
        f = gzip.open(apap_file)
    else:
        f = open(apap_file)
    foo = f.readline()
    while True:
        line = f.readline()
        if line[0:5] == "   -1":
            break
        lev = int(line[0:5])
        if lev <= max_levels:
            term1.append(line[25:26])
            term2.append(line[27:28])
            J.append(line[29:33])
            energy.append(line[34::])
    term = [t1+term2_dic[t2] for t1,t2 in zip(term1, term2)]
    
    lev_i = []
    lev_j = []
    As = []
    Omegas = []
    
    temps_str = f.readline().split()[2::]
    temp_array = np.array([conv_apap(temp) for temp in temps_str])
    N_temp = len(temp_array)
    while True:
        line = f.readline()
        if line[0:4] == "  -1":
            break
        i = int(line[1:4])
        j = int(line[4:8])
        if (i <= max_levels) and (j <= max_levels):
            lev_i.append(i)
            lev_j.append(j)
            strg5 = line[8:8+(N_temp+1)*8].split()
            As.append(conv_apap(strg5[0]))
            Omegas.append([conv_apap(omega) for omega in strg5[1::]])
    f.close()
    
    return(term, 
           np.array(energy, dtype='float'), 
           np.array(J, dtype='float'), 
           np.array(lev_i, dtype='int'), 
           np.array(lev_j, dtype='int'), 
           np.array(As, dtype='float'), 
           np.array(Omegas, dtype='float'), 
           temp_array)

def apap2ascii(apap_file, max_levels=10, source='http://www.apap-network.org/'):
    
    atom = apap_file.split('#')[1].split('.')[0]
    spectrum = int(re.findall(r'\d+', atom)[0]) + 1
    ioniz = int_to_roman(spectrum).lower()
    elem = re.findall('[a-zA-Z]+', atom)[0]
    ion_str = '{0}{1}'.format(elem.capitalize(), spectrum)
    print('Reading NIST data for {0}'.format(ion_str))
    NIST_levels = getLevelsNIST(ion_str)
    NIST_full_term = ['{0:2s}{1:.1f}'.format(l['term'][0:2],float(l['J'])) for l in  NIST_levels]
    
    term, energy, J, lev_i, lev_j, As, Omegas, temp_array = read_apap(apap_file, max_levels=max_levels)
    N_temp = len(temp_array)

    APAP_full_term = ['{0:2s}{1:.1f}'.format(t,float(j)) for t, j in zip(term, J)]
    
    N_levels = len(energy)
    stat_weight = J * 2 + 1
        
    Omega3d = np.zeros((N_levels, N_levels, N_temp))
    for i, j, Omega in zip(lev_i, lev_j, Omegas):
        try:
            i_NIST = NIST_full_term.index(APAP_full_term[i-1])
        except:
            i_NIST = None
        try:
            j_NIST = NIST_full_term.index(APAP_full_term[j-1])
        except:
            j_NIST = None
        if (i_NIST is not None) and (j_NIST is not None):
            Omega3d[i_NIST, j_NIST] = Omega
        else:
            pn.log_.warn('Levels not found in NIST: {0} {1}'.format(APAP_full_term[i-1], 
                                                                    APAP_full_term[j-1]))
            
    f_coll = open('{0}_{1}_coll_APAP14.dat'.format(elem, ioniz),'w')
    str_ = '0   0 '
    str_ += ''.join(' {0:10.4e}'.format(T) for T in temp_array)
    str_ += '\n'
    f_coll.write(str_)
    for i in 1+np.arange(N_levels):
        for j in i+1+np.arange(N_levels-i):
            str_ = '{0} {1} '.format(i, j)
            str_ += ''.join(' {0:10.8e}'.format(O) for O in Omega3d[j-1, i-1])
            str_ += '\n'
            f_coll.write(str_)
    try:
        f_coll.write("*** ATOM {0}\n".format(sym2name[elem.capitalize()]))
    except:
        pn.log_.warn('No full name for element {0}'.format(elem.capitalize()))

    f_coll.write("""*** SPECTRUM {0}
*** N_LEVELS {1}
*** T_UNIT 'K'
*** NOTE1 'Collision strengths'
*** SOURCE1 '{2}'
""".format(spectrum, N_levels, source))
    
    f_coll.close()


def make_all(dir_ = '.', max_levels=8):
    file_liste = os.listdir('.')
    for file_ in file_liste:
        if file_[-1] == 'z':
            atom = file_.split('#')[1].split('.')[0]
            if atom == 'fe':
                apap2ascii(file_, 34)
            else:
                apap2ascii(file_, max_levels)
            print('{0} done'.format(file_))