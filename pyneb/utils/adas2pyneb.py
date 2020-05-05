import numpy as np
import pandas as pd
import requests

def download_file(ion):
    """
    Download the file corresponding to the Mao, Badnell & Del Zanna 2020 paper on
    R-matrix electron-impact excitation data for the C-like iso-electronic sequence
    Rename them adas_ion.dat where ion runs on n1, o2, f3, ne4, na5, mg6, al7, and si8
    """
    url_base = "https://open.adas.ac.uk/download/adf04/copaw][c/clike_jm19]["
    url_file = ion + '.dat'
    url = url_base + url_file
    print('Download of', url, 'started')
    r = requests.get(url)
    with open('adas_' + url_file, 'wb') as f:
        f.write(r.content)
    print('Saving', 'adas_' + url_file, 'completed')
        
def make_pyneb(adas_filename, out_filename, ion_str, elem, ion, gsconf, NLevels=6):
    """
    Transform the ADAS adf04 data file into a pyneb file
    """
    use_data = False
    temp_str = None
    temp  = []
    I = []
    J = []
    Oms = []

    with open(adas_filename, 'r') as f:
        for l in f:
            if l == '   -1\n': # the data starts after this
                use_data = True
                continue
            if l == '  -1\n': # no more data after this
                break
            if use_data: # we read the data here
                if temp_str is None: # first data is temperature
                    temp_str = l.split()[2:]
                    for tt in temp_str:
                        temp.append(tt[0:-3] + 'e' + tt[-3:]) #add exponent (missing in original data)
                    continue
                dsplit = l.split() # now run on collision data
                I.append(int(dsplit[0])) # first field
                J.append(int(dsplit[1])) # secind field
                dsplit = dsplit[2:] # rest of the fields
                if len(dsplit[-1]) > 8: # the latest field can be 2 numbers together (because the 2nd one can be negative)
                    dd = dsplit[-1]
                    dsplit[-1] = dd[:7]
                    dsplit.append(dd[7:])
                dsplit2 = []
                for dd in dsplit:
                    dsplit2.append(dd[0:-3] + 'e' + dd[-3:]) # add exponent (missing in original data)
                Oms.append(dsplit2[1:-1])        
    temp = np.array(temp, dtype=float)
    Oms = pd.DataFrame(Oms, dtype=float)
    Oms.insert(loc=0, column='I', value=I)   # I and J are switched in PyNeb format
    Oms.insert(loc=0, column='J', value=J)   
    mask = Oms.I <= NLevels

    with open(out_filename, 'w') as f: # write the PyNeb data file
        f.write(' 0  0 ' + ' '.join(['{:.2e}'.format(float(tt)) for tt in temp]) + '\n')
        f.write(Oms[mask].to_string(index=False, header=False, float_format='{:.2e}'.format))
        f.write("""
*** {} collision strengths 
*** SOURCE1  J. Mao, N. R. Badnell, and G. Del Zanna 2020, A&A, 634, A7
*** NOTE1    Collision strengths
*** ATOM {}
*** SPECTRUM {}
*** N_LEVELS {}
*** GSCONFIG {} 
*** T_UNIT K
""".format(ion_str, elem, ion, NLevels, gsconf))
    print(out_filename, 'created from', adas_filename)

def download_all():
    """
    Download 8 files from ADAS
    """
    for ion in ('n1', 'o2', 'f3', 'ne4', 'na5', 'mg6', 'al7','si8'):
        download_file(ion)
        
def make_all(NLevels=9):
    """
    Write the 8 files in the PyNeb format
    """
    
    make_pyneb('adas_n1.dat','n_ii_coll_MBZ20.dat', 'N II', 'nitrogen', '2', 'p2', NLevels)
    make_pyneb('adas_o2.dat','o_iii_coll_MBZ20.dat', 'O III', 'oxygen', '3', 'p2', NLevels)
    make_pyneb('adas_f3.dat','f_iv_coll_MBZ20.dat', 'F IV', 'fluorine', '4', 'p2', NLevels)
    make_pyneb('adas_ne4.dat','ne_v_coll_MBZ20.dat', 'Ne V', 'neon', '5', 'p2', NLevels)
    make_pyneb('adas_na5.dat','na_vi_coll_MBZ20.dat', 'Na VI', 'sodium', '6', 'p2', NLevels)
    make_pyneb('adas_mg6.dat','mg_vii_coll_MBZ20.dat', 'Mg VII', 'magnesium', '7', 'p2', NLevels)
    make_pyneb('adas_al7.dat','al_viii_coll_MBZ20.dat', 'Al VIII', 'aluminium', '8', 'p2', NLevels)
    make_pyneb('adas_si8.dat','si_ix_coll_MBZ20.dat', 'Si IX', 'silicon', '9', 'p2', NLevels)
