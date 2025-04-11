"""
Read data from Del Zanna & Storey paper https://academic.oup.com/mnras/article/513/1/1198/6564189 
(MNRAS supplementary data)
Original data file is allem.d.gz
Will be stored in hdf5 format
"""
#%% Imports

import numpy as np
from astropy.table import Table
from astropy.io.misc.hdf5 import write_table_hdf5
import gzip
import pandas as pd
#%% Definition of the function that reads the data
def get_data():

    nd=9
    nt=20
    ntr=1032

    tab_log_dens = []
    tab_wl = []
    tab_tr = []
    tab_levs = []
    with gzip.open('allem.d.gz', 'r') as f: 
        first_iter = True
        for _i_d in range(nd):
            h1 = f.readline().split()
            tab_log_dens.append(float(h1[1].split(b'=')[1]))
            if first_iter:
                tab_temp = [10**float(t) for t in h1[4::]]
            for _i in range(8):
                foo = f.readline()
            for _i_tr in range(ntr):
                tr = f.readline().split()
                if first_iter:
                    tab_wl.append(tr[5].decode("utf-8"))
                    tab_levs.append([int(tr[0]), 
                                     int(tr[1]), 
                                     tr[2].decode("utf-8"), 
                                     int(tr[3]),
                                     tr[4].decode("utf-8"), 
                                     tr[5].decode("utf-8"), 
                                     tr[6].decode("utf-8")])
                tab_tr.extend([float(t) for t in tr[7::]])
            first_iter = False

    tab_tr = np.array(tab_tr).reshape(nd, ntr, nt)

    T = Table()
    T = Table(meta={'SOURCE': 'Del Zanna & Storey, 2022, MNRAS, 513, 1198'})
    T['TEMP'], T['DENS'] = np.meshgrid(tab_temp, tab_log_dens)

    for i_wl, wl in enumerate(tab_wl):
        T[wl] = tab_tr[:,i_wl,:]

    tab_levs = pd.DataFrame(tab_levs, 
                            columns=['S', 'n_low', 'L_low', 'n_up', 'L_up', 'wl', 'wl_air'])
    return T, tab_levs
#%% get the data   
T, tab_levs = get_data()
#%% Print out the table
#write_table_hdf5(T, '../he_i_rec_DZS22.hdf5', 'updated_data', overwrite=True)
#tab_levs.to_csv('../he_i_rec_DZS22_levels.csv', index=False)
# %%
