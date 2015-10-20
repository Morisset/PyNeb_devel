import numpy as np
import matplotlib.pyplot as plt
import pyCloudy as pc
import pyneb as pn
from pyCloudy.utils.misc import convert_label, sextract

emis_file = 'NGC5128_56_12b_emis.dat'

class In(object):

    def __init__(self, model_dir, name, r_in, dens, Teff, Q0, abund_dict, distance, grains=None):
        """
        Defining the parameters of the model
        """
        # combining dir and name
        self.model_name = '{0}/{1}'.format(model_dir, name)
        # set the input parameters to self. variables
        self.r_in = r_in
        self.dens = dens
        self.distance = distance
        self.abund_dict = abund_dict
        self.Q0 = Q0
        self.Teff = Teff
        self.grains = grains
        # define more options to Cloudy
        self.options = ('no molecules',
                        'no level2 lines',
                        'no fine opacities',
                        'atom h-like levels small',
                        'atom he-like levels small',
                        'COSMIC RAY BACKGROUND',
                        'element limit off -8',
                        )

        
    def print_model(self):
        """
        Preparing and printing the Cloudy input file
        """
        # define the name of the model
        model = pc.CloudyInput(self.model_name)
        # send the variables to the CloudyInput object to be printed 
        model.set_radius(self.r_in)
        model.set_cste_density(self.dens)
        model.set_distance(self.distance, unit='Mpc')
        model.set_abund(ab_dict = self.abund_dict)
        model.set_grains(self.grains)
        model.set_BB(Teff=self.Teff, lumi_unit='q(H)', lumi_value=self.Q0)
        # this is the file containing the list of emissivities we want
        model.read_emis_file(emis_file)
        model.set_iterate(0)
        model.set_sphere()
        model.set_other(self.options)
        # print the input file
        model.print_input(to_file = True, verbose = False)
        # store the model in a self variable to further interactions if needed
        self.model = model
        
    def run_model(self):
        # call the Cloudy runner
        self.model.run_cloudy()
        
class Outs(object):
    
    def __init__(self, model_dir, models):
        """
        if type(models) == type(''):
            models = [models]
        full_names = ['{0}/{1}'.format(model_dir, model) for model in models] 
        """
        self.Ms = pc.load_models('{0}/{1}'.format(model_dir, models))
        self.read_obs()
        
        
    def read_obs(self):
        self.obs_txt = np.genfromtxt(emis_file, dtype=["a11","float", "float"], 
                            delimiter=[11,8, 6], names = ['label', 'i_obs', 'e_obs'], usecols = (0, 1, 2))

        # redenning correction
        Hb = self.obs_txt['i_obs'][self.obs_txt['label'] == 'H  1  4861A']
        Ha = self.obs_txt['i_obs'][self.obs_txt['label'] == 'H  1  6563A']
        self.RC = pn.RedCorr(law = 'Fitz 99')
        self.RC.setCorr(Ha / Hb / 2.85, 6563, 4861)
        for line in self.obs_txt:
            lambda_ = np.float(line['label'][-5:-1])
            line['i_obs'] *= self.RC.getCorrHb(lambda_)
        
    def get_i_obs(self, label):
        i_label = (self.obs_txt['label'] == label)
        return self.obs_txt[i_label]['i_obs'][0]
        
    def pretty_print(self, str1, list1):
        if type(list1[0]) == type(''):
            print('{0:32s}'.format(str1) + ' '.join(['{0:>9}'.format(i) for i in list1]))
        else:
            print('{0:32s}'.format(str1) + ' '.join(['{0:>9.3f}'.format(i) for i in list1]))

    def print_ratio(self, label1, label2, title):
        ref_pycloudy1 = convert_label(label1)
        ref_pycloudy2 = convert_label(label2)
        obs_ratio = self.get_i_obs(label1) / self.get_i_obs(label2)
        mod_ratio = [M.get_emis_vol(ref_pycloudy1) / M.get_emis_vol(ref_pycloudy2) for M in self.Ms]
        str1 = '{0:12s} {1:>8.3f}'.format(title, obs_ratio)
        self.pretty_print(str1, mod_ratio)

    def print_res(self):
        
        names = [M.model_name_s for M in self.Ms]
        self.pretty_print('MODEL', names)
        r_in = [np.log10(M.r_in) for M in self.Ms]
        self.pretty_print('Inner radius', r_in)
        Teff = [np.float(sextract(M.out['Blackbody'], 'Blackbody ', ' '))/1e3 for M in self.Ms]
        self.pretty_print('Effective Temp kK', Teff)
        dens = [np.log10(M.nH[0]) for M in self.Ms]
        self.pretty_print('Hydrogen density', dens)
        Q0 = [np.log10(M.Q0) for M in self.Ms]
        self.pretty_print('Q0', Q0)
        logUmean = [M.log_U_mean for M in self.Ms]
        self.pretty_print('<logU>', logUmean)
        abHe = [M.abund['He'] for M in self.Ms]
        self.pretty_print('He/H', abHe)
        abC = [M.abund['C'] for M in self.Ms]
        self.pretty_print('C/H', abC)
        abN = [M.abund['N'] for M in self.Ms]
        self.pretty_print('N/H', abN)
        abO = [M.abund['O'] for M in self.Ms]
        self.pretty_print('O/H', abO)
        abNe = [M.abund['Ne'] for M in self.Ms]
        self.pretty_print('Ne/H', abNe)
        abS = [M.abund['S'] for M in self.Ms]
        self.pretty_print('S/H', abS)
        abAr = [M.abund['Ar'] for M in self.Ms]
        self.pretty_print('Ar/H', abAr)
        Hb = [np.log10(M.get_emis_vol('H__1__4861A')/self.RC.getCorr(4861)/
                       (4.*np.pi*(M.distance*pc.CST.KPC)**2)) for M in self.Ms]
        self.pretty_print('Hbeta         -16.170', Hb)
        
        for line in self.obs_txt:
            ref_pycloudy = convert_label(line['label'])
            try:
                mod_intens = [M.get_emis_vol(ref_pycloudy) / M.get_emis_vol('H__1__4861A') * 100 for M in self.Ms]
                str1 = '{0} {1:>9.3f} +/-{2:>6.3f} '.format(line['label'], line['i_obs'], line['e_obs'])
                self.pretty_print(str1, mod_intens)
            except:
                print('Something wrong with {0}'.format(line['label']))
        try:
            self.print_ratio('S II  6731A', 'S II  6716A', 'Dens(SII)')
        except:
            pass
        try:
            self.print_ratio('TOTL  3727A', 'O  3  5007A', 'OII/III')
        except:
            pass

if __name__ == '__main__':

    pc.log_.level = 2
    models_dir = '/Users/christophemorisset/DATAS/Choroni'
    # define the model name and properties
    model_name = 'M64_A'
    i = 5

    r_in = 16. 
    dens = np.log10(4e4) 
    Teff = 45000 
    Q0 = 47.5
    distance = 1.0
    ab_dict = {'He':-0.90, 'C':-3.60, 'N':-3.45, 'O':8.55-12 , 'Ne':-4.3, 'Mg':-4.95,
               'Si':-4.90, 'S':-5.35, 'Cl':-7.00, 'Ar':-6., 'Fe':-7.40}

    pc.log_.level = 3
    # create the object that generates the input files
    #for i, Teff in enumerate([120000, 120000, 130000]):
    Min = In(models_dir, '{0}_{1}'.format(model_name, i), r_in, dens, Teff, Q0,
             ab_dict, distance)
    Min.print_model()

    # run the models
    pc.run_cloudy(dir_=models_dir, n_proc=3, use_make=True,  model_name=model_name)

    # read the models
    pc.log_.level = 2
    Mouts = Outs(models_dir, model_name)
    # output the parameters and line intensities, with the observations in 1rst column
    Mouts.print_res()

