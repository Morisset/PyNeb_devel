import numpy as np
import matplotlib.pyplot as plt
import pyneb as pn
import pyCloudy as pc
from pyCloudy.utils.misc import make_mask, convert_label
from pyCloudy.utils.astro import conv_arc

class In(object):
    
    def __init__(self, name):
        """
        Defining the parameters of the model
        """
        
        self.model_name = name
        self.r_in = 16.09
        # Here we use a user-defined double gaussisn
        self.dlaw_params = [3, 2850. , 5256.6 , 16.28 , 15.96 , 9500.0 , 17.02 , 16.46]
        self.ff = 1.0
        self.distance = 1.25
        self.abunds = {'He':-0.92, 'C':-3.10, 'N':-4.00, 'O':-3.4, 'Ne':-4.0, 'Mg':-4.95,
                       'Si':-4.90, 'S':-5.35, 'Cl':-7.00, 'Ar':-5.8, 'Fe':-7.40}
        self.grains_type = ['orion graphite', 'orion silicate']
        self.grains = [0.2, 0.2]
        self.Q0 = 47.5
        self.SED = 'BB' # STAR or BB
        self.Teff = 39500.
        self.options = ('no molecules',
                        'no level2 lines',
                        'no fine opacities',
                        'atom h-like levels 10',
                        'atom he-like levels 10',
                        'COSMIC RAY BACKGROUND',
                        'element limit off -8',
                        )

        
    def print_model(self):
        """
        Preparing and printing the Cloudy input file
        """
        model = pc.CloudyInput(self.model_name)
        model.set_radius(self.r_in)
        model.set_dlaw(self.dlaw_params, self.ff)
        model.set_distance(self.distance)
        model.set_abund(ab_dict = self.abunds)
        for grains, grains_type in zip(self.grains, self.grains_type):
            model.set_grains('{0} {1}'.format(grains_type, grains))
        if self.SED == 'STAR':
            model.set_star('table star "mod103.mod"', 39390, 'q(H)', self.Q0)
        elif self.SED == 'BB':
            model.set_BB(Teff=self.Teff, lumi_unit='q(H)', lumi_value=self.Q0)
        else:
            pc.log_.error('unknown SED value: {0}'.format(self.SED))
        model.read_emis_file('ic418N.lines')
        model.set_iterate(0)
        model.set_sphere()
        model.set_other(self.options)
        model.print_input(to_file = True, verbose = False)
        self.model = model
        
    def run_model(self):
        self.model.run_cloudy()
        
        
class Out(object):
    
    def __init__(self, name):
        """
        Creating an object to manage the Cloudy output and compare to the observations
        """

        self.model = pc.CloudyModel(name)
        self.set_3D(use=False)
        self.m3d = None
        self.mask = None
        self.profile_defined = False
        self.dim_3D = 101
        self.mask_ap_center = (0., 3.5)
        self.mask_ap_size = (12., 1)
        # define a RedCorr object
        self.RC = pn.RedCorr(E_BV=0.26, R_V=3.6, law='Fitz 99')

    def read_obs(self):
        self.obs = np.genfromtxt('ic418N.lines', dtype=["a11","float","int"],
                            delimiter=[11,12,1], names=['label', 'i_obs', 'obs'])
        # Convert the Cloudy labels into pyCloudy labels
        self.obs['label'] = [convert_label(l) for l in self.obs['label']]
        self.IR_cross_calib_done = False
        self._cross_calibIR()
        self.UV_cross_calib_done = False
        self._cross_calibUV()

    def print_obs(self):
        for obs in self.obs:
            print('{0} : {1}'.format(obs['label'], obs['i_obs']))

    def _cross_calibIR(self):
        """
        Cross calibration of the IR line intensities, using the theoretical value of 4.05mu/100.Hbeta, as given by the model
        This method is not supposed to be used by the user (it begins with _).
        """
        if self.IR_cross_calib_done:
            pn.log_.warn('Are you trying to cross_calib data that are already calibrated???')
            return
        theo_ratio = self.get_emis('H__1_4051M') / self.get_emis('H__1__4861A') * 100
        i_4051 = (self.obs['label'] == 'H__1_4051M')
        norm = theo_ratio / self.obs['i_obs'][i_4051]
        i_IR = (self.obs['obs'] == 1)
        self.obs['i_obs'][i_IR] *= norm
        self.IR_cross_calib_done = True
        pc.log_.message('Correcting {0} IR lines so that 4.05mu/100.Hbeta = {1}, using norm={2}'.format(i_IR.sum(), theo_ratio, norm))

    def _cross_calibUV(self):
        """
        Cross calibration of the UV line intensities, using the theoretical value of [OII] 2471/7323, as given by the model.
        This method is not supposed to be used by the user (it begins with _).
        """        
        if self.UV_cross_calib_done:
            pn.log_.warn('Are you trying to cross_calib data that are already calibrated???')
            return
        theo_ratio = self.get_emis('O_II__2471A') /self.get_emis('O_II__7323A')
        i_2471 = (self.obs['label'] == 'O_II__2471A')
        i_7323 = (self.obs['label'] == 'O_II__7323A')
        norm = theo_ratio * self.obs['i_obs'][i_7323] / self.obs['i_obs'][i_2471]
        i_UV = (self.obs['obs'] == 2)
        self.obs['i_obs'][i_UV] *= norm
        self.UV_cross_calib_done = True
        pc.log_.message('Correcting {0} UV lines so that [OII] 2471/7323 = {1}, using norm={2}'.format(i_UV.sum(), theo_ratio, norm))
        
    def set_3D(self, use=False):
        if use:
            # already in 3D mode, nothing to do!
            if self.use_3D:
                return
            # create the 3D model
            # Very important to have center = True for the line profiles.
            self.m3d = pc.C3D(self.model, dims = self.dim_3D, center = True, n_dim = 1)
            self.use_3D = True
            # need to read the observations as the IR and UV cross-calibration depends on the apertures
            self.read_obs()
        else:
            self.use_3D = False
            self.read_obs()
            
    def get_emis(self, label):
        """
        Return the line emission. 
        If use_3D is set, the 3D model is used, otherwise the 1D
        """
        # test if we are in the 3D case or not
        if self.use_3D:
            # look for the observation type, according to the "label" parameter
            # in the table obs, look for the column 'obs'
            obs = self.obs['obs'][(self.obs['label'] == label)][0]
            if obs == 0: # this is the observation type for which a mask is required
                # the mask is already defined?
                if self.mask is None:
                    # if not, do it!
                    self.mask = self.get_mask()
                mask = self.mask
            else:
                # if obs != 1, the observation is of the full object, no need for a mask
                mask = 1.  
            return ((self.m3d.get_emis(label)).sum(axis=1) * mask).sum()
        else:
            return self.model.get_emis_vol(label)
        
    def get_mask(self, seeing = None):
        """
        Return a mask to apply to a 2D projection of a 3D models, to reproduce an aperture
        """
        self.set_3D(use=True)
        # convert the sizes in arcsec
        x_arc = pc.astro.conv_arc(dist_proj = self.m3d.cub_coord.x_vec, dist = self.model.distance)
        y_arc = pc.astro.conv_arc(dist_proj = self.m3d.cub_coord.y_vec, dist = self.model.distance)
        # produce two 2D tables for X and Y
        X, Y = np.meshgrid(x_arc, y_arc)
    
        mask = make_mask(X, Y, ap_center = self.mask_ap_center, ap_size = self.mask_ap_size, seeing = seeing)
        return mask

    def print_res(self):
        """
        Print the results
        """
        # Rout of the photoionization model
        r_out_arcsec = conv_arc(dist = self.model.distance, dist_proj = np.max(self.model.r_out))
        print('Outer radius in arcsec: {0:.2f}'.format(r_out_arcsec))
        # Hbeta total absolute intensity at earth (distance dependant)
        Hbeta_tot = self.model.get_emis_vol('H__1__4861A', at_earth=True)
        # print the reddened Hbeta intensity, to be compared to the observed value of -9.57
        print('Hbeta Abs   :  -9.57   {0:5.2f}'.format(np.log10(Hbeta_tot / self.RC.getCorr(4861.))))
        
        # the following prints line intensities, through an aperture or not, depending on the
        # value of the use_3D parameter.
        Hbeta = self.get_emis('H__1__4861A')
        for obs in self.obs:
            model = self.get_emis(obs['label']) / Hbeta * 100
            print('{0} : {1:6.2f} {2:6.2f} {3:6.2f}'.format(obs['label'], obs['i_obs'], model, 
                                                            (model-obs['i_obs'])/obs['i_obs']))
    
    def plot_im(self, label = 'H__1__4861A', axis=1, plot_mask=True):
        """
        Plot an image projected on an axis
        Parameters:
            - label: which emission line image to plot
            - axis: the projection axis
            - plot_mask: if True (default). overplot the mask
        """
        
        self.set_3D(use=True)
        # the size of the image in arcsec is determined from the distance of the object and 
        # the maximum of the x_vec variable
        max_x = conv_arc(dist = self.model.distance, dist_proj = np.max(self.m3d.cub_coord.x_vec))
        plt.imshow((self.m3d.get_emis(label)).sum(axis=axis), extent=[-max_x, max_x, -max_x, max_x])
        if plot_mask:
            if self.mask is None:
                self.mask = self.get_mask()
            plt.contour(self.mask,levels = [0.5], colors = 'black', linestyles='-', 
                        extent=[-max_x, max_x, -max_x, max_x])
            plt.contour(self.mask,levels = [0.1, 0.9], colors = 'black', linestyles='--', 
                        extent=[-max_x, max_x, -max_x, max_x])
        
    def plot_3col(self, list_emis = ['N__2__6548A', 'O__3__5007A', 'H__1__4861A']):
        """
        plot a 3-colors image, using the list_emis references for RGB channels.
        """
        self.set_3D(use=True)
        # the size of the image in arcsec is determined from the distance of the object and 
        # the maximum of the x_vec variable
        max_x = conv_arc(dist = self.model.distance, dist_proj = np.max(self.m3d.cub_coord.x_vec))
        plt.imshow(self.m3d.get_RGB(list_emis), extent=[-max_x, max_x, -max_x, max_x])
        
    def define_profiles(self, size_spectrum = 21, vel_max = 50., v_turb = 5., vel_params = [0., 10, 10]):
        """
        Define the line profiles by setting the velocity law and configuring the profiles.
        """
        self.set_3D(use=True)
        if self.mask is None:
            self.mask = self.get_mask()
        self.m3d.set_velocity(velocity_law = 'poly', params = vel_params)
        self.m3d.config_profile(size_spectrum = size_spectrum, vel_max = vel_max, 
                                v_turb = v_turb, profile_function = 'gaussian')
        self.profile_defined = True
        
    def plot_profile(self, x_pos = None, y_pos=None, legend=True, **kwargs):
        """
        plot the line profiles of Hbeta, [NII] and [OIII] at a given position (x_pos, y_pos).
        Arguments:
            - x_pos, y_pos: position at whci the line are axtracted
            - legend: If True (default), plot the legend.
            - Any other argument is passed to the 3 plots (e.g. linestyle)
        """
        if not self.profile_defined:
            pc.log_.error('Profile not defined!')
            return
        self.set_3D(use=True)
        if x_pos is None:
            x_pos = self.dim_3D/2
        if y_pos is None:
            y_pos = self.dim_3D/2
        plt.plot(self.m3d.vel_tab,self.m3d.get_profile('H__1__4861A', axis='x')[:,x_pos,y_pos] * 5, 
                 label = r'H$\beta$', **kwargs)
        plt.plot(self.m3d.vel_tab,self.m3d.get_profile('N__2__6584A', axis='x')[:,x_pos,y_pos] * 5, 
                 label = r'[NII]$\lambda$6584', **kwargs)
        plt.plot(self.m3d.vel_tab,self.m3d.get_profile('O__3__5007A', axis='x')[:,x_pos,y_pos], 
                 label = r'[OIII]$\lambda$5007', **kwargs)
        if legend:
            plt.legend()
        
    def plot_cont(self):
        
        plt.subplot(2,1,1)
        plt.plot(self.model.get_cont_x('eV'), self.model.get_cont_y(cont='incid', unit='esA'))
        plt.xlim(10, 54)
        plt.xlabel('Energy (eV)')
        plt.ylabel('erg.s-1.A-1')
        
        plt.subplot(2,1,2)
        plt.loglog(self.model.get_cont_x('Ang'), self.model.get_cont_y(cont='diffout', unit='Jy'))
        plt.xlim(1e3, 1e6)
        plt.xlabel('Wavelength (A)')
        plt.ylabel('Jy')
        