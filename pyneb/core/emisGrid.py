import os
import numpy as np

from pyneb import config, log_, atomicData
if config.INSTALLED['plt']:
    import matplotlib.pyplot as plt
from ..utils.misc import int_to_roman, parseAtom2
from ..core.pynebcore import Atom
from ..utils.saverestore import save, restore

class EmisGrid(object):
    """
    Instantiates an atom and compute the emissivities of all the atom lines for the (tem, den) 
    values of a regularly spaced grid.
    Each line is represented in a 2D array, and there are as many arrays as transitions in the atom.
    The results can be saved for a later use in a cPickle file.

    """

    def __init__(self, elem=None, spec=None, n_tem=100, n_den=100, tem_min=5000., tem_max=20000.,
                 den_min=10., den_max=1.e8, restore_file=None, atomObj=None, NLevels=None, **kwargs):

        """
        Container for the emission line intensity grids depending on Te and Ne.

        Usage:
            O3map = pn.EmisGrid('O', 3, n_tem=30, n_den=30)
            O3map = pn.EmisGrid(restore_file='plot/O3_10k.pypic') # Recovers the map from a previous computation

        Parameters:
            - elem               element (e.g. 'O')
            - spec               spectrum (e.g. 3)
            - n_tem              number of points of the electron temperature sample
            - n_den              number of points of the electron density sample
            - tem_min, tem_max   limits for the electron temperature
            - den_min, den_max   limits for the electron density
            - restore_file       if set, the emissivities are loaded from this file, 
                                 otherwise (None) it is computed
                                 If the tem or den table from the restored file does not match the given parameters,
                                 a new grid is computed and self.compute_new_grid is set to True. The same applies
                                 if the atomic data are not the same.
            - atomObj            an atom object. In this case, no need for elem and spec.

        """
        self.log_ = log_
        self.calling = 'EmisGrid'

        if restore_file is not None:
            old = restore(restore_file)
            self.elem = old['elem']
            self.spec = old['spec']
            if atomObj is None:
                self.atom = Atom(elem=self.elem, spec=self.spec, NLevels=NLevels, **kwargs)
            else:
                self.atom = atomObj
            if self.atom.atomFitsFile != old['atomFitsFile']:
                self.log_.error('You are using {0}, but restoring a file made with {1}'.format(self.atom.atomFitsFile, old['atomFitsFile']),
                                calling=self.calling)
            if self.atom.collFitsFile != old['collFitsFile']:
                self.log_.error('You are using {0}, but restoring a file made with {1}'.format(self.atom.collFitsFile, old['collFitsFile']),
                                calling=self.calling)
            
            self.compute_new_grid = False
            if n_tem != len(old['tem']):
                self.log_.warn('len(tem) does not match saved data. New grid is computed')
                self.compute_new_grid = True
            else:
                if not np.isclose(tem_min, min(old['tem'])):
                    self.log_.warn('Min(tem) does not match saved data. New grid is computed')
                    self.compute_new_grid = True
                if not np.isclose(tem_max, max(old['tem'])):
                    self.log_.warn('Max(tem) does not match saved data. New grid is computed')
                    self.compute_new_grid = True
            if n_den != len(old['den']):
                self.log_.warn('len(den) does not match saved data. New grid is computed')
                self.compute_new_grid = True
            else:
                if not np.isclose(den_min, min(old['den'])):
                    self.log_.warn('Min(den) does not match saved data. New grid is computed')
                    self.compute_new_grid = True
                if not np.isclose(den_max, max(old['den'])):
                    self.log_.warn('Max(den) does not match saved data. New grid is computed')
                    self.compute_new_grid = True
            if not self.compute_new_grid:
                self.tem = old['tem']
                self.den = old['den']
                self.n_tem = self.tem.size
                self.n_den = self.den.size
                self.tem_min = min(self.tem)
                self.tem_max = max(self.tem)
                self.den_min = min(self.den)
                self.den_max = max(self.den)
                self.emis_grid = old['emis_grid']
        else:
            self.compute_new_grid = True
            
        if self.compute_new_grid:
            if atomObj is None:
                self.elem = elem
                self.spec = spec
                self.atom = Atom(elem, spec, NLevels=NLevels, **kwargs)
            else:
                self.atom = atomObj
                self.elem = self.atom.elem
                self.spec = self.atom.spec
            self.n_tem = n_tem
            self.n_den = n_den
            self.tem_min = tem_min
            self.tem_max = tem_max
            self.den_min = den_min
            self.den_max = den_max
            self.tem = 10 ** np.linspace(np.log10(tem_min), np.log10(tem_max), n_tem)
            self.den = 10 ** np.linspace(np.log10(den_min), np.log10(den_max), n_den)
            try:
                self.atomFitsFile = self.atom.atomFitsFile
                self.collFitsFile = self.atom.collFitsFile
            except:
                self.atomFitsFile = None
                self.collFitsFile = None
            try:
                self.recFitsFile = self.atom.recFitsFile
            except:
                self.recFitsFile = None
            self.emis_grid = self.atom.getEmissivity(self.tem, self.den)
            if restore_file is not None:
                self.save(restore_file)
                self.log_.message('%s saved to %s' % ((self.atom), restore_file), calling=self.calling)
            #self.emis_grid = np.empty((self.atom.NLevels, self.atom.NLevels, n_tem, n_den))
        
        self.den2D, self.tem2D = np.meshgrid(self.den, self.tem)
        
#    def computeEmis(self):
#        """
#        This is where the emissivities are computed, by calling a function from Atom.
#        """
#        self.emis_grid = self.atom.getEmissivity(self.tem, self.den)

    def save(self, file_):
        """
        Save results for a later use.
        
        Usage:
            O3map.save('plot/O3_30by30.pypic')
        
        Parameter:
            file_  the name of the file in which the emissivity grids are stored

        """
        save(file_, emis_grid=self.emis_grid, tem=self.tem, den=self.den, elem=self.elem,
             spec=self.spec, atomFitsFile=self.atomFitsFile, collFitsFile=self.collFitsFile)


    def getGrid(self, lev_i=None, lev_j=None, wave= -1, to_eval=None, label=None):
        """
        2D array of a line emissivity for the (Te, Ne) values of a regularly spaced grid.
        The line is specified either as the wavelength or the levels. An expression can also be used, 
        for example to_eval = 'L(5007)/L(4959)'
        
        Parameters:
            - lev_i, lev_j   levels for the emission line
            - wave           wavelength
            - to_eval        algebraic expression of the line combination to evaluate
        
        """
        if wave != -1:
            lev_i, lev_j = self.atom.getTransition(wave)
        if label is not None:
            to_eval = 'S("{}")'.format(label)

        if to_eval is None:
            to_eval = 'I(' + str(lev_i) + ',' + str(lev_j) + ')'
            
        def I(lev_i, lev_j):
            return self.emis_grid[lev_i - 1, lev_j - 1, :, :]
        def L(wave):
            return self.getGrid(wave=wave)
        def S(label):
            return self.emis_grid[label]
        return eval(to_eval)
    
    
    def plotImage(self, to_eval=None, lev_i1=None, lev_j1=None, lev_i2=None, lev_j2=None,
                  wave1= -1, wave2= -1, cblabel='', **kwargs):
        """
        Draw a contour plot of a line ratio. The contour is obtained by a horizontal cut in a temden-emissivity array.
        More elaborated plots can be obtained with the Diagnostic object.
        
        Usage:
            O3map.plotImage('L(4363)/L(5007)')
            O3map.plotImage(wave1=4363, wave2=5007)
            
        Parameters:
            - to_eval                         algebraic expression of the line combination to draw. May combine L(lambda) and I(i,j).
            - lev_i1, lev_i2, lev_j1, lev_j2  levels of the a lines in case of simple line ratio plot.
            - wave1, wave2                    wavelengths of the two lines in case of simple line ratio plot.
            - cblabel                         a title for the colorbar
            - **kwargs                        any other parameter will be sent to contourf and pcolor

        """
        if not config.INSTALLED['plt']:
            log_.error('Matplotlib not available', calling=self.calling)
            return None
        L = lambda lam: self.getGrid(wave=lam)
        I = lambda i, j: self.getGrid(i, j)
        S = lambda label: self.getGrid(label)
        if to_eval is None:
            if wave1 != -1:
                lev_i1, lev_j1 = self.atom.getTransition(wave1)
            if wave2 != -1:
                lev_i2, lev_j2 = self.atom.getTransition(wave2)
            to_eval = 'I(' + str(lev_i1) + ',' + str(lev_j1) + ') / I(' + str(lev_i2) + ',' + str(lev_j2) + ')'
        X = np.log10(self.den2D)
        Y = self.tem2D / 1e4
        Z = self.getGrid(to_eval=to_eval)
        plt.contourf(X, Y, Z, **kwargs)
        emap = plt.pcolor(X, Y, Z, **kwargs)
        cbar = plt.colorbar(emap)
        cbar.set_label(cblabel)
        # plt.clabel(contour, inline=1)
        plt.xlabel(r'log Electron density (cm$^{-3}$)')
        plt.ylabel(r'Electron temperature ($10^4$K)')
        plt.show()


    def plotContours(self, to_eval, low_level=None, high_level=None, n_levels=20,
                      linestyles='-', clabels=True, log_levels=True, title=None, ax=None, **kwargs):
        """
        Plot a contour map of the diagnostic defined by to_eval.
        
        Usage:
            O3map.plotContours('L(4363)/L(5007)')
        
        Parameters:
            - to_eval                   algebraic definition of the line ratio to contour-plot
            - low_levels, high_levels   limit levels for the contour. If not set, the limit of the 
                                        ratio values are used.
            - n_levels                  number of levels (20 is the default)
            - linestyles                used for the contour lines (default: solid line)
            - clabels                   Boolean. Controls if line labels are printed (default: True)
            - log_levels                Boolean. If True (default), the log of the line ratio is used.
            - title                     plot title
            - **kwargs                  sent to plt.contour
            
        """ 
        if not config.INSTALLED['plt']:
            log_.error('Matplotlib not available', calling=self.calling)
            return None
        X = np.log10(self.den2D)
        Y = self.tem2D
        L = lambda lam: self.getGrid(wave=lam)
        I = lambda i, j: self.getGrid(i, j)
        S = lambda label: self.getGrid(label=label)
        try:
            diag_map = eval(to_eval)
        except:
            self.log_.warn('diag %s not found' % to_eval, calling=self.calling)
            return None
        if low_level is None:
            low_level = np.min(np.log10(diag_map))
        if high_level is None:
            high_level = np.max(np.log10(diag_map))
        if log_levels:
            Z = np.log10(diag_map)
            levels = np.linspace(low_level, high_level, n_levels)
            if title is None:
                title = '[%s%s]: log(%s)' % (self.elem, int_to_roman(int(self.spec)), to_eval)
        else:
            Z = diag_map
            levels = 10. ** (np.linspace(low_level, high_level, n_levels))
            if title is None:
                title = '[%s%s]: %s' % (self.elem, int_to_roman(int(self.spec)), to_eval)
        
        if ax is None:
            f, ax = plt.subplots()
        else:
            f = plt.gcf()
        CS = ax.contour(X, Y, Z, levels=levels, linestyles=linestyles, **kwargs)
        if clabels:
            ax.clabel(CS, CS.levels[::3], inline=True, fontsize=12, colors='black')
        ax.set_xlabel(r'log(n$_{\rm e}$) [cm$^{-3}$]')
        ax.set_ylabel(r'T$_{\rm e}$ [K]')
        ax.set_title(title)


    def plotLineRatio(self, to_eval, par=None, par_low=None, par_high=None, n_par=10,
                      linestyles='-', title=None, legend=True, loc=2, **kwargs):
        """
        Plot the diagnostic defined by to_eval as a function of either Ne or Te, with the other
        variable as a parameter.
        
        Usage:
            o3grid.plotLineRatio('L(4363)/L(5007)', par='den')
            s2grid.plotLineRatio('I(2,1)/I(3,1)', par='tem', par_low=5000, par_high=20000, n_par=4)
        
        Parameters:
            - to_eval      algebraic definition of the line ratio to contour-plot
            - par          quantity to be used as a parameter (either 'den' or 'tem')
            - par_low      lowest limit of the parameter range
            - par_high     highest limit of the parameter range
            - n_par        number of parameter's values (default: 10)
            - linestyles   used for the contour lines (default: solid)
            - title        plot title
            - legend       Boolean. If True, write legend title (default: True) 
            - **kwargs     sent to plt.contour
            
        """ 
        if not config.INSTALLED['plt']:
            log_.error('Matplotlib not available', calling=self.calling)
            return None

        L = lambda wav: self.getGrid(wave=wav)
        I = lambda i, j: self.getGrid(i, j)
        
        if par == 'tem':
            X = np.log10(self.den2D)
            Z = self.tem2D
            leg_format = "{0:.0f} K"
            leg_title = r'Temperature'
            xlabel = r'Log(N$_{\rm e}$) [cm$^{-3}$]'
        elif par == 'den':
            X = self.tem2D
            Z = np.log10(self.den2D)
            leg_format = "{0:.2f}"
            leg_title = r'Log(N$_{\rm e}$)'
            xlabel = r'T$_{\rm e}$ [K]'
        else:
            self.log_.error('The parameter must be either den or tem (%s given)' % par, calling=self.calling)

        try:
            Y = eval(to_eval)
        except:
            self.log_.warn('Diagnostic %s not found' % to_eval, calling=self.calling)            
            return None

        if par_low is None:
            par_low = np.min(Z)
        if par_high is None:
            par_high = np.max(Z)
        levels = np.linspace(par_low, par_high, n_par)

        CS = plt.contour(X, Y, Z, levels=levels, linestyles=linestyles, **kwargs)
        
        if legend:
            for i_level in range(len(levels)):
                CS.collections[i_level].set_label(leg_format.format(levels[i_level]))
                plt.legend(title=leg_title, loc=loc)

        if title is None:
            title = '[%s%s]' % (self.elem, int_to_roman(int(self.spec)))

        plt.xlabel(xlabel)
        plt.ylabel(to_eval)
        plt.title(title)
    
        plt.show()
        

def getEmisGridDict(elem_list=None, spec_list=None, atom_list=None, restore=True, pypic_path=None,
            n_tem=100, n_den=100, tem_min=5000., tem_max=20000., den_min=10., den_max=1e8, save=True,
            atomDict=None, createPypicsDir=True, computeIfAbsent=True):
    """
    Return a dictionary of EmisGrid objects referred to by their atom string (e.g. 'O3').
    
    Usage:
        emisgrids = pn.getEmisGridDict(['C', 'N', 'O'], [2, 3], save=True)
        emisgrids = pn.getEmisGridDict(atom_list=['N2', 'O2', 'O3'], pypic_path='/tmp/pypics/',
                        den_max = 1e6)
        emisgrids = pn.getEmisGridDict(atomDict=diags.atomDict)

    Parameters:
        - elem_list         list of elements
        - spec_list         list of spectrum values (integers)
        - atom_list         list of atoms (e.g. ['N2', 'O3'])
        - restore           Boolean. If True (default), the program will search for a previously computed grid.
                            If not found, it will compute the EmisGrid (unless computeIfAbsent is set to False)
        - pypic_path        directory where to look for the pypic emissivity files. It defaults to pn.config.pypic_path, 
                            which is defined when the PyNeb session begins and corresponds to $HOME/.pypics if
                            the ENVIRONMENT variable HOME is accessible or to /tmp/pypics otherwise
        - n_tem, n_den      number of points in the Te and Ne samples resp. 
        - tem_min, tem_max, den_min, den_max    limits for Te and Ne
        - save              Boolean. If True (default), the emissivity grid will be saved. It has no effect if the 
                            EmisGrid has just been restored.
        - atomDict          dictionary of Atom objects
        - computeIfAbsent   If the file to restore is absent and this parameter is True (default), the 
                            EmisGrid is computed
            
    """
    if pypic_path is None:
        pypic_path = config.pypic_path
    else:
        config.pypic_path = pypic_path
    calling = 'getEmisGridDict'
    if pypic_path is None:
        return None
    emis_grids = {}
    atoms = []
    if atomDict is not None:
        atoms = atomDict.keys()
    else:
        if atom_list is not None:
            atoms = atom_list
        else:
# VL Removed 13 March 2015 because it made the code ignore the elem keyword             
#            atoms = atomicData.getAllAtoms()
            if (elem_list is None) and (spec_list is None):
                atoms = atomicData.getAllAtoms()
            else:
                for elem in elem_list:
                    for spec in spec_list:
                        atoms.append(elem + str(spec)) 
    for atom in atoms:
        elem, spec, rec = parseAtom2(atom)
        if atomDict is not None:
            atomObj = atomDict[atom]
        else:
            atomObj = None
        file_ = '{0}/emis_{1}{2}{3}.pypic'.format(pypic_path, elem, spec, rec)
        if restore:
            if os.path.exists(file_):
                try:
                    emis_grids[elem + spec] = EmisGrid(restore_file=file_, n_tem=n_tem, n_den=n_den, 
                                                       tem_min=tem_min, tem_max=tem_max,
                                                       den_min=den_min, den_max=den_max, atomObj=atomObj)
                    log_.message('Read %s' % file_, calling=calling)
                    compute_this = False
                except:
                    if computeIfAbsent:
                        log_.warn('Wrong emission map: {0}, creating it'.format(file_), calling=calling)
                        compute_this = True
                    else:
                        log_.error('Wrong emission map: {0}'.format(file_), calling=calling)
                        compute_this = False
            else:
                log_.warn('Emission map not found: {0}'.format(file_), calling=calling)
                if computeIfAbsent:
                    compute_this = True
        else:
            compute_this = True
        if compute_this:
            try:
                emis_grids[atom] = EmisGrid(elem, spec, n_tem=n_tem, n_den=n_den, tem_min=tem_min, tem_max=tem_max,
                                              den_min=den_min, den_max=den_max, atomObj=atomObj)
                if save:
                    emis_grids[atom].save(file_)
                    log_.message('%s saved to %s' % ((atom), file_), calling=calling)
            except:
                log_.warn('No %s EmisGrid' % (atom), calling=calling)
    return emis_grids
