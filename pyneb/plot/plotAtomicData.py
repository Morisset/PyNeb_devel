"""@module plotAtomicData
Module to plot the atomic data
""" 

import pyneb as pn
import numpy as np
import string
if pn.config.INSTALLED['plt']:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.ticker
    from matplotlib.ticker import MaxNLocator

from pyneb.utils.misc import int_to_roman, parseAtom

class DataPlot(object):
    """
    Plot transition probabilities and collision strengths from different data sets
    
    """
    def __init__(self, elem=None, spec=None, all_data=[], atom=None, n_tem_points=10000, 
                 ref_tem=None, OmegaInterp='Cheb',NLevels=None):
        """
    Parameters:
        - elem         atomic elem 
        - spec       ionization stage in spectroscopic notation (I = 1, II = 2, etc.)
        - atom        e.g. 'O3'
        - all_data       dictionary of all_data to be compared (see above for format)
        - [n_tem_points] number of points in the fit (default=100; increase if fit is not smooth)
        - [ref_tem]      array of temperature values to be signaled in the plots
        - OmegaInterp    interpolating function between Omega values ('Cheb' [default], 'Linear')
    
    Example:
        dataplot = pn.DataPlot('O', 3) # initializes the plot
        dataplot.plotA() # transition probabilities plot 
        dataplot.plotRelA() # relative transition probabilities plot
        dataplot.plotOmega() # collision strength plot    

        
        """
        colors = np.array(['r', 'g', 'b', 'm', 'c', 'y'])
        self.calling = 'DataPlot'
        if atom is not None:
            self.atom = str.capitalize(atom)
            self.elem = parseAtom(self.atom)[0]
            self.spec = int(parseAtom(self.atom)[1])
        else:
            self.elem = str.capitalize(elem)
            self.spec = int(spec)
            self.atom = self.elem + str(self.spec)

        old_data = pn.atomicData.getDataFile(self.atom)
        # Check if matplotlib installed
        if not pn.config.INSTALLED['plt']:
            pn.log_.warn('Matplotlib not installed!', calling=self.calling)

        # Separate omega and A data sets    
        atom_data = []
        coll_data = []
        if all_data == []:
            all_data = pn.atomicData.getAllAvailableFiles(self.atom)
        i_colors = 0    
        for file_ in all_data:
            ID = file_.split('_')[-1]
            type_ = (file_.split('_')[2]).split('.')[0]
            if type_ in ['atom', 'coll']:
                data_list = type_ + '_data'
                vars()[data_list].append({'ID': ID, 'file_': file_, 'type': type_, 'color': colors[i_colors % len(colors)]})
                i_colors += 1
                
        self.atom_data = atom_data
        self.coll_data = coll_data

        # Temperature values to be signaled by vertical lines in omega plots. Can be changed by defining ref_tem
        if ref_tem is None:
            self.ref_tem = np.log10([5000., 10000., 20000.]) 
        else:
            self.ref_tem = ref_tem

        # If it's not standard effective collision strengths, we don't want it
        coll_data_temp=[]
        for item in self.coll_data: 
            coll_data_temp.append(item)
        for data in coll_data_temp:
            pn.atomicData.setDataFile(data['file_'])
            atom = pn.Atom(self.elem, self.spec, OmegaInterp=OmegaInterp, NLevels=NLevels)
            """
            try:
                if 'O_UNIT' in atom.CollData.comments.keys():
                    self.coll_data.remove(data)
            except:
                pass
            """
        # For each data set, an atom is built. 
        del(atom)
        for data in self.atom_data + self.coll_data:
            pn.atomicData.setDataFile(data['file_'])
            atom = pn.Atom(self.elem, self.spec, OmegaInterp=OmegaInterp, NLevels=NLevels)
            data['atom'] = atom
        for data in old_data:
            if data is not None:
                pn.atomicData.setDataFile(data)
        self.atom_rom = string.lower(self.elem + '_' + int_to_roman(int(self.spec)))
        self.n_tem_points = n_tem_points

        self.atom_n_max = 0
        for at in self.atom_data:
            if at['atom'].atomFileType != 'chianti':
                if at['atom'].atomNLevels > self.atom_n_max:
                    self.atom_n_max = at['atom'].atomNLevels
        self.coll_n_max = 0
        for at in self.coll_data:
            if at['atom'].collFileType != 'chianti':
                if at['atom'].collNLevels > self.coll_n_max:
                    self.coll_n_max = at['atom'].collNLevels
                                    
    def plotA(self, save=False, figsize=(18, 12), fignum=None, NLevels=None):
        """
        Plot the log of the A values of each data set 
        
        Parameters:
            - save     if True, saves the plot in a file
            - figsize  figure size (default: [18, 12])
            - fignum    figure Number

        """
        if NLevels is None:
            atom_n_max = self.atom_n_max
        else:
            atom_n_max = NLevels
        if not pn.config.INSTALLED['plt']:
            pn.log_.error('Matplotlib not installed!', calling=self.calling)
        plt.figure(fignum, figsize=figsize)
        plt.clf()
        ax = plt.subplot(111, axisbg='Ivory')
        x = np.arange(atom_n_max * (atom_n_max - 1) / 2.)
        ticks = x
        # Inventory of markers. Chosen to be distinguished even when overlapped
        mark = ['<', (5, 1), '>', 'o', '|']
        i_marker = 0
        # Background colors to distinguish lower levels 
        bg_color_sequence = ['g', 'r', 'b', 'y', 'm', 'c']
        # The following command replicates the string to allocate a sufficient number of colors
        bg_color = bg_color_sequence * (int(atom_n_max - 2) / len(bg_color_sequence) + 1)
        tick_label = []
        
        for i in range(atom_n_max - 1):
            # x0, x1 are start and end points of each level
            x0 = i * (atom_n_max - 1) - i * (i - 1) / 2 - 0.5
            width = atom_n_max - i - 1
            x1 = x0 + width
            plt.axvspan(x0, x1, facecolor=bg_color[i], alpha=0.05)
            # The x axis must stretch the maximum range (although some data set might have a lower n_level)
            for j in range(i + 1, atom_n_max):
                tick_label.append('(' + str(j + 1) + ', ' + str(i + 1) + ')')


        for data in self.atom_data:
            n_levels = data['atom'].atomNLevels
            Ay = []
            A = data['atom'].getA()
            color = data['color']
            try:
                for i in range(n_levels - 1):
                    for j in range(i + 1, n_levels):
                        if A[j, i] > 0:
                            Ay.append(np.log10(A[j, i]))
                        else:
                            Ay.append(np.NaN)
                plt.scatter(x, Ay, marker=mark[i_marker], s=300., c=color, alpha=0.35, linewidths=1, label='%s' % (data['ID']))
                ax.set_xticks(ticks)
            except:
                pn.log_.warn('Problem in plotting A', calling=self.calling + '.plotA')
            i_marker += 1
        ax.set_xticklabels(tick_label)

        # Plot features
        plt.xlabel('Transition')
        plt.ylabel('Log A(j, i)')
        plt.title('Transition probabilities for [%s %s]' % (self.elem, int_to_roman(int(self.spec))))
        plt.legend(loc='lower right', markerscale=1., scatterpoints=1, borderpad=1, labelspacing=1)
        plt.show()      
        if save:
            plt.figure(figsize=[18, 12])
            plt.savefig(self.atom_rom + '_' + data['ID'] + "-" + self.ref_data + "_A.pdf")



    def plotAllA(self, save=False, figsize=(18, 12), fignum=None, NLevels=None):
        """
        Plot the log of the A values of each data set 
        
        Parameters:
            - save     if True, saves the plot in a file
            - figsize  figure size (default: [18, 12])
            - fignum    figure Number

        """
        if NLevels is None:
            atom_n_max = self.atom_n_max
        else:
            atom_n_max = NLevels
        if not pn.config.INSTALLED['plt']:
            pn.log_.error('Matplotlib not installed!', calling=self.calling)
        fig = plt.figure(fignum, figsize=figsize)
        plt.clf()
        max_y_ticks = 6
        ticks = np.arange(len(self.atom_data)+1)[1:] 
        
        A = np.zeros([len(self.atom_data), atom_n_max, atom_n_max])
        color=[]
        for i, data in enumerate(self.atom_data):
            A_tmp = data['atom'].getA()
            if A_tmp.shape[0] > atom_n_max:
                A[i, :] = A_tmp[0:atom_n_max, 0:atom_n_max]
            else:
                A[i, 0:A_tmp.shape[0], 0:A_tmp.shape[1]] = A_tmp
            color.append(data['color'])
        
        A[np.where((A<=0))] = np.NaN
        lgA = np.log10(A)
        
        x = ticks
        color = np.array(color)
        for j in range(2, atom_n_max+1):
            for i in range(1, j):
                y = lgA[:, j-1, i-1]
                ax = plt.subplot(atom_n_max - 1, atom_n_max - 1, (atom_n_max-1)*(j-2) + i,
                                 axisbg='#FFFFE0')
                ax.set_xticks(x)
                plt.xlim((min(x)-0.5, max(x)+0.5))
                ax.yaxis.set_major_locator(MaxNLocator(max_y_ticks-1))
                try:
                    plt.scatter(x, y, c=color, label='_nolegend_', s=40)
                    lbl = '({0} -> {1})'.format(j,i)
                    ax.text(0.95, 0.95, lbl, fontsize=10, color="#660066", transform=ax.transAxes, ha="right", va="top")                                    
                except:
                    pn.log_.warn('Problem with plotting a subplot {} {}'.format(i,j), calling=self.calling)
                if (j==atom_n_max) & (i==1):  
                    plt.xlabel('Reference #', fontsize=8)
                    plt.ylabel('Log(A$_{ji}$)', fontsize=8)
        plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.12)
        plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
        ax.set_xticks(x)
        
        fig.text(.5, .05, "i (lower level)", fontsize=12, ha='center', color="#660066", va='top')
        fig.text(.05, .5, "j (upper level)", ha='left', va='center', fontsize=12, color="#660066", rotation=90)
        title = "Available transition probabilities for {0} {1}".format(self.elem, int_to_roman(int(self.spec)))
        fig.text(.5, .95, title, color="#191970", fontsize=14, ha='center')
        for i, data in enumerate(self.atom_data):
            x_txt = .95
            y_txt = .85 - .04*i 
            fig.text(x_txt, y_txt, "Ref. %d: %s" %(i+1, data['ID']), color=color[i], ha = 'right') 

        plt.tight_layout(pad=0.1)
        if save:
            plt.savefig(self.atom_rom + "_all_As.pdf")
    
        plt.show()      
        

    def plotRelA(self, ref_data=None, save=False, figsize=None, fignum=None, NLevels=None, fig=None):
        """
        Plot the relative difference of the A of each data set with respect to the reference one

        Parameters:
            - ref_data       reference data set for comparing transition probabilities (default=first data ID)
            - save           if True, save the plot in a file (default: False)
            - figsize        figure size (default: [18, 12])
            - fignum         figure number

        """
        if NLevels is None:
            atom_n_max = self.atom_n_max
        else:
            atom_n_max = NLevels
        if not pn.config.INSTALLED['plt']:
            pn.log_.error('Matplotlib not installed!', calling=self.calling)
        # Commented out because it produced an extra, empty frame - VL 30 Jul 2015
        #if fig is None:
        #    fig = plt.figure(fignum, figsize=figsize)
        #plt.clf()
        ticks = range(atom_n_max + 1)[1:]

        # The As of the reference data set are stored
        if ref_data is None:
            ref_data = self.atom_data[0]['ID']
        for data in self.atom_data:
            if (data['ID'] is ref_data):
                ref_A = data['atom'].getA()
                # Only non-zero values are taken into account to prevent dividing by zero
                nonzero_ref_A_indexes = np.nonzero(ref_A)

        for data in self.atom_data:
            if (data['ID'] is not ref_data):
                A = data['atom'].getA()
                try:
                    # The data to be compared might have fewer levels than the reference one. Only the ones in common are considered.
                    tmp_indexes = np.asarray(nonzero_ref_A_indexes) 
                    up_indexes = tmp_indexes[0][tmp_indexes[0] < data['atom'].atomNLevels]
                    lo_indexes = tmp_indexes[1][tmp_indexes[0] < data['atom'].atomNLevels]
                    # Indexes for which the reference data set is not zero
                    ref_indexes = (up_indexes, lo_indexes)
                    # Ratio of A/A_ref for non-zero A_ref values
                    A_ratio = A[ref_indexes] / ref_A[ref_indexes]
                    nonzero_A_ratio_indexes = np.nonzero(A_ratio)
                    rel_A = np.log10(A_ratio[nonzero_A_ratio_indexes])
                    # Plotting starts
                    #fig = plt.figure()
                    if fig is None:
                        fig = plt.figure(fignum, figsize=figsize)
                    plt.clf()
                    ax = plt.subplot(111, axisbg='#FFFFE0')
                    fig.subplots_adjust(top=0.9)
                    # Physical levels = array levels + 1
                    x = np.asarray(lo_indexes) + 1
                    y = np.asarray(up_indexes) + 1
                    # Size is proportional to difference between As
                    # Numerical value of size adjusted empirically, no particular meaning 
                    size = np.multiply(10000., np.abs(rel_A))
                    # A different color is assigned to data points depending on whether they are smaller or larger
                    # than the reference data.
                    # The ratio is not compared to exactly one to include only truly different values  
                    index_grt = np.where(rel_A > 0.0000)
                    color_grt = np.array(['#FF8C00']).repeat(x[index_grt].size)
                    index_sml = np.where(rel_A <= 0.0)
                    color_sml = np.array(['#FF1480']).repeat(x[index_sml].size)
                    # The range is adjusted to include only values for which there is a difference
                    xmin = np.max([np.min(x[index_sml]), np.min(x[index_grt])])
                    xmax = np.min([np.max(x[index_sml]), np.max(x[index_sml])])
                    ax.set_xticks(ticks)
                    ax.set_xlim((xmin - 0.25, xmax + 0.25)[:])
                    # The x-axis is plotted on the top to mimic the data array in the fits file
                    for tick in ax.xaxis.get_major_ticks():
                        tick.label1On = False
                        tick.label2On = True
                    ymin = np.min([np.min(y[index_sml]), np.min(y[index_grt])])
                    ymax = np.max([np.max(y[index_sml]), np.max(y[index_sml])])
                    ax.set_yticks(ticks)
                    # The y axis increases downward to mimic the data array in the fits file
                    ax.set_ylim((ymin - 0.25, ymax + 0.25)[::-1])
                    # Fake points, just to generate the labels (I don't know how to generate a reasonable label
                    # for data points with markers of different size)
                    plt.scatter(x[index_grt][0], y[index_grt][0], marker='o', s=.0002, c=color_grt[0],
                                label='Positive (max = %.2f)' % (np.max(rel_A)))
                    plt.scatter(x[index_grt][0], y[index_grt][0], marker='o', s=.0002, c=color_sml[0],
                                label='Negative (min = %.2f)' % (np.min(rel_A)))
                    # True data; size reflects divergence from reference data, color whether positive or negative
                    plt.scatter(x[index_grt], y[index_grt], marker='o', s=size[index_grt], c=np.asarray(color_grt), alpha=1.0)
                    plt.scatter(x[index_sml], y[index_sml], marker='o', s=size[index_sml], c=np.asarray(color_sml), alpha=1.0)
                    # X-axis label 
                    ax.text(.5, 1.05, '$i$ (lower level)', transform=ax.transAxes, ha='center', va='bottom')
                    #plt.xlabel('$i$ (lower level)')
                    # Y-axis label
                    plt.ylabel('$j$ (upper level)')
                    # Legend
                    plt.legend(loc='upper right', markerscale=1000., scatterpoints=1, borderpad=1,
                               labelspacing=1, prop=dict(size=12), title=u'Log (%s / %s)' % (data['ID'], ref_data))
                    # Plot title
                    ax.text(.5, -.10, "Relative difference between $A$s data for [%s %s]" 
                                % (self.elem, int_to_roman(int(self.spec))),
                                transform=ax.transAxes, ha='center', va='bottom', color='#191970', size=12)
                except:
                    pn.log_.warn('Problem in plotting relA', calling='DatasetPlot')

                if save:
                    plt.savefig(self.atom_rom + '_' + data['ID'] + "-" + ref_data + "_relA.pdf")
    
                plt.show()      


    def plotOmega(self, save=False, figsize=(18, 12), fignum=1, scan_orders=None, NLevels=None,
                  fig=None):
        """
        Plot the tabulated collision strengths of each data set and the fit that is performed by PyNeb
        
        Parameters:
            - save     Boolean. Determine if the plot is automatically saved in a file (default: False)
            - figsize  List. figure size in inches (default: [18, 12])
            - fignum   Figure Number
            - scan_orders = None or (min_order, max_order) or (min_order, -1) to go until the max. DEPRECATED!!!

        """
        if NLevels is None:
            coll_n_max = self.coll_n_max
        else:
            coll_n_max = NLevels
        # Plotting range somewhat larger than actual data range (less tight) 
        if not pn.config.INSTALLED['plt']:
            pn.log_.error('Matplotlib not installed!', calling=self.calling)
        if fig is None:
            fig = plt.figure(fignum, figsize=figsize)
        plt.autoscale(tight=False)
        # I need two 
        first = True
        first_done = False
        plt.clf()
        
        legend_text = []
        legend_lines = []
#        style_dic = ['-', '--', '-.', ':'] * 10
        for data in self.coll_data:
            # tem_points = tabulated temperature points (different for each data set)
            tem_points = self.tem_in_K(data['atom'].tem_units, data['atom'].getTemArray())
            tem_min = min(tem_points)
            tem_max = max(tem_points)
            # tem_funct = array of temperature values for which the fit is evaluated
            tem_funct = np.linspace(tem_min, tem_max, self.n_tem_points)
            x_dots = np.log10(tem_points)
            x_lines = np.log10(tem_funct)
            if 'COEFF' in data['atom'].CollData.comments.keys():
                coeff = float(data['atom'].CollData.comments['COEFF'])
            else:
                coeff = 1.0
            # Loops over all levels
            for i in range(1, coll_n_max+1):
                for j in range(i + 1, coll_n_max+1):                
                        # N levels require an N-1 x N-1 array of plots
                        # The subplots are arranged in rows and columns according to the upper and lower levels
 
                        ax = plt.subplot(coll_n_max-1, coll_n_max-1, 
                                         (coll_n_max-1) * (j - 2) + i,
                                         axisbg='#FFFFE0')
                        if (i <= data['atom'].collNLevels) and (j <= data['atom'].collNLevels):
                            y_dots = data['atom'].getOmega(tem_points, j, i)
                            try:
                                if 'O_UNIT' not in data['atom'].CollData.comments:
                                    y_dots = coeff * data['atom'].getOmegaArray(j, i)
                                    if y_dots.sum() > 0.0:
                                        ax.plot(x_dots, y_dots, color=data['color'],
                                                marker='*', linestyle='None', label='_nolegend_', markersize=10)
                                    
                                y_dots = data['atom'].getOmega(tem_points, j, i)
                                if y_dots.sum() > 0.0:
                                    ax.plot(x_dots, y_dots, color=data['color'],
                                            marker='o', linestyle='None', label='_nolegend_')
                                y_lines = data['atom'].getOmega(tem_funct, j, i)
                                if y_lines.sum() > 0.:
                                    ax.plot(x_lines, y_lines,
                                            color=data['color'], label='_nolegend_')
                            except:
                                pn.log_.warn('Problem with plotting a data set', calling=self.calling)
# Draws vertical lines at selected temperature values (only once for each subplot, hence "first") 
                        if (first and (data['atom'].collNLevels == coll_n_max)):
                            for tem in self.ref_tem:
                                plt.axvline(tem, c='blue', alpha=0.4, ls=':')
                            lbl = "$\Omega$" + "(" + str(j) + "," + str(i) + ")"
                            ax.text(0.95, 0.95, lbl, transform=ax.transAxes, ha="right", va="top")   
                            first_done = True                                 
                        # Maximum number of ticks to avoid overcrowding
                        ax.xaxis.set_major_locator(MaxNLocator(4))
                        ax.yaxis.set_major_locator(MaxNLocator(3))
            # The line type and identifier of each data set are stored for the legend
            legend_lines.append(mpl.lines.Line2D(tem_funct, tem_funct, color=data['color'], linestyle='-'))
            legend_text.append(data['ID'])
            if first_done:
                first = False
        # The label of the reference temperature values are only plotted above the last plot
        for tem in self.ref_tem:
            ax.text(tem, ax.get_ylim()[1], ' %0.0f K' % (10 ** tem), ha="left", va="bottom", 
                    rotation=65, fontsize=12, color='#FF0000', alpha=0.35)
        # Axis labels, title and legend    
        fig.text(.5, .05, "Log($T_e/K$)", fontsize=16, ha='center', va='top')
        fig.text(.05, .5, "$\Omega$", va='center', fontsize=20, rotation=90)
        fig.text(.5, .95, "[%s %s] collision strengths" % (self.elem, int_to_roman(int(self.spec))), 
                 color="#191970", fontsize=16, ha='center')
#        plt.legend(legend_lines, legend_text, loc='upper right', borderpad=1, 
#                   labelspacing=1, bbox_to_anchor=(1, 1 * coll_n_max))
        plt.legend(legend_lines, legend_text, loc='upper right', borderpad=1, 
                   labelspacing=1, bbox_to_anchor=(1, coll_n_max - 1))
        #plt.tight_layout()
        if save:
            plt.savefig(self.atom_rom + "_CS.pdf")
        
        plt.show()


    def tem_in_K(self, tem_units, tem):
        """
        Convert the temperature from the unit of the fits file into K

        Parameters:
            - tem_units     'log(K)' or 'K/1000'
            - tem           temperature
        """
        
        if (tem_units == "log(K)"):
            return np.power(10., tem)
        elif (tem_units == "K/10000"):
            return np.multiply(tem, 1.e4)
        else: #T in Kelvin in the fits file
            return tem
                
    
