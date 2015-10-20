import matplotlib.pyplot as plt
import pyneb as pn


def p1():
    # Set restore to True if the emission maps have already been generated, False otherwise
    ### General settings
    # Setting verbosity level. Enter pn.my_logging? for details
    pn.log_.level = 3
    # Adopt an extinction law
    #extinction_law = 'CCM 89'
    # Define the data file
    obs_data = 'IC2165.dat'
    
    ### Read observational data
    # define an Observation object and assign it to name 'obs'
    obs = pn.Observation()
    # fill obs with data read from file obs_data, with lines varying across rows and a default percent error on line intensities
    obs.readData(obs_data, fileFormat='lines_in_rows', err_default=0.05, corrected=True)
    
    ### Include the diagnostics of interest
    # instantiate the Diagnostics class
    diags = pn.Diagnostics()
    # include in diags the relevant line ratios
    diags.getAllDiagLabels()
    diags.getAllDiags()
    diags.addDiag(['[NI] 5198/5200',
                    '[NII] 5755/6584',
                    '[OI] 5579/6302', 
                    '[OII] 3726/3729', 
                    '[OII] 3727+/7325+',
                    '[OIII] 4363/5007+', 
                    '[SII] 6731/6716', 
                    '[SII] 4072+/6720+',
                    '[SIII] 6312/9069',
                    '[ArIII] 5192/7136',
                    '[ArIV] 4740/4711',
                    '[ClIII] 5538/5518',
                    ])
    
    diags.addDiag('[ClIV] 5323/7531',('Cl4', 'L(5323)/L(7531)', 'RMS([E(7531),E(5323)])'))
    
    # Create the emission maps to be compared to the observation data
    # To see the default parameters, do pn.getEmisGridDict? in ipython
    # The files go to ~/.pypics
    emisgrids = pn.getEmisGridDict(atomDict=diags.atomDict, den_max=1e6)
    
    ### Plot
    plt.figure(1)
    # Create the contour plot as the intersection of tem-den emission maps with dereddened line ratios
    diags.plot(emisgrids, obs)
    # Put the title
    plt.title('IC 2165 Optical diagnostics')
    # Display the plot
    plt.show()
    plt.savefig('IC2165-diag-opt.pdf')
    
def p2(restore=True):
    ### General settings
    # Setting verbosity level. Enter pn.my_logging? for details
    pn.log_.level = 3
    # Adopt an extinction law
    #extinction_law = 'CCM 89'
    # Define the data file
    obs_data = 'IC2165_UV.dat'
    
    ### Read observational data
    # define an Observation object and assign it to name 'obs'
    obs = pn.Observation()
    # fill obs with data read from file obs_data, with lines varying across rows and a default percent error on line intensities
    obs.readData(obs_data, fileFormat='lines_in_rows', err_default=0.05, corrected=True)
    
    ### Include the diagnostics of interest
    # instantiate the Diagnostics class
    diags = pn.Diagnostics()
    # include in diags the relevant line ratios
    diags.addDiag(['[NI] 5198/5200',
                    '[NII] 5755/6584',
                    '[OI] 5579/6302', 
                    '[OII] 3726/3729', 
                    '[OII] 3727+/7325+',
                    '[OIII] 4363/5007', 
                    '[SII] 6731/6716', 
                    '[SII] 4072+/6720+',
                    '[SIII] 6312/9069',
                    '[ArIII] 5192/7136',
                    '[ArIV] 4740/4711',
                    '[ClIII] 5538/5518',
                    '[CIII] 1909/1907', 
                    '[OIII] 1666/5007', 
                    '[NeV] 1575/3426'
                    ])
    
    diags.addDiag('[ClIV] 5323/7531',('Cl4', 'L(5323)/L(7531)', 'RMS([E(7531),E(5323)])'))
    
    # Create the emission maps to be compared to the observation data
    # To see the default parameters, do pn.getEmisGridDict? in ipython
    emisgrids = pn.getEmisGridDict(atomDict=diags.atomDict, den_max=1e6)
    
    ### Plot
    plt.figure(2)
    # Create the contour plot as the intersection of tem-den emission maps with dereddened line ratios
    diags.plot(emisgrids, obs)
    # Put the title
    plt.title('IC 2165 Optical + UV diagnostics')
    # Display the plot
    plt.show()
    plt.savefig('IC2165-diag-opt+UV.pdf')

def p3(restore = True, fignum=3, pypic_path=pn.config.pypic_path):
    # --- pregunta 3 Adding IR lines-------------------------------
    # Define where emission maps are to be stored (restore = False) or read from (restore = True)
    # Create the directory before running this program
    obs_data = 'IC2165_IR.dat'
    obs = pn.Observation()
    # fill obs with data read from file obs_data, with lines varying across rows and a default percent error on line intensities
    obs.readData(obs_data, fileFormat='lines_in_rows', err_default=0.05, corrected=True)
    temp = 14000
    dens = 10**3.5
    # Compute theoretical H5-4/Hbeta ratio from Hummer and Storey
    IR2Opt_theo = pn.getRecEmissivity(temp, dens, 5, 4) / pn.getRecEmissivity(temp, dens, 4, 2)
    IR2Opt_obs = obs.getLine(label='H1_4.1m').corrIntens / obs.getLine(label='H1_4861A').corrIntens
    for line in obs.lines: 
        if line.label[-1] == 'm':
            line.corrIntens *= IR2Opt_theo/IR2Opt_obs

    diags = pn.Diagnostics()

    diags.addDiag(['[NI] 5198/5200',
                    '[NII] 5755/6584',
                    '[OI] 5579/6302', 
                    '[OII] 3726/3729', 
                    '[OII] 3727+/7325+',
                    '[OIII] 4363/5007', 
                    '[SII] 6731/6716', 
                    '[SII] 4072+/6720+',
                    '[SIII] 6312/9069',
                    '[ArIII] 5192/7136',
                    '[ArIV] 4740/4711',
                    '[ClIII] 5538/5518',
                    '[CIII] 1909/1907', 
                    '[OIII] 1666/5007', 
                    '[NeV] 1575/3426',
                    '[OIII] 51m/88m',
                    '[NeIII] 15.6m/36.0m',
                    '[NeV] 14.3m/24.2m',
                    '[SIII] 18.7m/33.6m',
                    '[NeIII] 3869/15.6m',
                    '[OIII] 5007/88m',
                    '[ArIII] 7136/9m',
                    '[SIII] 6312/18.7m',
                    ])
    diags.addDiag('[ClIV] 5323/7531',('Cl4', 'L(5323)/L(7531)', 'RMS([E(7531),E(5323)])'))
    
    # Create the emission maps to be compared to the observation data
    # To see the default parameters, do pn.getEmisGridDict? in ipython
    emisgrids = pn.getEmisGridDict(atomDict=diags.atomDict, den_max=1e6, pypic_path=pypic_path)
    
    ### Plot
    plt.figure(fignum)
    # Create the contour plot as the intersection of tem-den emission maps with dereddened line ratios
    diags.plot(emisgrids, obs)
    # Put the title
    plt.title('IC 2165 Optical + UV + IR diagnostics')
    
    # Display the plot
    plt.show()
    plt.savefig('IC6165-diag-opt+UV+IR.pdf')
