import numpy as np
import pyneb as pn
from pyneb.utils.init import ELEM_LIST, SPEC_LIST
from pyneb.utils.misc import parseAtom
import random

class ICF(object):

    def __init__(self):
        """
        ICF tool.
        No parameters for the instantiation.
        """
        
        self.log_ = pn.log_ 
        self.calling = 'ICF'
                
        self._init_all_icfs()
        self._max_pass = 3
        
    def _init_all_icfs(self):
        # Dictionary of ICF recipes
        self.all_icfs = {'direct_He.23':{'elem': 'He',
                                       'atom': 'abun["He2"] + abun["He3"]',
                                       'icf': '1',
                                       'type': 'All',
                                       'comment': 'just summing visible ions'},
                         'direct_N.23':{'elem': 'N',
                                       'atom': 'abun["N2"] + abun["N3"]',
                                       'icf': '1',
                                       'type': 'All',
                                       'comment': 'just summing visible ions'},
                         'direct_N.234':{'elem': 'N',
                                       'atom': 'abun["N2"] + abun["N3"] + abun["N4"]',
                                       'icf': '1',
                                       'type': 'All',
                                       'comment': 'just summing visible ions'},
                         'direct_N.2345':{'elem': 'N',
                                       'atom': 'abun["N2"] + abun["N3"] + abun["N4"] + abun["N5"]',
                                       'icf': '1',
                                       'type': 'All',
                                       'comment': 'just summing visible ions'},
                         'direct_O.23':{'elem': 'O',
                                       'atom': 'abun["O2"] + abun["O3"]',
                                       'icf': '1',
                                       'type': 'All',
                                       'comment': 'just summing visible ions'},
                         'direct_O.234':{'elem': 'O',
                                       'atom': 'abun["O2"] + abun["O3"] + abun["O4"]',
                                       'icf': '1',
                                       'type': 'All',
                                       'comment': 'just summing visible ions'},
                         'direct_O.2345':{'elem': 'O',
                                       'atom': 'abun["O2"] + abun["O3"] + abun["O4"] + abun["O5"]',
                                       'icf': '1',
                                       'type': 'All',
                                       'comment': 'just summing visible ions'},
                         'direct_S.23':{'elem': 'S',
                                       'atom': 'abun["S2"] + abun["S3"]',
                                       'icf': '1',
                                       'type': 'All',
                                       'comment': 'just summing visible ions'},
                         'direct_S.234':{'elem': 'S',
                                       'atom': 'abun["S2"] + abun["S3"] + abun["S4"]',
                                       'icf': '1',
                                       'type': 'All',
                                       'comment': 'just summing visible ions'},
                         'direct_S.2345':{'elem': 'S',
                                       'atom': 'abun["S2"] + abun["S3"] + abun["S4"] + abun["S5"]',
                                       'icf': '1',
                                       'type': 'All',
                                       'comment': 'just summing visible ions'},
                         'direct_Ne.23':{'elem': 'Ne',
                                       'atom': 'abun["Ne2"] + abun["Ne3"]',
                                       'icf': '1',
                                       'type': 'All',
                                       'comment': 'just summing visible ions'},
                         'direct_Ne.235':{'elem': 'Ne',
                                       'atom': 'abun["Ne2"] + abun["Ne3"] + abun["Ne5"]',
                                       'icf': '1',
                                       'type': 'All',
                                       'comment': 'just summing visible ions'},
                         'direct_Ne.2345':{'elem': 'Ne',
                                       'atom': 'abun["Ne2"] + abun["Ne3"] + abun["Ne4"] + abun["Ne5"]',
                                       'icf': '1',
                                       'type': 'All',
                                       'comment': 'just summing visible ions'},
                         'direct_Ne.345':{'elem': 'Ne',
                                       'atom': 'abun["Ne3"] + abun["Ne4"] + abun["Ne5"]',
                                       'icf': '1',
                                       'type': 'All',
                                       'comment': 'just summing visible ions'},
                         'direct_Ne.2356':{'elem': 'Ne',
                                       'atom': 'abun["Ne2"] + abun["Ne3"] + abun["Ne5"] + abun["Ne6"]',
                                       'icf': '1',
                                       'type': 'All',
                                       'comment': 'just summing visible ions'},
                         'direct_Mg.45':{'elem': 'Mg',
                                       'atom': 'abun["Mg4"] + abun["Mg5"]',
                                       'icf': '1',
                                       'type': 'All',
                                       'comment': 'just summing visible ions'},
                         'direct_Ar.23':{'elem': 'Ar',
                                       'atom': 'abun["Ar2"] + abun["Ar3"]',
                                       'icf': '1',
                                       'type': 'All',
                                       'comment': 'just summing visible ions'},
                         'direct_Ar.234':{'elem': 'Ar',
                                       'atom': 'abun["Ar2"] + abun["Ar3"] + abun["Ar4"]',
                                       'icf': '1',
                                       'type': 'All',
                                       'comment': 'just summing visible ions'},
                         'direct_Ar.345':{'elem': 'Ar',
                                       'atom': 'abun["Ar3"] + abun["Ar4"] + abun["Ar5"]',
                                       'icf': '1',
                                       'type': 'All',
                                       'comment': 'just summing visible ions'},
                         'direct_Cl.23':{'elem': 'Cl',
                                       'atom': 'abun["Cl2"] + abun["Cl3"]',
                                       'icf': '1',
                                       'type': 'All',
                                       'comment': 'just summing visible ions'},
                         'direct_Cl.34':{'elem': 'Cl',
                                       'atom': 'abun["Cl3"] + abun["Cl4"]',
                                       'icf': '1',
                                       'type': 'All',
                                       'comment': 'just summing visible ions'},
                         'direct_Cl.234':{'elem': 'Cl',
                                       'atom': 'abun["Cl2"] + abun["Cl3"] + abun["Cl4"]',
                                       'icf': '1',
                                       'type': 'All',
                                       'comment': 'just summing visible ions'},
                         'TPP77_13': {'elem': 'O',
                                      'atom': 'abun["O2"] + abun["O3"]',
                                      'icf': '(abun["He2"] + abun["He3"]) / abun["He2"]',
                                      'type': 'PNe',
                                      'comment': ''},
                         'TPP77_14': {'elem': 'N',
                                      'atom': 'abun["N2"]',
                                      'icf': '(abun["O2"] + abun["O3"]) / abun["O2"]',
                                      'type': 'PNe',
                                      'comment': ''},
                         'TPP77_15': {'elem': 'Ne',
                                      'atom': 'abun["Ne3"]',
                                      'icf': '(abun["O2"] + abun["O3"]) / abun["O3"]',
                                      'type': 'PNe',
                                      'comment': ''},
                         'PHCD07_12': {'elem': 'Ne',
                                       'atom': 'abun["Ne3"]',
                                       'icf': '(abun["O2"] + abun["O3"]) / abun["O3"]',
                                       'type': 'HII',
                                       'comment': 'Based on a grid of photoionization models'},
                         'PHCD07_13': {'elem': 'Ne',
                                       'atom': 'abun["Ne3"]',
                                       'icf': '(0.753 + 0.142 * abun["O3"] / (abun["O2"] + abun["O3"]) + 0.171 / abun["O3"] * (abun["O2"] + abun["O3"]))',
                                       'type': 'HII',
                                       'comment': 'Based on a grid of photoionization models'},
# PHCD07_16 and 17 added on December 11, 2014
                         'PHCD07_16': {'elem': 'Ar',
                                       'atom': '(abun["Ar3"] + abun["Ar4"])',
                                       'icf': '(0.928 + 0.364 * (1 - abun["O3"] / (abun["O2"] + abun["O3"])) + 0.006 / (1 - abun["O3"] / (abun["O2"] + abun["O3"])))',
                                       'type': 'HII',
                                       'comment': 'Based on a grid of photoionization models'},
                         'PHCD07_17': {'elem': 'Ar',
                                       'atom': 'abun["Ar3"]',
                                       'icf': '(0.596 + 0.967 * (1 - abun["O3"] / (abun["O2"] + abun["O3"])) + 0.077 / (1 - abun["O3"] / (abun["O2"] + abun["O3"])))',
                                       'type': 'HII',
                                       'comment': 'Based on a grid of photoionization models'},
# Cited in PHCD07 but originally proposed by ITL94. Authorship changed in the code 3/Nov/2014. 
#                         'PHCD07_14': {'elem': 'Ar',
                         'ITL94_19': {'elem': 'Ar',
                                       'atom': '(abun["Ar3"] + abun["Ar4"])',
                                       'icf': '(0.99 + 0.091 * abun["O2"] / (abun["O2"] + abun["O3"]) - 1.14 * (abun["O2"] / (abun["O2"] + abun["O3"]))**2. + 0.077 * (abun["O2"] / (abun["O2"] + abun["O3"]))**3.)**-1.',
                                       'type': 'HII',
                                       'comment': 'Based on a grid of photoionization models'},
# Cited in PHCD07 but originally proposed by ITL94. Authorship changed in the code 3/Nov/2014. 
#                         'PHCD07_15': {'elem': 'Ar',
                         'ITL94_20': {'elem': 'Ar',
                                       'atom': 'abun["Ar3"]',
                                       'icf': '(0.15 + 2.39 * abun["O2"] / (abun["O2"] + abun["O3"]) - 2.64 * (abun["O2"] / (abun["O2"] + abun["O3"]))**2.)**-1.',
                                       'type': 'HII',   
                                       'comment': 'Based on a grid of photoionization models'},
                         'PTPR92_21':{'elem': 'He',
                                       'atom': 'abun["He2"]',
                                       'icf': '(1 + abun["S2"] / (elem_abun["KB94_A36.10"] -  abun["S2"]))',
                                       'type': 'HII',
                                       'comment': 'Based on ionization potentials'},
                         
                         'KB94_A0': {'elem': 'N',
                                    'atom': 'abun["N2"] + abun["N3"] + abun["N4"] + abun["N5"]',
                                    'icf': '1',
                                     'type': 'PNe',
                                     'comment': 'Sum of all N ions from N+ to N4+'},
                         
                         'KB94_A1.6': {'elem': 'N',
                                     'atom': 'abun["N2"]',
                                     'icf': 'elem_abun["KB94_A6"]  / abun["O2"]',
                                     'type': 'PNe',
# Comment clarified on 12 Dec 2014 
#                                      'comment': 'Based on a grid of photoionization models. To be used if both O4 and N5 detected'},
                                      'comment': 'Based on a grid of photoionization models. To be used if both O4 and N5 are detected, but O5 is not'},
                         'KB94_A1.8': {'elem': 'N',
                                     'atom': 'abun["N2"]',
                                     'icf': 'elem_abun["KB94_A8"]  / abun["O2"]',
                                     'type': 'PNe',
# wrong comment. Corrected 12 Dec 2014
#                                     'comment': 'Based on a grid of photoionization models. To be used if both O4 and N5 detected'},
                                     'comment': 'Based on a grid of photoionization models. To be used if both O3 and are N5 detected, but O4 is not'},
                         'KB94_A1.10': {'elem': 'N',
                                     'atom': 'abun["N2"]',
                                     'icf': 'elem_abun["KB94_A10"]  / abun["O2"]',
                                     'type': 'PNe',
# wrong comment. Corrected 12 Dec 2014
#                                      'comment': 'Based on a grid of photoionization models. To be used if both O4 and N5 detected'},
                                      'comment': 'Based on a grid of photoionization models. To be used if only O2 and O3 are detected'},
                         'KB94_A6': {'elem': 'O',
                                     'atom': 'abun["O2"] + abun["O3"] + abun["O4"]',
                                      'icf': '1 / (1 - 0.95 * abun["N5"] / (abun["N2"] + abun["N3"] + abun["N4"] + abun["N5"]))',
                                      'type': 'PNe',
# Comment clarified on 12 Dec 2014 
#                                      'comment': 'Based on a grid of photoionization models. To be used if both O4 and N5 detected'},
                                      'comment': 'Based on a grid of photoionization models. To be used if both O4 and N5 are detected, but O5 is not'},
                         'KB94_A8': {'elem': 'O',
                                     'atom': 'abun["O2"] + abun["O3"]',
                                     'icf': '(abun["N2"] + abun["N3"] + abun["N4"] + abun["N5"]) / (abun["N2"] + abun["N3"])',
                                     'type': 'PNe',
# Comment clarified on 12 Dec 2014 
#                                     'comment': 'Based on a grid of photoionization models. To be used if N5 detected and O4 not detected'},
                                     'comment': 'Based on a grid of photoionization models. To be used if N5 detected, but neither O3 nor O4 are'},
                         'KB94_A10': {'elem': 'O',
                                      'atom': 'abun["O2"] + abun["O3"]',
                                      'icf': '((abun["He2"] + abun["He3"]) / abun["He2"])**(2./3.)',
                                      'type': 'PNe',
                                      'comment': 'Based on a grid of photoionization models. To be used if only O2 and O3 are detected'},
                         'KB94_A10b': {'elem': 'O',
                                      'atom': 'abun["O2"] + abun["O3"]',
                                      'icf': '1',
                                      'type': 'PNe',
                                      'comment': 'Based on a grid of photoionization models. To be used if only O2 and O3 are detected and no HeII is seen'},
# Wrong  12 Dec 2014
#                         'KB94_A10.5': {'elem': 'C',
#                                        'atom': 'abun["C2"] + abun["C3"] + abun["C4"] + abun["C5"]',
#                                        'icf': '1',
#                                        'type': 'PNe',
#                                        'comment': 'Based on a grid of photoionization models. To be used if no He3 detected and C2, C3, and C4 detected'},
# Misleading comment clarified on 12 Dec 2014 
                        'KB94_A12': {'elem': 'C',
                                      'atom': 'abun["C3"] + abun["C4"]',
                                      'icf': '(abun["O2"] + abun["O3"]) / abun["O3"]',
                                      'type': 'PNe',
#                                      'comment': 'Based on a grid of photoionization models. To be used if no C II lines detected'},
                                      'comment': 'Based on a grid of photoionization models. To be used if no C2 lines available and C4 and He3 are not present'},
# Wrong because eq. A10 does not point to eq. A6 (no O4 present). 12 Dec 2014
#                         'KB94_A13.6': {'elem': 'C',
#                                      'atom': 'abun["C3"]',
#                                      'icf': 'elem_abun["KB94_A6"] / abun["O3"]',
#                                      'type': 'PNe',
#                                      'comment': 'Based on a grid of photoionization models. To be used if no C II lines detected'},
#
# Wrong because eq. A10 does not point to eq. A8 (no N5 present)
#                         'KB94_A13.8': {'elem': 'C',
#                                      'icf': 'elem_abun["KB94_A8"] / abun["O3"]',
#                                      'atom': 'abun["C3"]',
#                                      'type': 'PNe',
#                                      'comment': 'Based on a grid of photoionization models. To be used if no C II lines detected'},
# Misleading comment clarified on 12 Dec 2014 
                         'KB94_A13.10': {'elem': 'C',
                                      'atom': 'abun["C3"]',
                                      'icf': 'elem_abun["KB94_A10"] / abun["O3"]',
                                      'type': 'PNe',
#                                     'comment': 'Based on a grid of photoionization models. To be used if no C II lines detected'},
                                      'comment': 'Based on a grid of photoionization models. To be used if C5 is not present; only C3 lines are detected; O3 also detected; He2 and He3 also detected'},
                         'KB94_A13.10b': {'elem': 'C',
                                      'atom': 'abun["C3"]',
                                      'icf': 'elem_abun["KB94_A10b"] / abun["O3"]',
                                      'type': 'PNe',
#                                     'comment': 'Based on a grid of photoionization models. To be used if no C II lines detected'},
                                      'comment': 'Based on a grid of photoionization models. To be used if C5 is not present; only C3 lines are detected; O3 also detected; He2 also detected, He3 Not'},
# Clarified comment 12 Dec 2014
                         'KB94_A16': {'elem': 'C',
                                      'atom': 'abun["C2"] + abun["C3"] + abun["C4"]',
                                      'icf': '1 / (1 - 2.7 * abun["N5"] / (abun["N2"] + abun["N3"] + abun["N4"] + abun["N5"]))',
                                      'type': 'PNe',
#                                     'comment': 'Based on a grid of photoionization models. Valid for high excitation PNe if icf < 5'},
                                       'comment': 'Based on a grid of photoionization models. To be used in high excitation PNe where C5 is present (as inferred from He3 and N5). Valid if icf < 5'},
# Clarified comment 12 Dec 2014
                         'KB94_A19': {'elem': 'C',
                                      'atom': 'abun["C2"] + abun["C3"] + abun["C4"]',
                                      'icf': '(1 + abun["N5"] / (abun["N2"] + abun["N3"] + abun["N4"]))',
                                      'type': 'PNe',
#                                      'comment': 'Based on a grid of photoionization models. Valid for extremely high excitation PNe if icf > 5'},
                                       'comment': 'Based on a grid of photoionization models. To be used in very high excitation PNe where C5 is present (as inferred from He3 and N5). Valid if icf > 5'},
                         'KB94_A21': {'elem': 'C',
                                      'atom': 'abun["C2"] + abun["C3"] + abun["C4"]',
                                      'icf': '((abun["He2"] + abun["He3"]) / abun["He2"])**(1./3.)',
                                      'type': 'PNe',
                                      'comment': 'Based on a grid of photoionization models. Valid if He3 present and N5 absent'},
# Clarified comment 12 Dec 2014
                         'KB94_A26': {'elem': 'C',
                                      'atom': 'abun["C2"] + abun["C3"] + abun["C4"]',
                                      'icf': '(abun["O2"] + abun["O3"]) / abun["O3"] * ((abun["He2"] + abun["He3"]) / abun["He2"])**(1./3.)',
                                      'type': 'PNe',
#                                      'comment': 'Based on a grid of photoionization models. Valid if C5 is present but N4 or N5 not seen'},
                                      'comment': 'Based on a grid of photoionization models. Valid if C5 is present but only optical data are available, so N4 and N5 are not seen'},
                         'KB94_A27': {'elem': 'Ne',
                                      'atom': 'abun["Ne3"] + abun["Ne5"]',
                                      'icf': '1.5',
                                      'type': 'PNe',
# Wrong comment. Corrected 19 Dec 2014
#                                      'comment': 'Based on a grid of photoionization models. Valid if N4 or N5 not seen'},
                                      'comment': 'Based on a grid of photoionization models. Valid if Ne3 and Ne5 are seen and Ne4 is not'},
#
# Unsure about this icf. I don't know if eq. 28 can be mixed with eq. 6. In any case, I explicitly set out the conditions of eq. 6 in the comment.
                         'KB94_A28.6':{'elem': 'Ne',
                                       'atom': 'abun["Ne3"]',
                                       'icf': 'elem_abun["KB94_A6"]  / abun["O3"]',
                                      'type': 'PNe',
# Wrong comment. Corrected 19 Dec 2014
#                                      'comment': 'Based on a grid of photoionization models. Valid if N4 or N5 not seen'},
                                      'comment': 'Based on a grid of photoionization models. Valid if neither Ne4 nor Ne5 are seen, both O4 and N5 are detected, but O5 is not'},
#
# Unsure about this icf. I don't know if eq. 28 can be mixed with eq. 8. In any case, I explicitly set out the conditions of eq. 6 in the comment.
                         'KB94_A28.8':{'elem': 'Ne',
                                     'atom': 'abun["Ne3"]',
                                     'icf': 'elem_abun["KB94_A8"]  / abun["O3"]',
                                      'type': 'PNe',
# Wrong comment. Corrected 19 Dec 2014
#                                      'comment': 'Based on a grid of photoionization models. Valid if N4 or N5 not seen'},
                                      'comment': 'Based on a grid of photoionization models. Valid if neither Ne4 nor Ne5 are seen, N5 is detected, but neither O3 nor O4 are'},
#
# Unsure about this icf. I don't know if eq. 28 can be mixed with eq. 10. In any case, I explicitly set out the conditions of eq. 6 in the comment.
                         'KB94_A28.10':{'elem': 'Ne',
                                     'atom': 'abun["Ne3"]',
                                     'icf': 'elem_abun["KB94_A10"]  / abun["O3"]',
                                      'type': 'PNe',
# Wrong comment. Corrected 19 Dec 2014
#                                      'comment': 'Based on a grid of photoionization models. Valid if N4 or N5 not seen'},
                                      'comment': 'Based on a grid of photoionization models. Valid if neither Ne4 nor Ne5 are seen, if C5 is not present; only C3 lines are detected; O3 also detected; He2 and He3 also detected'},
                         'KB94_A28.10b':{'elem': 'Ne',
                                     'atom': 'abun["Ne3"]',
                                     'icf': 'elem_abun["KB94_A10b"]  / abun["O3"]',
                                      'type': 'PNe',
# Wrong comment. Corrected 19 Dec 2014
#                                      'comment': 'Based on a grid of photoionization models. Valid if N4 or N5 not seen'},
                                      'comment': 'Based on a grid of photoionization models. Valid if neither Ne4 nor Ne5 are seen, if C5 is not present; only C3 lines are detected; O3 also detected; He2 detected, He3 not.'},
# 26 Dec 2014 These icfs are temporarily commented out because they must be checked. The comments, at the very least, are wrong, and thee icfs may be physically wrong.  
                         'KB94_A30.0': {'elem': 'Ar',
                                     'atom': 'abun["Ar3"] + abun["Ar4"] + abun["Ar5"]',
                                     'type': 'PNe',
                                     'icf': '1./(1. - abun["N2"] / elem_abun["KB94_A0"])',
                                      'comment': 'Based on a grid of photoionization models. To be used if N+ to N4+ are seen'},
                         'KB94_A30.10': {'elem': 'Ar',
                                     'atom': 'abun["Ar3"] + abun["Ar4"] + abun["Ar5"]',
                                     'icf': '1./(1. - abun["N2"] / elem_abun["KB94_A1.10"])',
                                     'type': 'PNe',
                                      'comment': 'Based on a grid of photoionization models. To be used if He2 and He3 are detected.'},
                         'KB94_A30.10b': {'elem': 'Ar',
                                     'atom': 'abun["Ar3"] + abun["Ar4"] + abun["Ar5"]',
                                     'icf': '1./(1. - abun["N2"] / elem_abun["KB94_A1.10b"])',
                                     'type': 'PNe',
                                      'comment': 'Based on a grid of photoionization models. To be used if He2 is detected.and He3 not.'},
# end block of commented icfs
                         'KB94_A32': {'elem': 'Ar',
                                     'atom': 'abun["Ar3"]',
                                     'icf': '1.87',
                                     'type': 'PNe',
# Wrong comment. Corrected 26 Dec 2014
#                                      'comment': 'Based on a grid of photoionization models. To be used if both O4 and N5 detected'},
                                     'comment': 'Based on a grid of photoionization models. To be used when only Ar3 is observed'},
                         'KB94_A36.6':{'elem': 'S',
                                     'atom': 'abun["S2"] + abun["S3"]',
                                     'icf': ' (1 - (1 - abun["O2"]/elem_abun["KB94_A6"])**3)**(-1./3.)',
                                      'type': 'PNe',
# Wrong comment. Corrected 26 Dec 2014
#                                      'comment': 'Based on a grid of photoionization models. To be used if both O4 and N5 detected'},
                                      'comment': 'Based on a grid of photoionization models. Valid if both S2 and S3 detected'},
                         'KB94_A36.8':{'elem': 'S',
                                     'atom': 'abun["S2"] + abun["S3"]',
                                     'icf': ' (1-(1 - abun["O2"]/elem_abun["KB94_A8"])**3)**(-1./3.)',
                                      'type': 'PNe',
# Wrong comment. Corrected 26 Dec 2014
#                                      'comment': 'Based on a grid of photoionization models. To be used if both O4 and N5 detected'},
                                      'comment': 'Based on a grid of photoionization models. Valid if both S2 and S3 detected'},
                         'KB94_A36.10':{'elem': 'S',
                                     'atom': 'abun["S2"] + abun["S3"]',
                                     'icf': ' (1-(1 - abun["O2"]/elem_abun["KB94_A10"])**3)**(-1./3.)',
                                      'type': 'PNe',
# Wrong comment. Corrected 26 Dec 2014
#                                      'comment': 'Based on a grid of photoionization models. To be used if both O4 and N5 detected'},
                                      'comment': 'Based on a grid of photoionization models. Valid if both S2 and S3 detected. He2 and He3 detected.'},
                         'KB94_A36.10b':{'elem': 'S',
                                     'atom': 'abun["S2"] + abun["S3"]',
                                     'icf': ' (1-(1 - abun["O2"]/elem_abun["KB94_A10b"])**3)**(-1./3.)',
                                      'type': 'PNe',
# Wrong comment. Corrected 26 Dec 2014
#                                      'comment': 'Based on a grid of photoionization models. To be used if both O4 and N5 detected'},
                                      'comment': 'Based on a grid of photoionization models. Valid if both S2 and S3 detected. He2 detected, He3 not.'},
# Added 26 Dec 2014
# 26 Dec 2014 These icfs are temporarily commented out because they must be still be checked. 
#                         'KB94_A38.6':{'elem': 'S',
#                                     'atom': 'abun["S2"]',
#                                     'icf': '((1 - (1 - abun["O2"]/elem_abun["KB94_A6"])**3)**(-1./3.)*(5.677 + (abun["O3"]/abun["O2"])**(0.433))',
#                                      'type': 'PNe',
#                                      'comment': 'Based on a grid of photoionization models and the S3/S2 ratio of a sample of PNe. Valid if S2 is detected but S3 is not detected'},
# Added 26 Dec 2014
#                         'KB94_A38.8':{'elem': 'S',
#                                     'atom': 'abun["S2"]',
#                                     'icf': '((1 - (1 - abun["O2"]/elem_abun["KB94_A8"])**3)**(-1./3.)*(5.677 + (abun["O3"]/abun["O2"])**(0.433))',
#                                      'type': 'PNe',
#                                      'comment': 'BBased on a grid of photoionization models and the S3/S2 ratio of a sample of PNe. Valid if S2 is detected but S3 is not detected'},
# Added 26 Dec 2014
#                         'KB94_A38.10':{'elem': 'S',
#                                     'atom': 'abun["S2"]',
#                                     'icf': '((1 - (1 - abun["O2"]/elem_abun["KB94_A10"])**3)**(-1./3.)*(5.677 + (abun["O3"]/abun["O2"])**(0.433))',
#                                      'type': 'PNe',
#                                      'comment': 'Based on a grid of photoionization models and the S3/S2 ratio of a sample of PNe. Valid if S2 is detected but S3 is not detected'},
# end commented block
                         'KH01_4a': {'elem': 'He',
                                     'atom': 'abun["He2"] + abun["He3"]',
                                     'icf': '1',
                                     'type': 'PNe',
                                     'comment': 'Based on a grid of photoionization models'},
                         'KH01_4b': {'elem': 'O',
                                     'atom': 'abun["O2"] + abun["O3"]',
                                     'icf': '(abun["He2"] + abun["He3"]) / abun["He2"]',
                                     'type': 'PNe',
                                     'comment': 'Based on a grid of photoionization models'},
                         'KH01_4c': {'elem': 'N',
                                     'atom': 'abun["N2"]',
                                     'icf': '(abun["O2"] + abun["O3"]) / abun["O2"] * (abun["He2"] + abun["He3"]) / abun["He2"]',
                                     'type': 'PNe',
                                     'comment': 'Based on a grid of photoionization models'},
                         'KH01_4d': {'elem': 'Ne',
                                     'atom': 'abun["Ne3"]',
                                     'icf': '(abun["O2"] + abun["O3"]) / abun["O3"] * (abun["He2"] + abun["He3"]) / abun["He2"]',
                                     'type': 'PNe',
                                     'comment': 'Based on a grid of photoionization models'},
                         'KH01_4e': {'elem': 'S',
                                     'atom': 'abun["S2"] + abun["S3"]',
                                     'icf': '10**(-0.017+0.18*np.log10(elem_abun["KH01_4b"]/ abun["O2"])-0.11*np.log10(elem_abun["KH01_4b"]/ abun["O2"])**2+0.072*np.log10(elem_abun["KH01_4b"]/ abun["O2"])**3)',
                                     'type': 'PNe',
                                     'comment': 'Based on a grid of photoionization models'},
                         'KH01_4f': {'elem': 'Cl',
                                     'atom': 'abun["Cl3"] + abun["Cl4"]',
                                     'icf': '(abun["He2"] + abun["He3"]) / abun["He2"]',
                                     'type': 'PNe',
                                     'comment': 'Based on a grid of photoionization models'},
                         'KH01_4g': {'elem': 'Ar',
                                     'atom': 'abun["Ar3"] + abun["Ar4"]',
                                     'icf': '(abun["He2"] + abun["He3"]) / abun["He2"] / (1 -  abun["N2"] / (elem_abun["KH01_4c"]))',
                                     'type': 'PNe',
                                     'comment': 'Based on a grid of photoionization models'},
                         'KH01_4txt': {'elem': 'Ar',
                                     'atom': 'abun["Ar3"]',
                                     'icf': '1',
                                     'type': 'PNe',
                                     'comment': 'Based on a grid of photoionization models. To be used when log(O/O2) <= 0.6'},
                         'Ial06_16': {'elem': 'O',
                                      'atom': 'abun["O2"] + abun["O3"]',
                                      'icf': '1',
                                      'type': 'HII',
                                      'comment': 'Based on a grid of photoionization models. To be used if no 4686 detected'},
                         'Ial06_17': {'elem': 'O',
                                      'atom': 'abun["O2"] + abun["O3"]',
                                      'icf': '1 + 0.5 * abun["He3"] / (abun["He2"] + abun["He3"])',
                                      'type': 'HII',
                                      'comment': 'Based on a grid of photoionization models. To be used if 4686 detected'},
# Wrong. It is not a general ICF recipe but a particular case. Changed for PTPR92_21
#
#                         'GRal04_10a': {'elem': 'He',
#                                      'atom': 'abun["He2"]',
#                                      'type': 'HII',
#                                      'icf': '1.05',
#                                      'comment': 'Based on the ionization potential. To be used if t2 = 0.00'},
#                         'GRal04_10b': {'elem': 'He',
#                                      'atom': 'abun["He2"]',
#                                      'icf': '1.04',
#                                      'type': 'HII',
#                                      'comment': 'Based on the ionization potential. To be used if t2 > 0.00'},
                         'PTPR92_21': {'elem': 'He',
                                      'atom': 'abun["He2"]',
                                      'icf': '(1 + abun["S2"] / abun["S3"])',
                                      'type': 'HII',
                                      'comment': 'Based on the ionization potential, assumoing that S = S2 + S3'},
                         'PC69_40': {'elem': 'Ne',
                                    'atom': 'abun["Ne3"]',
                                   'icf': '(abun["O2"] + abun["O3"]) / abun["O3"]',
                                   'type': 'HII',
                                   'comment': 'High ionization degree'},
                         'S78_265b': {'elem': 'Ne',
                                      'atom': 'abun["Ne3"]',
                                      'icf': '(abun["O2"] + abun["O3"]) / (abun["O3"] - 0.2 * abun["O2"])',
                                      'type': 'HII',
                                      'comment': 'Based on photoionization models'},
# Superseded by RR05
#                         'RR04_1': {'elem': 'Fe',
#                                      'atom': 'abun["Fe3"]',
#                                      'icf': '(abun["O2"] / abun["O3"])**0.09 * (abun["O2"] + abun["O3"]) / abun["O2"]',
#                                      'type': 'HII',
#                                      'comment': 'Based on photoionization models and the assumption that O/H = (O+ + O++) / H+'},
                         'RR05_2': {'elem': 'Fe',
                                      'atom': 'abun["Fe3"]',
                                      'icf': '0.9 * (abun["O2"] / abun["O3"])**0.08 * (abun["O2"] + abun["O3"]) / abun["O2"]',
                                      'type': 'HII',
                                      'comment': 'Based on fit to data and the assumption that O/H = (O+ + O++) / H+'},
                         'RR05_3': {'elem': 'Fe',
                                      'atom': 'abun["Fe3"]',
                                      'icf': '1.1 * (abun["O2"] / abun["O3"])**0.58 * (abun["O2"] + abun["O3"]) / abun["O2"]',
                                      'type': 'All',
                                      'comment': 'Based on PNe and HII data with -1.35 < log(O2/O3) < -0.1. It is assumed that O/H = (O+ + O++) / H+'},
                         'RR05_4': {'elem': 'Fe',
                                      'atom': 'abun["Fe2"] + abun["Fe3"]',
                                      'icf': '(abun["O2"] + abun["O3"]) / abun["O2"]',
                                      'type': 'All',
                                      'comment': 'Based on PNe and HII data with log(O2/O3) > -0.1. It is assumed that O/H = (O+ + O++) / H+'},
                         'GKA07_1.p269': {'elem': 'Cl',
                                      'atom': 'abun["Cl3"]',
                                      'icf': '(elem_abun["direct_He.23"] / abun["He2"])**2',
                                      'type': 'PNe',
                                      'comment': 'Based on photoionization models'},
                         'mGKA07-PTPR92_p269': {'elem': 'Cl',
                                      'atom': 'abun["Cl3"]',
                                      'icf': '(elem_abun["PTPR92_21"] / abun["He2"])**2',
                                      'type': 'PNe',
                                      'comment': 'Based on photoionization models'},
                         'DIMS14_10': {'elem': 'He',
                                       'atom':'abun["He2"] + abun["He3"]', 
                                       'icf': '1.0',
                                       'type': 'PNe',
                                       'comment': 'Based on photoionization models'},
                         'DIMS14_12': {'elem': 'O',
                                       'atom':'abun["O2"] + abun["O3"]', 
                                       'icf': '10**((0.08 * abun["He3"] / (abun["He2"] + abun["He3"]) + 0.006 * (abun["He3"] / (abun["He2"] + abun["He3"]))**2) / (0.34 - 0.27 * abun["He3"] / (abun["He2"] + abun["He3"])))',
                                       'type': 'PNe',
                                       'comment': 'Based on photoionization models',
                                       'low_error': '0.03 + 0.5 * abun["He3"] / (abun["He2"] + abun["He3"]) - 0.2 * (abun["He3"] / (abun["He2"] + abun["He3"]))**2',
                                       'up_error': '0.03 + 0.5 * abun["He3"] / (abun["He2"] + abun["He3"]) - 0.2 * (abun["He3"] / (abun["He2"] + abun["He3"]))**2'},
                         'DIMS14_14': {'elem': 'N',
                                       'atom':'abun["N2"]', 
                                       'icf': '10**(-0.16 * abun["O3"] / (abun["O2"] + abun["O3"]) * (1.0 + np.log10(abun["He3"] / (abun["He2"] + abun["He3"])))) * elem_abun["DIMS14_12"] / abun["O2"]',
                                       'type': 'PNe',
                                       'comment': 'Based on photoionization models',
                                       'low_error': '0.32 * abun["O3"] / (abun["O2"] + abun["O3"])',
                                       'up_error': '0.50 * abun["O3"] / (abun["O2"] + abun["O3"])'},
                         'DIMS14_14b': {'elem': 'N',
                                       'atom':'abun["N2"]', 
                                       'icf': '10**(0.64 * abun["O3"] / (abun["O2"] + abun["O3"])) * elem_abun["DIMS14_12"] / abun["O2"]',
                                       'type': 'PNe',
                                       'comment': 'Based on photoionization models',
                                       'low_error': '0.32 * abun["O3"] / (abun["O2"] + abun["O3"])',
                                       'up_error': '0.50 * abun["O3"] / (abun["O2"] + abun["O3"])'},
                         'DIMS14_17a': {'elem': 'Ne',
                                       'atom':'abun["Ne3"]', 
                                       'icf': '(abun["O3"] / (abun["O2"] + abun["O3"]) + (0.014 / (abun["He3"] / (abun["He2"] + abun["He3"])) + 2 * (abun["He3"] / (abun["He2"] + abun["He3"]))**2.7)**3 * (0.7 + 0.2 * abun["O3"] / (abun["O2"] + abun["O3"]) - 0.8 * (abun["O3"] / (abun["O2"] + abun["O3"]))**2)) * elem_abun["DIMS14_12"] / abun["O3"]',
                                       'type': 'PNe',
                                       'comment': 'Based on photoionization models',
                                       'low_error': '0.17',
                                       'up_error': '0.12'},
                         'DIMS14_17b': {'elem': 'Ne',
                                       'atom':'abun["Ne3"]', 
                                       'icf': '(abun["O3"] / (abun["O2"] + abun["O3"]) + (0.014 / (0.01) + 2 * (0.01)**2.7)**3 * (0.7 + 0.2 * abun["O3"] / (abun["O2"] + abun["O3"]) - 0.8 * (abun["O3"] / (abun["O2"] + abun["O3"]))**2)) * elem_abun["DIMS14_12"] / abun["O3"]',
                                       'type': 'PNe',
                                       'comment': 'Based on photoionization models',
                                       'low_error': '0.17',
                                       'up_error': '0.12'},
                         'DIMS14_17c': {'elem': 'Ne',
                                       'atom':'abun["Ne3"]', 
                                       'icf': '(abun["O3"] / (abun["O2"] + abun["O3"]) + (0.014 / (0.015) + 2 * (0.015)**2.7)**3 * (0.7 + 0.2 * abun["O3"] / (abun["O2"] + abun["O3"]) - 0.8 * (abun["O3"] / (abun["O2"] + abun["O3"]))**2)) * elem_abun["DIMS14_12"] / abun["O3"]',
                                       'type': 'PNe',
                                       'comment': 'Based on photoionization models',
                                       'low_error': '0.17',
                                       'up_error': '0.12'},
                         'DIMS14_20': {'elem': 'Ne',
                                       'atom':'abun["Ne3"]+abun["Ne5"]', 
                                       'icf': '(1.31+12.68*(abun["He3"]/(abun["He2"]+abun["He3"]))**2.57)**0.27',
                                       'type': 'PNe',
                                       'comment': 'Based on photoionization models',
                                       'low_error': '0.20',
                                       'up_error': '0.17'},
                         'DIMS14_23': {'elem': 'S',
                                       'atom':'abun["S2"]', 
                                       'icf': '(10**(0.31 - 0.52 * abun["He3"] / (abun["He2"] + abun["He3"]))) * elem_abun["DIMS14_12"] / abun["O2"]',
                                       'type': 'PNe',
                                       'comment': 'Based on photoionization models',
                                       'low_error': '0.38',
                                       'up_error': '0.41'},
                         'DIMS14_26': {'elem': 'S',
                                       'atom':'abun["S2"] + abun["S3"]', 
                                       'icf': '(10**((-0.02 - 0.03*abun["O3"] / (abun["O2"] + abun["O3"]) - 2.31*(abun["O3"] / (abun["O2"] + abun["O3"]))**2 + 2.19 * (abun["O3"] / (abun["O2"] + abun["O3"]))**3) / (0.69 + 2.09 * abun["O3"] / (abun["O2"] + abun["O3"]) - 2.69 * (abun["O3"] / (abun["O2"] + abun["O3"]))**2))) * elem_abun["DIMS14_12"] / abun["O2"]',
                                       'type': 'PNe',
                                       'comment': 'Based on photoionization models',
                                       'low_error': '0.20',
                                       'up_error': '0.12'},
                         'DIMS14_29': {'elem': 'Cl',
                                       'atom':'abun["Cl3"]', 
                                       'icf': '((4.162 - 4.1622*(abun["O3"] / (abun["O2"] + abun["O3"]))**0.21)**.75) * elem_abun["DIMS14_12"] / abun["O2"]',
                                       'type': 'PNe',
                                       'comment': 'Based on photoionization models',
                                       'low_error': '0.27',
                                       'up_error': '0.15'},
                         'DIMS14_29b': {'elem': 'Cl',
                                       'atom':'abun["Cl2"] + abun["Cl3"]', 
                                       'icf': '(1.) * elem_abun["DIMS14_12"] / abun["O2"]',
                                       'type': 'PNe',
                                       'comment': 'Based on photoionization models',
                                       'low_error': '0.13',
                                       'up_error': '0.0'},
                         'DIMS14_32': {'elem': 'Cl',
                                       'atom':'abun["Cl2"] + abun["Cl3"] + abun["Cl4"]', 
                                       'icf': '0.98 + (0.56 - 0.57 * abun["He3"] / (abun["He2"] + abun["He3"]))**7.64',
                                       'type': 'PNe',
                                       'comment': 'Based on photoionization models',
                                       'low_error': '0.30',
                                       'up_error': '0.33'},
                         'DIMS14_35': {'elem': 'Ar',
                                       'atom':'abun["Ar3"]', 
                                       'icf': '(10**(0.05 / (0.06 + abun["O3"] / (abun["O2"] + abun["O3"])) - 0.07)) * elem_abun["DIMS14_12"] / (abun["O2"] + abun["O3"])',
                                       'type': 'PNe',
                                       'comment': 'Based on photoionization models',
                                       'low_error': '0.70',
                                       'up_error': '0.60'},
                         'DIMS14_36': {'elem': 'Ar',
                                       'atom':'abun["Ar3"]', 
                                       'icf': '(10**((0.03  * abun["O3"] / (abun["O2"] + abun["O3"])) / (0.4 - 0.3 * abun["O3"] / (abun["O2"] + abun["O3"])) - 0.05)) * elem_abun["DIMS14_12"] / (abun["O2"] + abun["O3"])',
                                       'type': 'PNe',
                                       'comment': 'Based on photoionization models',
                                       'low_error': '0.70',
                                       'up_error': '0.60'},
                         'DIMS14_39': {'elem': 'C',
                                       'atom':'abun["C3"]', 
                                       'icf': '(0.05 + 2.21 * abun["O3"] / (abun["O2"] + abun["O3"]) - 2.77 * (abun["O3"] / (abun["O2"] + abun["O3"]))**2 + 1.74 * (abun["O3"] / (abun["O2"] + abun["O3"]))**3) * elem_abun["DIMS14_12"] / abun["O3"]',
                                       'type': 'PNe',
                                       'comment': 'Based on photoionization models',
                                       'low_error': '0.19',
                                       'up_error': '(0.4 - 1.06 * abun["O3"] / (abun["O2"] + abun["O3"]) + 0.65 * (abun["O3"] / (abun["O2"] + abun["O3"]))**2 + 0.27 * (abun["O3"] / (abun["O2"] + abun["O3"]))**3) '},
# To be completed
#                         'Sal07_Kr1': {'elem': 'Kr',
#                                       'atom':'abun["Kr3"] + abun["Kr4"]', 
#                                       'icf': 'elem_abun["DIMS14_12"]/(abun["Ar3"] + abun["ar4"])',
#                                       'type': 'PNe',
#                                       'comment': 'Based on photoionization models'},}
                         
                         }
        
        Ial06 = [('18a', 'N', 'N2', (-0.825, 0.718, 0.853), 'v', 'low Z'),
                 ('18b', 'N', 'N2', (-0.809, 0.712, 0.852), 'v', 'int Z'),
                 ('18c', 'N', 'N2', (-1.476, 1.752, 0.688), 'v', 'high Z'),
                 ('19a', 'Ne', 'Ne3', (-0.385, 1.365, 0.022), 'w', 'low Z'),
                 ('19b', 'Ne', 'Ne3', (-0.405, 1.382, 0.021), 'w', 'int Z'),
                 ('19c', 'Ne', 'Ne3', (-0.591, 0.927, 0.546), 'w', 'high Z'),
                 ('21a', 'Cl', 'Cl3', (0.756, 0.648, 0.022), 'v', 'low Z'),
                 ('21b', 'Cl', 'Cl3', (0.814, 0.620, 0.128), 'v', 'int Z'),
                 ('21c', 'Cl', 'Cl3', (1.186, 0.357, 0.131), 'v', 'high Z'),
                 ('22a', 'Ar', 'Ar3', (0.278, 0.836, 0.121), 'v', 'low Z'),
                 ('22b', 'Ar', 'Ar3', (0.285, 0.833, 0.051), 'v', 'int Z'),
                 ('22c', 'Ar', 'Ar3', (0.517, 0.763, 0.042), 'v', 'high Z'),
                 ('24a', 'Fe', 'Fe3', (0.036, -0.146, 1.386), 'v', 'low Z'),
                 ('24b', 'Fe', 'Fe3', (0.301, -0.259, 1.367), 'v', 'int Z'),
                 ('24c', 'Fe', 'Fe3', (-1.377, 1.606, 1.045), 'v', 'high Z'),
                 ('20a', 'S', 'S2+S3', (0.121, 0.511, 0.161), 'v', 'low Z'),
                 ('20b', 'S', 'S2+S3', (0.155, 0.849, 0.062), 'v', 'int Z'),
                 ('20c', 'S', 'S2+S3', (0.178, 0.610, 0.153), 'v', 'high Z'),
                 ('23a', 'Ar', 'Ar3+Ar4', (0.158, 0.958, 0.004), 'v', 'low Z'),
                 ('23b', 'Ar', 'Ar3+Ar4', (0.104, 0.980, 0.001), 'v', 'int Z'),
                 ('23c', 'Ar', 'Ar3+Ar4', (0.238, 0.931, 0.004), 'v', 'high Z')]
        
        for item in Ial06:
            key = 'Ial06_' + item[0]
            atom = 'abun["' + item[2] + '"]'
            icf = self._calcIal06(item[3], item[4])
            self.all_icfs[key] = {'elem': item[1], 'atom': atom, 'icf': icf, 'type': 'HII', 'comment': item[5]}

        # Quickest way to account for the special case of Eqs. 20 and 23
        for key in ['Ial06_20a', 'Ial06_20b', 'Ial06_20c']:
            self.all_icfs[key]['atom'] = '(abun["S2"] + abun["S3"])'
        for key in ['Ial06_23a', 'Ial06_23b', 'Ial06_23c']:
            self.all_icfs[key]['atom'] = '(abun["Ar3"] + abun["Ar4"])'
        

        self.all_icf_refs = {'direct':{'ref': 'Direct determination by summing observed ions',
                                       'url':''},
                             'TPP77': {'ref': 'Torres-Peimbert and Peimbert 1977, RMAA, 2, 181',
                                       'url': 'http://adsabs.harvard.edu/abs/1977RMxAA...2..181T'},
                             'PHCD07': {'ref': 'Perez-Montero, Haegele, Contini, and Diaz 2007, MNRAS, 381, 125',
                                        'url': 'http://adsabs.harvard.edu/abs/2007MNRAS.381..125P'},
                             'KB94': {'ref': 'Kinsburgh & Barlow 1994, MNRAS, 271, 257',
                                      'url': 'http://adsabs.harvard.edu/abs/1994MNRAS.271..257K'},
                             'KH01': {'ref': 'Kwitter & Henry 2001, ApJ, 562, 804',
                                      'url': 'http://adsabs.harvard.edu/abs/2001ApJ...562..804K'},
                             'ITL94': {'ref': 'Izotov, Thuan & Lipovetsky 1994, ApJ, 435, 647',
                                       'url': 'http://http://adsabs.harvard.edu/abs/1994ApJ...435..647I'},
                             'Ial06': {'ref': 'Izotov et al 2006, A&A, 448, 955',
                                       'url': 'http://adsabs.harvard.edu/abs/2006A%26A...448..955I'},
                             'PC69': {'ref': 'Peimbert & Costero 1969, BOTT, 5, 3',
                                      'url': 'http://adsabs.harvard.edu/abs/1969BOTT....5....3P'},
                             'S78': {'ref': 'Stasinska 1978, A&A, 66, 257',
                                     'url': 'http://adsabs.harvard.edu/abs/1978A%26A....66..257S'},
# Superseded by RR05
#                             'RR04': {'ref': 'Rodriguez & Rubin 2004, IAUS, 217, 188',
#                                      'url': 'http://adsabs.harvard.edu/abs/2004IAUS..217..188R'},
                             'RR05': {'ref': 'Rodriguez & Rubin 2005, ApJ, 626, 900',
                                      'url': 'http://adsabs.harvard.edu/abs/2005ApJ...626..900R'},
# Corrected with PTPR92
#                             'GRal04': {'ref': 'Garcia-Rojas et al 2004, ApJS, 153, 501',
#                                      'url': 'http://adsabs.harvard.edu/abs/2004ApJS..153..501G'},
                             'GKA07': {'ref': 'Girard, Koeppen and Acker 2007, A&A, 463, 265',
                                      'url': 'http://adsabs.harvard.edu/abs/2007A%26A...463..265G'},
                             'PTPR92': {'ref': 'Peimbert, Torres-Peimbert and Ruiz 1992, RMAA, 24, 155',
                                        'url': 'http://adsabs.harvard.edu/abs/1992RMxAA..24..155P'},
                             'Sal07': {'ref': 'Sharpee et al 2007, MNRAS, 659, 1265',
                                        'url': 'http://adsabs.harvard.edu/abs/2007ApJ...659.1265S'},
                             'DIMS14': {'ref': 'Delgado-Inglada, Morisset and Stasinska, 2014, MNRAS',
                                        'url':'http://arxiv.org/abs/1402.4852'}
                             }

    def _calcIal06(self, coeff, degree):
        """
        Compute the analytical expressions of the Ial06 paper
        
        """
        if degree == 'v':
            k1 = 'abun["O2"] / (abun["O2"] + abun["O3"])'
            k2 = '(abun["O2"] + abun["O3"]) / abun["O2"]'
        else:
            k1 = 'abun["O3"] / (abun["O2"] + abun["O3"])'
            k2 = '(abun["O2"] + abun["O3"]) / abun["O3"]'
        return '(' + str(coeff[0]) + ' * ' + k1 + ' + ' + str(coeff[1]) + ' + ' + str(coeff[2]) + ' * ' + k2 + ')'


    def _getAtomExp(self, label):
        """
        Clean the ionic factor expressions for display purposes
        
        """
        return self.all_icfs[label]['atom'].replace('abun["', '').replace('"]', '')


    def _getIcfExp(self, label):
        """
        Clean the icf expressions for display purposes

        """        
        return self.all_icfs[label]['icf'].replace('elem_abun["', '').replace('abun["', '').replace('"]', '')

    
    def getAvailableICFs(self, elem_list=ELEM_LIST, type_=['HII', 'PNe', 'All']):
        """ 
        Get a list of all the available ICFs for the specified elements. 
        Details can be obtained by invoking getReference('label') and getExpression('label') or printAllICFs()
        
        Parameters:
            - elem_list    selected element or list of selected elements (default: all PyNeb elements)
            - type_         dbject class to which the icf is appliable (e.g. "PNe", "HII"; default: both)
        
        Usage:
        
        icf=pn.ICF()
        icf.getAvailableICFs()
        icf.getAvailableICFs('S')
        icf.getAvailableICFs('S', type_='HII')

        """        
        if type(type_) == type(''):
            type_ = [type_] 
        if elem_list.__class__ is str:
            elem_list = [elem_list]
        icf_dict = {}
        for elem in elem_list:
            # Each element is associated to an icf_list which holds all the icfs for that element 
            icf_list = []
            for key in self.all_icfs.keys():
                if (self.all_icfs[key]['elem'] == elem):
                    if self.all_icfs[key]['type'] in type_:
                        icf_list.append(key)
            if len(icf_list) > 0:
                icf_dict[elem] = icf_list        
        return icf_dict
   
    
    def printAllICFs(self, type_=['HII', 'PNe', 'All']):
        """ 
        Print a list of all the available ICFs. Details can be obtained by 
            invoking getReference('label') and getExpression('label')
        
        Parameters:
            - type_    object class to which the icf is appliable (e.g. "PNe", "HII"; default: both)

        """
        if type(type_) == type(''):
            type_ = [type_] 
        for label in self.all_icfs:
            if self.all_icfs[label]['type'] in type_:
                print(label + ': elem = ' + self.all_icfs[label]['elem'] + '; atom = ' + 
                      self._getAtomExp(label) + '; type = ' + self.all_icfs[label]['type'])
            
   
    def addICF(self, label, elem, atom, icf, type_='All', comment='', ref=None, url=None):
        """
        Add a new icf
        
        Parameters:
            - label    a label which identifies the new ICF. The suggested format for the label is "R_E.N", 
                        where R is a reference to the paper, E to the equation in the paper and the optional 
                        key N is an indication of the higher ion need to compute the icf.  
            - elem     element whose abundance is computed
            - atom     ion or ions whose abundance is multiplied by the icf to get the elemenmt abundance
            - icf      the correcting expression
            - type_     object class to which the icf is appliable (e.g. "PNe", "HII")
            - comment  additional comment, relevant to icf usage
            - ref      bibliographic reference, if any
            - url      URL of the source paper, if any
            
        Usage:
            icf=pn.ICF()
            icf.addICF(label='TEST_1',
                        elem='Ne',
                        atom ='abun["Ne2"] + abun["Ne3"]',
                        icf='(abun["O2"] + abun["O3"]) / abun["O2"]',
                        type_='PNe',
                        comment='',
                        ref='PyNeb international conference, 2025',
                        url='http://www.iac.es')

        """
        if label in self.all_icfs:
            pn.log_.warn('{0} is already an entry of the ICF label dictionary'.format(label), calling=self.calling)
            return None
        else:
            self.all_icfs[label] = {'elem': elem,
                                    'atom': atom,
                                    'icf': icf,
                                    'type': type_,
                                    'comment': comment}
            paper = label.split('_')[0]
            self.all_icf_refs[paper] = {'ref': ref, 'url': url}
    
    
    def delICF(self, label):
        """
        Delete an existing icf (only for current session)
        
        Parameters:
            - label    the label which identifies the new ICF. 
            
        Usage:
            icf=pn.ICF()
            icf.delICF(label='TEST_1')

        """
        del self.all_icfs[label]
    
    
    def printInfo(self, label, filter='licores'): 
        """ 
        Return complete information on the required icf or paper
        
        Parameters:
            label    label of selected ICF recipe or paper
            filter   string with initials of required fields (substring of 'licores'; all by default)   
            
        Usage:
        icf=pn.ICF()
        icf.printInfo('Ial06_20a')
                 
        """        
        for item in self.all_icfs:
            if (label in item):
                print('------------------------------------------------------------------------')
                if ('l' in filter): print('Label:', item)
                if ('e' in filter): print('Element:', self.all_icfs[item]['elem'])
                if ('r' in filter): print('Required ions:', self._getAtomExp(item))
                if ('i' in filter): print('ICF expression:', self._getIcfExp(item))
                if ('o' in filter): print('Object class:', self.all_icfs[item]['type'])
                if ('c' in filter): print('Comments:', self.all_icfs[item]['comment'])
                if ('s' in filter): print('Source:', self.getReference(item))
   
    
    def getReference(self, label): 
        """ 
        Return the reference of the selected ICF recipe 
        
        Parameters:
            - label    label of selected ICF recipe or paper

        """        
        paper = label.split('_')[0]
        if paper[0] == 'm':
            paper1 = paper.split('-')[0][1:]
            paper2 = paper.split('-')[1]
            mixed = self.all_icf_refs[paper1]['ref'] + ' using an abundance based on ' + self.all_icf_refs[paper2]['ref']
            return mixed
        else:
            return self.all_icf_refs[paper]['ref']
   
    
    def getURL(self, label): 
        """ 
        Return the ADS URL of the selected ICF recipe 
        
        Parameters:
            label    label of selected ICF recipe or paper

        """        
        paper = label.split('_')[0]
        if paper[0] == 'm':
            paper1 = paper.split('-')[0][1:]
            paper2 = paper.split('-')[1]
            mixed = self.all_icf_refs[paper1]['url'] + ' using an abundance based on ' + self.all_icf_refs[paper2]['url']
            return mixed
        else:
            return self.all_icf_refs[paper]['url']
   
    
    def getExpression(self, label): 
        """ 
        Return the analytical expression of the selected ICF 
        
        Parameters:
            label    label of selected ICF expression

        """
        
        elem = self.all_icfs[label]['elem']
        atom = self._getAtomExp(label)
        icf_factor = self._getIcfExp(label)
        if icf_factor != '1':
            return elem + " = (" + atom + ") * " + icf_factor
        else:
            return elem + " = " + atom
     
             
    def getType(self, label): 
        """ 
        Return the kind of object for which the selected ICF is suitable. 
        
        Parameters:
            label    label of selected ICF expression

        """
        return self.all_icfs[label]['type']
        

    def getComment(self, label): 
        """ 
        Return the comment associated to the selected ICF. 
        
        Parameters:
            label    label of selected ICF expression

        """
        return self.all_icfs[label]['comment']

    def getElemAbundance(self, atom_abun, icf_list=[], icf_family=None, absentIon=np.nan, use_coll=True,
                         use_MC=False):
        """
        Compute elemental abundances from ionic abundances. The complete iventory is printed through pn.ICF().printAllICFs().
        See the ICF class for more details.
        Store the result in ICF.elem_abund
        
        Parameters: 
           - atom_abun    a dictionary of ionic abundances
           - icf_list     a list of selected ICFs (default: all)
           - icf_family   a string common to all the icf one want to use, e.g. "DIMS14"
           - use_coll     In case both collision lines and recombination lines abundances exists:
                           if True (default), use collision line ionic abundance, if false use recombination ones
           - use_MC       If True and if the ICF definition includes error informations, the ICF will be multiplied
                           by some factor taking into account the error.
        Usage:
            atom_abun = {'O2': 0.001, 'O3': 0.002, 'Ne3': 1.2e-5}
            getElemAbundance(atom_abun)
         
        """
        if type(icf_list) == type(''):
            icf_list = [icf_list]
        # List of all existing atoms. Necessary to initialize the ionic abundance dictionary
        if icf_family is not None:
            icf_list = [icf for icf in self.all_icfs.keys() if icf_family in icf]
        atom_list = []
        for elem in ELEM_LIST:
            for spec in SPEC_LIST:
                atom_list.append(elem + str(spec)) 
        atom_list.extend(['He2', 'He3', 'He1r', 'He2r', 'O2r'])
    
        # Initialize the ionic abundances so that the code does not crash when a specific abundance is invoked
        # TODO : Check that the same ion is not present as collisional and recombination !!
        abun = {}
        for atom in atom_list:
            # Set those abundances which are different from 0. These determine which ICFs can be computed        
            if atom in atom_abun:
                if atom[-1] == 'r':
                    elem, spec = parseAtom(atom[:-1])
                    atomc = elem + str(int(spec)+1)
                    if (not use_coll) or (atomc not in atom_abun):
                        abun[atomc] = atom_abun[atom]
                else:
                    elem, spec = parseAtom(atom)
                    atomr = elem + str(int(spec)-1) + 'r'
                    if use_coll or (atomr not in atom_abun):
                        abun[atom] = atom_abun[atom]
            else:
                abun[atom] = absentIon
        self.abun = abun
        # Initialize the lists of elemental abundances
        self.icf_value = {}
        elem_abun = {}
        # Compute either all the available ICFs or just a selected subset of them
        if len(icf_list) == 0:
            icf_list = self.all_icfs.keys()

        do_one_more_pass = True
        i_pass = 0
        while (do_one_more_pass and (i_pass < self._max_pass)):
            do_one_more_pass = False
            i_pass += 1
            for icf_label in icf_list:
                atom = self.all_icfs[icf_label]['atom']
                try:
                    atom_value = eval(atom)
                except:
                    pn.log_.warn('Unable to eval {}'.format(atom))
                if self.all_icfs[icf_label]['elem'] is not None:
                    pn.log_.debug('Doing {}'.format(icf_label), calling='ICF')
                    this_icf = self.all_icfs[icf_label]
                    elem = this_icf['elem']
                    try:
                        icf_value = eval(this_icf['icf'])
                        self.icf_value[icf_label] = icf_value
                        if use_MC and 'low_error' in this_icf:
                            try:
                                low_error = eval(this_icf['low_error'])
                            except:
                                pn.log_.warn('low error {0} cannot be evaluated for {1}'.format(this_icf['low_error'], icf_label))
                            try:
                                up_error = eval(this_icf['up_error'])
                            except:
                                pn.log_.warn('up error {0} cannot be evaluated for {1}'.format(this_icf['up_error'], icf_label))
                            try:
                                if np.size(low_error) == 1:
                                    D_icf_MC = np.random.uniform(-low_error, up_error, np.size(atom_value))
                                else:
                                    D_icf_MC = np.random.uniform(-low_error, up_error)
                                D_icf_MC[0] = 0.0
                            except:
                                pn.log_.warn('MC not working for ICF {}'.format(this_icf['icf']), calling=self.calling)
#                             print '--------IN---------'
#                             print this_icf['icf']
#                             print this_icf['low_error']
#                             print 'low_error'
#                             print low_error
#                             print this_icf['up_error']
#                             print 'up_error'
#                             print up_error
#                             if icf_label == 'DIMS14_36':
#                                 print 'icf_MC for {}'.format(icf_label)
#                                 print Dicf_MC
#                                 print 'icf_value'                            
#                                 print icf_value
                            icf_value = icf_value + D_icf_MC
                            pn.log_.message('Using MC for ICF {}'.format(icf_label))
#                             print 'icf_value2'                            
#                             print icf_value
#                             print '---------OUT---------'
                        elem_abun[icf_label] = atom_value * icf_value
                    except:
                        do_one_more_pass = True
                        if i_pass == self._max_pass:
                            pn.log_.warn('{0} cannot be evaluated for {1}'.format(this_icf['icf'], icf_label),
                                         calling=self.calling)
                else:
                    pn.log_.debug('{} not present'.format(icf_label), calling='ICF')
        self.elem_abun = elem_abun
        return self.elem_abun
#            elem_abun[icf.all_icfs[icf_ref]['elem']] = np.log10(elem_abun_dic[icf_ref])+12

    def __getElemAbundanceFromStategy(self, atom_abun, strategy):
        """
        Return the elemental abundances in a dictionary, computed by applying a given strategy.
        
        Parameters:
            atom_abund: dictionary of ionic abundances
            strategy: dictionary of the ICF to apply in order to obtain elemental abundances.
                example: strategy = {'He':['direct_He.2', 'PTPR92_21'], 
                                   'N': 'KB94_A1.10'}
            final_abund = icf.getElemAbundanceFromStrategy(ion_ab_dic, icf_strategies)
        """
        
        # the following is great, as it only computes the abundances from a subset, but
        # it does not work when an ICF requires an ab undance not in the list!
        icf_list = []
        for elem in strategy:
            if type(strategy[elem]) == type(''):
                icf_list.append(strategy[elem])
            else:
                for meth in strategy[elem]:
                    icf_list.append(meth)
                
        allAbund = self.getElemAbundance(atom_abun, icf_list, absentIon=np.nan)
        elemAbund = {}
        for elem in strategy:
            if type(strategy[elem]) == type(''):
                elemAbund[elem] = allAbund[strategy[elem]]
            else:
                elemAbund[elem] = allAbund[strategy[elem][0]]
                for meth in strategy[elem][1::]:
                    tt = np.isnan(elemAbund[elem])
                    pn.log_.message('Computing {0} abundance using {1} for {2} observations.'.format(elem,
                                                                                                     meth,
                                                                                                     tt.sum()),
                                    calling=self.calling)
                    elemAbund[elem][tt] = allAbund[meth][tt]
        return elemAbund


    def test(self):
        """ 
        Test function for the ICF class 

        """
        label1 = 'Ial06_18a'
        label2 = 'TPP77_14'
        label3 = 'KH01_4b'
        label4 = 'PHCD07_13'
        print('All the available ICFs, listed by element:')
        print(self.getAvailableICFs())
        print('\nAll the available ICFs, with some detail:')
        print(self.printAllICFs())
        print('\nAll the available ICFs for S and Ne:' )  
        print(self.getAvailableICFs(['S', 'Ne']))
        print('\nAnalytical expression for ' + label1 + ':')
        print(self.getExpression(label1))
        print('\nBibliographic reference for ' + label2 + ':')
        print(self.getReference(label2))
        print('\nURL for ' + label3 + ':')
        print(self.getURL(label3))
        print('\nType of object to which ' + label4 + ' can be applied:')
        print(self.getType(label4))
        print('\nAdditional details of ' + label4 + ':')
        print(self.getComment(label4) )

        atom_abun = {}
        atom_abun['O2'] = 0.001
        atom_abun['O3'] = 0.002
        atom_abun['Ne2'] = 1.0e-4
        atom_abun['Ne3'] = 1.2e-5
        atom_abun['Ar3'] = 4.e-6
        atom_abun['Ar4'] = 1.e-6
        atom_abun['N2'] = 1.e-4
        atom_abun['He2'] = 1.e-1
        atom_abun['He3'] = 1.e-2
        atom_abun['He2'] = 1.e-1
        atom_abun['Cl3'] = 4.e-6
        atom_abun['Cl4'] = 1.e-6
        atom_abun['Ar3'] = 6.e-5
        atom_abun['Ar4'] = 5.e-7
        atom_abun['S3'] = 1.e-5

        self.getElemAbundance(atom_abun, icf_list=[label4])

