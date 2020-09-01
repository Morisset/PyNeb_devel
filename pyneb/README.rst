Introduction
============

PyNeb (by V. Luridiana, C. Morisset, and R. A. Shaw) is the last in a lineage of tools dedicated to the analysis of emission lines, which includes FIVEL and nebular (see also here). PyNeb is completely written in python and is designed to be more user-friendly and powerful than its predecessors, its functionality representing a giant leap forward with respect to them in terms of speed, easiness of use, graphic visualization, and accessible information. 

In PyNeb, the atom is represented as an n-level atom. For given density and temperature, PyNeb's machinery solves the equilibrium equations and determines the level populations. These are some of the things it can do:

* compute physical conditions from suitable diagnostic line ratios;
* compute level populations, critical densities and line emissivities;
* compute and display emissivity grids as a function of Te and Ne;
* deredden line intensities;
* read and manage observational data;
* plot and compare atomic data from different publications;
* compute ionic abundances from line intensities and physical conditions;
* compute elemental abundances from ionic abundances and icfs;
* compute electron temperature from Balmer/Pashen jump, at any wavelength and normalized to any HI line.

PyNeb also incorporates emissivity tables or emissivities fits of recombination lines. 
The interpolated emissivities can be used by any of the modules that rely on the n-level line emissivities to carry out the actions listed above.

PyNeb uses part of Chiantipy:

Utility functions, many for reading the CHIANTI database files.

Copyright 2009, 2010 Kenneth P. Dere

This software is distributed under the terms of the GNU General Public License
that is found in the LICENSE file

and
 
FortranFormat  

Written by Konrad Hinsen <hinsen@cnrs-orleans.fr>
With contributions from Andreas Prlic <andreas@came.sbg.ac.at>
last revision: 2006-6-23
