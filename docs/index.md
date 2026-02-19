# PyNeb's Home Page

## Introduction

PyNeb (Luridiana V., Morisset C. and Shaw, R. A 2013) is a modern python tool to compute emission line emissivities (recombination and collisionally excited lines).

In PyNeb, the atom is represented as an n-level atom. For given density and temperature, PyNeb's machinery solves the equilibrium equations and determines the level populations. These are some of the things it can do:

- compute physical conditions from suitable diagnostic line ratios.
- compute level populations, critical densities and line emissivities 
- compute and display emissivity grids as a function of Te and Ne
- deredden line intensities
- read and manage observational data
- plot and compare atomic data from different publications
- compute ionic abundances from line intensities and physical conditions
- compute elemental abundances from ionic abundances and icfs.

PyNeb also incorporates emissivity tables of recombination lines for a few atoms. The interpolated emissivities can be used by any of the module that rely on the n-level line emissivities to carry out the actions listed above.

Emission line ratios are used to self consistently determine electron temperature and density and ionic abundances  
Diagnostic diagrams can easily be plotted.  
Various ionization correction factors (ICFs) from the literarure are available to obtain total elemental abundances from the ionic abundances.  
Atomic data can easily be changed and updated.  
Additional tools are provided, like reddening determination and correction procedures, Balmer/Pashen jump  
temperature determinations.

## Citation

If you use PyNeb in your research, please cite the following paper:

- Luridiana, V., Morisset, C. and Shaw, R. A. 2013, A&A, 558, A57  
  http://adsabs.harvard.edu/abs/2015A%26A...573A..42L

- Morisset, C., Luridiana, V., García-Rojas, J., Gómez-Llanos, V., Bautista, M., & Mendoza, C. 2020, Atoms, 8, 66,  
  «Atomic Data Assessment with PyNeb»  
  https://ui.adsabs.harvard.edu/abs/2020Atoms...8...66M

- Mendoza, C., Méndez-Delgado, J. E., Bautista, M., García-Rojas, J., & Morisset, C. 2023, Atoms, 11, 63,  
  «Atomic Data Assessment with PyNeb: Radiative and Electron Impact Excitation Rates for [Fe ii] and [Fe iii]»  
  https://ui.adsabs.harvard.edu/abs/2023Atoms..11...63M

## Requirements

PyNeb uses numpy, matplotlib, pyfits, scipy and other standard python libraries.

## Installation

You may find useful to download, install and upgrade PyNeb using [pip](http://www.pip-installer.org/en/latest/index.html).

For example:

- `pip install -U PyNeb`

Note: you MAY need `--user` if you installed python without Anaconda or Canopy.

Updates use the same command.

You can also install from the github repository:

- `pip install -U git+https://github.com/Morisset/PyNeb_devel.git`

To use the development branch (at your own risks!!!):

- `pip install -U git+https://github.com/Morisset/PyNeb_devel.git@devel`

## Warranty

PyNeb is provided as it is. No warranty at all.

## Manual

- The manuals are here: <https://github.com/Morisset/PyNeb_devel/tree/master/docs>
- The reference manual is accessible from <http://morisset.github.io/PyNeb_devel/>

## Discussion Groups

- https://groups.google.com/forum/#!forum/pyneb
- Send a mail to the group: pyneb@googlegroups.com

## Acknowledgements

This project is partly supported by grants DGAPA/PAPIIT-107215 and CONACyT-CB2015-254132.

PyNeb uses part of Chiantipy:

- Utility functions, many for reading the CHIANTI database files:

Copyright 2009, 2010 Kenneth P. Dere  
This software is distributed under the terms of the GNU General Public License that is found in the LICENSE file

- FortranFormat: Written by Konrad Hinsen <hinsen@cnrs-orleans.fr> With contributions from Andreas Prlic <andreas@came.sbg.ac.at> last revision: 2006-6-23

- Many thanks to Gloria, Jorge, Manu, Adal, Grazyna, Rubén, Enrique, César, Arturo and many other early adopters of the code.

## PyNeb's song

If [Cloudy](http://www.nublado.org) has a favorite [song](ftp://gradj.pa.uky.edu//gary//cloudy_gold_old//c94//Cloudy_Simon_and_Garfunkle.mp3), why shouldn't PyNeb have [one](http://research.iac.es/proyecto/PyNeb//files/elements.mp3) too?

