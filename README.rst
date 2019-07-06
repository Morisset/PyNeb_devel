PyNeb (Luridiana V., Morisset C. and Shaw, R. A 2013) is a modern python tool to compute emission line emissivities (recombination and collisionally excited lines).

In PyNeb, the atom is represented as an n-level atom. For given density and temperature, PyNeb's machinery solves the equilibrium equations and determines the level populations. These are some of the things it can do:

* compute physical conditions from suitable diagnostic line ratios.
* compute level populations, critical densities and line emissivities 
* compute and display emissivity grids as a function of Te and Ne
* deredden line intensities
* read and manage observational data
* plot and compare atomic data from different publications
* compute ionic abundances from line intensities and physical conditions
* compute elemental abundances from ionic abundances and icfs.

PyNeb also incorporates emissivity tables of recombination lines for a few atoms. The interpolated emissivities can be used by any of the module that rely on the n-level line emissivities to carry out the actions listed above.

Emission line ratios are used to self consistently determine electron temperature and density and ionic abundances
Diagnostic diagrams can easily be plotted.
Various ionization correction factors (ICFs) from the literarure are available to obtain total elemental abundances from the ionic abundances.
Atomic data can easily be changed and updated.
Additional tools are provided, like reddening determination and correction procedures.

Requirements
============

PyNeb uses numpy, matplotlib, pyfits, scipy and other standard python libraries.

Installation
============

You may find useful to download, install and upgrade PyNeb using `pip <http://www.pip-installer.org/en/latest/index.html>`_.

For example:

* pip install -U PyNeb

Note: you MAY need --user if you installed python without Anaconda or Canopy

Updates use the same command.

You can also install from the github repository:

* pip install -U git+https://github.com/Morisset/PyNeb_devel.git

To use the development branch (at your own risks!!!):

* pip install -U git+https://github.com/Morisset/PyNeb_devel.git@devel

Warranty
========

PyNeb is provided as it is. No warranty at all.

Manual
======

* The manuals are here: `<https://github.com/Morisset/PyNeb_devel/tree/master/docs>`_.

* The reference manual is accessible from `<https://morisset.github.io/PyNeb_Manual/html/index.html>`_.

Discussion Groups
=================
* https://groups.google.com/forum/#!forum/pyneb
* Send a mail to the group: pyneb@googlegroups.com

Acknowledgements
================

This project is partly supported by grants DGAPA/PAPIIT-107215 and CONACyT-CB2015-254132.
