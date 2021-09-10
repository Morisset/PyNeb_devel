# PyNeb's Home Page

## Introduction

PyNeb (by V. Luridiana, C. Morisset, and R. A. Shaw) is the last in a lineage of tools dedicated to the analysis of emission lines, which includes [FIVEL](https://ui.adsabs.harvard.edu/abs/1987JRASC..81..195D/abstract) and [nebular](http://research.iac.es/proyecto/PyNeb/) (see also [here](https://ui.adsabs.harvard.edu/abs/1998ASPC..145..192S/abstract)). PyNeb is completely written in [python](https://www.python.org) and is designed to be more user-friendly and powerful than its predecessors, its functionality representing a giant leap forward with respect to them in terms of speed, easiness of use, graphic visualization, and accessible information. 

In PyNeb, the atom is represented as an n-level atom. For given density and temperature, PyNeb's machinery solves the equilibrium equations and determines the level populations. These are some of the things it can do:

compute physical conditions from suitable diagnostic line ratios;
compute level populations, critical densities and line emissivities;
compute and display emissivity grids as a function of Te and Ne;
deredden line intensities;
read and manage observational data;
plot and compare atomic data from different publications;
compute ionic abundances from line intensities and physical conditions;
compute elemental abundances from ionic abundances and icfs.
PyNeb also incorporates emissivity tables of recombination lines for a few atoms (currently, H and He ions and O II). The interpolated emissivities can be used by any of the modules that rely on the n-level line emissivities to carry out the actions listed above.

## Download and install

To install or upgrade PyNeb, you need pip. Get it [here](https://pypi.org/project/pip/) and install it using:

* `sudo easy_install pip`

Once pip is installed, enter from the command line:

* `pip install --user PyNeb`

to install PyNeb, and:

* `pip install --upgrade --user PyNeb`

to upgrade previous versions of PyNeb (but omit the option --user if your python has been installed with Ureka or Canopy).

Uninstalling PyNeb is easy as well

* `pip uninstall PyNeb`

Problems? Goto [Troubleshooting](http://research.iac.es/proyecto/PyNeb//troubleshooting.html)

## Manual

Download PyNeb's manual [here](https://github.com/Morisset/PyNeb_devel/tree/master/docs).

## Support

There is a [discussion](https://groups.google.com/accounts/SetOSID) group to post requests, help other users, and share your problems.

To get a quick answer to a problem, send an email to the [group account](mailto:pyneb@googlegroups.com). 

## Sample scripts

Several useful scripts to get started with PyNeb can be found in the final section of [PyNeb's manual](https://github.com/Morisset/PyNeb_devel/tree/master/docs).

## Acknowledgments

Many thanks to Gloria, Jorge, Manu, Adal, Grazyna, Rubén, Enrique, César, Arturo and many other early adopters of the code.

## PyNeb's song

If [Cloudy](http://www.nublado.org) has a favorite [song](ftp://gradj.pa.uky.edu//gary//cloudy_gold_old//c94//Cloudy_Simon_and_Garfunkle.mp3), why shouldn't PyNeb have [one](http://research.iac.es/proyecto/PyNeb//files/elements.mp3) too?

