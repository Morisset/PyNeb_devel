#!/usr/bin/env python

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))
# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

from pyneb.version import __version__

setup(name='PyNeb', 
      version=__version__,
      description='Nebular tools',
      author='Christophe Morisset, Valentina Luridiana',
      long_description=long_description,
      author_email='chris.morisset@gmail.com',
      url='http://www.iac.es/proyecto/PyNeb/',
      py_modules=['pyneb.test.test_pyneb'],
      packages=['pyneb',
		'pyneb.core',
		'pyneb.plot',
		'pyneb.utils',
		'pyneb.extinction'],
      package_data={'pyneb':['atomic_data_fits/*.fits',
			     'atomic_data_fits/*.hdf5',
			     'atomic_data_fits/old_fits/*.fits',
			     'atomic_data/*.dat',
			     'atomic_data/*.func',
			     'atomic_data/*.pickle',
			     'atomic_data/levels/*.dat'],
                    'pyneb.extinction':['*.txt'],
                    'pyneb.plot':['level_diagrams/*.png']
                    }
     )

