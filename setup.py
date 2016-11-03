#!/usr/bin/env python
'''
Distutils setup.py file for the pymeanfield python module.

TODO:
- set up compilation of Fortran extension upon install

'''
import os
from setuptools import setup

with open('README.md') as file:
    long_description = file.read()


setup(
    name = 'pymeanfield',
    version = '0.1.0',
    maintainer = 'Jannis Schuecker',
    maintainer_email = 'j.schuecker@fz-juelich.de',
    url = 'https://github.com/INM-6/neural_network_meanfield',
    download_url = 'https://github.com/INM-6/neural_network_meanfield/tarball/v0.1.0',
    packages = ['pymeanfield'],
    provides = ['pymeanfield'],
    package_data = {},
    description = 'population rate and population spectra from mean-field network models',
    long_description = long_description,
    license='LICENSE',
        classifiers=[
            'License :: OSI Approved :: GNU General Public License (GPL)',
            'Programming Language :: Python',
            'Programming Language :: Python :: 2.7',
            'Operating System :: OS Independent',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Physics',
            'Intended Audience :: Developers',
            'Intended Audience :: Science/Research',
            'Development Status :: 4 - Beta',
            ],
    install_requires = [
        'numpy', 'scipy', 'matplotlib',
        ],

)
