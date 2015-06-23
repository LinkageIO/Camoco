#!/usr/bin/env python3

from setuptools import setup, find_packages, Extension
import os

setup(
    name = 'camoco',
    version = '0.1.8',
    packages = find_packages(),
    scripts = [
        'bin/BootstrapLocality.py',
        'interfaces/CamocoWeb.py'
    ],

    install_requires = [
        'cython>=0.16',    
        'igraph>=0.1.5',
        'matplotlib>=1.4.3',
        'numpy>=1.9.1',
        'pandas>=0.16',
        'scipy>=0.15',
        'termcolor>=1.1.0'
    ],

    package_data = {
        '':['*.cyx']    
    },
    include_package_data=True,

    author = 'Rob Schaefer',
    author_email = 'schae234@gmail.com',
    description = 'Library for Co-Analysis of Molecular Componenets.',
    license = "Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License",
    url = 'https://github.com/schae234/camoco'
)
