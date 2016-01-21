#!/usr/bin/env python3

from setuptools import setup, find_packages, Extension
from Cython.Distutils import build_ext
import os

pccup = Extension(
    'PCCUP',
    sources=['camoco/PCCUP.pyx']
)
refgendist = Extension(
    'RefGenDist',
    sources=['camoco/RefGenDist.pyx']
)

setup(
    name = 'camoco',
    version = '0.1.8',
    packages = find_packages(),
    scripts = [
        'camoco/cli/camoco'
    ],
    ext_modules = [pccup,refgendist],
    cmdclass = {'build_ext': build_ext},

    package_data = {
        '':['*.cyx']    
    },
    install_requires = [		
        'cython>=0.16',    		
        'igraph>=0.1.5',		
        'matplotlib>=1.4.3',		
        'numpy>=1.9.1',		
        'pandas>=0.16',		
        'scipy>=0.15',		
        'termcolor>=1.1.0'		
    ],
    include_package_data=True,

    author = 'Rob Schaefer',
    author_email = 'schae234@gmail.com',
    description = 'Library for Co-Analysis of Molecular Componenets.',
    license = "Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License",
    url = 'https://github.com/schae234/camoco'
)
