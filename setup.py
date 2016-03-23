#!/usr/bin/env python3

from setuptools import setup, find_packages, Extension
from Cython.Distutils import build_ext
import os
import numpy

pccup = Extension(
    'PCCUP',
    sources=['camoco/PCCUP.pyx'],
    extra_compile_args=['-ffast-math'],
    inlcude_dirs=[numpy.get_include()]
)
refgendist = Extension(
    'RefGenDist',
    sources=['camoco/RefGenDist.pyx'],
    extra_compile_args=['-ffast-math'],
    inlcude_dirs=[numpy.get_include()]
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
        'termcolor>=1.1.0',
        'powerlaw==1.3.5'
    ],
    include_package_data=True,

    author = 'Rob Schaefer',
    author_email = 'schae234@gmail.com',
    description = 'Library for Co-Analysis of Molecular Componenets.',
    license = "Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License",
    url = 'https://github.com/schae234/camoco'
)
