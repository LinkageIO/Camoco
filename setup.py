#!/usr/bin/env python3

from setuptools import setup, find_packages, Extension
from Cython.Distutils import build_ext
import os

import io
import re
import numpy

def read(*names, **kwargs):
    with io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")
    ) as fp:
        return fp.read()

def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

pccup = Extension(
    'camoco.PCCUP',
    sources=['camoco/PCCUP.pyx'],
#    extra_compile_args=['-ffast-math'],
#    include_dirs=[numpy.get_include()]
)
refgendist = Extension(
    'camoco.RefGenDist',
    sources=['camoco/RefGenDist.pyx'],
#    extra_compile_args=['-ffast-math'],
#    include_dirs=[numpy.get_include()]
)

setup(
    name = 'camoco',
    version = find_version('camoco','__init__.py'),
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
