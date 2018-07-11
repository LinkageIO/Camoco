#!/usr/bin/env python3

from setuptools import setup, find_packages, Extension
from Cython.Distutils import build_ext
from setuptools.command.develop import develop
from setuptools.command.install import install

from subprocess import check_call

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
    include_dirs=[numpy.get_include()]
)
refgendist = Extension(
    'camoco.RefGenDist',
    sources=['camoco/RefGenDist.pyx'],
#    extra_compile_args=['-ffast-math'],
    include_dirs=[numpy.get_include()]
)

setup(
    name = 'camoco',
    version = find_version('camoco','__init__.py'),
    packages = find_packages(),
    scripts = [
        'camoco/cli/camoco'
    ],
    ext_modules = [pccup,refgendist],
    cmdclass = {
        'build_ext': build_ext,
    },

    package_data = {
        '':['*.cyx']    
    },
    install_requires = [		
        'minus80==0.1.3',
        'flask==0.12.2',
        'cython==0.16',    		
        'igraph==0.1.5',		
        'matplotlib==2.2.2',		
        'numpy==1.14.5',		
        'pandas==0.22.0',		
        'scipy==1.1.0',		
        'termcolor==1.1.0',
        'scikit-learn==0.19.1',
        'powerlaw==1.3.5',
        'statsmodels==0.9.0'
    ],
    include_package_data=True,

    author = 'Rob Schaefer',
    author_email = 'rob@linkgae.io',
    description = 'Library for Co-Analysis of Molecular Componenets.',
    license = "Copyright Linkage Analytics 2017, Available under the MIT License",
    url = 'https://github.com/schae234/camoco'
)
