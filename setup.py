#!/usr/bin/env python3

from setuptools import setup, find_packages, Extension

import os
import io
import re
import sys
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
    include_dirs=[numpy.get_include()]
)

# Check if running Windows
if os.name == 'nt':
    raise ValueError("Camoco isn't supported for windows! Camoco runs best on Linux!")
    sys.exit(1)

setup(
    name = 'camoco',
    version = find_version('camoco','__init__.py'),
    packages = find_packages(),
    scripts = [
    ],
    ext_modules = [pccup],
    cmdclass = {
    },
    package_data = {
        '':['*.cyx']    
    },
    python_requires='>=3.8',
    setup_requires = [
        # Setuptools 18.0 properly handles Cython extensions.
        'setuptools>=45.1.0',
        'numpy>=1.19.4',
        'cython>=0.29.13',
    ],
    install_requires = [		
        'minus80>=1.0.0',
        'locuspocus>=1.0.2',
        'tinymongo>=0.2.0',
        'fastcluster>=1.1.25',
        'scipy>=1.4.1',		
        'pandas>=1.0.1',		
        'fa2>=0.3.5',
        'statsmodels>=0.8.0',
       #'powerlaw>=1.4.6',
        'markov-clustering>=0.0.6.dev0',
        'psutil>=5.6.3',
    ],
    extras_require={
        'docs' : [
            'ipython>=6.5.0',
            'matplotlib>=2.2.3',
            'numpydoc',
            'sphinx_materialdesign_theme',
            'sphinx_rtd_theme',
            'sphinxcontrib-programoutput'
        ]
    },
    dependency_links = [
    ],
    include_package_data=True,
    author = 'Rob Schaefer',
    author_email = 'rob@linkage.io',
    description = 'Library for Co-Analysis of Molecular Componenets.',
    license = "Copyright Linkage Analytics 2017, Available under the MIT License",
    url = 'https://github.com/schae234/camoco'
)
