#!/usr/bin/env python3

from setuptools import setup, find_packages, Extension

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
#   extra_compile_args=['-ffast-math'],
    include_dirs=[numpy.get_include()]
)
refgendist = Extension(
    'camoco.RefGenDist',
    sources=['camoco/RefGenDist.pyx'],
#   extra_compile_args=['-ffast-math'],
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
        'camoco/cli/camoco'
    ],
    ext_modules = [pccup,refgendist],
    cmdclass = {
    },

    package_data = {
        '':['*.cyx']    
    },
    python_requires='>=3.6, <3.7',
    setup_requires = [
        # Setuptools 18.0 properly handles Cython extensions.
        'setuptools>=18.0',
        #'numpy==1.14.5',
        'cython',
    ],
#   include_dirs=['camoco/include'],
    install_requires = [		
        'minus80==0.3.3',
        'cython==0.29.10',    		
        'igraph==0.1.11',		
        'pyyaml==5.1.1',
        'matplotlib==3.1.0',		
        'numpy==1.16.4',		
        'scipy==1.2.1',		
        'pandas==0.23.4',		

        'scikit-learn==0.21.2',
        'statsmodels==0.8.0',
        'termcolor==1.1.0',
        'powerlaw==1.4.6',
        'flask==1.0.3',
        'markov-clustering==0.0.6.dev0',
        'fastcluster==1.1.25'
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
        #'git+https://github.com/LinkageIO/Minus80'
        #'git+https://github.com/rogerbinns/apsw' 
    ],
    include_package_data=True,

    author = 'Rob Schaefer',
    author_email = 'rob@linkgae.io',
    description = 'Library for Co-Analysis of Molecular Componenets.',
    license = "Copyright Linkage Analytics 2017, Available under the MIT License",
    url = 'https://github.com/schae234/camoco'
)
