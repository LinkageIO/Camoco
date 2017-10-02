"""

Camoco Library - CoAnalysis of Molecular Components

CacheMoneyCorn

"""

__license__ = """
    Creative Commons Non-Commercial 4.0 Generic
    http://creativecommons.org/licenses/by-nc/4.0/
"""

__version__ = '0.3.1.dev'

import sys
import os
import numpy

import pyximport
pyximport.install(setup_args={
    "include_dirs":numpy.get_include() 
})

import matplotlib
matplotlib.use('Agg')
from .Config import cf
from .Camoco import Camoco
from .Expr import Expr
from .COB import COB
from .RefGen import RefGen
from .RefGenDist import *
from .PCCUP import *
from .Ontology import Ontology,Term
from .GWAS import GWAS
from .Locus import Locus
from .Tools import available_datasets,del_dataset
from .Tools import mv_dataset,redescribe_dataset
from .GEO import Family
from .GOnt import GOnt
from .Overlap import Overlap
#from .Annotation import GWASData

# Create yourself
Camoco.create('Camoco','Mother Database')
