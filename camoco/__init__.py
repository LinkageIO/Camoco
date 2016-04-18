"""

Camoco Library - CoAnalysis of Molecular Components

CacheMoneyCorn

"""

__license__ = """
    Creative Commons Non-Commercial 4.0 Generic
    http://creativecommons.org/licenses/by-nc/4.0/
"""

__version__ = '0.1.8'

import sys
import os

#sys.path.insert(0,os.path.expanduser('~/.camoco/.conda/camoco/lib/python3.4/site-packages'))

from .Config import cf
from .Camoco import Camoco
from .Expr import Expr
from .COB import COB
from .RefGen import RefGen
from .RefGenDist import *
from .PCCUP import *
from .Ontology import Ontology,Term
from .GWAS import GWAS
from .HapMap import HapMap
from .Locus import Locus
from .Tools import available_datasets,del_dataset
from .Tools import mv_dataset,redescribe_dataset
from .GEO import Family
from .GOnt import GOnt 
