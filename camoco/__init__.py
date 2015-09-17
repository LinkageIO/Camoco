"""

Camoco Library - CoAnalysis of Molecular Components

CacheMoneyCorn

"""

__license__ = """
    Creative Commons Non-Commercial 4.0 Generic
    http://creativecommons.org/licenses/by-nc/4.0/
"""

import pyximport
pyximport.install()

from .Config import cf
from .Camoco import Camoco
from .Expr import Expr
from .COB import COB
from .RefGen import RefGen
from .Ontology import Ontology,Term
from .GWAS import GWAS
from .HapMap import HapMap
from .Locus import Locus
from .Tools import available_datasets,del_dataset
from .Tools import mv_dataset,redescribe_dataset
from .GEO import Family
from .GOnt import GOnt 
