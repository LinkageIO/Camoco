""" 

Camoco Library - CoAnalysis of Molecular Components 

CacheMoneyCorn

"""

__license__ = """ 
    Creative Commons Non-Commercial 4.0 Generic
    http://creativecommons.org/licenses/by-nc/4.0/
"""

import pyximport; pyximport.install() 

from .Camoco import Camoco
from .Expr import Expr
from .COB import COB
from .RefGen import RefGen
from .Ontology import Ontology,Term
from .HapMap import HapMap
from .Locus import Locus
from .Tools import available_datasets,del_dataset,\
                   mv_dataset,redescribe_dataset
from .Config import cf
from .GEO import Family
