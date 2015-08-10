""" 

Camoco Library - CoAnalysis of Molecular Components 

CacheMoneyCorn

"""

__license__ = """ 
    Creative Commons Non-Commercial 4.0 Generic
    http://creativecommons.org/licenses/by-nc/4.0/
"""

import pyximport; pyximport.install() 

from camoco.Camoco import Camoco
from camoco.Expr import Expr
from camoco.COB import COB
from camoco.RefGen import RefGen
from camoco.Ontology import Ontology,Term
from camoco.HapMap import HapMap
from camoco.Locus import Locus
from camoco.Tools import available_datasets,del_dataset,\
                         mv_dataset,redescribe_dataset
from camoco.Config import cf
from camoco.GEO import Family
