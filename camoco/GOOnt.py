#!/usr/bin/python3

from camoco.Ontology import Term,Ontology
from camoco.Camoco import Camoco
from camoco.RefGen import RefGen
from camoco.Locus import Locus
from camoco.Tools import log

from collections import defaultdict
import networkx as nx
import pandas as pd
import numpy as np

class GOOnt(Ontology):
    '''Ontology extension for GO'''
