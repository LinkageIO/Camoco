#!/usr/bin/python3

from camoco.Ontology import Term,Ontology
from camoco.Locus import Locus,Gene
from camoco.Camoco import Camoco
from camoco.RefGen import RefGen
from camoco.Tools import log
from camoco.COB import COB

from pandas import DataFrame
import pandas as pd

class CompCOB(Camoco):
    def __init__(COB,Ont,min_members=3,max_members=300):
        super().__init__(name=(COB.name + '_v_' + Ont.name))
        self.min = min_members
        self.max = max_members
        self.cob = cob
        self.ont = Ont

    def _find_term_groups(self,use_ont_terms=True):
        ont_list = self.ont.terms
        ont_terms = dict()
        if use_ont_terms:
            for id in ont_list:
                term = self.ont[id]
                if len(term.locus_list) >= self.min and len(term.locus_list) <= self.max:
                    ont_terms[id] = term
            self.ont_terms = terms
        else:
            pass
        return
