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

class GWASOnt(Ontology):
    '''Ontology extension for GWAS'''
    def __init__():





    @classmethod
    def from_DataFrame(cls, df, name, description, refgen,
            term_col='Term', chr_col=None, pos_col=None,
            start_col=None, end_col=None, id_col=None
            ):
        '''
            Import an ontology from a pandas dataframe.
            Groups by term_col, then iterates over the rows in the group.
            It adds each row as a locus to the term.
        '''
        self = cls.create(name, description, refgen)
        # group each trait by its name
        for term_id, df in df.groupby(term_col):
            term = Term(term_id)
            # we have a SNP
            if pos_col is not None:
                for i, row in df.iterrows():
                    # make sure there are no collisions with the SNP init func names
                    # this is hackey and I dont like it
                    kwargs = {key:val for key, val in dict(row).items() if key not in Locus.__init__.__code__.co_varnames}
                    snp = Locus(row[chr_col], int(row[pos_col]), int(row[pos_col]), gene_build=self.refgen.build, **kwargs)
                    term.locus_list.add(snp)
            self.log("Importing {}", term)
            self.add_term(term)
        return self
