#!/usr/bin/env python3

import numpy as np
import pandas as pd
import minus80 as m80
import locuspocus as lp

from enum import Enum

class Method(Enum):
    DENSITY
    LOCALITY

class Overlap(m80.Freezable):
    """
    The Overlap class represents the statistical enrichment of co-expression
    among a set of loci.

    """

    def __init__(self, Coex):
        self._create_tables()


    def _create_tables(self):
        cur = self.db.cursor()
        cur.execute("""
            CREATE TABLE IF NOT EXISTS overlap (
              OID INTEGER PRIMARY KEY AUTOINCREMENT,
              ontology TEXT,
              term TEXT,
              method TEXT,
              score REAL,
            );
        """)


    def overlap(self, loci, bootstrap=False, iter_name=None):
        """
        Calculate Network-Term Overlap based on the method in CLI args
        """
        # generate emirical and bootstrap values
        if self.args.method == "density":
            return self.cob.trans_locus_density(
                loci,
                flank_limit=self.args.candidate_flank_limit,
                by_gene=True,
                bootstrap=bootstrap,
                iter_name=iter_name,
            )
        elif self.args.method == "locality":
            return self.cob.trans_locus_locality(
                loci,
                flank_limit=self.args.candidate_flank_limit,
                by_gene=True,
                bootstrap=bootstrap,
                iter_name=iter_name,
                include_regression=True,
            ).rename(columns={"resid": "score"})
