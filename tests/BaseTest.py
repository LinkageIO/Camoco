#!/usr/bin/python3

import unittest
import os
import numpy as np

import camoco as co
import pandas as pd
import itertools
from camoco.Config import cf

class LocusBase(unittest.TestCase):
    def test_locus_initialization(self):
        # numeric chromosomes
        a = co.Locus(1,500)
        self.assertIsInstance(a,co.Locus)

class BaseRefGen(unittest.TestCase):
    @unittest.skipUnless(
            co.available_datasets('RefGen','Zm5bFGS') and 
            co.available_datasets('Ontology','ZmIonome'),'Not Adequate Datasets')

    def test_candidate_vs_bootstrap_length(self):
        ZM = co.RefGen('Zm5bFGS')
        snps = co.Ontology("ZmIonome")['Fe57'].effective_snps(window_size=50000)
        self.assertEqual(
            len(ZM.candidate_genes(snps)),
            len(ZM.bootstrap_candidate_genes(snps)),
            'Bootstrap vs emirical candidate gene sizes do not match'
        )

    @unittest.skipUnless(co.available_datasets('RefGen','Zm5bFGS'),'Not Adequate Datasets')
    def test_get_attr_and_generating_from_ids(self):
        ZM = co.RefGen('Zm5bFGS')
        self.assertIsInstance(ZM['GRMZM2G000014'],co.Locus)
        self.assertIsInstance(ZM.from_ids(['GRMZM2G000014'])[0],co.Locus)

class BaseCOB(unittest.TestCase):

    @unittest.skipUnless(co.available_datasets('Expr','ZmRoot'),'ZmRoot not defined')
    def setUp(self):
        self.cob = co.COB('ZmRoot')

    def ZmRoot(self):
        self.assertEqual(len(self.cob.coex), comb(self.cob.num_genes(),2),
                "coex dimentions mismatch")

    def test_coexpress_concordance(self):
        for a,b in itertools.combinations([self.cob.refgen.random_gene() for x in range(50)],2):
            self.assertTrue(abs(self.cob.coexpression(a,b).score - self.cob._coex_concordance(a,b)) < 0.001)
            dis_dif = abs(self.cob.coexpression(a,b).distance - abs(a-b))
            self.assertTrue( np.isnan(dis_dif) or dis_dif < 0.001)

    def test_num_neighbors_equals_degree(self):
        random_gene = self.cob.refgen.random_gene()
        self.assertTrue(len(self.cob.neighbors(random_gene) == self.cob.global_degree(random_gene)))

if __name__ == '__main__':
    unittest.main()
