#!/usr/bin/python3

import unittest
import os
import numpy as np

import camoco as co
import pandas as pd
import itertools
from camoco.Config import cf

class Locus(unittest.TestCase):
    def test_locus_initialization(self):
        # numeric chromosomes
        a = co.Locus(1,500)
        self.assertIsInstance(a,co.Locus)

class RefGen(unittest.TestCase):

    def setUp(self):
        self.Fe = co.Ontology('ZmIonome')['Fe57']
        self.ZM = co.RefGen('Zm5bFGS')

    def test_candidate_vs_bootstrap_length(self):
        snps = self.Fe.effective_snps(window_size=50000)
        self.assertEqual(
            len(self.ZM.candidate_genes(snps)),
            len(self.ZM.bootstrap_candidate_genes(snps)),
            'Bootstrap vs emirical candidate gene sizes do not match'
        )

    def test_get_attr_and_generating_from_ids(self):
        self.assertIsInstance(ZM['GRMZM2G000014'],co.Locus)
        self.assertIsInstance(ZM.from_ids(['GRMZM2G000014'])[0],co.Locus)

class COB(unittest.TestCase):

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
