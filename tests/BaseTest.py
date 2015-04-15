#!/usr/bin/python3

import unittest
import os
import numpy as np

import camoco as co
import pandas as pd
import itertools
from camoco.Config import cf

# Set the basedir to the testdir
cf.set('options','basedir', cf.get('options','testdir'))

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

class BaseCOB(unittest.TestCase):

    @unittest.skipUnless(co.available_datasets('Expr','ZmRoot'),'ZmRoot not defined')
    def setUp(self):
        self.cob = co.COB('ZmRoot')

    @unittest.skipUnless(co.available_datasets('Expr','ZmRoot'),'ZmRoot not defined')
    def ZmRoot(self):
        self.assertEqual(len(self.cob.coex), comb(self.cob.num_genes(),2),
                "coex dimentions mismatch")

    @unittest.skipUnless(co.available_datasets('Expr','ZmRoot'),'ZmRoot not defined')
    def test_coexpress_concordance(self):
        for a,b in itertools.combinations([self.cob.refgen.random_gene() for x in range(50)],2):
            self.assertTrue(abs(self.cob.coexpression(a,b).score - self.cob._coex_concordance(a,b)) < 0.001)
            dis_dif = abs(self.cob.coexpression(a,b).distance - abs(a-b))
            self.assertTrue( np.isnan(dis_dif) or dis_dif < 0.001)

if __name__ == '__main__':
    unittest.main()
