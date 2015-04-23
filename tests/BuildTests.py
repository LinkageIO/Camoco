#!/usr/bin/python3

import unittest
import os
import numpy as np

import camoco as co
import pandas as pd
import itertools
from camoco.Config import cf


# write test case to import refgen from GFF

class RefGen(unittest.TestCase):
    def Tair10(self):
        gff = os.path.join(cf.get('options','testdir'),'raw','TAIR10_GFF3_genes.gff')
        co.del_dataset('RefGen','T10',safe=False)
        T10 = co.RefGen.from_gff(gff,'T10','Tair 10','10','Arabidopsis')
        self.assertIsInstance(T10,co.RefGen)

    def Zm5bFGS(self):
        gff = os.path.join(cf.get('options','testdir'),'raw','ZmB73_5b_FGS.gff')
        co.del_dataset('RefGen','Zm5bFGS',safe=False)
        ZM = co.RefGen.from_gff(gff,'Zm5bFGS','Maize 5b Filtered Gene Set','5b','Zea Mays')
        self.assertIsInstance(ZM,co.RefGen)

class Ontology(unittest.TestCase):
    def Ionome(self):
        csv = os.path.join(cf.get('options','testdir'),'raw','sigGWASsnpsCombinedIterations.longhorn.allLoc.csv')
        ZM = co.RefGen('Zm5bFGS')
        df = pd.DataFrame.from_csv(csv,index_col=None)
        co.del_dataset('Ontology','ZmIonome')
        IONS  = co.Ontology.from_DataFrame(df,'ZmIonome','Maize Ionome',ZM,term_col='el',chr_col='chr',pos_col='pos');
        self.assertIsInstance(IONS,co.Ontology)

class COB(unittest.TestCase):

    def AtLeaf(self):
        Leaf = ['GSE14578','GSE5630','GSE13739', #'GSE26199',
                'GSE5686','GSE5615','GSE5620','GSE5628','GSE5624','GSE5626','GSE5621','GSE5622',
                'GSE5623','GSE5625','GSE5688']
        LeafFam = sum([co.Family.from_file("raw/GSE/{}_family.soft".format(x)) for x in Leaf ])
        # Generate the LeafKeep file
        #LeafFam.to_keepfile("LeafKeep.tsv",keep_hint="lea")
        AtLeaf = co.COB.from_DataFrame(LeafFam.series_matrix(keepfile="raw/GSE/LeafKeep.tsv"),'AtLeaf','Arabidopsis Leaf',co.RefGen('Tair10'),rawtype='MICROARRAY')
        self.assertIsInstance(AtLeaf,co.COB)

    def AtSeed(self):
        Seed = ['GSE12404',#'GSE30223',
                'GSE1051','GSE11852','GSE5634']
        SeedFam = sum([co.Family.from_file("raw/GSE/{}_family.soft".format(x)) for x in Seed ])
        #SeedFam.to_keepfile("SeedKeep.tsv",keep_hint='seed')
        AtSeed = co.COB.from_DataFrame(SeedFam.series_matrix(keepfile="raw/GSE/SeedKeep.tsv"),'AtSeed','Arabidopsis Seed',co.RefGen('Tair10'),rawtype='MICROARRAY')

    def AtRoot(self):
        Root = ['GSE14578','GSE46205','GSE7631','GSE10576','GSE42007','GSE34130','GSE21611','GSE22966','GSE7641','GSE5620',
                'GSE8934','GSE5628','GSE30095','GSE30097','GSE5624','GSE5626','GSE5749','GSE5621','GSE5622','GSE5623','GSE5625','GSE5688']
        RootFam = sum([co.Family.from_file("raw/GSE/{}_family.soft".format(x)) for x in Root ])
        #RootFam.to_keepfile("RootKeep.tsv",keep_hint='root')
        AtRoot = co.COB.from_DataFrame(RootFam.series_matrix(keepfile="raw/GSE/RootKeep.tsv"),'AtRoot','Arab Root',co.RefGen('Tair10'),rawtype='MICROARRAY')

    def AtGen(self):
        General = ['GSE18975','GSE39384','GSE19271','GSE5632','GSE39385','GSE5630','GSE15617','GSE5617','GSE5686','GSE2473',
                   'GSE5633','GSE5620','GSE5628','GSE5624','GSE5626','GSE5621','GSE5622','GSE5623','GSE5625','GSE5688']
        GenFam = sum([co.Family.from_file("raw/GSE/{}_family.soft".format(x)) for x in General ])
        #GenFam.to_keepfile("GenKeep.tsv")
        AtGen = co.COB.from_DataFrame(GenFam.series_matrix(keepfile="raw/GSE/GenKeep.tsv"),'AtGen','Arab General',co.RefGen('Tair10'),rawtype='MICROARRAY')

    def ZmSAM(self):
        co.del_dataset('Expr','ZmSAM',safe=False)
        ZM = co.RefGen('Zm5bFGS')
        ZmSAM = co.COB.from_table(
            os.path.join(
                cf.get('options','testdir'),'raw','Expression',
                'TranscriptomeProfiling_B73_Atlas_SAM_FGS_LiLin_20140316.txt'
            ),
            'ZmSAM',
            'Maize Root Network',
            ZM,
            rawtype='RNASEQ',
            max_gene_missing_data=0.4,
            min_expr=0.1,
            quantile=False,
            dry_run=False,
            max_val=300
        )

    def ZmPAN(self):
        co.del_dataset('Expr','ZmPAN',safe=False)
        ZM = co.RefGen('Zm5bFGS')
        ZmPAN = co.COB.from_table(
            os.path.join(
                cf.get('options','testdir'),'raw','Expression',
                'PANGenomeFPKM.txt'
            ),
            'ZmPAN',
            'Maize Root Network',
            ZM,
            rawtype='RNASEQ',
            max_gene_missing_data=0.4,
            min_expr=0.1,
            quantile=False,
            dry_run=False,
            sep=',',
            max_val=300
        )

    def ZmRoot(self):
        co.del_dataset('Expr','ZmRoot',safe=False)
        ZM = co.RefGen('Zm5bFGS')
        ZmRoot = co.COB.from_table(
            os.path.join(cf.get('options','testdir'),'raw','Expression','ROOTFPKM.tsv'),
            'ZmRoot',
            'Maize Root Network',
            ZM,
            rawtype='RNASEQ',
            max_gene_missing_data=0.4,
            min_expr=0.1,
            quantile=False,
            max_val=300
        )

if __name__ == '__main__':
    unittest.main()
