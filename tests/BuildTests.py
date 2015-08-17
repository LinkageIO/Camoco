#!/usr/bin/python3

import itertools
import os
import glob

import numpy as np
import camoco as co
import pandas as pd

from camoco.Config import cf
from os.path import join as pjoin 

import unittest

# write test case to import refgen from GFF

class RefGen(unittest.TestCase):

    def Tair10(self):
        gff = os.path.expanduser(
            os.path.join(
                cf['options']['testdir'],
                'raw','RefGen','TAIR10_GFF3_genes.gff'
            )
        )
        co.del_dataset('RefGen','T10',safe=False)
        T10 = co.RefGen.from_gff(gff,'T10','Tair 10','10','Arabidopsis')
        self.assertIsInstance(T10,co.RefGen)

    def Zm5bFGS(self):
        gff = os.path.expanduser(
            os.path.join(
                cf['options']['testdir'],
                'raw','RefGen','ZmB73_5b_FGS.gff.gz'
            )
        )
        co.del_dataset('RefGen','Zm5bFGS',safe=False)
        ZM = co.RefGen.from_gff(
            gff,'Zm5bFGS','Maize 5b Filtered Gene Set','5b','Zea Mays'
        )
        self.assertIsInstance(ZM,co.RefGen)

def RefGenSuite():
    suite = unittest.TestSuite()
    suite.addTest(RefGen('Tair10'))
    suite.addTest(RefGen('Zm5bFGS'))
    return suite

class Ontology(unittest.TestCase):

    def AtSeedIonome(self):
        co.del_dataset('Ontology','AtSeedIonome',safe=False)
        # glob glob is god
        csvs = glob.glob(os.path.join(
            cf.get('options','testdir'),
            'raw','GWAS','AtLeaf',
            '*.sigsnps.csv'
        ))
        # Read in each table individually then concat for GIANT table
        df = pd.concat([pd.read_table(x,sep=',') for x in csvs])
        # Add 'Chr' to chromosome column
        df.CHR = df.CHR.apply(lambda x: 'Chr'+str(x))
        # Chase dat refgen
        T10 = co.RefGen('T10')
        # Import class from dataframe
        AtSeedIonome = co.Ontology.from_DataFrame(
            df,'AtSeedIonome','Arabidopsis 1.6M EmmaX GWAS',
            T10,term_col='Trait',chr_col='CHR',pos_col='BP'
        )
        self.assertIsInstance(AtSeedIonome,co.Ontology)


    def AtLeafIonome(self):
        co.del_dataset('Ontology','AtLeafIonome',safe=False)
        # glob glob is god
        csvs = glob.glob(os.path.join(
            cf.get('options','testdir'),
            'raw','GWAS','AtLeaf',
            '*.sigsnps.csv'
        ))
        # Read in each table individually then concat for GIANT table
        df = pd.concat([pd.read_table(x,sep=',') for x in csvs])
        # Add 'Chr' to chromosome column
        df.CHR = df.CHR.apply(lambda x: 'Chr'+str(x))
        # Chase dat refgen
        T10 = co.RefGen('T10')
        # Import class from dataframe
        AtLeaf = co.Ontology.from_DataFrame(
            df,'AtLeafIonome','Arabidopsis 1.6M EmmaX GWAS',
            T10,term_col='Trait',chr_col='CHR',pos_col='BP'
        )
        self.assertIsInstance(AtLeaf,co.Ontology)

    def AtLeafHydroIonome(self):
        co.del_dataset('Ontology','AtLeafHydroIonome',safe=False)
        # glob glob is god
        csvs = glob.glob(os.path.join(
            cf.get('options','testdir'),
            'raw','GWAS','AtLeafHydro',
            '*.sigsnps.csv'
        ))
        # Read in each table individually then concat for GIANT table
        df = pd.concat([pd.read_table(x,sep=',') for x in csvs])
        # Add 'Chr' to chromosome column
        df.CHR = df.CHR.apply(lambda x: 'Chr'+str(x))
        # Chase dat refgen
        T10 = co.RefGen('T10')
        # Import class from dataframe
        AtLeafHydroIonome = co.Ontology.from_DataFrame(
            df,'AtLeafHydroIonome','Arabidopsis 1.6M EmmaX GWAS',
            T10,term_col='Trait',chr_col='CHR',pos_col='BP'
        )
        self.assertIsInstance(AtLeafHydroIonome,co.Ontology)


    def AtRootHydroIonome(self):
        co.del_dataset('Ontology','AtRootHydroIonome',safe=False)
        # glob glob is god
        csvs = glob.glob(os.path.join(
            cf.get('options','testdir'),
            'raw','GWAS','AtRootHydro',
            '*.sigsnps.csv'
        ))
        # Read in each table individually then concat for GIANT table
        df = pd.concat([pd.read_table(x,sep=',') for x in csvs])
        # Shorten the term name 
        df.Trait = df.Trait.apply(lambda x: x.replace('RootHydro.',''))
        # Add 'Chr' to chromosome column
        df.CHR = df.CHR.apply(lambda x: 'Chr'+str(x))
        # Chase dat refgen
        T10 = co.RefGen('T10')
        # Import class from dataframe
        AtRootHydroIonome = co.Ontology.from_DataFrame(
            df,'AtRootHydroIonome','Arabidopsis 1.6M EmmaX GWAS',
            T10,term_col='Trait',chr_col='CHR',pos_col='BP'
        )
        self.assertIsInstance(AtRootHydroIonome,co.Ontology)

    def ZmIonome(self):
        # Delete the old dataset
        co.del_dataset('Ontology','ZmIonome',safe=False)
        # Grab path the csv
        csv = os.path.join(
            cf.get('options','testdir'),
            'raw','GWAS','Ionome',
            'sigGWASsnpsCombinedIterations.longhorn.allLoc.csv'
        )
        # Define our reference geneome
        ZM = co.RefGen('Zm5bFGS')
        df = pd.DataFrame.from_csv(csv,index_col=None)
        # Import class from dataframe
        IONS  = co.Ontology.from_DataFrame(
            df,'ZmIonome','Maize Ionome',
            ZM,term_col='el',chr_col='chr',pos_col='pos'
        )
        IONS.del_term('Co59')
        # I guess we need a test in here too
        self.assertIsInstance(IONS,co.Ontology)

def OntologySuite():
    suite = unittest.TestSuite()
    suite.addTest(Ontology('AtSeedIonome'))
    suite.addTest(Ontology('AtLeafIonome'))
    suite.addTest(Ontology('AtLeafHydroIonome'))
    suite.addTest(Ontology('AtRootHydroIonome'))
    suite.addTest(Ontology('ZmIonome'))
    return suite


class COB(unittest.TestCase):

    def setUp(self):
        self.tdir = pjoin(cf['options']['testdir'])
        self.rawdir = pjoin(cf['options']['testdir'],'raw')

    def AtLeaf(self):
        co.del_dataset('Expr','AtLeaf',safe=False)
        Leaf = ['GSE14578','GSE5630','GSE13739', #'GSE26199',
                'GSE5686','GSE5615','GSE5620','GSE5628',
                'GSE5624','GSE5626','GSE5621','GSE5622',
                'GSE5623','GSE5625','GSE5688']
        LeafFam = sum(
            [co.Family.from_file(
                pjoin(
                    self.rawdir,'GSE','{}_family.soft'.format(x)
                )
            )
            for x in Leaf ]
        )
        # Generate the LeafKeep file
        #LeafFam.to_keepfile("LeafKeep.tsv",keep_hint="lea")
        AtLeaf = co.COB.from_DataFrame(
            LeafFam.series_matrix(
                keepfile=pjoin(self.rawdir,'GSE','LeafKeep.tsv')
            ),
            'AtLeaf',
            'Arabidopsis Leaf',
            co.RefGen('T10'),
            rawtype='MICROARRAY',
            max_gene_missing_data=0.3,
            min_expr=0.01,
            quantile=True,
        )
        self.assertIsInstance(AtLeaf,co.COB)

    def AtSeed(self):
        co.del_dataset('Expr','AtSeed',safe=False)
        Seed = ['GSE12404',#'GSE30223',
                'GSE1051','GSE11852','GSE5634']
        SeedFam = sum(
            [co.Family.from_file(
                pjoin(
                    self.rawdir,'GSE','{}_family.soft'.format(x)
                )
            )
            for x in Seed ]
        )
        #SeedFam.to_keepfile("SeedKeep.tsv",keep_hint='seed')
        AtSeed = co.COB.from_DataFrame(
            SeedFam.series_matrix(
                keepfile=pjoin(self.rawdir,'GSE','SeedKeep.tsv')
            ),
            'AtSeed','Arabidopsis Seed',
            co.RefGen('T10'),rawtype='MICROARRAY',
            quantile=True

        )
        self.assertIsInstance(AtSeed,co.COB)

    def AtRoot(self):
        co.del_dataset('Expr','AtRoot',safe=False)
        Root = ['GSE14578','GSE46205','GSE7631','GSE10576','GSE42007',
                'GSE34130','GSE21611','GSE22966','GSE7641','GSE5620',
                'GSE8934','GSE5628','GSE30095','GSE30097','GSE5624',
                'GSE5626','GSE5749','GSE5621','GSE5622',
                'GSE5623','GSE5625','GSE5688']
        RootFam = sum(
            [co.Family.from_file(
                os.path.join(
                    cf['options']['testdir'],
                    'raw','GSE','{}_family.soft'.format(x)
                )
            )
            for x in Root ]
        )
        #RootFam.to_keepfile("RootKeep.tsv",keep_hint='root')
        AtRoot = co.COB.from_DataFrame(
            RootFam.series_matrix(
                keepfile=pjoin(self.rawdir,'GSE','RootKeep.tsv')
            ),
            'AtRoot','Arab Root',
            co.RefGen('T10'),
            rawtype='MICROARRAY',
            quantile=True
        )
        self.assertIsInstance(AtRoot,co.COB)


    def AtGen(self):
        General = ['GSE18975','GSE39384','GSE19271','GSE5632','GSE39385',
                   'GSE5630','GSE15617','GSE5617','GSE5686','GSE2473',
                   'GSE5633','GSE5620','GSE5628','GSE5624',
                   'GSE5626','GSE5621','GSE5622','GSE5623','GSE5625','GSE5688']
        GenFam = sum(
            [co.Family.from_file(
                pjoin(
                    self.rawdir,'GSE','{}_family.soft'.format(x)
                )
            )
            for x in General ]
        )
        #GenFam.to_keepfile("GenKeep.tsv")
        AtGen = co.COB.from_DataFrame(
            GenFam.series_matrix(
                keepfile=pjoin(self.rawdir,'GSE','GenKeep.tsv')
            ),
            'AtGen','Arab General',
            co.RefGen('T10'),
            rawtype='MICROARRAY'
        )
        self.assertIsInstance(AtGen,co.COB)

    def ZmSAM(self):
        co.del_dataset('Expr','ZmSAM',safe=False)
        ZM = co.RefGen('Zm5bFGS')
        ZmSAM = co.COB.from_table(
            os.path.join(
                cf.get('options','testdir'),
                'raw','Expression','RNASEQ',
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
            max_val=250
        )
        self.assertIsInstance(ZmSAM,co.COB)

    def ZmPAN(self):
        co.del_dataset('Expr','ZmPAN',safe=False)
        ZM = co.RefGen('Zm5bFGS')
        ZmPAN = co.COB.from_table(
            os.path.join(
                cf.get('options','testdir'),'raw','Expression','RNASEQ',
                'PANGenomeFPKM.txt'
            ),
            'ZmPAN',
            'Maize Root Network',
            ZM,
            rawtype='RNASEQ',
            max_gene_missing_data=0.4,
            min_expr=1,
            quantile=False,
            dry_run=False,
            sep=',',
            max_val=300
        )
        self.assertIsInstance(ZmPAN,co.COB)

    def ZmRoot(self):
        co.del_dataset('Expr','ZmRoot',safe=False)
        ZM = co.RefGen('Zm5bFGS')
        ZmRoot = co.COB.from_table(
            os.path.join(cf.get('options','testdir'),'raw','Expression',
                'RNASEQ','ROOTFPKM.tsv'),
            'ZmRoot',
            'Maize Root Network',
            ZM,
            rawtype='RNASEQ',
            max_gene_missing_data=0.3,
            max_accession_missing_data=0.08,
            min_single_sample_expr=1,
            min_expr=0.001,
            quantile=False,
            max_val=300
        )
        self.assertIsInstance(ZmRoot,co.COB)

def COBSuite():
    suite = unittest.TestSuite()
    suite.addTest(COB('AtLeaf'))
    suite.addTest(COB('AtRoot'))
    suite.addTest(COB('AtGen'))
    suite.addTest(COB('ZmSAM'))
    suite.addTest(COB('ZmPAN'))
    suite.addTest(COB('ZmRoot'))
    return suite

def COBAtSuite():
    suite = unittest.TestSuite()
    suite.addTest(COB('AtLeaf'))
    suite.addTest(COB('AtRoot'))
    suite.addTest(COB('AtGen'))
    return suite


if __name__ == '__main__':
    unittest.main()
