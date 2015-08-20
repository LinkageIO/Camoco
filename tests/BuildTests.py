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
    suite.addTest(COB('AtSeed'))
    suite.addTest(COB('AtGen'))
    return suite


if __name__ == '__main__':
    unittest.main()
