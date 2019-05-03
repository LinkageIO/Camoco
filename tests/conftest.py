import pytest
import os
import glob


from camoco.Config import cf
from camoco import Tools as tools

import camoco as co
import pandas as pd

''' -------------------------------------------------------------------------
            'Test' Fixtures
'''

@pytest.fixture(scope='module')
def testRefGen(Zm5bFGS):
    # This was a mistake
    return Zm5bFGS

@pytest.fixture(scope='module')
def testGWAS(ZmWallace):
    return ZmWallace

@pytest.fixture(scope='module')
def testCOB(ZmRNASeqTissueAtlas):
    return ZmRNASeqTissueAtlas

''' -------------------------------------------------------------------------
            RefGen Fixtures
'''

@pytest.fixture(scope="module")
def Zm5bFGS():
    if cf.test.force.RefGen:
        tools.del_dataset('RefGen', 'Zm5bFGS', force=True)
    if not tools.available_datasets('RefGen', 'Zm5bFGS'):
        # We have to build it
        gff = os.path.expanduser(
            os.path.join(
                cf.options.testdir,
                'raw', 'RefGen', 'ZmB73_5b_FGS.gff.gz'
            )
        )
        # This is stupid and necessary because pytables wont let me open
        # more than one table
        co.RefGen.from_gff(
            gff, 'Zm5bFGS', 'Maize 5b Filtered Gene Set', '5b', 'Zea Mays'
        )
    return co.RefGen('Zm5bFGS')

@pytest.fixture(scope="module")
def Zm5bFGS_duplicate(Zm5bFGS):
    # Build a refgen over an already built refgen 
    # We have to build it
    gff = os.path.expanduser(
        os.path.join(
            cf.options.testdir,
            'raw', 'RefGen', 'ZmB73_5b_FGS.gff.gz'
        )
    )
    # This is stupid and necessary because pytables wont let me open
    # more than one table
    co.RefGen.from_gff(
        gff, 'Zm5bFGS', 'Maize 5b Filtered Gene Set', '5b', 'Zea Mays'
    )
    return co.RefGen('Zm5bFGS')




''' -------------------------------------------------------------------------
            COB Fixtures
'''

@pytest.fixture(scope="module")
def ZmRNASeqTissueAtlas(Zm5bFGS):
    if cf.test.force.COB:
        print('Rebuilding ZmRNASeqTissueAtlas')
        tools.del_dataset('COB', 'ZmRNASeqTissueAtlas', force=True)
        tools.del_dataset('Expr', 'ZmRNASeqTissueAtlas', force=True)
    if not tools.available_datasets('Expr', 'ZmRNASeqTissueAtlas'):
        # Build it
        return co.COB.from_table(
            os.path.join(cf.options.testdir,
                'raw', 'Expr', 'RNASEQ',
                'MaizeRNASeqTissue.tsv.bz2',
            ),
            'ZmRNASeqTissueAtlas',
            'Maize RNASeq Tissue Atlas Network, Sekhon 2013, PLoS ONE',
            Zm5bFGS,
            rawtype='RNASEQ',
            max_gene_missing_data=0.3,
            max_accession_missing_data=0.08,
            min_single_sample_expr=1,
            min_expr=0.001,
            quantile=False,
            max_val=300,
            dry_run=False
        )
    else:
        return co.COB('ZmRNASeqTissueAtlas')

@pytest.fixture(scope="module")
def ZmRoot(Zm5bFGS):
    if cf.test.force.COB:
        tools.del_dataset('Expr','ZmRoot',force=True)
    if not tools.available_datasets('Expr','ZmRoot'):
        return co.COB.from_table(
            os.path.join(
                cf.options.testdir,
                'raw','Expr',
                'RNASEQ','ROOTFPKM.tsv.gz'
            ),
            'ZmRoot',
            'Maize Root Network',
            Zm5bFGS,
            rawtype='RNASEQ',
            max_gene_missing_data=0.3,
            max_accession_missing_data=0.08,
            min_single_sample_expr=1,
            min_expr=0.001,
            quantile=False,
            max_val=300
        )
    else:
        return co.COB('ZmRoot')

@pytest.fixture(scope="module")
def ZmSAM(Zm5bFGS):
    if cf.test.force.COB:
        tools.del_dataset('Expr','ZmSAM',force=True)
    if not tools.available_datasets('Expr','ZmSAM'):
        return co.COB.from_table(
            os.path.join(
                cf.options.testdir,
                'raw','Expr','RNASEQ',
                'TranscriptomeProfiling_B73_Atlas_SAM_FGS_LiLin_20140316.txt.gz'
            ),
            'ZmSAM',
            'Maize Root Network',
            Zm5bFGS,
            rawtype='RNASEQ',
            max_gene_missing_data=0.4,
            min_expr=0.1,
            quantile=False,
            dry_run=False,
            max_val=250
        )
    else:
        return co.COB('ZmSAM')

@pytest.fixture(scope="module")
def ZmSAM2(Zm5bFGS):
    if cf.test.force.COB:
        tools.del_dataset('Expr','ZmSAM2',force=True)
    if not tools.available_datasets('Expr','ZmSAM2'):
        return co.COB.from_table(
            os.path.join(
                cf.options.testdir,
                'raw','Expr','RNASEQ',
                'TranscriptomeProfiling_B73_Atlas_SAM_FGS_LiLin_20140316.txt.gz'
            ),
            'ZmSAM2',
            'Maize Root Network, but loose',
            Zm5bFGS,
            rawtype='RNASEQ',
            max_gene_missing_data=0.4,
            min_single_sample_expr=1,
            min_expr=0.01,
            quantile=False,
            dry_run=False,
            max_val=250
        )
    else:
        return co.COB('ZmSAM2')


@pytest.fixture(scope="module")
def ZmPAN2(Zm5bFGS):
    if cf.test.force.COB:
        tools.del_dataset('Expr','ZmPAN2',force=True)
    if not tools.available_datasets('Expr','ZmPAN2'):
        return co.COB.from_table(
            os.path.join(
                cf.options.testdir,
                'raw','Expr','RNASEQ',
                'PANGenomeFPKM.txt.gz'
            ),
            'ZmPAN2',
            'Maize Root Network but extra loose',
            Zm5bFGS,
            rawtype='RNASEQ',
            max_gene_missing_data=0.4,
            min_single_sample_expr=1,
            min_expr=0.01,
            quantile=False,
            dry_run=False,
            sep=',',
            max_val=300
        )
    else:
        return co.COB('ZmPAN2')


@pytest.fixture(scope="module")
def ZmPAN(Zm5bFGS):
    if cf.test.force.COB:
        tools.del_dataset('Expr','ZmPAN',force=True)
    if not tools.available_datasets('Expr','ZmPAN'):
        return co.COB.from_table(
            os.path.join(
                cf.options.testdir,
                'raw','Expr','RNASEQ',
                'PANGenomeFPKM.txt.gz'
            ),
            'ZmPAN',
            'Maize Root Network',
            Zm5bFGS,
            rawtype='RNASEQ',
            max_gene_missing_data=0.4,
            min_expr=1,
            quantile=False,
            dry_run=False,
            sep=',',
            max_val=300
        )
    else:
        return co.COB('ZmPAN')


''' -------------------------------------------------------------------------
            GWAS Fixtures
'''
@pytest.fixture(scope='module')
def testGWAS(testRefGen):
    if cf.test.force.Ontology:
        tools.del_dataset('GWAS','testGWAS',force=True)
    df = pd.DataFrame({
        'Trait' : ['a','a','b','b'],
        'CHR' : ['chr1','chr2','chr3','chr4'],
        'POS' : [100,200,300,400],
        'Start' : [100,200,300,400],
        'End' : [1000,20000,3000,4000],
        'id' : ['snp1','snp2','snp3','snp4'],
        'pval' : [0.05,0.05,0.01,0.01]
    }) 
    gwas = co.GWAS.from_DataFrame(
        df,
        'testGWAS',
        'Test GWAS Dataset',
        testRefGen,
        chr_col='CHR',
        pos_col='POS',
        id_col='id',
        term_col='Trait'
    ) 
    return gwas

@pytest.fixture(scope='module')
def ZmWallace(Zm5bFGS):
    if cf.test.force.Ontology:
        tools.del_dataset('GWAS','ZmWallace',force=True)
    if not tools.available_datasets('GWAS','ZmWallace'):
        # Grab path the csv
        csv = os.path.join(
            cf.options.testdir,
            'raw','GWAS','WallacePLoSGenet',
            'Wallace_etal_2014_PLoSGenet_GWAS_hits-150112.txt.gz'
        )
        # Define our reference geneome
        df = pd.DataFrame.from_csv(csv,index_col=None,sep='\t')
        # Import class from dataframe
        gwas  = co.GWAS.from_DataFrame(
            df, 'ZmWallace', 'Wallace PLoS ONE Dataset.',
            Zm5bFGS,
            term_col='trait', chr_col='chr', pos_col='pos'
        )
        return gwas
    else:
        return co.GWAS('ZmWallace')


@pytest.fixture(scope="module")
def ZmIonome(Zm5bFGS):
        # Delete the old dataset
    if cf.test.force.Ontology:
        tools.del_dataset('GWAS','ZmIonome',force=True)
    if not tools.available_datasets('GWAS','ZmIonome'):
        # Grab path the csv
        csv = os.path.join(
            cf.options.testdir,
            'raw','GWAS','Ionome',
            'sigGWASsnpsCombinedIterations.longhorn.allLoc.csv.gz'
        )
        # Define our reference geneome
        df = pd.DataFrame.from_csv(csv,index_col=None)
        # Import class from dataframe
        IONS  = co.GWAS.from_DataFrame(
            df,'ZmIonome','Maize Ionome',
            Zm5bFGS,
            term_col='el',chr_col='chr',pos_col='pos'
        )
        # Get rid of pesky Cobalt
        IONS.del_term('Co59')
        # I guess we need a test in here too
        return IONS
    else:
        return co.GWAS('ZmIonome')


'''----------------------------------------------------------------------------
    GOnt Fixtures
----------------------------------------------------------------------------'''

@pytest.fixture(scope='module')
def testGO(Zm5bFGS):
    if cf.test.force.Ontology:
        tools.del_dataset('GOnt','TestGO',force=True)
    if not tools.available_datasets('GOnt','TestGO'):
        obo = os.path.join(
            cf.options.testdir,
            'raw','GOnt','go.test.obo'
        )
        gene_map_file = os.path.join(
            cf.options.testdir,
            'raw','GOnt','go.test.tsv'
        )
        return co.GOnt.from_obo(
           obo, gene_map_file, 'TestGO',
           'Test GO', Zm5bFGS
        )
    else:
        return co.GOnt('testGO')



@pytest.fixture(scope="module")
def testZmGO(Zm5bFGS):
    if cf.test.force.Ontology:
        tools.del_dataset('GOnt','ZmGO',force=True)
    if not tools.available_datasets('GOnt','ZmGO'):
        obo = os.path.join(
            cf.options.testdir,
            'raw','GOnt','go.obo.gz'
        )
        gene_map_file = os.path.join(
            cf.options.testdir,
            'raw','GOnt','zm_go.tsv.gz'
        )
        return co.GOnt.from_obo(
           obo, gene_map_file, 'ZmGO',
           'Maize Gene Ontology', Zm5bFGS
        )
    else:
        return co.GOnt('ZmGO')

