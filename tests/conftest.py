import pytest
import os

from camoco.Config import cf

import camoco as co

@pytest.fixture()
def Zm5bFGS():
    if cf['test'].getboolean('force'):
        co.del_dataset('RefGen','Zm5bFGS',safe=False)
    if not co.available_datasets('RefGen','Zm5bFGS'):
        # We have to build it
        gff = os.path.join(
            cf['options']['testdir'],
            'raw','RefGen','ZmB73_5b_FGS.gff.gz'
        )
        co.RefGen.from_gff(
            gff,'Zm5bFGS','Maize 5b Filtered Gene Set','5b','Zea Mays'
        )
    return co.RefGen('Zm5bFGS')

@pytest.fixture()
def AtTair10():
    if cf['test'].getboolean('force'):
        co.del_dataset('RefGen','AtTair10',safe=False)
    if not co.available_datasets('RefGen','AtTair10'):
        gff = os.path.expanduser(
            os.path.join(
                cf['options']['testdir'],
                'raw','RefGen','TAIR10_GFF3_genes.gff.gz'
            )
        )
        co.RefGen.from_gff(
            gff,'AtTair10','Tair 10','10','Arabidopsis'
        )
    return co.RefGen('AtTair10')


@pytest.fixture()
def ZmRNASeqTissueAtlas(Zm5bFGS):
    if cf['test'].getboolean('force'):
        co.del_dataset('COB','ZmRNASeqTissueAtlas',safe=False)
        co.del_dataset('Expr','ZmRNASeqTissueAtlas',safe=False)
    if not co.available_datasets('COB','ZmRNASeqTissueAtlas'):
        # Build it 
        ZmRoot = co.COB.from_table(
            os.path.join(cf.get('options','testdir'),
                'raw', 'Expr', 'MaizeRNASeqTissue.tsv.gz',
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
            max_val=300
        )
    return co.COB('ZmTissueAtlas')

@pytest.fixture()
def AtSeed(AtTair10):
    if cf['test'].getboolean('force'):
        co.del_dataset('Expr','AtSeed',safe=False)
    if not co.available_datasets('COB','AtSeed'):
        Seed = ['GSE12404', #'GSE30223',
                'GSE1051', 'GSE11852', 'GSE5634']
        SeedFam = sum(
            [co.Family.from_file(
                os.path.join(
                    cf.get('options','testdir'),
                    'raw','GSE','{}_family.soft'.format(x)
                )
            )
            for x in Seed ]
        )
        #SeedFam.to_keepfile("SeedKeep.tsv",keep_hint='seed')
        AtSeed = co.COB.from_DataFrame(
            SeedFam.series_matrix(
                keepfile=os.path.join(
                    cf.get('options','testdir'),
                    'raw','GSE','SeedKeep.tsv'
                )
            ),
            'AtSeed','Arabidopsis Seed',
            co.RefGen('T10'),rawtype='MICROARRAY',
            quantile=True
    
        )
    return AtSeed


