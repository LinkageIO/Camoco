import pytest
import os

from camoco.Config import cf

import camoco as co

@pytest.fixture
def Zm5bFGS():
    if cf['test'].getboolean('force'):
        co.del_dataset('RefGen','Zm5bFGS',safe=False)
    if not co.available_datasets('RefGen','Zm5bFGS'):
        # We have to build it
        gff = os.path.expanduser(
            os.path.join(
                cf['options']['testdir'],
                'raw','RefGen','ZmB73_5b_FGS.gff.gz'
            )
        )
        # This is stupid and necessary because pytables wont let me open
        # more than one table
        return co.RefGen.from_gff(
            gff,'Zm5bFGS','Maize 5b Filtered Gene Set','5b','Zea Mays'
        )
    return co.RefGen('Zm5bFGS')

@pytest.fixture
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
        return co.RefGen.from_gff(
            gff,'AtTair10','Tair 10','10','Arabidopsis'
        )
    else:
        return co.RefGen('AtTair10')

@pytest.fixture
def ZmRNASeqTissueAtlas(Zm5bFGS):
    if cf['test'].getboolean('force'):
        co.del_dataset('COB','ZmRNASeqTissueAtlas',safe=False)
        co.del_dataset('Expr','ZmRNASeqTissueAtlas',safe=False)
    if not co.available_datasets('COB','ZmRNASeqTissueAtlas'):
        # Build it 
        return co.COB.from_table(
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
            max_val=300,
            dry_run=True
        )
    else:
        return co.COB('ZmTissueAtlas')

@pytest.fixture
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
        return co.COB.from_DataFrame(
            SeedFam.series_matrix(
                keepfile=os.path.join(
                    cf.get('options','testdir'),
                    'raw','GSE','SeedKeep.tsv'
                )
            ),
            'AtSeed','Arabidopsis Seed',
            AtTair10,
            rawtype='MICROARRAY',
            quantile=True
    
        )
    else:
        return co.COB('AtSeed')


@pytest.fixture
def AtGen(AtTair10):
    if cf['test'].getboolean('force'):
        co.del_dataset('Expr','AtGen',safe=False)
    if not co.available_datasets('COB','AtGen'):
        General = ['GSE18975','GSE39384','GSE19271','GSE5632','GSE39385',
                'GSE5630','GSE15617','GSE5617','GSE5686','GSE2473',
                'GSE5633','GSE5620','GSE5628','GSE5624',
                'GSE5626','GSE5621','GSE5622','GSE5623','GSE5625','GSE5688']
        GenFam = sum(
            [co.Family.from_file(
                os.path.join(
                    cf.get('options','testdir'),
                    'raw','GSE','{}_family.soft'.format(x)
                )
            )
            for x in General ]
        )
        #GenFam.to_keepfile("GenKeep.tsv")
        return co.COB.from_DataFrame(
            GenFam.series_matrix(
                keepfile=os.path.join(
                    cf.get('options','testdir'),
                    'raw','GSE','GenKeep.tsv'
                )
            ),
            'AtGen','Arab General',
            AtTair10,
            rawtype='MICROARRAY',
            quantile=True
        )
    else:
        return co.COB('AtGen')

@pytest.fixture
def AtLeaf(AtTair10):
    if cf['test'].getboolean('force'):
        co.del_dataset('Expr','AtLeaf',safe=False)
    if not co.available_datasets('COB','AtLeaf'):
        Leaf = ['GSE14578','GSE5630','GSE13739', #'GSE26199',
                'GSE5686','GSE5615','GSE5620','GSE5628',
                'GSE5624','GSE5626','GSE5621','GSE5622',
                'GSE5623','GSE5625','GSE5688']
        LeafFam = sum(
            [co.Family.from_file(
                os.path.join(
                    cf.get('options','testdir'),
                    'raw','GSE','{}_family.soft'.format(x)
                )
            )
            for x in Leaf ]
        )
        #LeafFam.to_keepfile("LeafKeep.tsv",keep_hint="lea")
        return co.COB.from_DataFrame(
            LeafFam.series_matrix(
                keepfile=os.path.join(
                    cf.get('options','testdir'),
                    'raw','GSE','LeafKeep.tsv'
                )
            ),
            'AtLeaf','Arabidopsis Leaf',
            AtTair10,
            rawtype='MICROARRAY',
            max_gene_missing_data=0.3,
            min_expr=0.01,
            quantile=True,
        )
    else:
        return co.COB('AtLeaf')

@pytest.fixture
def AtRoot(AtTair10):
    if cf['test'].getboolean('force'):
        co.del_dataset('Expr','AtRoot',safe=False)
    if not co.available_datasets('COB','AtRoot'):
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
        return co.COB.from_DataFrame(
            RootFam.series_matrix(
                keepfile=os.path.join(
                    cf.get('options','testdir'),
                    'raw','GSE','RootKeep.tsv')
            ),
            'AtRoot','Arab Root',
            AtTair10,
            rawtype='MICROARRAY',
            quantile=True
        )
    else:
        return co.COB('AtRoot')




