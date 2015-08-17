import pytest
import os

from camoco.Config import cf

import camoco as co

@pytest.fixture()
def Zm5bFGS():
    if not co.available_datasets('RefGen','Zm5bFGS'):
        # We have to build it
        gff = os.path.join(
            cf['options']['testdir'],
            'raw','ZmB73_5b_FGS.gff'
        )
        return co.RefGen.from_gff(
            gff,'Zm5bFGS','Maize 5b Filtered Gene Set','5b','Zea Mays'
        )
    return co.RefGen('Zm5bFGS')

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

