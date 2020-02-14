
import pytest

import camoco as co
import pandas as pd
import minus80 as m80
import locuspocus as lp

from pathlib import Path

# Create a Loci object from the Zm5b Refgen

@pytest.fixture(scope='module')
def Zm5bFGS():
    if not m80.Tools.available('Loci','Zm5bFGS'):
        # Create the refgen object
        loci = lp.Loci.from_gff(
            'Zm5bFGS',
            'raw/RefGen/ZmB73_5b_FGS.gff.gz'
        )
    else:
        # Load it
        loci = lp.Loci('Zm5bFGS')
    return loci


@pytest.fixture(scope='module')
def ZmTissue(Zm5bFGs):
    if not m80.Tools.available('CoexNet','ZmTissue'):
        # Get the pd DataFrame
        df = pd.read_table(
            'raw/Expr/RNASEQ/MaizeRNASeqTissue.tsv.bz2'
        )
        # Create the refgen object
        ZmTissue = co.Coex.from_DataFrame(
            'ZmTissue',
            'Maize RNASeq Tissue Atlas Network, Sekhon 2013, PLoS ONE',
            df,
            Zm5bFGS,
            min_expr=0.001,
            max_locus_missing_data=0.3,
            max_accession_missing_data=0.08,
            min_single_accession_expr=1,
        )
    else:
        # Load it
        loci = lp.Loci('Zm5bFGS')
    return loci

