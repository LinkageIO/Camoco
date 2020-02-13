
import pytest

import minus80 as m80
import locuspocus as lp

from pathlib import Path

# Create a Loci object from the Zm5b Refgen

@pytest.fixture(scope='module')
def Zm5bFGS():
    if not m80.Tools.available('Loci','Zm5bFGS'):
        # Create the refgen object
        loci = lp.Loci.from_gff('Zm5bFGS','raw/RefGen/ZmB73_5b_FGS.gff.gz')
    else:
        # Load it
        loci = lp.Loci('Zm5bFGS')
    return loci


@pytest.fixture(scope='module')
def ZmTissue(Zm5bFGs):
    if not m80.Tools.available('CoexNet','ZmTissue'):
        # Create the refgen object
        loci = lp.Loci.from_DataFrame(
            'Zm5bFGS',
            'raw/RefGen/ZmB73_5b_FGS.gff.gz'
        )
    else:
        # Load it
        loci = lp.Loci('Zm5bFGS')
    return loci

