'''
    COB Tests
'''

import pytest
import camoco as co

from scipy.misc import comb

''' ------------------------------------------------------------------------
          Maize
'''

def test_build(ZmRNASeqTissueAtlas):
    assert ZmRNASeqTissueAtlas

def test_ZmRoot(ZmRoot):
    assert ZmRoot

def test_ZmSAM(ZmSAM):
    assert ZmSAM

def test_ZmPAN(ZmPAN):
    assert ZmPAN

''' ------------------------------------------------------------------------
        Arabidopsis
'''

def test_AtSeed(AtSeed):
    assert AtSeed

def test_AtGen(AtGen):
    assert AtGen

def test_AtRoot(AtRoot):
    assert AtRoot

def test_AtLeaf(AtLeaf):
    assert AtLeaf


def test_shape(ZmRNASeqTissueAtlas):
    cob = ZmRNASeqTissueAtlas
    assert len(self.cob.coex) == comb(self.cob.num_genes(), 2)

def test_coexpress_concordance(ZmRNASeqTissueAtlas):
    cob = ZmRNASeqTissueAtlas
    for a, b in itertools.combinations(
            [self.cob.refgen.random_gene() for x in range(50)], 2
        ):
        assert abs(
            self.cob.coexpression(a, b).score \
            - self.cob._coex_concordance(a, b)
        ) < 0.001
        dis_dif = abs(self.cob.coexpression(a, b).distance - abs(a-b))
        assert np.isnan(dis_dif) or dis_dif < 0.001

def test_num_neighbors_equals_degree(ZmRNASeqTissueAtlas):
    cob = ZmRNASeqTissueAtlas
    random_gene = self.cob.refgen.random_gene()
    assert len(
        self.cob.neighbors(random_gene)
        == self.cob.global_degree(random_gene)
    )

