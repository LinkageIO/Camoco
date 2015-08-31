'''
    COB Tests
'''

import pytest
import camoco as co
import itertools

import numpy as np

from scipy.misc import comb



def test_shape(ZmRNASeqTissueAtlas):
    cob = ZmRNASeqTissueAtlas
    assert len(cob.coex) == comb(cob.num_genes(), 2)

def test_coexpress_concordance(ZmRNASeqTissueAtlas):
    cob = ZmRNASeqTissueAtlas
    for a, b in itertools.combinations(
            [cob.refgen.random_gene() for x in range(50)], 2
        ):
        assert (
            abs(
                cob.coexpression(a, b).score \
                - cob._coex_concordance(a, b)) < 0.001
            ) or ( np.isnan(cob.coexpression(a, b).score) 
                and np.isnan(cob._coex_concordance(a, b))
            )
        dis_dif = abs(cob.coexpression(a, b).distance - abs(a-b))
        assert np.isnan(dis_dif) or dis_dif < 0.001

def test_num_neighbors_equals_degree(ZmRNASeqTissueAtlas):
    cob = ZmRNASeqTissueAtlas
    random_gene = cob.refgen.random_gene()
    assert len(cob.neighbors(random_gene)) \
        == cob.global_degree(random_gene)
