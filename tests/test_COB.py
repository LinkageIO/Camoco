'''
    COB Tests
'''

import pytest
import camoco as co
import itertools

import numpy as np

from scipy.misc import comb
from camoco import cf

def test_coordination_between_expr_and_expr_index(ZmRNASeqTissueAtlas):
    cob = ZmRNASeqTissueAtlas
    for i,x in enumerate(cob._expr.index):
        assert i == cob._expr_index[x]

def test_coordination_between_expr_index_and_coex_index(ZmRNASeqTissueAtlas):
    cob = ZmRNASeqTissueAtlas
    assert set(itertools.chain(*cob.coex.index.values)) == set(cob._expr.index.values)

def test_expr_nans_in_same_place(ZmRNASeqTissueAtlas):
    pass

def test_shape(ZmRNASeqTissueAtlas):
    cob = ZmRNASeqTissueAtlas
    assert len(cob.coex) == comb(cob.num_genes(), 2)

def test_coex_score_concordance(ZmRNASeqTissueAtlas):
    cob = ZmRNASeqTissueAtlas
    for a, b in itertools.combinations(
            [cob.refgen.random_gene() for x in range(cf.test.num)], 2
        ):
        assert (
            abs(
                cob.coexpression(a, b).score \
                - cob._coex_concordance(a, b)) < 0.001
            ) or ( np.isnan(cob.coexpression(a, b).score) 
                and np.isnan(cob._coex_concordance(a, b))
            )

def test_coex_distance_concordance(ZmRNASeqTissueAtlas):
    cob = ZmRNASeqTissueAtlas
    for a, b in itertools.combinations(
            [cob.refgen.random_gene() for x in range(cf.test.num)], 2
        ):
        dis_dif = abs(cob.coexpression(a, b).distance - abs(a-b))
        assert np.isnan(dis_dif) or dis_dif < 0.001

def test_coex_id_concordance(ZmRNASeqTissueAtlas):
    cob = ZmRNASeqTissueAtlas
    for a, b in itertools.combinations(
            [cob.refgen.random_gene() for x in range(cf.test.num)], 2
        ):
        assert cob.coexpression(a,b).name == tuple(sorted([a.id,b.id]))


def test_num_neighbors_equals_degree(ZmRNASeqTissueAtlas):
    cob = ZmRNASeqTissueAtlas
    random_gene = cob.refgen.random_gene()
    assert len(cob.neighbors(random_gene)) \
        == cob.global_degree(random_gene)

def test_subnetwork_contains_only_input_genes(ZmRoot):
    cob = ZmRoot
    random_genes = set(cob.refgen.random_genes(n=cf.test.num))
    subnet = cob.subnetwork(random_genes,sig_only=False)
    assert set(itertools.chain(*subnet.index.values)) == set([x.id for x in random_genes])
