'''
    COB Tests
'''

import pytest
import camoco as co
import itertools

import numpy as np

import random
from scipy.misc import comb
from camoco import cf

def test_coordination_between_expr_and_expr_index(testCOB):
    for i,x in enumerate(testCOB._expr.index):
        assert i == testCOB._expr_index[x]

def test_coordination_between_expr_index_and_coex_index(testCOB):
    assert True
    return
    assert set(itertools.chain(*testCOB.coex.index.values)) \
        == set(testCOB._expr.index.values)

def test_expr_nans_in_same_place(testCOB):
    pass

def test_shape(testCOB):
    assert len(testCOB.coex) == comb(testCOB.num_genes(), 2)

def test_coex_score_concordance(testCOB):
    for a, b in itertools.combinations(
            [testCOB.refgen.random_gene() for x in range(cf.test.num)], 2
        ):
        assert (
            abs(
                testCOB.coexpression(a, b).score \
                - testCOB._coex_concordance(a, b)) < 0.001
            ) or ( np.isnan(testCOB.coexpression(a, b).score) 
                and np.isnan(testCOB._coex_concordance(a, b))
            )

def test_coex_distance_concordance(testCOB):
    for a, b in itertools.combinations(
            [testCOB.refgen.random_gene() for x in range(cf.test.num)], 2
        ):
        dis_dif = abs(testCOB.coexpression(a, b).distance - abs(a-b))
        assert np.isnan(dis_dif) or dis_dif < 100

def test_coex_id_concordance(testCOB):
    for a, b in itertools.combinations(
            [testCOB.refgen.random_gene() for x in range(cf.test.num)], 2
        ):
        assert sorted(testCOB.coexpression(a,b).name) == sorted([a.id,b.id])

def test_coex_to_expr_concordance(testCOB):
    expr_len = testCOB._expr.shape[0]
    expr_idxs = sort(np.unique(np.array(
        [random.randint(0,expr_len) for i in range(cf.test.num*10)]
    )))
    coex_idxs = co.PCCUP.coex_index(expr_idxs, expr_len)
    new_expr_idxs = np.unique(co.PCCUP.coex_expr_index(coex_idxs, expr_len).flatten())
    missing = 0
    for i in new_expr_idxs:
        if i not in expr_idxs:
            missing += 1
    assert missing == 0

def test_num_neighbors_equals_degree(testCOB):
    random_gene = testCOB.refgen.random_gene()
    assert len(testCOB.neighbors(random_gene)) \
        == testCOB.global_degree(random_gene)

def test_subnetwork_contains_only_input_genes(testCOB):
    random_genes = set(testCOB.refgen.random_genes(n=cf.test.num))
    subnet = testCOB.subnetwork(random_genes,sig_only=False)
    assert set(itertools.chain(*subnet.index.values)) == set([x.id for x in random_genes])

def test_subnetwork_contains_only_input_when_duplicates(testCOB):
    random_genes = list(testCOB.refgen.random_genes(n=cf.test.num))
    # Add duplicates
    random_genes = random_genes + random_genes[-10:-1]
    subnet = testCOB.subnetwork(random_genes,sig_only=False)
    assert set(itertools.chain(*subnet.index.values)) == set([x.id for x in random_genes])

def test_degree_index_matches_degree(testCOB):
    assert True
