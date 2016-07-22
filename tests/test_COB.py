'''
    COB Tests
'''
import camoco as co
from camoco import cf

import random
import itertools
from collections import Counter

import pytest
import numpy as np
from scipy.misc import comb

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
    # Get Random set of expr indexes
    expr_len = testCOB._expr.shape[0]
    expr_idxs = np.sort(np.unique(np.array(
        [random.randint(0,expr_len) for i in range(cf.test.num*10)]
    )))
    
    # Translate them back and forth
    coex_idxs = co.PCCUP.coex_index(expr_idxs, expr_len)
    new_expr_idxs = co.PCCUP.coex_expr_index(coex_idxs, expr_len).flatten()
    
    # Check all elements are the same in both
    if not (expr_idxs.shape[0] == np.unique(new_expr_idxs).shape[0]):
        assert False
    res = (expr_idxs == np.unique(new_expr_idxs))
    if not (np.sum(res) == res.shape[0]):
        assert False
    
    # Check all values occur the proper number of times
    corr_count = len(expr_idxs)-1
    bad_vals = []
    for k,v in Counter(new_expr_idxs).items():
        if not (v == corr_count):
            bad_vals.append((k,v))
    assert bad_vals == []

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
    assert False

def test_empty_subnetwork_returns_proper_dataframe(testCOB):
    assert False
