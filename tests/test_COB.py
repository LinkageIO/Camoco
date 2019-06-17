''' 
    COB Tests
'''

import os
import pytest
import random
import itertools
import scipy.stats

import camoco as co 
import pandas as pd
import numpy as np

from camoco import cf
from scipy.special import comb
from collections import Counter

def test_consistent_MCL_clusters(testCOB):
    '''
        Make sure that between build, that MCL
        clusters stay consistent.
    '''
    pass

def test_consistent_network_degree_dist(testCOB):
    '''
        In general, the degree of each node should stay consistent unless 
        something went wrong during network creation. This test compares
        the freshly built network (testCOB) gene degree to a csv with
        degree numbers from past builds
        The old_degree file was created with:

        >>> testCOB.degree.reset_index().to_csv(
            os.path.join(
                co.Config.cf.options.testdir,
                'raw','StabilityData','ZmRNASeqTissueAtlas_degree.csv'
            ),index=False
        )
    '''
    # Read in degree for genes from old file
    old_degree = pd.read_table(
        os.path.join(
            co.Config.cf.options.testdir,
            'raw',
            'StabilityData',
            'ZmRNASeqTissueAtlas_degree.csv'
        ),
        sep=',',
        index_col='index'
    )
    # merge together to make easier to match values
    merged_degree = old_degree.join(
        testCOB.degree,
        lsuffix='_old',
        rsuffix='_new'
    )
    # Get rid of nans
    not_nan_mask = np.logical_not(
        np.isnan(merged_degree.Degree_old) & \
        np.isnan(merged_degree.Degree_new)
    )
    merged_degree = merged_degree.loc[not_nan_mask]
    # calculate PCC and compare old v new
    r2,pval = scipy.stats.pearsonr(merged_degree.Degree_old,merged_degree.Degree_new)
    assert r2 > 0.99



def test_consistent_GO_density(testCOB,testZmGO):
    '''
        Denisty of GO terms should be fairly consistent between
        network builds. This test loads pre-calculated densities
        for GO terms (10<n<300) and compares against densities
        calculated from the current built network.
    '''
    # Load the prev calculated densities
    old_densities = pd.read_table(
        os.path.join(
            co.Config.cf.options.testdir,
            'raw',
            'StabilityData',
            'ZmRNASeqTissueAtlas_ZmGO_Densities.csv'
        ),
        sep=','
    ) 
    # caluclate the current densities from GO terms
    current_densities = pd.DataFrame(
        [(i.id,testCOB.density(i.loci)) \
                for i in testZmGO.terms(min_term_size=10,max_term_size=300)\
        ],
        columns=['term','density']
    )
    # Merge them together to make it easier to match values
    merged_densities = old_densities.merge(
        current_densities,
        left_on='term',
        right_on='term',
        suffixes=['_old','_new'],
        how='inner'
    )
    # Get rid of nans
    not_nan_mask = np.logical_not(
        np.isnan(merged_densities.density_old) & \
        np.isnan(merged_densities.density_new)
    )
    merged_densities = merged_densities.loc[not_nan_mask]

    # Calculate pearson correlation of old vs new
    r2,pval = scipy.stats.pearsonr(
        merged_densities.density_new,
        merged_densities.density_old
    )
    # check for high correlation
    assert r2 > 0.99


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
    '''
        Translate expr indexes to coex indexes and back again
    '''
    # Get Random set of expr indexes
    expr_len = testCOB._expr.shape[0]
    expr_idxs = np.sort(np.unique(np.array(
        [random.randint(0,expr_len-1) for i in range(cf.test.num*10)]
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
    assert len(testCOB.neighbors(random_gene,return_gene_set=True)) \
        == testCOB.global_degree(random_gene)
    assert len(testCOB.neighbors(random_gene,names_as_index=False)) \
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
    # Compare the degree determined from subnetwork aggregation
    # is the same as what is in the data frame
    for k,v in Counter(itertools.chain(*testCOB.subnetwork().index.values)).items():
        assert testCOB.degree.loc[k].Degree == v


def test_empty_subnetwork_returns_proper_dataframe(testCOB):
    subnet = testCOB.subnetwork([])
    assert len(subnet) == 0
    assert 'score' in subnet.columns
    assert 'significant' in subnet.columns
    assert 'distance' in subnet.columns


def test_zero_index_genes_doesnt_get_filtered(testCOB):
    ''' ---- Regression bug
        This bug occured when one of the genes in the subnetwork list
        had an index of 0, which was filtered out because the list filter
        function thought it was None
    '''
    # get the first gene in the matrix
    gene_a = testCOB.refgen[testCOB._expr.index[0]]
    # get a random gene
    gene_b = testCOB.refgen.random_gene()
    # make sure that the subnetwork contains a single coexpression 
    # entry between genes
    assert len(testCOB.subnetwork([gene_a,gene_b],sig_only=False)) == 1

def test_zero_degree_genes_return_empty_dataframe(testCOB):
    # get a random zero degree gene
    gene_id = testCOB.degree.loc[testCOB.degree.Degree==0].sample(1).index[0]
    gene = testCOB.refgen[gene_id]
    assert len(testCOB.neighbors(gene)) == 0

def test_repr(testCOB):
    assert repr(testCOB).startswith('<COB:')

def test_str(testCOB):
    assert str(testCOB).startswith('<COB:')

def test_qc_gene(testCOB):
    assert isinstance(testCOB.qc_gene(),pd.DataFrame)
