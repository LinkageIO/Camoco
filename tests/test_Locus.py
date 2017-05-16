import pytest

from itertools import chain
from camoco import Locus
from camoco.Config import cf

@pytest.fixture
def simple_Locus():
    return Locus(1,100,200) 

@pytest.fixture
def LocusX():
    return Locus(1,100,200) 

@pytest.fixture
def LocusY():
    return Locus(1,300,400) 

def test_locus_initialization(simple_Locus):
    # numeric chromosomes
    assert simple_Locus.chrom == '1'
    assert simple_Locus.start == 100
    assert simple_Locus.end == 200
    assert len(simple_Locus) == 101

def test_distance_between_loci():
    x = Locus(1,100,200) 
    y = Locus(1,300,400) 
    assert x - y == 99

def test_combine_loci(LocusX,LocusY):
    z = LocusX + LocusY
    assert len(z) == 301

def test_candidate_vs_bootstrap_length(testRefGen,testGWAS):
    Term = next(testGWAS.iter_terms())
    snps = Term.effective_loci(window_size=50000)
    candidates = testRefGen.candidate_genes(snps,chain=False)
    bootstraps = testRefGen.bootstrap_candidate_genes(snps,chain=False)
    # Make sure we are pulling out the same number of random genes for
    # Each locus
    for c,b in zip(candidates,bootstraps):
        assert len(c) == len(b)
    assert len(set(chain(*candidates))) == len(set(chain(*bootstraps)))

def test_generate_from_id(Zm5bFGS):
   random_gene = Zm5bFGS.random_gene()
   assert random_gene == Zm5bFGS[random_gene.id]


