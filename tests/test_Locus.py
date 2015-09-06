import pytest
from camoco import Locus
from camoco.Config import cf

@pytest.fixture
def simple_Locus():
    return Locus(1,100,200) 

def test_locus_initialization(simple_Locus):
    # numeric chromosomes
    assert simple_Locus.chrom is 1
    assert simple_Locus.start is 100
    assert simple_Locus.end is 200
    assert len(simple_Locus) == 100

def test_candidate_vs_bootstrap_length(testRefGen,testGWAS):
    Term = testGWAS[cf.test.term]
    snps = Term.effective_loci(window_size=50000)
    assert \
          len(testRefGen.candidate_genes(snps)) \
       == len(testRefGen.bootstrap_candidate_genes(snps))

def test_generate_from_id(Zm5bFGS):
   random_gene = Zm5bFGS.random_gene()
   assert random_gene == Zm5bFGS[random_gene.id]


