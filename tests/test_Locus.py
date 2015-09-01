import pytest
from camoco import Locus
from camoco.Config import cf

@pytest.fixture
def simple_Locus():
    return Locus(1,100) 


def test_locus_initialization(simple_Locus):
    # numeric chromosomes
    assert simple_Locus


def test_candidate_vs_bootstrap_length(Zm5bFGS,ZmIonome):
    Term = ZmIonome[cf.test.term]
    snps = Term.effective_loci(window_size=50000)
    assert \
       len(Zm5bFGS.candidate_genes(snps)) \
       == len(Zm5bFGS.bootstrap_candidate_genes(snps))


def test_generate_from_id(Zm5bFGS):
   assert Zm5bFGS[cf.test.gene]
