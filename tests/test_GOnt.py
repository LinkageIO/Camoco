import pytest

# GO Term: Biological Process
bp_term_id = 'GO:0008150'


def test_GO(TestGO):
    assert TestGO

def test_len_of_GO(TestGO):
    assert len(TestGO) == 10

def test_go_term_propagation(TestGO):
    assert len(TestGO['GO:0000001'].loci) == 15

def test_get_item(ZmGO):
    assert ZmGO[bp_term_id]

def test_biological_process_has_genes(ZmGO):
    assert len(ZmGO['GO:0008150'].loci) != 0

