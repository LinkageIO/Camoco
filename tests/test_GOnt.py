import pytest

# GO Term: Biological Process
bp_term_id = 'GO:0008150'


def test_GO(TestGO):
    assert TestGO

def test_len_of_GO(TestGO):
    assert len(TestGO) == 10

def test_go_term_propagation(TestGO):
    assert len(TestGO['GO:0000001'].loci) == 15

def test_is_a_propagation(TestGO):
    # Looks like:
    #       5
    #      /
    # 6   1--4
    #  \ / \  \ 
    #   2   3--8-10
    #  /       \
    # 7         9
    assert 'GO:0000001' in TestGO['GO:0000002'].is_a
    assert 'GO:0000001' in TestGO['GO:0000003'].is_a
    assert 'GO:0000001' in TestGO['GO:0000004'].is_a
    assert 'GO:0000001' in TestGO['GO:0000005'].is_a
    assert 'GO:0000002' in TestGO['GO:0000006'].is_a
    assert 'GO:0000002' in TestGO['GO:0000007'].is_a
    assert 'GO:0000003' in TestGO['GO:0000008'].is_a
    assert 'GO:0000004' in TestGO['GO:0000008'].is_a
    assert 'GO:0000008' in TestGO['GO:0000009'].is_a
    assert 'GO:0000008' in TestGO['GO:0000010'].is_a

def test_get_item(ZmGO):
    assert ZmGO[bp_term_id]

def test_biological_process_has_genes(ZmGO):
    assert len(ZmGO['GO:0008150'].loci) != 0

