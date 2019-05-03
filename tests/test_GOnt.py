import pytest

# GO Term: Biological Process
bp_term_id = 'GO:0008150'


def test_GO(testGO):
    assert testGO

def test_len_of_GO(testGO):
    assert len(testGO) == 10

def test_go_term_propagation(testGO):
    assert len(testGO['GO:0000001'].loci) == 15

def test_go_namespace_in_attrs(testGO):
    assert testGO['GO:0000001'].attrs['namespace'] == 'biological_process'
    assert testGO['GO:0000002'].attrs['namespace'] == 'biological_process'
    assert testGO['GO:0000003'].attrs['namespace'] == 'biological_process'
    assert testGO['GO:0000004'].attrs['namespace'] == 'molecular_function'
    assert testGO['GO:0000005'].attrs['namespace'] == 'molecular_function' 
    assert testGO['GO:0000006'].attrs['namespace'] == 'molecular_function'
    assert testGO['GO:0000007'].attrs['namespace'] == 'molecular_function'
    assert testGO['GO:0000008'].attrs['namespace'] == 'molecular_function'
    assert testGO['GO:0000009'].attrs['namespace'] == 'molecular_function'
    assert testGO['GO:0000010'].attrs['namespace'] == 'biological_process'

def test_is_a_propagation(testGO):
    # Looks like:
    #       5
    #      /
    # 6   1--4
    #  \ / \  \ 
    #   2   3--8-10
    #  /       \
    # 7         9
    assert 'GO:0000001' in testGO['GO:0000002'].is_a
    assert 'GO:0000001' in testGO['GO:0000003'].is_a
    assert 'GO:0000001' in testGO['GO:0000004'].is_a
    assert 'GO:0000001' in testGO['GO:0000005'].is_a
    assert 'GO:0000002' in testGO['GO:0000006'].is_a
    assert 'GO:0000002' in testGO['GO:0000007'].is_a
    assert 'GO:0000003' in testGO['GO:0000008'].is_a
    assert 'GO:0000004' in testGO['GO:0000008'].is_a
    assert 'GO:0000008' in testGO['GO:0000009'].is_a
    assert 'GO:0000008' in testGO['GO:0000010'].is_a

def test_get_item(testZmGO):
    assert testZmGO[bp_term_id]

def test_biological_process_has_genes(testZmGO):
    assert len(testZmGO['GO:0008150'].loci) != 0

