#!/usr/bin/python3

import pytest
from camoco import Term
from camoco import Locus

@pytest.fixture
def testTerm():
    loci = [
        # Overlapping Loci, No windows
        Locus(1,100,500,score=0), Locus(1,400,700,score=5),
        # Loci with Overlapping windows
        Locus(2,100,200,window=100,score=0), 
        Locus(2,300,500,window=100,score=5),
        # SNPs with overlapping windows
        Locus(3,100,window=50,score=5), 
        Locus(3,200,window=50,score=0),
        # Three overlapping loci, one not
        Locus(4,100,window=80,score=1), 
        Locus(4,200,window=80,score=2),
        Locus(4,300,window=80,score=3), 
        Locus(4,400,window=10,score=4) # <- one not
    ]
    return Term('test',desc='hello',loci=loci,attr1=True,attr2=False)


def test_term_init(testTerm):
    assert testTerm.id == 'test'
    assert testTerm.desc == 'hello'
    assert testTerm['attr1'] == True
    assert testTerm['attr2'] == False

def test_add_Locus(testTerm):
    new_locus = Locus(6,100)
    testTerm.add_locus(new_locus)
    assert new_locus in testTerm.loci
    testTerm.loci.remove(new_locus)

def test_term_len(testTerm):
    assert len(testTerm) == len(testTerm.loci)

def test_effective_loci(testTerm):
    assert len(testTerm.effective_loci()) == 5 

def test_effective_loci_custom_windoe(testTerm):
    assert len(testTerm.effective_loci(window_size=150)) == 4

def test_effective_loci_lens(testTerm):
    assert list(map(len,testTerm.effective_loci())) == [600,400,100,200,1]

def test_strongest_loci(testTerm):
    assert list(
        map(lambda x:x.start, testTerm.strongest_loci('score',lowest=False))
    ) == [400,300,100,300,400]
