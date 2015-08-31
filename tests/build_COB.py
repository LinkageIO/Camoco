import pytest

''' ------------------------------------------------------------------------
          Maize
'''

def test_build(ZmRNASeqTissueAtlas):
    assert ZmRNASeqTissueAtlas

def test_ZmRoot(ZmRoot):
    assert ZmRoot

def test_ZmSAM(ZmSAM):
    assert ZmSAM

def test_ZmPAN(ZmPAN):
    assert ZmPAN

''' ------------------------------------------------------------------------
        Arabidopsis
'''

def test_AtSeed(AtSeed):
    assert AtSeed

def test_AtGen(AtGen):
    assert AtGen

def test_AtRoot(AtRoot):
    assert AtRoot

def test_AtLeaf(AtLeaf):
    assert AtLeaf


