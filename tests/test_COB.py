'''
    COB Tests
'''

import pytest
import camoco as co

def test_build(ZmRNASeqTissueAtlas):
    assert ZmRNASeqTissueAtlas

def test_AtSeed(AtSeed):
    assert AtSeed

def test_AtGen(AtGen):
    assert AtGen

def test_AtRoot(AtRoot):
    assert AtRoot

def test_AtLeaf(AtLeaf):
    assert AtLeaf
