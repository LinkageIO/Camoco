import pytest

'''
    Init tests
'''

def build_ZmWallace(ZmWallace):
    assert ZmWallace 

def build_ZmIonome(ZmIonome):
    assert ZmIonome

def build_AtLeafIonome(AtLeafIonome):
    assert AtLeafIonome
    
def build_AtSeedIonome(AtSeedIonome):
    assert AtSeedIonome

def build_AtLeafHydroIonome(AtLeafHydroIonome):
    assert AtLeafHydroIonome

def build_AtRootHydroIonome(AtRootHydroIonome):
    assert AtRootHydroIonome

def build_AtIonome(AtLeafIonome,AtSeedIonome,
    AtLeafHydroIonome,AtRootHydroIonome):
    assert AtRootHydroIonome
    assert AtLeafHydroIonome
    assert AtSeedIonome
    assert AtLeafIonome
