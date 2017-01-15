#!/usr/bin/python3

import pytest
import pandas as pd
import camoco as co

def test_init(testGWAS):
    assert len(testGWAS) == 2

def test_getitem(testGWAS):
    assert testGWAS['a'].id == 'a'

def test_attrs(testGWAS):
    for snp in testGWAS['a'].loci:
        assert snp['pval'] == '0.05'  
    for snp in testGWAS['b'].loci:
        assert snp['pval'] == '0.01'  


def test_fromTerms(testRefGen):
    '''
        Test GWAS creation from terms
    '''
    pass

def test_fromDataFrame(testRefGen):
    '''
        Test GWAS creation from DataFrame
    '''
    co.del_dataset('GWAS','testGWAS',safe=False)
    df = pd.DataFrame({
        'Trait' : ['a','a','b','b'],
        'CHR' : ['chr1','chr2','chr3','chr4'],
        'POS' : [100,200,300,400],
        'Start' : [100,200,300,400],
        'End' : [1000,20000,3000,4000],
        'id' : ['snp1','snp2','snp3','snp4'],
        'pval' : [0.05,0.05,0.01,0.01]
    }) 

    gwas = co.GWAS.from_DataFrame(
        df,
        'testGWAS',
        'Test GWAS Dataset',
        testRefGen,
        chr_col='CHR',
        pos_col='POS',
        id_col='id',
        term_col='Trait'
    ) 
    assert len(gwas) == 2
