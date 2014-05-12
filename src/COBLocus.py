#!/usr/bin/python

from COBDatabase import COBDatabase
import pandas as pd
import numpy as np

class COBLocus(COBDatabase):
    def __init__(self,chromosome,chromo_start,chromo_end,gene_build='5b'):
        # Make sure we arent going nutso
        self.chromosome   = chromosome
        self.chromo_start = chromo_start
        self.chromo_end   = chromo_end
        self.gene_build   = gene_build
        super(COBLocus,self).__init__()

    def genes(self,gene_build=None):
        ''' returns genes within the locus '''
        if not gene_build:
           gene_build = self.gene_build 
        genes = self.query('''SELECT * FROM genes WHERE chromosome = {}
            AND chromo_start >= {}
            AND chromo_end <={}
            AND build = '{}'
            ORDER BY chromo_start''', self.chromosome, self.chromo_start, self.chromo_end, gene_build)
        return genes
    def upstream_genes(self,limit=100):
        ''' returns the X amount of upstream genes '''
        genes = self.query('''SELECT * from genes WHERE chromosome = {}
            AND chromo_end < {}
            AND build = '{}'
            ORDER BY chromo_start DESC
            LIMIT {}''',self.chromosome,self.chromo_start,self.gene_build,int(limit))      
        return genes   
    
    def downstream_genes(self,limit=100):
        ''' returns the X amount of upstream genes '''
        genes = self.query(''' SELECT * FROM genes WHERE chromosome = {}
            AND chromo_start > {} 
            AND build = '{}'
            ORDER BY chromo_start
            LIMIT {}''',self.chromosome,self.chromo_end,self.gene_build,int(limit))
        return genes
    def flanking_genes(self,limit=100):
        return pd.concat([
            self.upstream_genes(limit=np.ceil(limit/2)),
            self.downstream_genes(limit=np.ceil(limit/2))
        ])
    def __contains__(self,locus):
        if (locus.chromosome == self.chromosome and
               (( locus.chromo_start >= self.chromo_start and locus.chromo_start <= self.chromo_end)
               or(locus.chromo_end   <= self.chromo_end   and locus.chromo_end   >= self.chromo_start)
            )):
            return True
        else:
            return False
    def __len__(self):
        return self.chromo_end - self.chromo_start

    def __str__(self):
        return "Locus Object: {".format(str(self.id))
    def __repr__(self):
        return str(self.id)

    
