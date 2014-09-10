#!/usr/bin/python
import numpy as np

class Genome(object):   
    def __init__(self,id,chroms=list()):
        self.id = id
        self.chroms = chroms 
    def add_chromosone(self,chrom):
        ''' adds a chromosome to itself '''
        self.chroms.append(chrom)
    def rChrom(self):
        ''' returns a random chromosome object from genome '''
        rindex = np.random.randint(0,len(self.chroms))
        return self.chroms[rindex]
    def rLocus(self,length):
        ''' returns a random QTL of specified lengths '''
        return self.rChrom().rLocus(length)
    def rSNP(self):
        ''' returns a random 'SNP' from the genome '''
        return self.rChrom().rSNP()
    def __repr__(self):
        return "\n".join(map(str,self.chroms))
 
