#!/usr/bin/python

import numpy as np
from camoco.Locus import Locus

class Chrom(object):
    def __init__(self,id,length):
        self.id = id
        self.length = length

    def rLocus(self,length):
        ''' returns a random QTL within the chromosome '''
        start = np.random.randint(0,self.length)
        while start+length > self.length:
            start = np.random.randint(0,self.length)
        return Locus(chrom=self.id,start=start,end=start+length,id="rLocus-{}".format(length))

    def rSNP(self):
        ''' returns a random single nucleotide (position) from within the chromosome '''
        pos = np.random.randint(0,self.length)
        return Locus(chrom=self.id,start=pos,end=pos,id='rSNP-chr{}:{}'.format(self.id,pos))

    def __len__(self):
        return self.length

    def __str__(self):
        return "Chromosome {}: {} bp".format(self.id,self.length)

    def __repr__(self):
        return str(self)


