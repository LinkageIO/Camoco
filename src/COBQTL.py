#!/usr/bin/python

from COBLocus import *

class COBQTL(COBLocus):
    def __init__(self,chrom,start,end,id=None):
        if id == None:
            self.id = "QTL-chr{}:{}-{}".format(chrom,start,end)
        else:
            self.id = id
        super(COBQTL,self).__init__(chrom,start,end,self.id)
    def __str__(self):
        return "{}".format(self.id)


