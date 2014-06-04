#!/usr/bin/python

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
    def rQTL(self,length):
        ''' returns a random QTL of specified lengths '''
        return self.rChrom().rQTL(length)
    def rLocus(self):
        ''' returns a random 'SNP' from the genome '''
        return self.rChrom().rLocus()

    def __repr__(self):
        return "\n".join(map(str,self.chroms))
 
