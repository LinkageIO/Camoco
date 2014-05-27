#!/usr/bin/python3


class Locus(object):
    def __init__(self, chrom, start, end, id=None ,gene_build='5b', organism='Zea'):
        self.chrom = int(chrom)
        self.start = int(start)
        self.end = int(end)
        self.id = id
        self.gene_build = gene_build
        self.organism = organism

    @property
    def stop(self):
        ''' because im an idiot ''' 
        return self.end

    def __contains__(self,locus):
        if (locus.chrom == self.chrom and
               (( locus.start >= self.start and locus.start <= self.end)
               or(locus.end   <= self.end   and locus.end   >= self.start)
            )):
            return True
        else:
            return False
    def __len__(self):
        ''' inclusive length of locus '''
        return self.stop - self.start + 1

    def __cmp__(self,snp):
        return self.start - snp-start

    def __lt__(self,snp):
        return self.start < snp.start    
    def __gt__(self,snp):
        return self.start > snp.start
    def __sub__(self,other):
        return self.start - other.start
    def __str__(self):
        return '''
            organism:{} 
            type: {}
            chromosome: {}
            start: {}
            end: {}'''.format(self.organism,self.__class__,self.chrom,self.start,self.stop)
    def __repr__(self):
        return str(self)

class SNP(Locus):
    def __init__(self, chrom, pos, id=None ,gene_build='5b', organism='Zea'):
        self.pos = int(pos)
        super(self.__class__,self).__init__(chrom,pos,pos,id,gene_build,organism)
       

class Gene(Locus):
    def __init__(self,chrom=None,start=None,end=None,strand=None,id=None,build='5b',organism='Zea'):
        super().__init__(chrom,start,end,id,build,organism) 
