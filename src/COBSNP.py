#!/usr/bin/python

from COBLocus import *

class COBSNP(COBLocus):
    def __init__(self,chrom,pos,alt=None,ref=None,id=None):
        super(COBSNP,self).__init__(chrom,pos,pos,id)
        self.chromosome = chrom
        self.pos = pos
        self.id = id
        self.alt = alt
        self.ref = ref

    def nearby_genes(self,window=100000):
        downstream = np.floor(self.pos+(window/2))
        upstream   = np.floor(self.pos-(window/2))
        genes = self.query('''SELECT * FROM mzn_gene_loci WHERE chromosome = {}
            AND chromo_start > {}
            AND chromo_end < {}
            AND gene_build = '{}'
            ORDER BY chromo_start''', self.chromosome, upstream, downstream, self.gene_build)
        return genes

    def nearest_gene(self,window=100000):
        genes = self.query('''SELECT * FROM mzn_gene_loci WHERE chromosome = {}
            AND chromo_start < {}
            AND chromo_end > {}
            AND gene_build = '{}'
            ORDER BY chromo_start''', self.chromosome, self.pos, self.pos, self.gene_build)
        if genes.empty:
            downstream = np.floor(self.pos+(window/2))
            upstream   = np.floor(self.pos-(window/2))
            genes = self.query('''SELECT * FROM mzn_gene_loci WHERE chromosome = {}
                AND chromo_start > {}
                AND chromo_end < {}
                AND gene_build = '{}'
                ORDER BY chromo_start''', self.chromosome, upstream, downstream, self.gene_build)
            if len(genes) > 2:
                if sum(genes.chromo_end < self.pos) == len(genes.chromo_end):
                    # All the genes are before the SNP, return the last
                    return genes.tail(1)
                elif sum(genes.chromo_start > self.pos) == len(genes.chromo_start):
                    # All genes are before SNP, return the first
                    return genes.head(1)
                else:
                    downstream_index = ((genes.chromo_start > self.pos) * 1).idxmax()
                    upstream_index = downstream_index -1
                    return genes.ix[[upstream_index,downstream_index]]


