#!/usr/bin/python

import pandas as pd
from COBLocus import COBLocus

class COBGene(COBLocus):
    def __init__(self,GrameneID,gene_build='5b'):
        self.GrameneID = GrameneID
        # Create a dummy locus
        super(COBGene,self).__init__(0,0,0,gene_build)
        try:
            (self.ID,
            self.chromosome,
            self.chromo_start,
            self.chromo_end) = self.fetchone('''SELECT ID, chromosome, chromo_start, chromo_end
                FROM genes WHERE GrameneID = '{}' and build = '{}' '''.format(self.GrameneID,self.gene_build))
        except (TypeError,ValueError) as e: 
            self.log("{} Not found for build {}".format(GrameneID,gene_build))
            self.ID = None

    @property
    def common_name(self):
        pass

    @property
    def arab_ortho(self):
        pass

    @property
    def go_terms(self):
        pass

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "COBGene: {}".format(self.GrameneID)

def import_gene_list(filename,delim="\n",gene_build='5a'):
    ''' This reads a file in line by line and creates a list of gene object
        the filter call is to get rid of genes not in the database, since the
        COBGene class returns None for genes it does not find in database '''
    return filter(lambda x:x, [COBGene(x,gene_build) for x in open(filename,'r').read().strip().split(delim)])

def gene_list(gene_list,gene_build='5a'):
    ''' This creates a list of COBGene objects from grmzm names '''
    return filter(lambda x:x, [COBGene(x,gene_build) for x in gene_list])
