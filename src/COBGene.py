#!/usr/bin/python

from COBLocus import COBLocus

class COBGene(COBLocus):
    def __init__(self,GrameneID,gene_build='5b'):
        self.GrameneID = GrameneID
        # Create a dummy locus
        super(COBGene,self).__init__(0,0,0,gene_build)
        (self.ID,
         self.chromosome,
         self.chromo_start,
         self.chromo_end) = self.query('''SELECT ID, chromosome, chromo_start, chromo_end
            FROM genes WHERE GrameneID = '{}' and build = '{}' '''.format(self.GrameneID,self.gene_build)).iloc[0]
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


def from_GrameneID(GrameneID):
    ''' Creates a COBGene object from a gramene ID '''
    
