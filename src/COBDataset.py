#!/usr/bin/python
from __future__ import division,print_function

from COBDatabase import *
from COBGene import COBGene
from numpy import matrix,arcsinh,zeros,arange
from scipy.stats import pearsonr
from pandas import Series
import time
import math as math
import multiprocessing as multiprocessing
import itertools
import os as os

def _PCCUp(tpl):
    ''' Returns a tuple containing indices i,j and their Pearson Correlation Coef '''
    i,m = (tpl)
    return [pearsonr(m[i,:],m[j,:])[0] for j in arange(i+1,len(m))] 

class COBDataset(COBDatabase):
    def __init__(self,name,description, FPKM=True,gene_build='5b'):
        super(COBDataset,self).__init__()
        self.name = name
        self.description = description
        # rows of matrix are genes and columns are accessions
        self.expr = matrix([0])
        self.genes = Series()
        self.accessions = Series()
        self.FPKM = FPKM
        self.gene_build = gene_build
        # Database specific stuff
        self.id = None

    def from_csv(self,filename,sep="\t"):
        ''' Returns a COBDataset from a CSV '''
        df = pd.read_table(filename,sep=sep)
        self.from_DataFrame(df)

    def from_DataFrame(self,df):
        all_genes = Series([COBGene(GrameneID,gene_build=self.gene_build) for GrameneID in df.index])
        # Filter out genes which are not in Database
        valid_gene_indices = [x.ID is not None for x in all_genes]
        self.genes = all_genes[valid_gene_indices]
        self.accessions = Series(df.columns)
        self.expr = df[valid_gene_indices].as_matrix()
        # Log how many genes we dropped
        
        self.log("{} out of {} genes were kept for {}",len(self.expr),len(df),self.name)
        if self.FPKM:
            # FPKM values get the inverse hyperbolic sine transform
            self.expr = arcsinh(self.expr)

    @property
    def num_genes(self):
        return len(self.genes)

    def coex(self,multi=True,UseGramene=False,DatasetID=0):
        if not multi: 
            scores = itertools.imap(_PCCUp,( (i,self.genes,self.expr,UseGramene) for i in arange(self.num_genes)))
        else:
            pool = multiprocessing.Pool()
            scores = np.array(list(itertools.chain(*pool.imap(_PCCUp,[ (i,self.expr) for i in arange(self.num_genes)]))))
            pool.close()
            pool.join()
            return scores

    def dat(self):
        tbl = pd.DataFrame(list(itertools.combinations([gene.GrameneID for gene in self.genes],2)),columns=['gene_a','gene_b'])
        scores = self.coex() 
        assert len(scores) == len(tbl)
        tbl['scores' ] = scores
        return tbl
    

    def to_Dat(self,filename,UseGramene=True,**kw):
        with open(filename,'w') as FOUT:
            for gene in enumerate([x.ID for x in self.genes]):
                print("")
                for i,j,r in chunk:
                    print("{}\t{}\t{}".format(i,j,r),file=FOUT)

    def test(self):
        old_per = 0
        total_genes = (((len(self.genes)**2)/2)-len(self.genes))
        for i,j in enumerate(self.coex()):
            cur_per = math.floor((i / total_genes) * 100)
            if cur_per > old_per:
                old_per = cur_per
                print("{} - {}% complete ({} edges)".format(time.ctime(),cur_per,i))

    def __str__(self):
        return '''
            COB Dataset:
            {} Genes
            {} Accessions
        '''.format(len(self.genes),len(self.accessions))

    def __repr__(self):
        return str(self)
