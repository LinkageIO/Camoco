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
    i,m,l = (tpl) # values are index, exprMatrix, and lock
    total = ((m.shape[0]**2)-m.shape[0])/2 # total number of calculations we need to do
    left = ((m.shape[0]-i)**2-(m.shape[0]-i))/2 # How many we have left
    percent_done = math.floor((1-left/total)*100) # percent complete based on i and m
    if percent_done % 10 == 0: # report at 10% intervals
        l.acquire() # We need to lock because we are running in parallel
        print("\t\t{} {}% Complete".format(time.ctime(),percent_done))
        l.release()
    return [pearsonr(m[i,:],m[j,:])[0] for j in arange(i+1,len(m))] 

class COBDatasetBuilder(COBDatabaseBuilder):
    def __init__(self,name,description='',FPKM=True,gene_build='5b'):
        super(COBDatasetBuilder,self).__init__()
        # Warn us if there is already a Dataset named this
        assert self.fetchone("SELECT COUNT(*) FROM datasets WHERE name = '{}'",name)[0] == 0,\
            self.log("There is already a dataset named '{}'. Please choose something else.",name) 
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
            manager = multiprocessing.Manager()
            lock = manager.Lock()
            pool = multiprocessing.Pool()
            scores = np.array(list(itertools.chain(*pool.imap(_PCCUp,[ (i,self.expr,lock) for i in arange(self.num_genes)]))))
            pool.close()
            pool.join()
            return scores

    def dat(self,score_cutoff=None,transform=np.arctanh,normalize=True):
        tbl = pd.DataFrame(list(itertools.combinations([gene.GrameneID for gene in self.genes],2)),columns=['gene_a','gene_b'])
        scores = self.coex() 
        assert len(scores) == len(tbl)
        if transform:
            scores[scores == 1] = 0.999999
            scores = transform(scores)
        if normalize:
            scores = (scores -scores.mean()) / scores.std()
        tbl['scores' ] = scores
        # If there is a cutoff, filter it out
        if score_cutoff:
            tbl = tbl[tbl.scores > score_cutoff]
        return tbl
   
    def compare_to_dat(self,filename,sep="\t",score_cutoff=3):
        ''' Compare the number of significant edges with a DAT file '''
        ref_dat = pd.read_table(filename,sep=sep,names=['gene_a','gene_b','scores'])
        dat = self.dat(score_cutoff=score_cutoff) 
        # Get a list of genes in the dataset
        genes = set(dat.gene_a).union(set(dat.gene_b))
        for gene in genes:
            # Compare the degree of each gene with the reference dat
            ref_degree = len(ref_dat[(ref_dat.gene_a == gene) | (ref_dat.gene_b == gene)])
            degree     = len(dat[(dat.gene_a == gene)|(dat.gene_b == gene)])
            if ref_degree != degree:
                self.log("{} degree didn't match! ref: {} vs {}",gene,ref_degree,degree)
                print(ref_dat[(ref_dat.gene_a == gene) | (ref_dat.gene_b == gene)])
                print("----------------------------------------------------")
                print(dat[(dat.gene_a == gene)|(dat.gene_b == gene)])
            else:
                self.log("{} matches!".format(gene))

    def to_Dat(self,filename,UseGramene=True,**kw):
        with open(filename,'w') as FOUT:
            for gene in enumerate([x.ID for x in self.genes]):
                print("")
                for i,j,r in chunk:
                    print("{}\t{}\t{}".format(i,j,r),file=FOUT)

    def __str__(self):
        return '''
            COB Dataset Builder: {}
            desc: '{}'
            {} Genes
            {} Accessions
        '''.format(self.name,self.description,len(self.genes),len(self.accessions))

    def __repr__(self):
        return str(self)
