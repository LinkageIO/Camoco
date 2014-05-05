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

class progress_bar(object):
    def __init__(self,manager,interval=10): 
        self.interval = interval
        self.lock = manager.Lock()
        self.old_per = manager.Value(float,0.0)
    def update(self,cur_per):
        cur_per = math.floor(cur_per)
        if cur_per > self.old_per.value and cur_per % self.interval == 0:
            self.lock.acquire()
            print("\t\t{} {}% Complete".format(time.ctime(),cur_per))
            self.old_per.set(cur_per)
            self.lock.release()

def _PCCUp(tpl):
    ''' Returns a tuple containing indices i,j and their Pearson Correlation Coef '''
    i,m,pb = (tpl) # values are index, exprMatrix, and lock
    total = ((m.shape[0]**2)-m.shape[0])/2 # total number of calculations we need to do
    left = ((m.shape[0]-i)**2-(m.shape[0]-i))/2 # How many we have left
    percent_done = (1-left/total)*100 # percent complete based on i and m
    pb.update(percent_done)
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

    @property
    def num_genes(self):
        return len(self.genes)

    def coex(self,multi=True,UseGramene=False,DatasetID=0):
        if not multi: 
            scores = itertools.imap(_PCCUp,( (i,self.genes,self.expr,UseGramene) for i in arange(self.num_genes)))
        else:
            progress = progress_bar(multiprocessing.Manager())
            pool = multiprocessing.Pool()
            scores = np.array(list(itertools.chain(*pool.imap(_PCCUp,[(i,self.expr,progress) for i in arange(self.num_genes)]))))
            pool.close()
            pool.join()
            return scores

    def profile(self,gene_names):
        '''return the expression profile based on gene name '''
        matches = [ i for i,x in enumerate(self.genes) if x.GrameneID in gene_names ]
        return self.expr[matches,:]
       
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
        ref_genes = set(ref_dat.gene_a).union(set(ref_dat.gene_b))
        assert False
        for gene in genes:
            # Compare the degree of each gene with the reference dat
            ref_edges = ref_dat[(ref_dat.gene_a == gene) | (ref_dat.gene_b == gene)]
            edges     = dat[(dat.gene_a == gene)|(dat.gene_b == gene)]
            ref_neighbors = set(ref_edges.gene_a).union(set(ref_edges.gene_b))
            neighbors     = set(edges.gene_a).union(set(edges.gene_b))
            if len(ref_edges) > len(edges): # the number of edges is the degree
                self.log("{} reg degree bigger! ref: {} vs {}",gene,len(ref_edges),len(edges))
                print(ref_neighbors - neighbors)
                print("----------------------------------------------------")
            if len(ref_edges) < len(edges): # the number of edges is the degree
                self.log("{} degree bigger! ref: {} vs {}",gene,len(ref_edges),len(edges))
                print(neighbors - ref_neighbors)
                print("----------------------------------------------------")
            else:
                self.log("{} matches!".format(gene))

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
