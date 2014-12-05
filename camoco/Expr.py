#! /usr/bin/env python3

from camoco.Camoco import Camoco
from camoco.RefGen import RefGen

import pandas as pd
import numpy as np
import matplotlib.pylab as plt

class Expr(Camoco):
    ''' A representation of gene expression. '''
    def __init__(self,name,description=None,basedir='~/.camoco'):
        super().__init__(name=name,description=description,type='Expr',basedir=basedir) 
        self._create_tables()
        if self.refgen:
            self.refgen = RefGen(self.refgen)

    def __repr__(self):
        pass
    def __str__(self):
        pass
    def summary(self):
        pass
   
    @classmethod
    def from_table(cls,filename,name,description,refgen,rawtype=None,basedir='~/.camoco',sep='\t',normalize=True,quality_control=True,**kwargs):
        # we are all pandas on the inside O.O
        tbl = pd.read_table(filename,sep=sep)
        self = cls(name=name,description=description,basedir=basedir)
        self._global("refgen",refgen.name)
        if rawtype is None:
            self.log('WARNING: not passing in a rawtype makes downstream normalization hard...')
        self._global('rawtype',rawtype)
        self.refgen = refgen 
        # put raw values into the database
        self.log('Importing Raw Expression Values')
        cur = self.db.cursor()
        try:
            # lets be safe about this
            cur.execute('BEGIN TRANSACTION')
            self.log('Adding accessions...')
            cur.executemany("INSERT OR IGNORE INTO accessions (name) VALUES (?)",
                [(a,) for a in tbl.columns]) 
            self.log('... imported {} accessions',len(tbl.columns))
            self.log('Importing genes ...')
            cur.executemany('''
                INSERT OR REPLACE INTO raw_expression (gene, accession, value) VALUES (?,?,?)''',
                [(accession,gene.upper(),float(value)) for (gene,accession),value in tbl.unstack().iteritems()] 
                #              ^ All genes are uppercase in database 
            )
            self.log('... imported {} raw gene values',len(tbl.index))
            cur.execute('END TRANSACTION')
        except Exception as e:
            cur.execute('ROLLBACK') 
        if quality_control:
            self.log('Performing Quality Control on genes')
            self.quality_control(**kwargs)
        if normalize:
            self.log('Performing Normalization')
            self.normalize()
        return self

    def normalize(self,method=None):
        ''' evaluates qc expression data and re-enters 
            normaized data into database '''
        df = self.expr(raw=False,long=True)
        if method is not None:
            method = method
        elif self.rawtype.upper() == 'RNASEQ':
            method = np.arctanh
        elif self.rawtype.upper() == 'MICROARRAY':
            method = np.log2
        else:
            raise ValueError('Could not guess correct normalization for {}, pass in function through method argument.'.format(self.rawtype))
        self._global('normalization',method.__name__)
        df.value = method(df.value)
        cur = self.db.cursor()
        try:
            # lets be safe about this
            cur.execute('BEGIN TRANSACTION')
            self.log('Importing normalized values ...')
            cur.executemany('''
                INSERT OR REPLACE INTO expression (gene, accession, value) VALUES (?,?,?)''',
                [(accession,gene.upper(),float(value)) for (gene,accession),value in df.unstack().iteritems()] 
                #              ^ All genes are uppercase in database 
            )
            self.log('Imported {} values',df.shape[0]*df.shape[1])
            cur.execute('END TRANSACTION')
        except Exception as e:
            cur.execute('ROLLBACK') 

    def quality_control(self,min_expr=1,max_gene_missing_data=0.2,min_single_sample_expr=5, \
        max_accession_missing_data=0.5,membership=None):
        ''' Sets quality control flag for all values in expression table '''        
        # remember how we set the flags
        self._global('qc_min_expr',min_expr)
        self._global('qc_max_gene_missing_data',max_gene_missing_data)
        self._global('qc_min_single_sample_expr',min_single_sample_expr)
        self._global('qc_max_accession_missing_data',max_accession_missing_data)
        self._global('qc_membership',str(membership))
        # retrieve raw data as a data frame
        df = self.expr(raw=True,long=False)
        self.log('Starting set: {} genes {} accessions'.format(len(df.index),len(df.columns)))
        # ------ Membership test
        if not membership:
            membership = self.refgen
        self.log("Filtering out genes not in {}",membership)
        df = df[[x in membership for x in df.index]]
        self.log('Kept: {} genes {} accessions'.format(len(df.index),len(df.columns)))
        # -----------------------------------------
        # Set minimum FPKM threshold
        self.log("Filtering out expression values lower than {}",min_expr)
        df_flt = df.copy()
        df_flt[df < min_expr] = np.nan
        df = df_flt
        self.log('Kept: {} genes {} accessions'.format(len(df.index),len(df.columns)))
        # -----------------------------------------
        # filter out genes with too much missing data
        self.log("Filtering out genes with > {} missing data",max_gene_missing_data)
        df = df.loc[df.apply(lambda x : ((sum(np.isnan(x))) < len(x)*max_gene_missing_data),axis=1),:] 
        self.log('Kept: {} genes {} accessions'.format(len(df.index),len(df.columns)))
        # -----------------------------------------
        # Filter out accession with too much missing data
        self.log("Filtering out accessions with > {} missing data",max_accession_missing_data)
        accession_mask = df.apply(lambda x : ((sum(np.isnan(x))) / len(x)),axis=0)
        for i,percent in enumerate(accession_mask):
            if percent > max_accession_missing_data:
                self.log("\tRemoved: {}: missing {} data",df.columns[i],percent*100)
        df = df.loc[:,np.logical_not(accession_mask > max_accession_missing_data)] 
        self.log('Kept: {} genes {} accessions'.format(len(df.index),len(df.columns)))
        # -----------------------------------------
        # filter out genes which do not meet a minimum expr threshold in at least one sample
        self.log("Filtering out genes which do not have one sample above {}",min_single_sample_expr)
        df = df[df.apply(lambda x: any(x >= min_single_sample_expr),axis=1)]
        self.log('Kept: {} genes {} accessions'.format(len(df.index),len(df.columns)))
        cur = self.db.cursor()
        try:
            # lets be safe about this
            cur.execute('BEGIN TRANSACTION')
            self.log('Importing qc genes ...')
            cur.executemany('''
                INSERT OR REPLACE INTO expression (gene, accession, value) VALUES (?,?,?)''',
                [(accession,gene.upper(),float(value)) for (gene,accession),value in df.unstack().iteritems()] 
                #              ^ All genes are uppercase in database 
            )
            cur.execute('END TRANSACTION')
        except Exception as e:
            cur.execute('ROLLBACK') 

    def expr(self,genes=None,accessions=None,long=False,raw=False):
        ''' returns expression data '''
        # Set up dynamic query
        tbl = 'raw_expression' if raw else 'expression'
        gene_filter = " WHERE gene in ('{}')".format("','".join(genes)) if genes is not None else ""
        accession_filter = " AND accession in ('{}')".format("','".join(accessions)) if accessions is not None else ""
        query = 'SELECT * FROM {} {} {};'.format(tbl,gene_filter,accession_filter)
        # pull and create data frame
        df = pd.DataFrame(
            self.db.cursor().execute(query).fetchall(),
            columns = ['gene','accession','value']
        )
        if long is False:
            df = df.pivot('gene','accession','value')
        return df

    def quantile(self):
        pass

    def plot_value_hist(self,groupby='accession',raw=False,bins=50,figsize=(16,16),title='',log=False):
        ''' Plots Value histograms on one of the expression matrix axis'''
        for group,df in self.expr(long=True,raw=raw).groupby(groupby):
            self.log('Plotting values for {}',group)
            plt.clf()
            plt.hist(
                list(filter(lambda x: not np.isnan(x),df.value)),
                bins=bins,
                log=log
            )
            plt.title(group+title)
            plt.xlabel('Expression')
            plt.ylabel('Frequency')
            plt.savefig(
                filename="{}{}_VALUES.png".format(group,title),
                figsize=figsize
            )

    def _create_tables(self):
        cur = self.db.cursor()
        cur.execute('PRAGMA page_size = 1024;') 
        cur.execute('PRAGMA cache_size = 100000;') 
        cur.execute(''' 
            CREATE TABLE IF NOT EXISTS accessions (
                name TEXT,
                type TEXT,
                description TEXT
            );
            CREATE TABLE IF NOT EXISTS raw_expression (
                gene TEXT,
                accession TEXT,
                value REAL
            );
            CREATE TABLE IF NOT EXISTS expression (
                gene TEXT,
                accession TEXT,
                value REAL
            );
        ''')
 
    def _build_indices(self):
        pass
