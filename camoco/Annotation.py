#!/usr/bin/python3

import os
import glob
import camoco as co
import pandas as pd
from .Camoco import Camoco

class RefGenFunc(Camoco):
    def __init__(self, refgen, description=None):
        if not isinstance(refgen, str):
            refgen = refgen.name
        super().__init__(refgen,type='RefGenFunc')
        self._global('refgen',refgen)
        self.refgen = co.RefGen(refgen)
        self._create_tables() 
    
    def __getitem__(self,item):
        # Build the query from all the genes provided
        if isinstance(item,(set,list)):
            ls = "{}".format("','".join([str(x) for x in item]))
            single = False
        else:
            ls = item
            single = True
        query = "SELECT * FROM func WHERE id IN ('{}');".format(ls)
        
        # Run the query and turn the result into a list of tuples
        cur = self.db.cursor()
        cur.execute(query)
        annotes = cur.fetchall()
        
        # If a list of genes was passed in, return a dictionary of lists
        if not single:
            res = {}
            for id,desc in annotes:
                if id in res:
                    res[id].append(desc)
                else:
                    res[id] = [desc]
        
        # Otherwise just return the list annotations
        else:
            res = []
            for id,desc in annotes:
                res.append(desc)
        return res
    
    def to_csv(self, filename=None, sep="\t"):
        '''
            Make a table of all functional annotations.
        '''
        # Find the default filename
        if filename == None:
            filename = self.name + '_func.tsv'
        
        # Pull them all from sqlite
        cur = self.db.cursor()
        cur.execute("SELECT * FROM func;")
        
        # Used pandas to save it
        df = pd.DataFrame(cur.fetchall(),columns=['gene','desc']).set_index('gene')
        df.to_csv(filename,sep=sep)
    
    def add_table(self, filename, sep="\t", gene_col=0, skip_cols=None):
        ''' 
            Imports Annotation relationships from a csv file. By default will
            assume gene names are first column
        '''
        # import from file, assume right now that in correct order
        tbl = pd.read_table(filename,sep=sep,dtype=object)
        idx_name = tbl.columns[gene_col]
        tbl[idx_name] = tbl[idx_name].str.upper()
        
        # Drop columns if we need to
        if skip_cols is not None:
            # removing certain columns
            tbl.drop(tbl.columns[skip_cols],axis=1,inplace=True)
        
        # Get rid of any genes not in the refence genome
        refcur = self.refgen.db.cursor()
        refcur.execute('SELECT id FROM genes;')
        rm = set(tbl[idx_name]) - set([id[0] for id in refcur.fetchall()])
        tbl.drop(rm,axis=0,inplace=True)
        del rm, refcur
        
        # One Annotation per row, drop the nulls and duplicates
        tbl = pd.melt(tbl,id_vars=idx_name,var_name='col',value_name='desc')
        tbl.drop('col',axis=1,inplace=True)
        tbl.dropna(axis=0,inplace=True)
        tbl.drop_duplicates(inplace=True)
        
        # Run the transaction to throw them in there
        cur = self.db.cursor()
        try:
            cur.execute('BEGIN TRANSACTION')
            cur.executemany(
                'INSERT INTO func VALUES (?,?)'
                ,tbl.itertuples(index=False))
            cur.execute('END TRANSACTION')
        
        except Exception as e:
            self.log("import failed: {}",e)
            cur.execute('ROLLBACK')
        
        # Make sure the indices are built
        self._build_indices()
    
    @classmethod
    def create(cls, refgen, description=None):
        if not isinstance(refgen, str):
            refgen = refgen.name
        self = super().create(refgen, description, type='RefGenFunc')
        return self
    
    @classmethod
    def from_table(cls, filename, refgen, description=None, 
        sep="\t", gene_col=0, skip_cols=None):
        self = cls.create(refgen, description)
        self.add_table(filename, sep=sep, gene_col=gene_col, skip_cols=skip_cols)
        
    # We also have the groundwork for a table of ortholog annotations
    # Currently there is no interface to either add or access it
    # We may add it in the future though
    def _create_tables(self):
        cur = self.db.cursor()
        cur.execute('''
            CREATE TABLE IF NOT EXISTS func (
                id TEXT,
                desc TEXT,
                UNIQUE(id,desc) ON CONFLICT IGNORE
            );
            CREATE TABLE IF NOT EXISTS ortho_func (
                id TEXT,
                desc TEXT,
                UNIQUE(id,desc) ON CONFLICT IGNORE
            );
        ''')
    
    def _build_indices(self):
        cur = self.db.cursor()
        cur.execute('CREATE INDEX IF NOT EXISTS id ON func(id)')
        cur.execute('CREATE INDEX IF NOT EXISTS id ON ortho_func(id)')

class GWASData(Camoco):
    def __init__(self, gwas):
        if not isinstance(gwas, str):
            gwas = gwas.name
        super().__init__(gwas,type='GWASData')
        self._global('gwas',gwas)
        self._create_tables()
    
    # Function to get data using any of the parameters it is normaly queried by
    def get_data(self, gene=None, cob=None, term=None,
        windowSize=None, flankLimit=None):
        # Base section of query
        query = "SELECT * FROM gwas_data WHERE"
        
        # Throw arguments into a dictionary for effective looping
        args = {
            'gene':gene, 'cob':cob, 'term':term,
            'window':windowSize, 'flank':flankLimit}
        
        # For each argument, add a clase to the SQL query
        for k,v in args.items():
            if not(v is None):
                if isinstance(v,(set,list)):
                    ls = "{}".format("','".join([str(x) for x in v]))
                else:
                    ls = v
                query += " {} IN ('{}') AND".format(k,ls)
        
        # Peel off unneeded things at the end of the query, depending args
        if(query[-4:] == ' AND'):
            query = query[:-4]
        else:
            query = query[:-6]
        query += ';'
        
        # Run the query
        cur = self.db.cursor()
        cur.execute(query)
        
        # Throw the results into a DataFrame for conviniece
        return pd.DataFrame(cur.fetchall(),columns=[
            'gene','COB','Term','WindowSize','FlankLimit','score',
            'zscore','fdr','num_real','num_random','bs_mean','bs_std',
            'NumSNPs','NumBootstraps'
        ]).set_index('gene')
        
    def add_table(self, filename, sep='\t'):
        ''' 
            Imports GWAS Data from the output file.
        '''
        # Import from file keep and correctly order the needed columns
        tbl = pd.read_table(filename,sep=sep,dtype=object)
        tbl = tbl[['gene','COB','Term','WindowSize','FlankLimit','score',
            'zscore','fdr','num_real','num_random','bs_mean','bs_std',
            'NumSNPs','NumBootstraps']]
        
        # Run the transaction to put them in the db
        cur = self.db.cursor()
        cur.execute('BEGIN TRANSACTION')
        try:
            cur.executemany(
            'INSERT INTO gwas_data VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)'
            ,tbl.itertuples(index=False))
        except Exception as e:
            self.log("Import of '{}' failed: {}",filename, e)
            cur.execute('ROLLBACK')
        
        # Cleanup
        self._build_indices()
        cur.execute('END TRANSACTION')
    
    def add_dir(self, dir, sep='\t'):
        '''
            Wrapper to import all GWAS output from a directory.
        '''
        # Setup db connection
        cur = self.db.cursor()
        cur.execute('BEGIN TRANSACTION')
        
        # Find all the files in the folder related to our GWAS
        fileList = glob.glob(os.path.join(dir, ('*' + self.name +'*.[tc]sv')))

        # For each file, import to the db
        for fn in fileList:
            # Import from file keep and correctly order the needed columns
            tbl = pd.read_table(fn,sep=sep,dtype=object)
            tbl = tbl[['gene','COB','Term','WindowSize','FlankLimit','score',
                'zscore','fdr','num_real','num_random','bs_mean','bs_std',
                'NumSNPs','NumBootstraps']]
            try:
                cur.executemany('''
                    INSERT INTO gwas_data VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)
                ''',tbl.itertuples(index=False))
                self.log('Imported file: {}',fn)
            except Exception as e:
                self.log("Import of '{}' failed: {}",fn, e)
                cur.execute('ROLLBACK')
        
        # Cleanup
        self._build_indices()
        cur.execute('END TRANSACTION')
    
    @classmethod
    def create(cls, gwas, description=None):
        if not isinstance(gwas, str):
            gwas = gwas.name
        self = super().create(gwas, description, type='GWASData')
        return self
    
    @classmethod
    def from_dir(cls, dir, gwas, description=None):
        self = cls.create(gwas, description)
        self.add_dir(dir)

    def _create_tables(self):
        cur = self.db.cursor()
        cur.execute('''
            CREATE TABLE IF NOT EXISTS gwas_data (
                gene TEXT,
                cob TEXT,
                term TEXT,
                window INTEGER,
                flank INTEGER,
                score REAL,
                zscore REAL,
                fdr REAL,
                num_real REAL,
                num_random REAL,
                bs_mean REAL,
                bs_std REAL,
                num_snps INTEGER,
                num_bootstraps INTEGER,
                UNIQUE(gene,cob,term,window,flank) ON CONFLICT REPLACE
            );
        ''')
    
    def _build_indices(self):
        cur = self.db.cursor()
        cur.execute('CREATE INDEX IF NOT EXISTS gene ON gwas_data(gene)')
