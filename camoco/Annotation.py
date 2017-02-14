#!/usr/bin/python3

import os
import glob
import camoco as co
import pandas as pd
from .Camoco import Camoco

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
        fileList = glob.glob(os.path.join(dir, ('*.[tc]sv')))

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
