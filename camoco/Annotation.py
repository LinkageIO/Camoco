#!/usr/bin/python

from camoco import Camoco
import camoco as co
import pandas as pd

class Annotation(Camoco):
    def __init__(self, name, description=None, refgen=None, basedir='~/.camoco'):
        super().__init__(name,description,type='Annotation',basedir=basedir)
        self._create_tables()
        try:
            self.refgen = co.RefGen(self.refgen)
        except NameError as e:
            # no refgen has been assigned yet, check if it was passed in  
            self._global('refgen',refgen.name)
            self.refgen = co.RefGen(self.refgen)

    def __getitem__(self,item):
        if isinstance(item,(set,list)):
            cur = self.db.cursor().execute('''
                SELECT * 
                from annotations WHERE gene IN ('{}')
            '''.format("','".join([gene.id for gene in item])))
            # need to check if the query was empty, if so return an empty data frame
            return pd.DataFrame(cur.fetchall(),columns=['gene','name','desc']).pivot('gene','name','desc')
        else:
            cur = self.db.cursor().execute('''
                SELECT * 
                from annotations WHERE gene = ?
            ''',(item.id,))
            return pd.Series(
                index = [x[0] for x in cur.getdescription()],
                data=cur.fetchone(),
                dtype=object
            )

    def add_csv(self, filename, sep="\t", skip_cols=None, gene_col=0):
        ''' Imports Annotation relationships from a csv file. By default will assume gene names 
            are first column'''
        # import from file, assume right now that in correct order
        table = pd.read_table(filename,sep=sep,dtype=object).fillna('-') 
        if skip_cols is not None:
            # removing certain columns
                table.drop(table.columns[skip_cols],axis=1,inplace=True)
        cur = self.db.cursor()
        try:
            cur.execute('BEGIN TRANSACTION')
            cur.executemany(''' 
                INSERT INTO annotations VALUES (?,?,?)
            ''',pd.melt(table,id_vars=table.columns[gene_col]).itertuples(index=False))
            cur.execute('END TRANSACTION')
        except Exception as e:
            self.log("import failed: {}",e)
            cur.execute('ROLLBACK')     
        self._build_indices()
    
    def _build_indices(self):
        cur = self.db.cursor()
        cur.execute('CREATE INDEX IF NOT EXISTS gene ON annotations(gene)')

    def _create_tables(self):
        cur = self.db.cursor()
        cur.execute('''
            CREATE TABLE IF NOT EXISTS annotations (
                gene TEXT,
                name TEXT,
                desc TEXT,
                UNIQUE(gene,name) ON CONFLICT REPLACE
            );
        ''')
        self._build_indices()


