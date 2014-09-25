#!/usr/bin/python

from camoco import Camoco
import camoco as co
import pandas as pd

class Annotation(Camoco):
    def __init__(self,name,description=None,basedir='~/.camoco'):
        super().__init__(name,description,type='Annotation',basedir=basedir)
        if self.refgen:
            self.refgen = co.RefGen(self.refgen)
        else:
            self.log('Refgen not assigned')

    def __getitem__(self,item):
        if isinstance(item,list):
            return self.db.cursor().execute('''
                SELECT * from annotations WHERE gene IN ('{}')
            '''.format("','".join(item))).fetchall()
        else:
            # return a tuple of annotations
            return self.db.cursor().execute('''
                SELECT * from annotations WHERE gene = ?''',(item,)).fetchone()

    @classmethod
    def from_csv(cls,name,description,refgen,filename,sep="\t"):
        self = cls(name,description)
        self._global('refgen',refgen.name)
        self._create_tables()
        # import from file, assume right now that in correct order
        table = pd.read_table(filename,sep=sep,dtype=object).fillna('-') 
        #table[['Start Position','End Position']] = table[['Start Position','End Position']].astype(int)
        cur = self.db.cursor()
        try:
            cur.execute('BEGIN TRANSACTION')
            cur.executemany(''' 
                INSERT INTO annotations VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
            ''',table.itertuples(index=False))
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
        cur.execute("PRAGMA page_size = 1024;")
        cur.execute("PRAGMA cache_size = 100000;")
        cur.execute('''
            CREATE TABLE IF NOT EXISTS annotations (
                gene TEXT UNIQUE,
                model TEXT,
                chrom TEXT,
                start INTEGER,
                end INTEGER,
                strand TEXT,
                maizeseq TEXT,
                bfgr TEXT,
                pfam TEXT,
                rice_orth TEXT,
                rice_annot TEXT,
                sorghum_orth TEXT,
                sorghum_annot TEXT,
                panicum_orth TEXT,
                panicum_annot TEXT,
                grassius_tf TEXT,
                mapman TEXT,
                classical TEXT,
                efp TEXT,
                functional_annot TEXT
            )
        ''')

