#!/usr/bin/python3

from .Ontology import Ontology
from .Term import Term
from .Locus import Locus


class GWAS(Ontology):
    '''
        Ontology extension for GWAS. This class implements when you
        do not have only loci found in a reference genome.
    '''
    def __init__(self,name):
        super().__init__(name,type='GWAS')


    def __getitem__(self, id):
        ''' retrieve a term by id '''
        try:
            (id, desc) = self.db.cursor().execute(
                'SELECT * from terms WHERE id = ?', (id, )
            ).fetchone()
            term_loci = [
                Locus(chrom,start,end,id=id,) \
                for term, id, chrom, start, end, window  in \
                self.db.cursor().execute(''' 
                    SELECT * from term_loci WHERE 
                    term = ?
                ''',(id,))
            ]
            for locus in term_loci:
                for key,val in self.db.cursor().execute('''
                        SELECT key,val FROM loci_attr
                        WHERE term = ? AND id = ?
                        ''',(id,locus.id)).fetchall():
                    locus.attr[key] = val
            return Term(id, desc=desc, loci=term_loci)
        except TypeError as e: #Not in database
            raise e

    def add_term(self, term, overwrite=True):
        ''' 
            This will add a single term to the ontology.
            Extends the functionality of the Ontology class.
        '''
        cur = self.db.cursor()
        if overwrite:
            self.del_term(term.id)
        cur.execute('BEGIN TRANSACTION')
        # Add the term name and description
        cur.execute('''
            INSERT OR REPLACE INTO terms (id, desc)
            VALUES (?, ?)''', (term.id, term.desc)
        )
        # Add the term loci
        for locus in term.loci:
            cur.execute('''
                INSERT OR REPLACE INTO term_loci 
                (term, id, chrom, start, end, window)
                VALUES (?, ?, ?, ? ,? ,?);
                ''', (term.id, locus.id, locus.chrom,
                      locus.start, locus.end, locus.window)
            )
            cur.executemany('''
                INSERT OR REPLACE INTO loci_attr
                (term,id,key,val) VALUES (?,?,?,?);
            ''',[(term.id,locus.id,key,val) for key,val in locus.attr.items()])
        cur.execute('END TRANSACTION')

    def del_term(self,id):
        super().del_term(id)
        # Get rid of the loci_attr also
        self.db.cursor().execute(''' 
            DELETE FROM loci_attr WHERE term = ?
        ''',(id,))

    def _create_tables(self):
        super()._create_tables()
        # Add the loci table so it works with SNPs
        cur = self.db.cursor()
        cur.execute('''
            -- More concise to drop table and create new one.
            DROP TABLE term_loci;
            CREATE TABLE IF NOT EXISTS term_loci (
                term TEXT,
                -- Locus information
                id TEXT,
                chrom TEXT,
                start INT,
                end INT,
                window INT,
                -- Relationship here is b/w term and locus-id
                PRIMARY KEY(term,id)
            );
            -- Keep track of 
            CREATE TABLE loci_attr (
                -- Some loci may be associated with multiple terms, so term
                -- key is necessary
                term TEXT,
                id TEXT,
                key TEXT,
                val TEXT,
                PRIMARY KEY(term,id,key)
            );
        ''')

    
    ''' -----------------------------------------------------------------------
            Class Methods -- Factory Methods
    '''

    @classmethod
    def create(cls, name, description, refgen, type='GWAS'):
        return super().create(name, description, refgen, type='GWAS')

    @classmethod
    def from_DataFrame(cls, df, name, description, refgen,
            term_col='Term', chr_col=None, pos_col=None,
            start_col=None, end_col=None, id_col=None
            ):
        '''
            Import an ontology from a pandas dataframe.
            Groups by term_col, then iterates over the rows in the group.
            It adds each row as a locus to the term.
        '''
        self = cls.create(name, description, refgen)
        # group each trait by its name
        for term_id, df in df.groupby(term_col):
            term = Term(term_id)
            # we have a SNP
            if pos_col is not None:
                for i, row in df.iterrows():
                    # make sure there are no collisions with the Locus instance
                    # function names. This is hackey and I dont like it
                    kwargs = {
                        key:val for key, val in dict(row).items() \
                        if key not in Locus.__init__.__code__.co_varnames \
                        and key not in [chr_col,pos_col,start_col,end_col,id_col]
                    }
                    snp = Locus(
                        row[chr_col], int(row[pos_col]), int(row[pos_col]), 
                        gene_build=self.refgen.build, **kwargs
                    )
                    term.loci.add(snp)
            self.log("Importing {}", term)
            self.add_term(term)
        return self
