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
            # Iterate through loci and get attrs
            for locus in term_loci:
                for key,val in self.db.cursor().execute('''
                        SELECT key,val FROM loci_attr
                        WHERE term = ? AND id = ?
                        ''',(id,locus.id)).fetchall():
                    locus.attr[key] = val
            return Term(id, desc=desc, loci=term_loci)
        except TypeError as e: #Not in database
            raise e

    def add_term(self, term, cursor=None, overwrite=True):
        ''' 
            This will add a single term to the ontology.
            Extends the functionality of the Ontology class.
        '''
        try:
            if not cursor:
                # create a new cursor and initiate a transaction
                cur = self.db.cursor()
                cur.execute('BEGIN TRANSACTION')
            else:
                # Otherwise, assume that another transaction was initiated
                # perhaps by self.add_terms (notice the plurality)
                cur = cursor
            if overwrite:
                self.del_term(term.id,cursor=cur)
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
            if not cursor: 
                # Still assume that 
                cur.execute('END TRANSACTION')
        except Exception as e:
            cur.execute('ROLLBACK')
            raise e

    def del_term(self,term,cursor=None):
        super().del_term(term,cursor=cursor)

        if not isinstance(term, str):
            id = term.id
        else:
            id = term
        # Get rid of the loci_attr also
        if not cursor:
            cur = self.db.cursor()
        else:
            cur = cursor
        cur.execute(''' 
            DELETE FROM loci_attr WHERE term = ?
        ''',(id,))


    ''' -----------------------------------------------------------------------
            Internal Methods -- Factory Methods
    '''

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
            -- Keep track of loci attributes
            CREATE TABLE IF NOT EXISTS loci_attr (
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
        '''
            Create an empty GWAS dataset.

            Parameters
            ----------
            name : str
                Name of the GWAS dataset
            description : str
                Short description of the GWAS dataset
            refgen : camoco.RefGen object
                The corresponding RefGen object that is 
                affiliated with the GWAS dataset. (For 
                extracting candidate genes, etc)

            Returns
            -------
            A camoco.GWAS object

            Note
            ----
            See classmethods: from_DataFrame and from_terms
            for help building GWAS datasets from other common
            data types.

        '''
        return super().create(name, description, refgen, type='GWAS')

    @classmethod
    def from_DataFrame(cls, df, name, description, refgen,
            term_col='Term', chr_col='CHR', pos_col=None,
            start_col=None, end_col=None, id_col=None, 
            strongest_attr='pval', strongest_higher=True
            ):
        '''
            Import an GWAS dataset from a pandas dataframe.
            Groups by term_col, then iterates over the rows in the group.
            It adds each row as a locus to the term.


            Parameters
            ----------
            df : Pandas.DataFrame
                The data frame containing information for GWAS. Terms will
                be created by grouping the "term_col" column. Loci for that
                term will be created from the 'chr_col' and either the 'pos_col'
                OR the 'start_col' and 'end_col' depending on what is specified. 
                See docs for those parameters for more info.
            name : str
                The name of the GWAS dataset
            description : str
                The description of the GWAS dataset
            refgen : camoco.RefGen object
                The corresponding reference genome for the loci that are within
                the terms for the GWAS dataset

            Keyword Parameters
            ------------------
            term_col : str (default: 'Term')
                The rows in *df* will be grouped by when adding loci to the terms.
                Each different value in the term column represents a different GWAS
                experiment, i.e. that loci in those rows will be added 
            chr_col : str (default: 'CHR')
                This column designates the chromosome that the row is affiliated with.
           |pos_col : str (default: None)
           |    If the pos_column is designated, the GWAS is assumed to be SNPs (i.e. 
           |    the start end end position of the loci are the same). There is no need
           |    to include the start_col or end_col
           |start_col : str (default: None)
           |    Represents the start position of the loci thats associated with the trait.
           |    If the start_col is designated, the end_col must also be designated.
           |    When included, the loci are are assumed to be QTL or genomic spans
           |end_col : str (default: None)
               Must be included if the start_col is inlcuded. Represents the end
                position of the locus associated with the trait.
            id_col : str (default: None)
                Assign an id to the locus
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
        self.set_strongest(attr=strongest_attr,higher=strongest_higher)
        return self
