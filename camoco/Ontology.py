#!/usr/bin/python3

from .Camoco import Camoco
from .RefGen import RefGen
from .Locus import Locus
from .Term import Term

from pandas import DataFrame
from scipy.stats import hypergeom
from itertools import chain
from functools import lru_cache
from collections import OrderedDict

import sys
import numpy

class Ontology(Camoco):
    '''
        An Ontology is just a collection of terms. Each term is just a
        collection of genes. Sometimes terms are related or nested
        within each other, sometimes not. Simple enough.
        
        Parameters
        ----------
        name : unique identifier

        Returns
        -------
        An Ontology Object

    '''
    def __init__(self, name, type='Ontology'):
        super().__init__(name, type=type)
        if self.refgen:
            self.refgen = RefGen(self.refgen)

    def __len__(self):
        '''
            Return the number of non-empty terms
        '''
        return self.num_terms(min_term_size=1)

    def num_terms(self,min_term_size=0,max_term_size=10e10):
        '''
            Returns the number of terms in the Ontology
            within the min_term_size and max_term_size

            Parameters
            ----------
            min_term_size (default:0)
                The minimum number of loci associated with the term 
            max_term_size (default: 10e10)
                The maximum number of loci associated with the term

            Returns
            -------
            the number of terms that fit the criteria

        '''
        return self.db.cursor().execute(
            '''SELECT COUNT(*) FROM (
                SELECT DISTINCT(term) FROM term_loci 
                GROUP BY term 
                HAVING COUNT(term) >= ? 
                    AND  COUNT(term) <= ?
            );''',
            (min_term_size, max_term_size)
        ).fetchone()[0]

    @lru_cache(maxsize=131072)
    def __getitem__(self, id):
        ''' retrieve a term by id '''
        try:
            (id, desc) = self.db.cursor().execute(
                'SELECT * from terms WHERE id = ?', (id, )
            ).fetchone()
            term_loci = [
                self.refgen[gene_id] for gene_id, in self.db.cursor().execute(
                ''' SELECT id FROM term_loci WHERE term = ?''', (id, )
            ).fetchall()]
            term_attrs = {k:v for k,v in self.db.cursor().execute(
                ''' SELECT key,val FROM term_attrs WHERE term = ?''',(id,)         
                )
            }
            return Term(id, desc=desc, loci=term_loci,**term_attrs)
        except TypeError as e: # Not in database
            raise e

    def num_distinct_loci(self):
        return self.db.cursor().execute(
            'SELECT COUNT(DISTINCT(id)) FROM term_loci;'
        ).fetchone()[0]

    def iter_terms(self,min_term_size=0,max_term_size=10000000):
        '''
            Return a generator that iterates over each term in the ontology.
        '''
        terms = self.db.cursor().execute('''
            SELECT term from term_loci
            GROUP BY term
            HAVING COUNT(term) >= ?
                AND COUNT(term) <= ?
        ''',(min_term_size,max_term_size))
        for id, in terms:
            yield self[id]

    def terms(self):
        return list(self.iter_terms())

    def summary(self):
        return "Ontology:{} - desc: {} - contains {} terms for {}".format(
            self.name, self.description, len(self), self.refgen)

    def rand(self, n=1, min_term_size=1, max_term_size=100000):
        '''
            Return a random Term from the Ontology

            Parameters
            ----------
            n : int (default: 1)
                The number of random terms to return
            min_term_size : int (default: 1)
                The smallest acceptable term size
                i.e. the number of genes annotated to the term
            max_term_size : int (default: 100000)
                The largest acceptable term size
        '''
        cur = self.db.cursor()
        ids = cur.execute(''' 
            SELECT term FROM term_loci 
            GROUP BY term 
            HAVING COUNT(term) >= ?
                AND COUNT(term) <= ?
            ORDER BY RANDOM() 
            LIMIT ?;
        ''',(min_term_size,max_term_size,n)).fetchall()
        if len(ids) == 0:
            raise ValueError(
                'No Terms exists with this criteria '
                '{} < len(term) < {}:'.format(min_term_size,max_term_size)
            )
        terms = [self[id[0]] for id in ids]
        if len(terms) == 1:
            return terms[0]
        else:
            return terms

    def add_term(self, term, cursor=None, overwrite=False):
        ''' This will add a single term to the ontology

        Parameters
        ----------
        term : Term object
            The term object you wish to add.
        cursor : apsw cursor object
            A initialized cursor object, for batch operation. This will
            allow for adding many terms in one transaction as long as the 
            passed in cursor has executed the "BEGIN TRANSACTION" command.
        overwrite : bool
            Indication to delete any existing entry before writing'''
        if overwrite:
            self.del_term(term.id)
        if not cursor:
            cur = self.db.cursor()
            cur.execute('BEGIN TRANSACTION')
        else:
            cur = cursor

        # Add the term id and description
        cur.execute('''
            INSERT OR ABORT INTO terms (id, desc)
            VALUES (?, ?)''', (term.id, term.desc))

        # Add the term loci
        if term.loci:
            for locus in term.loci:
                cur.execute('''
                    INSERT OR ABORT INTO term_loci (term, id)
                    VALUES (?, ?)
                    ''', (term.id, locus.id))

        if not cursor:
            cur.execute('END TRANSACTION')

    def del_term(self, term, cursor=None):
        ''' This will delete a single term to the ontology

        Parameters
        ----------
        term : Term object or str
            The term object or id you wish to remove.
        cursor : apsw cursor object
            A initialized cursor object, for batch operation.'''

        if not cursor:
            cur = self.db.cursor()
            cur.execute('BEGIN TRANSACTION')
        else:
            cur = cursor

        if not isinstance(term, str):
            id = term.id
        else:
            id = term

        cur.execute('''
            DELETE FROM term_loci WHERE term = ?;
            DELETE FROM terms WHERE id = ?;
            ''', (id, id))
        if not cursor:
            cur.execute('END TRANSACTION')

    def add_terms(self, terms, overwrite=True):
        '''
            A Convenience function to add terms from an iterable.

            Parameters
            ----------
            terms : iterable of camoco.Term objects
        '''
        if overwrite:
            self.del_terms(terms)
        cur = self.db.cursor()
        cur.execute('BEGIN TRANSACTION')
        for term in terms:
            self.add_term(term, cursor=cur, overwrite=False)
        cur.execute('END TRANSACTION')

    def del_terms(self, terms):
        '''
            A Convenience function to delete many term object

            Parameters
            ----------
            terms : iterable of camoco.Term objects.
        '''
        cur = self.db.cursor()
        cur.execute('BEGIN TRANSACTION')
        for term in terms:
            self.del_term(term, cursor=cur)
        cur.execute('END TRANSACTION')

    @classmethod
    def create(cls, name, description, refgen, type='Ontology'):
        '''
            This method creates a fresh Ontology with nothing it it.
        '''
        # run the inherited create method from Camoco
        self = super().create(name, description, type=type)
        # set the refgen for the current instance
        self.refgen = refgen
        # add the global refgen
        self._global('refgen', refgen.name)
        # build the tables
        self._create_tables()
        return self

    @classmethod
    def from_terms(cls, terms, name, description, refgen):
        '''
            Convenience function to create a Ontology from an iterable
            terms object. 

            Parameters
            ----------
            terms : iterable of camoco.GOTerm objects
                Items to add to the ontology. The key being the name
                of the term and the items being the loci.
            name : str
                The name of the camoco object to be stored in the database.
            description : str
                A short message describing the dataset.
            refgen : camoco.RefGen
                A RefGen object describing the genes in the dataset
        '''
        self = cls.create(name,description,refgen)
        self.log('Adding {} terms to the database.',len(terms))
        self.add_terms(terms, overwrite=True)
        # Build the indices
        self.log('Building the indices.')
        self._build_indices()

        self.log('Your gene ontology is built.')
        return self

    def _create_tables(self):
        cur = self.db.cursor()
        cur.execute('''
            CREATE TABLE IF NOT EXISTS terms (
                id TEXT UNIQUE,
                desc TEXT
            )'''
        )
        cur.execute('''
            CREATE TABLE IF NOT EXISTS term_loci (
                term TEXT, 
                id TEXT
            );'''
        )
        cur.execute('''
            CREATE TABLE IF NOT EXISTS term_attrs (
                term TEXT,
                key TEXT,
                val TEXT
            );
        ''')

    def _clear_tables(self):
        cur = self.db.cursor()
        cur.execute('DELETE FROM terms; DELETE FROM term_loci;')

    def _build_indices(self):
        cursor = self.db.cursor()
        cursor.execute('CREATE INDEX IF NOT EXISTS termIND ON terms (id)')
        cursor.execute('CREATE INDEX IF NOT EXISTS lociIND ON term_loci (term,id)')

    def _drop_indices(self):
        cursor = self.db.cursor()
        cursor.execute('DROP INDEX IF EXISTS termIND; DROP INDEX IF EXISTS lociIND;')

    def enrichment(self, locus_list, pval_cutoff=0.05, max_term_size=300,
                   min_term_size=2, num_universe=None, return_table=False,
                   include_genes=False):
        '''
            Evaluates enrichment of loci within the locus list in terms within
            the ontology. NOTE: this only tests terms that have at least one
            locus that exists in locus_list.

            Parameters
            ----------
            locus_list : list of co.Locus
                A list of loci for which to test enrichment. i.e. is there
                an over-representation of these loci within and the terms in
                the Ontology.
            pval_cutoff : float (default: 0.05)
                Report terms with a pval lower than this value
            max_term_size : int (default: 300)
                The maximum term size for which to test enrichment. Useful
                for filtering out large terms that would otherwise be 
                uninformative (e.g. top level GO terms)
            min_term_size : int (default: 5)
                The minimum term size for which to test enrichment. Useful
                for filtering out very small terms that would be uninformative
                (e.g. single gene terms)
            num_universe : int (default: None)
                Use a custom universe size for the hypergeometric calculation, 
                for instance if you have a reduced number of genes in a reference
                co-expression network. If None, the value will be calculated as
                the total number of distinct genes that are observed in the 
                ontology.
            include_genes : boo (default: False)
                Include comma delimited genes as a field
        '''
        terms = self.db.cursor().execute('''SELECT DISTINCT term 
        FROM term_loci WHERE id IN ('{}')
        '''.format(
            "','".join([x.id for x in locus_list])
        )).fetchall()
        
        self.log('Getting GOTerm Objects')
        terms = list(
            filter(
                lambda t: (len(t) >= min_term_size) and (len(t) <= max_term_size),
                [self[name] for name, in terms]
            )
        )
        if num_universe is None:
            #num_universe = len(set(chain(*[x.loci for x in terms])))
            num_universe = self.num_distinct_loci() 
        self.log(
            'Loci occur in {} terms, containing {} genes'.format(
                len(terms), num_universe
            )
        )
        
        significant_terms = []
        locus_list = set(locus_list)
        for term in terms:
            term_genes = set(term.loci)
            if len(term_genes) > max_term_size:
                continue
            num_common = len(term_genes.intersection(locus_list))
            num_in_term = len(term_genes)
            num_sampled = len(locus_list)
            # the reason this is num_common - 1 is because we are looking for 1 - cdf
            # and we need to greater than OR equal to
            pval = hypergeom.sf(num_common-1,num_universe,num_in_term,num_sampled)
            if pval <= pval_cutoff:
                term.attrs['pval'] = pval
                term.attrs['hyper'] = OrderedDict([
                    ('pval'        , pval),
                    ('term_tested' , len(terms)),
                    ('num_common'  , num_common),
                    ('num_universe', num_universe),
                    ('term_size'   , num_in_term),
                    ('num_terms'   , len(self)),
                    ('sum_sampled' , num_sampled)
                ])
                if include_genes == True:
                    term.attrs['hyper']['genes'] = ",".join(
                        [x.id for x in term_genes.intersection(locus_list)]
                    )
                significant_terms.append(term)
        if return_table == True:
            tbl = []
            for x in significant_terms:
                val = OrderedDict([
                    ('term', x.name),
                    ('id'  , x.id)
                ])
                val.update(x.attrs['hyper'])
                tbl.append(val)
            return DataFrame.from_records(tbl).sort_values(by='pval')
        else:
            return sorted(significant_terms,key=lambda x: x.attrs['pval'])



