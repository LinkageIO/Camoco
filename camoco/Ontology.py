#!/usr/bin/python3

from .Camoco import Camoco
from .RefGen import RefGen
from .Locus import Locus
from .Term import Term

from pandas import DataFrame
from scipy.stats import hypergeom

import sys

class Ontology(Camoco):
    '''
        An Ontology is just a collection of terms. Each term is just a
        collection of genes. Sometimes terms are related or nested
        within each other, sometimes not. Simple enough.
    '''
    def __init__(self, name, type='Ontology'):
        super().__init__(name, type=type)
        if self.refgen:
            self.refgen = RefGen(self.refgen)

    def __len__(self):
        return self.db.cursor().execute(
                "SELECT COUNT(*) FROM terms;"
        ).fetchone()[0]

    def __getitem__(self, id):
        ''' retrieve a term by id '''
        try:
            (id, desc) = self.db.cursor().execute(
                'SELECT * from terms WHERE id = ?', (id, )
            ).fetchone()
            term_loci = [
                self.refgen[gene_id] for gene_id in self.db.cursor().execute(
                ''' SELECT id FROM term_loci WHERE term = ?''', (id, )
            ).fetchall()]
            return Term(id, desc=desc, loci=term_loci)
        except TypeError as e: # Not in database
            raise e

    def iter_terms(self):
        '''
            Return a generator that iterates over each term in the ontology.
        '''
        for id, in self.db.cursor().execute("SELECT id FROM terms"):
            yield self[id]

    def terms(self):
        return list(self.iter_terms())

    def summary(self):
        return "Ontology:{} - desc: {} - contains {} terms for {}".format(
            self.name, self.description, len(self), self.refgen)

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
            )'''
        )

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

    def enrichment(self, locus_list, pval_cutoff=0.05, gene_filter=None,
            label=None, max_term_size=300):
        raise NotImplementedError('This is broken')
        
