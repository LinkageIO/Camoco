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
            return Term(id, desc=desc, locus_list=term_loci)
        except TypeError as e: # Not in database
            raise e

    def iter_terms(self):
        for id, in self.db.cursor().execute("SELECT id FROM terms"):
            yield self[id]

    def terms(self):
        return list(self.iter_terms())

    def print_term_stats(
        self, cob_list, filename=None, window_size=100000,
        gene_limit=4, num_bootstrap=50, bootstrap_density=2):
        '''
            Print various statistics about a term
        '''
        for term in self.iter_terms():
            term.print_stats(
                cob_list, filename, window_size=window_size,
                gene_limit=gene_limit, num_bootstraps=num_bootstraps,
                bootstrap_density=boostrap_density
            )

    def summary(self):
        return "Ontology:{} - desc: {} - contains {} terms for {}".format(
            self.name, self.description, len(self), self.refgen
        )

    def del_term(self, id):
        '''
        Remove a term from the dataset.

        Parameters
        ----------
        id : string
            The term id you wish to remove.

        Returns
        -------
            bool indicating success
        '''
        cur = self.db.cursor()
        cur.execute('''
            BEGIN TRANSACTION;
            DELETE FROM term_loci WHERE term = ?;
            DELETE FROM terms WHERE id = ?;
            END TRANSACTION;
        ''', (id, id))

    def add_term(self, term, overwrite=True):
        ''' This will add a single term to the ontology '''
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
        for locus in term.locus_list:
            cur.execute('''
                INSERT OR REPLACE INTO term_loci (term, id)
                VALUES (?, ?)
                ''', (term.id, locus.id)
            )
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
                desc TEXT,
                PRIMARY KEY(id)
            );
        ''')
        cur.execute('''
            CREATE TABLE IF NOT EXISTS term_loci (
                term TEXT,
                id TEXT,
                PRIMARY KEY(term, id)
            );
        ''')
        self._build_indices()

    # Unimplemented

    def _build_indices(self):
        pass

    def enrichment(self, gene_list, pval_cutoff=0.05, gene_filter=None, label=None, max_term_size=300):
        raise NotImplementedError('This is broken')
        # extract possible terms for genes
        if label:
            self.log("Caculating Enrichemnt for {}", label)
        cur = self.db.cursor()
        terms = [ x[0] for x in cur.execute(
            '''SELECT DISTINCT(term) FROM gene_terms
            WHERE gene IN ('{}');'''.format(
                "','".join([x.id for x in gene_list])
            )
        )]
        # compute hypergeometric for each term
        enrichment = []
        for id in terms:
            # Generate a list of records (tuples) and append to enrichment
            # list, which will eventually be a DataFrame.
            try:
                (id, name, type, desc) = cur.execute(
                    "SELECT * FROM terms WHERE id = ?", (id, )
                ).fetchone()
            except TypeError as e:
                self.log("No information for ontology term {}", id)
            genes_in_term = [x[0] for x in cur.execute(
                '''SELECT gene FROM gene_terms WHERE term = ?''', (id, ))
            ]
            if len(genes_in_term) > max_term_size:
                self.log(
                    "Skipping {} due to size ({})", 
                    name, len(genes_in_term)
                )
                continue
            if gene_filter:
                genes_in_term = [
                    gene for gene in genes_in_term if gene in gene_filter
                ]
            num_genes_in_term = len(genes_in_term)
            overlap = set(genes_in_term).intersection(
                set([x.id for x in gene_list])
            )
            num_genes_total, = cur.execute(
                'SELECT COUNT(DISTINCT(gene)) FROM gene_terms;'
            ).fetchone()
            pval = hypergeom.sf(
                len(overlap)-1, num_genes_total,
                num_genes_in_term, len(gene_list)
            )
            term_genes = ", ".join(overlap)
            enrichment.append(
                (id, name, pval, num_genes_in_term, len(overlap), 
                 len(gene_list), num_genes_total, type, term_genes, desc)
            )
        try:
            enrichment = DataFrame(enrichment,
                columns = [
                    'TermID', 'Name', 'pval', 'LenTerm', 'LenOverlap',
                    'LenList', 'LenTotal', 'Type', 'TermGenes', 'Desc'
                ]
            ).sort('pval', ascending=True)
            enrichment.index = enrichment.TermID
        except ValueError as e:
            self.log(
                "No enrichment for {}",
                ",".join([x.id for x in gene_list])
            )
            return DataFrame()
        if label:
            enrichment['Label'] = label
        return enrichment[enrichment.pval <= pval_cutoff]


