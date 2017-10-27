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
import copy
import numpy
import camoco as co
import pandas as pd

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

    def __iter__(self):
        return self.iter_terms()

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

    def terms_containing(self,locus_list,max_term_size=10e10,min_term_size=0):
        '''
            Retrurns the set of terms which contains the 
            specified loci.

            Parameters
            ----------
            locus_list : iterable of type Locus
                The list of loci for which to retrieve 
                corresponding terms.
            max_term_size : int (default: 10e10)
                The maximum term size for which to test enrichment. Useful
                for filtering out large terms that would otherwise be 
                uninformative (e.g. top level GO terms)
            min_term_size : int (default: 0)
                The minimum term size for which to test enrichment. Useful
                for filtering out very small terms that would be uninformative
                (e.g. single gene terms)

            Returns
            -------
            list of terms which contain provided loci
        '''
        # Filter to unique set
        locus_list = set(locus_list)
        # query the database
        terms = self.db.cursor().execute('''SELECT DISTINCT term 
        FROM term_loci WHERE id IN ('{}')
        '''.format(
            "','".join([x.id for x in locus_list])
        )).fetchall()
        # Fetch the terms with the proper size
        terms = list(
            filter(
                lambda t: (len(t) >= min_term_size) and (len(t) <= max_term_size),
                [self[name] for name, in terms]
            )
        )
        return terms


    def num_distinct_loci(self):
        return self.db.cursor().execute(
            'SELECT COUNT(DISTINCT(id)) FROM term_loci;'
        ).fetchone()[0]

    def iter_terms(self,min_term_size=0,max_term_size=10e10):
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

    def terms(self,min_term_size=0,max_term_size=10e10):
        return list(self.iter_terms(min_term_size=min_term_size,max_term_size=max_term_size))

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
        ''' 
        This will add a single term to the ontology

        Parameters
        ----------
        term : Term object
            The term object you wish to add.
        cursor : apsw cursor object
            A initialized cursor object, for batch operation. This will
            allow for adding many terms in one transaction as long as the 
            passed in cursor has executed the "BEGIN TRANSACTION" command.
        overwrite : bool
            Indication to delete any existing entry before writing
        '''

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

        # Add the term attrs
        if term.attrs:
            for key,val in term.attrs.items():
                cur.execute('''
                    INSERT OR ABORT INTO term_attrs (term,key,val)
                    VALUES (?,?)
                ''',(term.id,key,val))

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

        try:
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
        except Exception as e:
            cur.execute('ROLLBACK')
            raise e


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

    def set_strongest(self,attr=None,higher=None):
        '''
            Convinience function that allows you to set default values for
            strongest SNP2Gene mapping tasks.

            Parameters
            ----------
            attr:     The locus attr used to determine which locus is the 
                      strongest locus.
            
            higher:   Flag indicating whether the value in --strongest-attr
                      is stronger if it is higher. Default behavior is to
                      treatlower values as stronger (i.e. p-vals)
        '''
        if not(attr is None):
            self._global('strongest_attr',attr)
        if not(higher is None):
            self._global('strongest_higher',higher)

    def get_strongest_attr(self):
        '''
            Convinience function that allows you to get the default value for
            the locus attr used to determine which locus is the strongest locus
            strongest SNP2Gene mapping.
        '''
        return self._global('strongest_attr')
    
    def get_strongest_higher(self):
        '''
            Convinience function that allows you to get default values for
            the flag indicating whether the value in `strongest-attr` is
            is stronger if higher for strongest SNP2Gene mapping tasks.
        '''
        return self._global('strongest_higher')


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
    def from_DataFrame(cls, dataframe, name, description, refgen,
                      gene_col='gene',term_col='Term'):
        '''
            Convenience function to create a Ontology from an iterable
            terms object. 

            Parameters
            ----------
            dataframe : pandas.DataFrame
                A pandas dataframe containing the mapping betweeen gene ids
                and 
            name : str
                The name of the camoco object to be stored in the database.
            description : str
                A short message describing the dataset.
            refgen : camoco.RefGen
                A RefGen object describing the genes in the dataset

            Optional Parameters
            -------------------
            gene_col : str (default: gene)
                The string designating the column in the dataframe containing
                gene names (ids)
            term_col : str (default: Term)
                The string designating the column in the dataframe containing
                the term name.

        '''
        self = cls.create(name,description,refgen)
        # create terms from 
        terms = [
            Term(id,loci=refgen[set(df[gene_col])]) \
            for id,df in dataframe.groupby(term_col)
        ]
        self.log('Adding {} terms to the database.',len(terms))
        self.add_terms(terms, overwrite=True)
        # Build the indices
        self.log('Building the indices.')
        self._build_indices()
        self.log('Your gene ontology is built.')
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
                   label=None,include_genes=False,bonferroni_correction=True,
                   min_overlap=1):
        '''
            Evaluates enrichment of loci within the locus list for terms within
            the ontology. NOTE: this only tests terms that have at least one
            locus that exists in locus_list.

            Parameters
            ----------
            locus_list : list of co.Locus *of* instance of co.Ontology
                A list of loci for which to test enrichment. i.e. is there
                an over-representation of these loci within and the terms in
                the Ontology. If an ontology is passed, each term in the ontology
                will be iterated over and tested as if they were a locus_list.
            pval_cutoff : float (default: 0.05)
                Report terms with a pval lower than this value
            bonferroni_correction : bool (default: True)
                correct for testing multiple terms using Bonferroni correction
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
            include_genes : bool (default: False)
                Include comma delimited genes as a field
            return_table : bool (default: False)
                If True, return results as a data frame
            label: str (default: None)
                If a label is specified, it will be inlcuded in the results
            min_overlap : int (default: 1)
                The minimum overlap between genes in the term and genes in
                the locus list. Increasing this value can minimize spurious
                or uninformative terms
        '''
        if isinstance(locus_list,co.Ontology):
            ontology = locus_list
            self.log('Calculating enrichment for an  Ontology: {}',ontology.name)
            enrich = []
            for term in ontology:
                e = self.enrichment(
                    term.loci,
                    pval_cutoff=pval_cutoff,
                    max_term_size=max_term_size,
                    min_term_size=min_term_size,
                    num_universe=num_universe,
                    return_table=return_table,
                    label=term.id,
                    include_genes=include_genes,
                    bonferroni_correction=bonferroni_correction,
                    min_overlap=min_overlap,
                ) 
                enrich.append(e)
            if return_table:
                return pd.concat(enrich)
            else:
                return enrich
        # return a new copy of each 

        terms = [copy.copy(term) for term in self.terms_containing(
            locus_list,
            min_term_size=min_term_size,
            max_term_size=max_term_size
        )]
        # Calculate the size of the Universe
        if num_universe is None:
            num_universe = self.num_distinct_loci() 
        self.log(
            '{}: Loci occur in {} terms, containing {} genes'.format(
                label,len(terms), num_universe
            )
        )
        significant_terms = []
        for term in terms:
            term_genes = set(term.loci)
            if len(term_genes) > max_term_size:
                continue
            num_common = len(term_genes.intersection(locus_list))
            num_in_term = len(term_genes)
            num_sampled = len(locus_list)
            # the reason this is num_common - 1 is because we are looking for 1 - cdf
            # and we need to greater than OR EQUAL TO num_common
            # Look. Do this in ipython:
            '''
                In [99]: probs = [hypergeom.pmf(x,100,5,10) for x in range(0,6)]
                In [100]: probs
                Out[100]: 
                [0.58375236692612187,
                 0.33939091100357333,
                 0.070218809173150043,
                 0.006383528106649855,
                 0.00025103762217164457,
                 3.3471682956218215e-06]
                In [103]: 1-sum(probs[0:3]) 
                # Get the probs of drawing 3 or more
                Out[103]: 0.006637912897154763
                # Remember slicing is exclusive for the end value
                In [105]: hypergeom.sf(3,100,5,10)
                # That aint right
                Out[105]: 0.00025438479046726637
                In [106]: hypergeom.sf(3-1,100,5,10)
                # See Dog? You wnat num_common - 1
                Out[106]: 0.0066379128971171221
                # can we go back to drinking coffee now?
            '''
            pval = hypergeom.sf(num_common-1,num_universe,num_in_term,num_sampled)
            if pval <= pval_cutoff and num_common >= min_overlap:
                if label != None:
                    term.attrs['label'] = label
                term.attrs['hyper'] = OrderedDict([
                    ('pval'        , pval),
                    ('term_tested' , len(terms)),
                    ('num_common'  , num_common),
                    ('num_universe', num_universe),
                    ('term_size'   , num_in_term),
                    ('num_terms'   , len(self)),
                    ('num_sampled' , num_sampled)
                ])
                if bonferroni_correction == True:
                    # Right now this isn't true bonferroni, its only correcting for
                    # the number of terms that had term genes in it
                    if pval > pval_cutoff / len(terms):
                        term.attrs['hyper']['bonferroni'] = False
                    else:
                        term.attrs['hyper']['bonferroni'] = True
                term.attrs['pval'] = pval
                if include_genes == True:
                    term.attrs['hyper']['genes'] = ",".join(
                        [x.id for x in term_genes.intersection(locus_list)]
                    )
                significant_terms.append(term)
        self.log('\t\tFound {} was significant for {} terms',label,len(significant_terms))
        if return_table == True:
            tbl = []
            for x in significant_terms:
                val = OrderedDict([
                    ('name', x.name),
                    ('id'  , x.id)
                ])
                val.update(x.attrs['hyper'])
                val.update(x.attrs)
                del val['hyper']
                tbl.append(val)
            tbl = DataFrame.from_records(tbl)
            if label != None:
                tbl['label'] = label
            if len(tbl) > 0:
                tbl = tbl.sort_values(by='pval')
            return tbl
        else:
            return sorted(significant_terms,key=lambda x: x.attrs['pval'])



