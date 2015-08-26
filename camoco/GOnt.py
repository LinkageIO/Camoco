#!/usr/bin/python3

from .Ontology import Ontology
from .Term import Term
from .Camoco import Camoco
from .RefGen import RefGen
from .Locus import Locus
from .Tools import log

from collections import defaultdict
import networkx as nx
import pandas as pd
import apsw as lite
import numpy as np
import re

class GOTerm(Term):
    '''
        Subclass to handle the intricacies of GO Terms
    '''
    def __init__(self, id, name='', desc='', alt_id=None, is_a=None, locus_list=None, **kwargs):
        super().__init__(id, desc=desc, locus_list=locus_list, **kwargs)
        self.name = name
        self.is_a = set(is_a) if is_a else set()
        self.alt_id = set(alt_id) if alt_id else set()

    def __str__(self):
        return "Term: {}, Name: {}, Desc: {}, {} Loci".format(self.id, self.name, self.desc, len(self))

    def __repr__(self):
        return str(self)

    def add_parent(self,termname):
        'Add parent to term'
        self.is_a.add(termname)

    def add_alt(self, altID):
        'Add alternate name to term'
        self.alt_id.add(altID)

class GOnt(Ontology):
    '''Ontology extension for GO'''
    def __init__(self, name, type='GOnt'):
        super().__init__(name, type=type)

    def __getitem__(self, id):
        ''' retrieve a term by id '''
        main_id = self.db.cursor().execute('SELECT main FROM alts WHERE alt = ?', (id, )).fetchone()
        if main_id:
            id = main_id

        try:
            # Get the term information
            (id, desc, name) = self.db.cursor().execute(
                'SELECT * from terms WHERE id = ?', (id, )).fetchone()
        except TypeError: # Not in database
            raise ValueError('This term is not in the database.')

        # Get the loci associated with the term
        term_loci = [self.refgen[gene_id] for gene_id in self.db.cursor().execute(
            'SELECT id FROM term_loci WHERE term = ?', (id, )).fetchall()]

        # Get the isa relationships
        is_a = set(termID for termID in self.db.cursor().execute(
            'SELECT parent FROM rels WHERE child = ?', (id, )).fetchall())

        # Get the alternate ids of the term
        alts = set(termID for termID in self.db.cursor().execute(
            'SELECT alt FROM alts WHERE main = ?', (id, )).fetchall())

        return GOTerm(id, name=name, desc=desc, alt_id=alts, is_a=is_a, locus_list=term_loci)

    def _build_indices(self):
        cursor = self.db.cursor()
        cursor.execute('''
            BEGIN TRANSACTION;
            CREATE INDEX IF NOT EXISTS termind ON terms (id);
            CREATE INDEX IF NOT EXISTS lociind ON term_loci (term,id);
            CREATE INDEX IF NOT EXISTS relsind ON rels (parent,child);
            CREATE INDEX IF NOT EXISTS altsind ON alts (alt,main);
            END TRANSACTION;
            ''')

    def _drop_indices(self):
        cursor = self.db.cursor()
        cursor.execute('''
            BEGIN TRANSACTION;
            DROP INDEX IF EXISTS termind;
            DROP INDEX IF EXISTS lociind;
            DROP INDEX IF EXISTS relsind;
            DROP INDEX IF EXISTS altsind;
            END TRANSACTION;
            ''')

    def _add_term(self, term, cursor):
        ''' This will add a single term to the ontology '''

        # Add the term name and description
        cursor.execute('INSERT OR ABORT INTO terms (id, desc, name) VALUES (?, ?, ?)', (term.id, term.desc, term.name))

        if term.locus_list:
            cursor.executemany('INSERT OR ABORT INTO term_loci (term, id) VALUES (?, ?)', [(term.id, locus.id) for locus in term.locus_list])

        if term.is_a:
            cursor.executemany('INSERT OR ABORT INTO rels (parent, child) VALUES (?, ?)', [(parent, term.id) for parent in term.is_a])

        if term.alt_id:
            cursor.executemany('INSERT OR ABORT INTO alts (alt, main) VALUES (?,?)', [(altID, term.id) for altID in term.alt_id])

    def _del_term(self, term, cursor):
        '''This will remove a single term from the ontology.'''
        main_id = cursor.execute('SELECT main FROM alts WHERE alt = ?', (term.id, )).fetchone()
        if main_id:
            id = main_id
        else:
            id = term.id
        cursor.execute('''
            DELETE FROM terms WHERE id = ?;
            DELETE FROM term_loci WHERE term = ?;
            DELETE FROM rels WHERE child = ?;
            DELETE FROM rels WHERE parent = ?;
            DELETE FROM alts WHERE main = ?;
            ''', (id, id, id, id, id))

    def add_terms(self, terms, overwrite=False):
        if overwrite:
            self.del_terms(terms)

        cursor = self.db.cursor()
        cursor.execute('BEGIN TRANSACTION')
        for term in terms:
            self._add_term(term, cursor)
        cursor.execute('END TRANSACTION')

    def del_terms(self, terms):
        cursor = self.db.cursor()
        cursor.execute('BEGIN TRANSACTION')
        for term in terms:
            self._del_term(term, cursor)
        cursor.execute('END TRANSACTION')

    @classmethod
    def create(cls, name, description, refgen, overwrite=True, type='GOnt'):
        '''This method creates a fresh GO Ontology with nothing it it.'''

        # run the inherited create method from Camoco
        self = super().create(name, description, refgen, type=type)

        # Alter the tables as needed
        cur = self.db.cursor()
        try:
            cur.execute('ALTER TABLE terms ADD COLUMN name TEXT;')
        except lite.SQLError:
            pass
        cur.execute('CREATE TABLE IF NOT EXISTS rels (parent TEXT, child TEXT);')
        cur.execute('CREATE TABLE IF NOT EXISTS alts (alt TEXT UNIQUE, main TEXT);')

        if overwrite:
            self._drop_indices()
            cur.execute('''
                DELETE FROM terms;
                DELETE FROM term_loci;
                DELETE FROM rels;
                DELETE FROM alts;
            ''')

        return self

    @classmethod
    def from_obo(cls, obo_file, gene_map_file ,name, description, refgen, go_col=1, overwrite=True):
        ''' Convenience function for importing GO obo files '''

        # Internal function to handle propagating is_a relationships
        def getISA(terms, term):
            tot = set()
            if not terms[term]['is_a']:
                return tot
            else:
                for item in terms[term]['is_a']:
                    tot.add(item)
                    tot = tot | getISA(terms, item)
                return tot

        self = cls.create(name, description, refgen, overwrite=overwrite)

        # Importing the obo information
        self.log('Importing OBO: {}', obo_file)
        terms = defaultdict(dict)
        alt_map = dict()
        cur_term = ''
        isa_re = re.compile('is_a: (.*) !.*')
        with open(obo_file, 'r') as INOBO:
            for line in INOBO:
                line = line.strip()
                if line.startswith('id: '):
                    cur_term = line.replace('id: ', '')
                elif line.startswith('name: '):
                    terms[cur_term]['name'] = line.replace('name: ', '')
                    terms[cur_term]['desc'] = ''
                    terms[cur_term]['is_a'] = set()
                    terms[cur_term]['alt_id'] = set()
                    terms[cur_term]['genes'] = set()
                elif line.startswith('namespace: '):
                    terms[cur_term]['type'] = line.replace('namespace: ', '')
                elif line.startswith('alt_id: '):
                    the_id = line.replace('alt_id: ', '')
                    terms[cur_term]['alt_id'].add(the_id)
                    alt_map[the_id] = cur_term
                elif line.startswith('def: '):
                    terms[cur_term]['desc'] += line.replace('def: ', '')
                elif line.startswith('comment: '):
                    terms[cur_term]['desc'] += line.replace('comment: ', '')
                elif line.startswith('is_a: '):
                    terms[cur_term]['is_a'].add(isa_re.match(line).group(1))

        # Propagating the relationships using the embeded function
        self.log('Propagating is_a relationships.')
        for cur_term in terms:
            terms[cur_term]['is_a'] = getISA(terms,cur_term)

        # Importing gene map information, and cross referencing with obo information
        self.log('Importing Gene Map: {}', gene_map_file)
        genes = dict()
        gene = ''
        cur_term = ''
        INMAP = open(gene_map_file)
        garb = INMAP.readline()
        for line in INMAP.readlines():
            row = line.strip().split('\t')
            gene = row[0].split('_')[0].strip()
            cur_term = row[go_col]
            if gene not in genes:
                genes[gene] = set([cur_term])
            else:
                genes[gene].add(cur_term)
        del cur_term
        INMAP.close()

        # Get the requisite gene objects
        self.log('Mixing genes and term data.')
        for (cur_gene, cur_terms) in genes.items():
            for cur_term in cur_terms:
                if terms[cur_term] == {}:
                    cur_term = alt_map[cur_term]
                terms[cur_term]['genes'].add(cur_gene)
                for parent in terms[cur_term]['is_a']:
                    terms[parent]['genes'].add(cur_gene)
        del genes

        # Making the Term objects
        self.log('Making the objects to insert into the database.')
        termObjs = []
        for (term, info) in terms.items():
            if terms[term] == {}:
                continue
            # Make the Gene Objects from the refgen
            geneObjs = refgen.from_ids(terms[term]['genes'])
            # Make the Term object and add it to the list
            termObjs.append(GOTerm(term, name=terms[term]['name'], alt_id=terms[term]['alt_id'], is_a=terms[term]['is_a'],
            desc=terms[term]['desc'], locus_list=geneObjs))
        del terms

        # Add them all to the database
        self.log('Adding {} terms to the database.',len(termObjs))
        self.add_terms(termObjs, overwrite=False)
        del termObjs

        # Build the indices
        self.log('Building the indices.')
        self._build_indices()

        self.log('Ontology is built!')
        return self
