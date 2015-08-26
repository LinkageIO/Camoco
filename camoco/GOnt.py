#!/usr/bin/python3

from camoco.Ontology import Term,Ontology
from camoco.Camoco import Camoco
from camoco.RefGen import RefGen
from camoco.Locus import Locus
from camoco.Tools import log

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

    def add_term(self, term, overwrite=True):
        ''' This will add a single term to the ontology '''
        cur = self.db.cursor()
        if overwrite:
            self.del_term(term.name)
        cur.execute('BEGIN TRANSACTION')

        # Add the term name and description
        cur.execute('INSERT OR REPLACE INTO terms (id, desc, name) VALUES (?, ?, ?)',
            (term.id, term.desc, term.name))

        # Add the term loci
        if term.locus_list:
            cur.executemany('INSERT OR REPLACE INTO term_loci (term, id) VALUES (?, ?)',
                [(term.id, locus.id) for locus in term.locus_list])

        # Add the is_a relationships if they are there
        if term.is_a:
            cur.executemany('INSERT OR REPLACE INTO rels (parent, child) VALUES (?, ?)',
                [(parent, term.id) for parent in term.is_a])

        # Add the alternate ids to their database
        if term.alt_id:
            cur.executemany('INSERT OR REPLACE INTO alts (alt, main) VALUES (?,?)',
                [(altID, term.id) for altID in term.alt_id])

        cur.execute('END TRANSACTION')

        def del_term(self, id):
            '''This will remove a single term from the ontology.'''
            cur = self.db.cursor()
            cur.execute('BEGIN TRANSACTION')
            main_id = self.db.cursor().execute('SELECT main FROM alts WHERE alt = ?', (id, )).fetchone()
            if main_id:
                id = main_id
            cur.execute('''
                DELETE FROM terms WHERE id = ?;
                DELETE FROM term_loci WHERE term = ?;
                DELETE FROM rels WHERE child = ?;
                DELETE FROM alts WHERE main = ?;
                END TRANSACTION;
            ''', (id, id, id, id))

    @classmethod
    def create(cls, name, description, refgen, type='GOnt'):
        '''This method creates a fresh GO Ontology with nothing it it.'''

        # run the inherited create method from Camoco
        self = super().create(name, description, refgen, type=type)

        # Alter the tables as needed
        cur = self.db.cursor()
        try:
            cur.execute('ALTER TABLE terms ADD COLUMN name TEXT;')
        except lite.SQLError:
            pass
        cur.execute('CREATE TABLE IF NOT EXISTS rels (parent TEXT, child TEXT, PRIMARY KEY(parent, child));')
        cur.execute('CREATE TABLE IF NOT EXISTS alts (alt TEXT, main TEXT, PRIMARY KEY(alt, main));')
        return self

    @classmethod
    def from_obo(cls, obo_file, gene_map_file ,name, description, refgen, go_col=1, debug=False):
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

        self = cls.create(name, description, refgen)
        if debug:
            self.log.warn('Debugging on, will not record all terms.')

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
                term = alt_map[term]
            # Make the Gene Objects from the refgen
            geneObjs = refgen.from_ids(terms[term]['genes'])

            # Make the Term object and add it to the list
            termObjs.append(GOTerm(term, name=terms[term]['name'], alt_id=terms[term]['alt_id'], is_a=terms[term]['is_a'],
            desc=terms[term]['desc'], locus_list=geneObjs))
        del terms

        if debug:
            termObjs = termObjs[:100]

        # Add them all to the database
        total = len(termObjs)
        n = 0
        for term in termObjs:
            if n % 1000 == 0:
                self.log('Added {}/{} terms.', n, total)
            self.add_term(term)
            print(term.id)
            n += 1
        self.log('Added all {} records to the databse.', n)
        del termObjs
        return self
