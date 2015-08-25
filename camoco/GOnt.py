#!/usr/bin/python3

from camoco.Ontology import Term,Ontology
from camoco.Camoco import Camoco
from camoco.RefGen import RefGen
from camoco.Locus import Locus
from camoco.Tools import log

from collections import defaultdict
import networkx as nx
import pandas as pd
import numpy as np

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
        if self.refgen:
            self.refgen = RefGen(self.refgen)

    def __getitem__(self, id):
        ''' retrieve a term by id '''
        try:
            id = self.db.cursor().execute(
                'SELECT primary FROM alts WHERE alt = ?', (id, )).fetchone()
        except TypeError as e:
            id = id

        try:
            # Get the term information
            (id, desc, name) = self.db.cursor().execute(
                'SELECT * from terms WHERE id = ?', (id, )).fetchone()
        except TypeError as e: # Not in database
            raise

        # Get the loci associated with the term
        term_loci = [self.refgen[gene_id] for gene_id in self.db.cursor().execute(
            'SELECT id FROM term_loci WHERE term = ?', (id, )).fetchall()]

        # Get the isa relationships
        is_a = set(termID for termID in self.db.cursor().execute(
            'SELECT parent FROM rels WHERE child = ?', (id, )).fetchall())

        # Get the alternate ids of the term
        alts = set(termID for termID in self.db.cursor().execute(
            'SELECT alt FROM alts WHERE primary = ?', (id, )).fetchall())

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
            cur.executemany('INSERT OR REPLACE INTO alts (alt, primary) VALUES (?,?)',
                [(altID, term.id) for altID in term.alt_id])

        cur.execute('END TRANSACTION')

    @classmethod
    def create(cls, name, description, refgen, type='GOnt'):
        '''
            This method creates a fresh GO Ontology with nothing it it.
        '''
        # run the inherited create method from Camoco
        self = super().create(name, description, refgen, type=type)

        # Alter the tables as needed
        cur = self.db.cursor()
        cur.execute('ALTER TABLE terms ADD COLUMN name TEXT;')
        cur.execute('CREATE TABLE IF NOT EXISTS rels (parent TEXT, child TEXT, PRIMARY KEY(parent, child));')
        cur.execute('CREATE TABLE IF NOT EXISTS alts (alt TEXT UNIQUE, primary TEXT, PRIMARY KEY(alt, primary));')
        return self

    @classmethod
    def from_obo(cls,obo_file,gene_map_file,name,description,refgen,go_col=1):
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

        # Importing the obo information
        self.log('Importing OBO: {}', obo_file)
        terms = defaultdict(dict)
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
                    terms[cur_term]['alt_id'].add(line.replace('alt_id: ', ''))
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
        self.log('Mixing genes and obo data.')
        for (cur_gene, cur_terms) in genes.items():
            for cur_term in cur_terms:
                terms[cur_term]['genes'].add(cur_gene)
                for parent in terms[cur_term]['is_a']:
                    terms[parent]['genes'].add(cur_gene)
        del genes

        # Making the Term objects
        self.log('Making the objects to insert into the database.')
        termObjs = []
        for (term, info) in terms.items():
            # Make the Gene Objects from the refgen
            geneObjs = refgen.from_ids(info['genes'])

            # Make the Term object and add it to the list
            termObjs.append(Term(term, name=info['name'], type=info['type'],
            desc=info['desc'], locus_list=geneObjs))
        del terms

        # Add them all to the database
        self.log('Adding the records to the databse.')
        n = 1
        for term in termObjs:
            if n % 10000 == 0:
                self.log(str(n)+' terms added.')
            self.add_term(term)
            n += 1
        del termObjs
        return self
