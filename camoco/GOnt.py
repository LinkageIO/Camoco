#!/usr/bin/python3

from .Ontology import Ontology
from .Term import Term
from .Camoco import Camoco
from .RefGen import RefGen
from .Locus import Locus
from .Tools import log,rawFile

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
    def __init__(self, id, name='', desc='', alt_id=None, 
        is_a=None, loci=None, **kwargs):
        '''
            Initialize a GOTerm
        '''
        super().__init__(id, desc=desc, loci=loci, **kwargs)
        self.name = name
        self.is_a = set(is_a) if is_a else set()
        self.alt_id = set(alt_id) if alt_id else set()

    def __str__(self):
        return "Term: {}, Name: {}, Desc: {}, {} Loci".format(
            self.id, self.name, self.desc, len(self)
        )

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
        main_id = self.db.cursor().execute(
            'SELECT main FROM alts WHERE alt = ?',
            (id, )
        ).fetchone()
        if main_id:
            (id,) = main_id
        try:
            # Get the term information
            (id, desc, name) = self.db.cursor().execute(
                'SELECT * from terms WHERE id = ?', (id, )).fetchone()
        except TypeError: # Not in database
            raise KeyError('This term is not in the database.')

        # Get the loci associated with the term
        term_loci = self.refgen.from_ids([
            gene_id[0] for gene_id in self.db.cursor().execute(
            'SELECT id FROM term_loci WHERE term = ?', (id, )).fetchall()
        ])

        # Get the isa relationships
        is_a = set(termID for termID in self.db.cursor().execute(
            'SELECT parent FROM rels WHERE child = ?', (id, )).fetchall())

        # Get the alternate ids of the term
        alts = set(termID for termID in self.db.cursor().execute(
            'SELECT alt FROM alts WHERE main = ?', (id, )).fetchall())

        return GOTerm(
            id, name=name, desc=desc, alt_id=alts, 
            is_a=is_a, loci=term_loci
        )


    def add_term(self, term, cursor=None, overwrite=False):
        ''' This will add a single term to the ontology '''
        if overwrite:
            self.del_term(term.id)
        if not cursor:
            cur = self.db.cursor()
            cur.execute('BEGIN TRANSACTION')
        else:
            cur = cursor

        # Add the term name and description
        cur.execute(
            'INSERT OR ABORT INTO terms (id, desc, name) VALUES (?, ?, ?)', 
            (term.id, term.desc, term.name)
        )
        if term.loci:
            cur.executemany(
                'INSERT OR REPLACE INTO term_loci (term, id) VALUES (?, ?)', 
                [(term.id, locus.id) for locus in term.loci]
            )
        if term.is_a:
            cur.executemany(
                'INSERT OR REPLACE INTO rels (parent, child) VALUES (?, ?)',
                [(parent, term.id) for parent in term.is_a]
            )
        if term.alt_id:
            cur.executemany(
                'INSERT OR REPLACE INTO alts (alt, main) VALUES (?,?)',
                [(altID, term.id) for altID in term.alt_id]
            )
        if not cursor:
            cur.execute('END TRANSACTION')

    def del_term(self, term, cursor=None):
        '''This will remove a single term from the ontology.'''
        if not cursor:
            cur = self.db.cursor()
            cur.execute('BEGIN TRANSACTION')
        else:
            cur = cursor
        if not isinstance(term, str):
            id = term.id
        else:
            id = term
        main_id = cur.execute('SELECT main FROM alts WHERE alt = ?', (id, )).fetchone()
        if main_id:
            id = main_id
        else:
            id = term.id
        cur.execute('''
            DELETE FROM terms WHERE id = ?;
            DELETE FROM term_loci WHERE term = ?;
            DELETE FROM rels WHERE child = ?;
            DELETE FROM rels WHERE parent = ?;
            DELETE FROM alts WHERE main = ?;
            ''', (id, id, id, id, id))
        if not cursor:
            cur.execute('END TRANSACTION')

    @staticmethod
    def getISA(terms, term):
        '''
            Internal function to handle propagating is_a relationships
        '''
        tot = set()
        if not terms[term]['is_a']:
            return tot
        else:
            for item in terms[term]['is_a']:
                tot.add(item)
                tot = tot | GOnt.getISA(terms, item)
            return tot

    @classmethod
    def create(cls, name, description, refgen, overwrite=True, type='GOnt'):
        '''This method creates a fresh GO Ontology with nothing it it.'''

        # run the inherited create method from Camoco
        self = super().create(name, description, refgen, type=type)

        if overwrite:
            self._drop_indices()
            self._clear_tables()

        return self

    def _parse_obo(self,obo_file):
        '''
            Parses a GO obo file
            
            Returns
            -------
            A dictionary containing id:fields for the tersm  
        '''
        # Importing the obo information
        self.log('Importing OBO: {}', obo_file)
        terms = defaultdict(dict)
        cur_term = ''
        isa_re = re.compile('is_a: (.*) !.*')
        with rawFile(obo_file) as INOBO:
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
                    terms[the_id] = terms[cur_term]
                elif line.startswith('def: '):
                    terms[cur_term]['desc'] += line.replace('def: ', '')
                elif line.startswith('comment: '):
                    terms[cur_term]['desc'] += line.replace('comment: ', '')
                elif line.startswith('is_a: '):
                    terms[cur_term]['is_a'].add(isa_re.match(line).group(1))

        # Propagating the relationships using the embeded function
        self.log('Propagating is_a relationships.')
        for cur_term in terms:
            terms[cur_term]['is_a'] = self.getISA(terms,cur_term)
        return terms

    def _parse_gene_term_map(self,gene_map_file,headers=True,
            go_col=1,id_col=0):
        # Importing gene map information, and cross referencing with obo information
        self.log('Importing Gene Map: {}', gene_map_file)
        genes = dict()
        gene = ''
        cur_term = ''
        with rawFile(gene_map_file) as INMAP:
            if headers:
                garb = INMAP.readline()
            for line in INMAP.readlines():
                row = line.strip().split('\t')
                gene = row[id_col].split('_')[0].strip()
                cur_term = row[go_col]
                # Make a map between genes and associated GO terms
                if gene not in genes:
                    genes[gene] = set([cur_term])
                else:
                    genes[gene].add(cur_term)
        return genes

    @classmethod
    def from_obo(cls, obo_file, gene_map_file ,name, description, 
            refgen, go_col=1, id_col=0, headers=True, overwrite=True):
        ''' 
            Convenience function for importing GOnt from obo files 

            Parameters
            ----------
            obo_file : str
                Path to the obo file
            gene_map_file : str
                Path to the file which specifies what GO term each 
                gene is a part of.
            name : str
                The name of the camoco object to be stored in the database.
            description : str
                A short message describing the dataset.
            refgen : camoco.RefGen
                A RefGen object describing the genes in the dataset
            go_col : int (default: 1)
                The index column for GO term in the gene_map_file
            id_col : int (default: 0)
                The index column for gene id in the gene_map_file
            headers : bool (default: True)
                A flag indicating whether or not there is a header line
                in the gene_map_file
            overwrite : bool (default: True)
                Kill old instances.
        '''
        self = cls.create(name, description, refgen, overwrite=overwrite)
        # Parse the input files
        terms = self._parse_obo(obo_file)
        genes = self._parse_gene_term_map(gene_map_file,
            headers=headers,go_col=go_col,id_col=id_col
        )

        # Get the requisite gene objects
        self.log('Mixing genes and term data.')
        for (cur_gene, cur_terms) in genes.items():
            for cur_term in cur_terms:
                if terms[cur_term] == {}:
                    #self.log('Could not find term for {}',cur_term)
                    continue
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
            termObjs.append(
                GOTerm(term, name=terms[term]['name'], 
                    alt_id=terms[term]['alt_id'], is_a=terms[term]['is_a'],
                    desc=terms[term]['desc'], loci=geneObjs
                )
            )
        del terms

        # Add them all to the database
        self.log('Adding {} terms to the database.',len(termObjs))
        self.add_terms(termObjs, overwrite=False)
        del termObjs

        # Build the indices
        self.log('Building the indices.')
        self._build_indices()

        self.log('Your gene ontology is built.')
        return self

    def _create_tables(self):
        # Alter the tables as needed
        super()._create_tables()
        cur = self.db.cursor()
        try:
            cur.execute('ALTER TABLE terms ADD COLUMN name TEXT;')
        except lite.SQLError:
            pass
        cur.execute('CREATE TABLE IF NOT EXISTS rels (parent TEXT, child TEXT);')
        cur.execute('CREATE TABLE IF NOT EXISTS alts (alt TEXT UNIQUE, main TEXT);')

    def _clear_tables(self):
        super()._clear_tables()
        cur = self.db.cursor()
        cur.execute('DELETE FROM rels; DELETE FROM alts;')

    def _build_indices(self):
        super()._build_indices()
        cursor = self.db.cursor()
        cursor.execute('''
            CREATE INDEX IF NOT EXISTS relsIND ON rels (child);
            CREATE INDEX IF NOT EXISTS altsIND ON alts (main);
            CREATE INDEX IF NOT EXISTS term_loci_ID 
                ON term_loci (term);
        ''')

    def _drop_indices(self):
        super()._drop_indices()
        cursor = self.db.cursor()
        cursor.execute('''
            DROP INDEX IF EXISTS relsIND; 
            DROP INDEX IF EXISTS altsIND;
        ''')
