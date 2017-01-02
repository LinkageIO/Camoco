#!/usr/bin/python3

from .Ontology import Ontology
from .Term import Term
from .Camoco import Camoco
from .RefGen import RefGen
from .Locus import Locus
from .Tools import log,rawFile

from collections import defaultdict
from itertools import chain
from functools import lru_cache

import pandas as pd
import apsw as lite
import numpy as np
import networkx as nx
import re

class GOTerm(Term):
    '''
        Subclass to handle the intricacies of GO Terms

        GO term are groups of genes with an evidence based function.
        They are curated by the GeneOntology consortium:
        (http://geneontology.org/page/download-ontology)

        GO Terms are just special cases of the camoco.Term class.
        Unlike normal Terms, GO Terms are related to one another. 
        They form a Directed Acyclic Graph, i.e. they form a tree
        where terms near the leaves are more functionally specific.

        Normally, ontologies are split into three clades: Cellular Component,
        Biological Process, and Molecular Function. 

        Each GO Terms has these properties:
        - GO Id
        - Name 
        - Namespace (BP, CC, pr MF)
        - Definition
        - Relationship to other terms
        - A Set of genes it describes
        Terms Optionally have:
        - Alt/Secondary IDs
        - Synonyms
        - Database cross references
        - Comments
        - Tags


    '''
    def __init__(self, id, name='', desc='', alt_id=None, 
        is_a=None, loci=None, **kwargs):
        '''
            Initialize a GOTerm

            Parameters
        '''
        super().__init__(id, desc=desc, loci=loci, **kwargs)
        self.name = name
        # GO terms have parents as well as alternate IDs
        self.is_a = set(is_a) if is_a else set()
        self.alt_id = set(alt_id) if alt_id else set()

    @property
    def namespace(self):
        return self.attrs['namespace']

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

    @lru_cache(maxsize=2**17)
    def __getitem__(self, id):
        ''' retrieve a term by id '''
        cur = self.db.cursor()
        main_id = cur.execute(
            'SELECT main FROM alts WHERE alt = ?',
            (id, )
        ).fetchone()
        if main_id:
            (id,) = main_id
        try:
            # Get the term information
            (id, desc, name) = cur.execute(
                'SELECT * from terms WHERE id = ?', (id, )).fetchone()
        except TypeError: # Not in database
            raise KeyError('This term is not in the database.')

        # Get the loci associated with the term
        term_loci = self.refgen.from_ids([
            gene_id[0] for gene_id in cur.execute(
            'SELECT id FROM term_loci WHERE term = ?', (id, )).fetchall()
        ])

        # Get the isa relationships
        is_a = set(chain(*cur.execute(
            'SELECT parent FROM rels WHERE child = ?', (id, )).fetchall())
        )

        # Get the alternate ids of the term
        alts = set(chain(*cur.execute(
            'SELECT alt FROM alts WHERE main = ?', (id, )).fetchall())
        )
        
        # Causing issues for me, but I have to rebuild and see if that fixes
        #term_attrs = {k:v for k,v in self.db.cursor().execute(
        #    ''' SELECT key,val FROM term_attrs WHERE term = ?''',(id,)         
        #    )
        #}

        return GOTerm(
            id, name=name, desc=desc, alt_id=alts, 
            is_a=is_a, loci=term_loci
        )


    def add_term(self, term, cursor=None, overwrite=False):
        ''' 
            Add a single term to the ontology 
        
        '''
        self.__getitem__.cache_clear()
        if overwrite:
            self.del_term(term.id)
        if not cursor:
            cur = self.db.cursor()
            cur.execute('BEGIN TRANSACTION')
        else:
            cur = cursor

        # Add the term name and description
        cur.execute(
            'INSERT OR REPLACE INTO terms (id, desc, name) VALUES (?, ?, ?)', 
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

    def parents(self, term):
        '''
            Return an iterable containing the parents of a term.
            Parents are determined via the is_a property of the term.

            Parameters
            ----------
            term : GOTerm

            Returns
            -------
            An iterable containing the parents

            NOTE: note guaranteed to be in order NOR guaranteed 
            to not contain duplicates
        '''
        for parent_term in term.is_a:
            yield self[parent_term]
            for grand_parent in self.parents(self[parent_term]):
                yield grand_parent

    def graph(self, terms=None):
        '''
            Create a NetworkX graph from terms
        '''
        # Create the graph object
        if terms == None:
            terms = self.iter_terms()
        else:
            # terms include 
            terms = list(chain(*[self.parents(term) for term in terms]))

        # 
        G = nx.Graph()
        G.add_nodes_from(set(terms))
        # build a list of all the edges
        edges = []
        for term in terms:
            for parent in term.is_a:
                edges.append((term.id,parent))
        G.add_edges_from(edges)
        return G


    @classmethod
    def create(cls, name, description, refgen, overwrite=True, type='GOnt'):
        '''This method creates a fresh GO Ontology with nothing it it.'''

        # run the inherited create method from Camoco
        self = super().create(name, description, refgen, type=type)

        if overwrite:
            self._drop_indices()
            self._clear_tables()

        return self

    ''' --------------------------------------------------------------------
        Internal Methods
     --------------------------------------------------------------------'''

    def _parse_obo(self,obo_file):
        '''
            Parses a GO obo file

            Paramters
            ---------
            obo_file : filename
                The path the the obo file. You can download
                it here: http://geneontology.org/page/download-ontology
            
            Returns
            -------
            A list of Empty (i.e. no genes) GO Terms 
        '''
        # Importing the obo information
        self.log('Importing OBO: {}', obo_file)
        cur_term = None
        terms = []
        isa_re = re.compile('is_a: (.*) !.*')
        with rawFile(obo_file) as INOBO:
            for line in INOBO:
                line = line.strip()
                if line.startswith('id: '):
                    if cur_term is not None:
                        terms.append(cur_term)
                    GOid = line.replace('id: ', '')
                    cur_term = GOTerm(GOid)
                elif line.startswith('name: '):
                    GOname = line.replace('name: ', '')
                    cur_term.name = GOname
                elif line.startswith('namespace: '):
                    cur_term.attrs['namespace'] = line.replace('namespace: ', '')
                elif line.startswith('alt_id: '):
                    alt_id = line.replace('alt_id: ', '')
                    cur_term.alt_id.add(alt_id)
                elif line.startswith('def: '):
                    cur_term.desc += line.replace('def: ', '')
                elif line.startswith('comment: '):
                    cur_term.desc += line.replace('comment: ', '')
                elif line.startswith('is_a: '):
                    cur_term.is_a.add(isa_re.match(line).group(1))
        terms.append(cur_term)
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
                row = line.strip('\n').split('\t')
                gene = row[id_col].split('_')[0].strip()
                cur_term = row[go_col]
                # Make a map between genes and associated GO terms
                if gene not in genes:
                    genes[gene] = set([cur_term])
                else:
                    genes[gene].add(cur_term)
        return genes

    @classmethod
    def from_terms(cls, terms, name, description, refgen):
        '''
            Convenience function to create a GOnt from an iterable
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
        self.add_terms(terms, overwrite=False)
        # Build the indices
        self.log('Building the indices.')
        self._build_indices()

        self.log('Your gene ontology is built.')
        return self

         

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
        # Parse the OBO file and get a list of empty terms
        terms = {term.id : term for term in self._parse_obo(obo_file)}
        # Add terms to the database so that the parent method works
        self.add_terms(terms.values(), overwrite=False)
        # Build the indices
        self.log('Building the indices.')
        self._build_indices()
        # Parse the Gene/Term mapping 
        genes = self._parse_gene_term_map(gene_map_file,
            headers=headers,go_col=go_col,id_col=id_col
        )
        # Assign genes to terms AND to parent terms
        self.log("Adding GO-gene assignments")
        for gene_id,term_ids in genes.items():
            # Get a gene object from the refgen
            try:
                gene = self.refgen[gene_id]
                for term_id in term_ids:
                    # Add gene to each term its annotated to
                    if term_id not in terms:
                        self.log("{} not in Ontology",term_id)
                        continue
                    terms[term_id].loci.add(gene)
                    # Propogate gene to each parental term
                    for parent in self.parents(terms[term_id]):
                        terms[parent.id].loci.add(gene)
            except ValueError as e:
                pass
        self.add_terms(terms.values(), overwrite=False)
        self.log('Build Sucessful.')
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
