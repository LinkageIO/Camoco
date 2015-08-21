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
    def __init__(self, id, name='', desc='', alt_id=None, is_a=None,locus_list=None, **kwargs):
        super().__init__(id, desc=desc, locus_list=locus_list, **kwargs)
        self.name = name
        self.is_a = set(is_a) if is_a else set()
        self.alt_id = set(alt_id) if alt_id else set()

    def __str__(self):
        return "Term: {}, Name: {}, Desc: {}, {} Loci".format(self.id, self.name, self.desc, len(self))

    def add_parent(self,termname):
        self.is_a.add(termname)

class GOOnt(Ontology):
    '''Ontology extension for GO'''
    def __init__(self, name, type='GOOnt'):
        super().__init__(name, type=type)
        if self.refgen:
            self.refgen = RefGen(self.refgen)

    @classmethod
    def create(cls, name, description, refgen, type='GOOnt'):
        '''
            This method creates a fresh GO Ontology with nothing it it.
        '''
        # run the inherited create method from Camoco
        self = super().create(name, description, refgen, type=type)

        # Alter the tables as needed
        cur = self.db.cursor()
        cur.execute('ALTER TABLE terms ADD COLUMN name TEXT;')
        cur.execute('CREATE TABLE IF NOT EXISTS rels (parent TEXT, child TEXT);')
        cur.execute('CREATE TABLE IF NOT EXISTS alts (primary TEXT, alt TEXT UNIQUE);')
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
            parents = getISA(terms,cur_term)
            terms[cur_term]['is_a'] = parents

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
