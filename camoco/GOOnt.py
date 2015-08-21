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

class GOOnt(Ontology):
    '''Ontology extension for GO'''
    def __init__():


    @classmethod
    def create(cls, name, description, refgen, type='GOOnt'):
        '''
            This method creates a fresh Ontology with nothing it it.
        '''
        # run the inherited create method from Camoco
        self = super().create(name, description, refgen, type=type)
        # build the tables
        self._create_tables()
        return self



    @classmethod
    def from_obo(cls,obo_file,gene_map_file,name,description,refgen,go_col=1):
        ''' Convenience function for importing GO obo files '''

        # Internal function to handle propagating is_a relationships
        def getISA(terms, term):
            tot = set()
            if terms[term]['is_a'] == []:
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
        alt_terms = []
        isa_re = re.compile('is_a: (.*) !.*')
        with open(obo_file, 'r') as INOBO:
            for line in INOBO:
                line = line.strip()
                if line.startswith('id: '):
                    for alt in alt_terms:
                        terms[alt] = terms[cur_term].copy()
                    alt_terms = []
                    cur_term = line.replace('id: ', '')
                elif line.startswith('name: '):
                    terms[cur_term]['name'] = line.replace('name: ', '')
                    terms[cur_term]['desc'] = ''
                    terms[cur_term]['is_a'] = []
                    terms[cur_term]['genes'] = set()
                elif line.startswith('namespace: '):
                    terms[cur_term]['type'] = line.replace('namespace: ', '')
                elif line.startswith('alt_id: '):
                    alt_terms.append(line.replace('alt_id: ', ''))
                elif line.startswith('def: '):
                    terms[cur_term]['desc'] += line.replace('def: ', '')
                elif line.startswith('comment: '):
                    terms[cur_term]['desc'] += line.replace('comment: ', '')
                elif line.startswith('is_a: '):
                    terms[cur_term]['is_a'].append(isa_re.match(line).group(1))

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

    def _create_tables(self):
        super()._create_tables():
        cur = self.db.cursor()
        cur.execute('''
            CREATE TABLE IF NOT EXISTS rels (
                parent TEXT,
                child TEXT,
            );
        ''')

def addGenes(terms, term, gene):
    x = getISA(terms, term)
    for item in x:
        terms[item]['genes'].add(gene)
    return

def getISA(terms, term):
    tot = set([term])
    if terms[term]['is_a'] == []:
        return tot
    else:
        for item in terms[term]['is_a']:
            tot.add(item)
            tot = tot | getISA(terms, item)
        return tot
