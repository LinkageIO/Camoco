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
import numpy as np

class MapOnt(Ontology):
    '''Ontology extension for MapMan'''


    @classmethod
    def from_mapman(self, filename):
        '''
            Convenience function for files provided by MapMan, columns are
            CODE, NAME, Gene, DESC, TYPE seperated by space and
            enclosed in single quotes
        '''
        self.log('Importing MAPMAN text file: {}', filename)
        terms = dict()
        is_a = dict()
        gene_terms = list()
        transcript_strip = re.compile("_T\d+$")
        is_a_pattern = re.compile('\.\d+$')
        with open(filename, 'r') as INMM:
            headers = INMM.readline()
            for line in INMM:
                # the map just takes out leading/trailing single quotes
                (term, name, gene, desc, *type) = [x.strip("'") for x in  line.strip().split("\t")]
                # strip transcript out of gene name
                gene = transcript_strip.sub('', gene.upper())
                terms[term] = (term, name, '', '')
                gene_terms.append((gene, term))
                # add if there is a relationship there
                if is_a_pattern.match(term):
                    is_a[term] = is_a_pattern.sub('', term)
        self.log("Dumping {} terms and {} gene-terms", len(terms), len(gene_terms))
        cur = self.db.cursor()
        cur.execute('BEGIN TRANSACTION')
        cur.executemany('''INSERT INTO terms (id, name, type, desc) VALUES (?, ?, ?, ?)''', terms.values())
        cur.executemany('''INSERT INTO relationships (term, is_a) VALUES (?, ?) ''', is_a.items())
        cur.executemany('''INSERT INTO gene_terms (gene, term) VALUES (?, ?)''', gene_terms)
        cur.execute('END TRANSACTION')
