#!/usr/bin/python

from Camoco import Camoco
from collections import defaultdict
from pandas import DataFrame
from scipy.stats import hypergeom
import re

class Ontology(Camoco):
    def __init__(self,name,basedir="~/.camoco"):
        super(self.__class__,self).__init__(name,type='Ontology',basedir=basedir)

    def term_genes(self,term):
        return [x[0] for x in self.db.cursor().execute("SELECT gene from gene_terms WHERE term = ?",(term,)).fetchall()]

    def terms(self):
        return self.db.cursor().execute("SELECT * from terms").fetchall()

    def enrichment(self,gene_list,pval_cutoff=0.05):
        # extract possible terms for genes
        cur = self.db.cursor()
        terms = [ x[0] for x in cur.execute(
            '''SELECT DISTINCT(term) FROM gene_terms 
            WHERE gene IN ('{}');'''.format("','".join([x.id for x in gene_list]))
        )]
        # computer hypergeometric for each term
        enrichment = []
        for id in terms:
            try:
                (id,name,type,desc) = cur.execute("SELECT * FROM terms WHERE id = ?",(id,)).fetchone()
            except TypeError as e:
                self.log("No information for ontology term {}",id)
            genes_in_term = [x[0] for x in cur.execute(
                '''SELECT gene FROM gene_terms WHERE term = ?''',(id,))
            ]
            num_genes_in_term = len(genes_in_term)
            overlap = set(genes_in_term).intersection(set([x.id for x in gene_list]))
            num_genes_total, = cur.execute('SELECT COUNT(DISTINCT(gene)) FROM gene_terms;').fetchone()
            pval = hypergeom.sf(len(overlap)-1,num_genes_total,num_genes_in_term,len(gene_list))
            enrichment.append(
                (id,name,pval,num_genes_in_term,len(overlap),len(gene_list),num_genes_total,type,desc)
            )
        try:
            enrichment = DataFrame(enrichment,
                columns = ['TermID','Name','pval','LenTerm','LenOverlap','LenList','LenTotal','Type','Desc']
            ).sort('pval',ascending=True)
            enrichment.index = enrichment.TermID
        except ValueError as e:
            self.log("No enrichment for {}",",".join([x.id for x in gene_list]))
            return DataFrame()
        return enrichment[enrichment.pval <= pval_cutoff]

    def term(self,id):
        return self.db.cursor().execute('SELECT * FROM terms WHERE id = ?',(id,)).fetchall()
            

class OntologyBuilder(Camoco):
    def __init__(self,name,description='',basedir="~/.camoco"):
        super().__init__(name,description=description,type="Ontology",basedir=basedir)
        self._create_tables()

    def import_gene_terms(self, filename, gene_col=1, term_col=2, term_filter=".*",
        gene_filter=".*", skip=0, sep="\t"):
        ''' import tool for gene terms ''' 
        term_filter = re.compile(term_filter)
        gene_filter = re.compile(gene_filter)
        gene_terms = []
        self.log("Reading in term file {}",filename)
        with open(filename,'r') as IN:
            for x in range(0,skip):
                header = IN.readline()
            for line in IN:
                term = ''
                gene = ''
                cols = line.strip().split(sep)
                tmatch = term_filter.match(cols[term_col - 1])
                if tmatch is None:
                    continue
                elif len(tmatch.groups()) == 0:
                    term =  tmatch.string
                else:
                    term = tmatch.group(1)
                gmatch = gene_filter.match(cols[gene_col - 1])
                if gmatch is None:
                    continue
                elif len(gmatch.groups()) == 0:
                    gene =  gmatch.string
                else:
                    gene = gmatch.group(1)
                gene_terms.append((gene,term))
        self.log("Inserting {} gene term pairs",len(gene_terms))
        cur = self.db.cursor()
        cur.execute('BEGIN TRANSACTION')
        cur.executemany(''' 
            INSERT INTO gene_terms VALUES (?,?)''', gene_terms
        )
        cur.execute('END TRANSACTION')

    def import_obo(self,filename):
        ''' Convenience function for importing GO obo files '''
        self.log('importing OBO: {}',filename)
        terms= defaultdict(dict)
        is_a = list()
        cur_term = ''
        isa_re = re.compile('is_a: (.*) !.*')
        with open(filename,'r') as INOBO:
            for line in INOBO:
                line = line.strip()
                if line.startswith('id: '):
                    cur_term = line.replace('id: ','') 
                elif line.startswith('name: '):
                    terms[cur_term]['name'] = line.replace('name: ','')
                    terms[cur_term]['desc'] = ''
                elif line.startswith('namespace: '):
                    terms[cur_term]['type'] = line.replace('namespace: ','')
                elif line.startswith('def: '):
                    terms[cur_term]['desc'] += line.replace('def: ','')
                elif line.startswith('comment: '):
                    terms[cur_term]['desc'] += line.replace('comment: ','')
                elif line.startswith('is_a: '):
                    is_a.append((cur_term,isa_re.match(line).group(1)))
        self.log("Dumping {} annotations and {} relationships",len(terms),len(is_a))
        cur = self.db.cursor()
        cur.execute('BEGIN TRANSACTION')
        cur.executemany('''
            INSERT INTO terms (id,name,type,desc) VALUES(?,?,?,?)''', 
            [ (key,val['name'],val['type'],val['desc']) for key,val in terms.items()]
        )
        self.log('Done inserting terms')
        cur.executemany(''' 
            INSERT INTO relationships (term,is_a) VALUES (?,?)''',
            is_a
        )
        cur.execute('END TRANSACTION')
        self._build_indices()

    def import_mapman(self,filename):
        ''' Convenience function for files provided by MapMan, columns are 
            CODE,NAME,Gene,DESC,TYPE seperated by space and enclosed in single quotes'''
        self.log('Importing MAPMAN text file: {}',filename)
        terms = dict()
        is_a = dict()
        gene_terms = list()
        transcript_strip = re.compile("_T\d+$")
        is_a_pattern = re.compile('\.\d+$')
        with open(filename,'r') as INMM:
            headers = INMM.readline()
            for line in INMM:
                # the map just takes out leading/trailing single quotes
                (term,name,gene,desc,*type) = [x.strip("'") for x in  line.strip().split("\t")]
                # strip transcript out of gene name
                gene = transcript_strip.sub('',gene.upper())
                terms[term] = (term,name,'','') 
                gene_terms.append((gene,term))
                # add if there is a relationship there
                if is_a_pattern.match(term):
                    is_a[term] = is_a_pattern.sub('',term)
        self.log("Dumping {} terms and {} gene-terms",len(terms),len(gene_terms))
        cur = self.db.cursor()
        cur.execute('BEGIN TRANSACTION')
        cur.executemany('''INSERT INTO terms (id,name,type,desc) VALUES (?,?,?,?)''',terms.values())
        cur.executemany('''INSERT INTO relationships (term,is_a) VALUES (?,?) ''',is_a.items()) 
        cur.executemany('''INSERT INTO gene_terms (gene,term) VALUES (?,?)''',gene_terms)
        cur.execute('END TRANSACTION')
        

    def _build_indices(self):
        cur = self.db.cursor()
        cur.execute('''CREATE INDEX IF NOT EXISTS termid ON terms (id)''')
        cur.execute('''CREATE INDEX IF NOT EXISTS termtype ON terms (type)''')
        cur.execute('''CREATE INDEX IF NOT EXISTS relsource ON relationships (term)''')
        cur.execute('''CREATE INDEX IF NOT EXISTS reltarget ON relationships (is_a)''')
        cur.execute('''CREATE INDEX IF NOT EXISTS gene_terms_gene ON gene_terms (gene)''')
        cur.execute('''CREATE INDEX IF NOT EXISTS gene_terms_term ON gene_terms (term)''')

    def _create_tables(self):
        cur = self.db.cursor()
        cur.execute("PRAGMA page_size = 1024;")
        cur.execute("PRAGMA cache_size = 100000;")
        cur.execute(''' 
            CREATE TABLE IF NOT EXISTS terms (
                id TEXT UNIQUE,
                name TEXT,
                type TEXT,
                desc TEXT   
            ); 
        ''')
        cur.execute(''' 
            CREATE TABLE IF NOT EXISTS relationships (
                term TEXT,
                is_a TEXT
            ) 
        ''')
        cur.execute(''' 
            CREATE TABLE IF NOT EXISTS gene_terms (
                gene TEXT,
                term TEXT
            );            
        ''')
