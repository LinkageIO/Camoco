#!/usr/bin/python

from camoco.Camoco import Camoco
from camoco.RefGen import RefGen
from camoco.Tools import log
from camoco.Locus import Locus


from collections import defaultdict
from pandas import DataFrame
from scipy.stats import hypergeom
import matplotlib.pylab as plt
import itertools
import pandas as pd
import re
import sys

class Term(object):
    '''
        A Term is a just group of loci that are related.
    '''
    def __init__(self,name,desc='',locus_list=None,**kwargs):
        self.name = name
        self.desc = desc
        self.locus_list = set()
        self.attr = {}
        if locus_list:
            self.locus_list = set(locus_list)
        for key,val in kwargs.items():
            self.attrs[key] = val

    def __len__(self):
        ''' 
            Returns the number of loci in the term.
        '''
        return len(self.locus_list)

    def add_locus(self,locus):
        '''
            Adds a locus to the Term.
        '''
        self.locus_list.add(locus) 
    
    def flanking_snps(self,gene,window_size=100000):
        ''' 
            returns any nearby Term SNPs to a gene 
        '''
        return [locus for locus in self.locus_list if abs(gene-locus) <= window_size]

    def effective_snps(self,window_size=None,max_genes_between=1):
        ''' 
            Collapse down loci that have overlapping windows.
            Also collapses down snps that have 
        '''
        locus_list = sorted(self.locus_list)
        if window_size is not None:
            for locus in locus_list:
                locus.window = window_size
        collapsed = [locus_list.pop(0)]
        for locus in locus_list:
            # if they have overlapping windows, collapse
            if locus in collapsed[-1]:
                # Collapse if the windows overlap
                collapsed[-1] = collapsed[-1] + locus
            else:
                collapsed.append(locus)
        log('{}: Found {} SNPs -> {} effective SNPs',self.name,len(self.locus_list),len(collapsed))
        return collapsed
        
    def __str__(self):
        return "Term: {}, Desc: {}, {} Loci".format(self.name,self.desc,len(self))

    def __repr__(self):
        return str(self.name)

class Ontology(Camoco):
    ''' An Ontology is just a collection of terms. Each term is just a collection of genes. 
        Sometimes terms are related or nested within each other, sometimes not. Simple enough.  
    '''
    def __init__(self,name):
        super().__init__(name,type='Ontology')
        if self.refgen:
            self.refgen = RefGen(self.refgen)

    def __len__(self):
        return self.db.cursor().execute("SELECT COUNT(*) FROM terms;").fetchone()[0]

    def __getitem__(self,name):
        ''' retrieve a term by name '''
        try:
            desc = self.db.cursor().execute(
                'SELECT desc from terms WHERE name = ?',(name,)
            ).fetchone()
            term_loci = [Locus.from_record(record) for record in self.db.cursor().execute(
                'SELECT chrom,start,end,name,window,id FROM term_loci WHERE term = ?',(name,)
            ).fetchall()]
            return Term(name,desc,locus_list=term_loci)
        except TypeError as e: # Not in database
            raise

    def term_like(self,like="%"):
        return [self.term(x[0]) for x in self.db.cursor().execute('SELECT id FROM terms WHERE name LIKE ?',(like,)).fetchall()]

    def iter_terms(self):
        for id, in self.db.cursor().execute("SELECT name FROM terms"):
            yield self[id]

    def terms(self):
        return list(self.iter_terms())

    def enrichment(self,gene_list,pval_cutoff=0.05,gene_filter=None,label=None,max_term_size=300):
        # extract possible terms for genes
        if label:
            self.log("Caculating Enrichemnt for {}",label)
        cur = self.db.cursor()
        terms = [ x[0] for x in cur.execute(
            '''SELECT DISTINCT(term) FROM gene_terms 
            WHERE gene IN ('{}');'''.format("','".join([x.id for x in gene_list]))
        )]
        # compute hypergeometric for each term
        enrichment = []
        for id in terms:
            try:
                (id,name,type,desc) = cur.execute("SELECT * FROM terms WHERE id = ?",(id,)).fetchone()
            except TypeError as e:
                self.log("No information for ontology term {}",id)
            genes_in_term = [x[0] for x in cur.execute(
                '''SELECT gene FROM gene_terms WHERE term = ?''',(id,))
            ]
            if len(genes_in_term) > max_term_size:
                self.log("Skipping {} due to size ({})",name,len(genes_in_term))
                continue
            if gene_filter:
                genes_in_term = [gene for gene in genes_in_term if gene in gene_filter]
            num_genes_in_term = len(genes_in_term)
            overlap = set(genes_in_term).intersection(set([x.id for x in gene_list]))
            num_genes_total, = cur.execute('SELECT COUNT(DISTINCT(gene)) FROM gene_terms;').fetchone()
            pval = hypergeom.sf(len(overlap)-1,num_genes_total,num_genes_in_term,len(gene_list))
            term_genes = ",".join(overlap)
            enrichment.append(
                (id,name,pval,num_genes_in_term,len(overlap),len(gene_list),num_genes_total,type,term_genes,desc)
            )
        try:
            enrichment = DataFrame(enrichment,
                columns = ['TermID','Name','pval','LenTerm','LenOverlap','LenList','LenTotal','Type','TermGenes','Desc']
            ).sort('pval',ascending=True)
            enrichment.index = enrichment.TermID
        except ValueError as e:
            self.log("No enrichment for {}",",".join([x.id for x in gene_list]))
            return DataFrame()
        if label:
            enrichment['Label'] = label
        return enrichment[enrichment.pval <= pval_cutoff]

    def print_term_stats(self, cob_list, filename=None, window_size=100000, gene_limit=4,num_bootstrap=50,bootstrap_density=2):
        for term in self.iter_terms():
            term.print_stats(cob_list,filename,window_size=window_size,gene_limit=gene_limit,num_bootstraps=num_bootstraps,bootstrap_density=boostrap_density)

    def summary(self):
        return "Ontology:{} - desc: {} - contains {} terms for {}".format(self.name,self.description,len(self),self.refgen)

    def del_term(self,name):
        '''
        Remove a term from the dataset.

        Parameters
        ----------
        name : string
            The term name you wish to remove.

        Returns
        -------
            bool indicating success
        '''
        cur = self.db.cursor()
        cur.execute('''
            BEGIN TRANSACTION;
            DELETE FROM loci_attr WHERE loci_id IN (
                SELECT id FROM term_loci WHERE term = ?
            );
            DELETE FROM term_loci WHERE term = ?;
            DELETE FROM terms WHERE name = ?;
            END TRANSACTION;
        ''',(name,name,name))

    def add_term(self,term,overwrite=True):
        ''' This will add a single term to the ontology '''
        cur = self.db.cursor()
        if overwrite:
            self.del_term(term.name)
        cur.execute('BEGIN TRANSACTION')
        # Add the term name and description
        cur.execute('''
            INSERT OR REPLACE INTO terms (name,desc) 
            VALUES (?,?)''',(term.name,term.desc)
        )
        # Add the term loci
        for locus in term.locus_list:
            cur.execute('''
                INSERT OR REPLACE INTO term_loci (term,chrom,start,end,name,window,id)
                VALUES (?,?,?,?,?,?,?)
                ''',(term.name,) + locus.as_record()
            )
            # Add the loci attrs
            cur.executemany('''
                INSERT OR REPLACE INTO loci_attr (term,loci_id,key,val)
                VALUES (?,?,?,?)
            ''',[(term.name,locus.id,key,val) for key,val in locus.attr.items()])
        cur.execute('END TRANSACTION')

    @classmethod
    def create(cls,name,description,refgen):
        ''' 
            This method creates a fresh Ontology with nothing it it. 
        '''
        # run the inherited create method from Camoco
        self = super(Ontology,cls).create(name,description,type='Ontology')
        # add the global refgen
        self._global('refgen',refgen.name)
        # set the refgen for the current instance
        self.refgen = refgen
        # build the tables
        self._create_tables()
        # create indices
        self._build_indices()
        return self

    @classmethod
    def from_DataFrame(cls,df,name,description,refgen,
            term_col='Term', chr_col=None, pos_col=None,
            start_col=None, end_col=None, id_col=None
            ):
        '''
            Import an ontology from a pandas dataframe. 
            Groups by term_col, then iterates over the rows in the group.
            It adds each row as a locus to the term.
        '''
        self = cls.create(name,description,refgen)
        # group each trait by its name
        for term_name,df in df.groupby(term_col):
            term = Term(term_name)
            # we have a SNP
            if pos_col is not None:
                for i,row in df.iterrows():
                    # make sure there are no collisions with the SNP init func names
                    # this is hackey and I dont like it
                    kwargs = {key:val for key,val in dict(row).items() if key not in Locus.__init__.__code__.co_varnames}
                    snp = Locus(row[chr_col],int(row[pos_col]),int(row[pos_col]),gene_build=self.refgen.build,**kwargs)
                    term.locus_list.add(snp) 
            self.log("Importing {}",term)
            self.add_term(term)
        return self
            

    @classmethod
    def from_obo(self,filename,name,description,refgen):
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

    @classmethod
    def from_mapman(self,filename):
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
        pass    

    def _create_tables(self):
        cur = self.db.cursor()
        cur.execute(''' 
            CREATE TABLE IF NOT EXISTS terms (
                name TEXT UNIQUE,
                desc TEXT   
            ); 
        ''')
        cur.execute(''' 
            CREATE TABLE IF NOT EXISTS loci_attr (
                term TEXT,
                loci_id TEXT,
                key TEXT,
                val TEXT,
                PRIMARY KEY(loci_id,key)
            );
        ''')
        cur.execute('''
            CREATE TABLE IF NOT EXISTS term_loci (
                term TEXT,
                chrom TEXT,
                start INT,
                end INT,
                name TEXT,
                window INT,
                id TEXT,
                PRIMARY KEY(term,id)
            );
        ''')
        self._build_indices()
    
