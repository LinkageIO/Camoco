#!/usr/bin/python

from camoco.Camoco import Camoco
from camoco.RefGen import RefGen
from camoco.Locus import SNP
from camoco.Tools import log

from collections import defaultdict
from pandas import DataFrame
from scipy.stats import hypergeom
import matplotlib.pylab as plt
import itertools
import pandas as pd
import re
import sys



class Term(object):
    def __init__(self,id,name='',type='',desc='',gene_list=None,snp_list=None):
        self.id = id
        self.name = name
        self.type = type
        self.desc = desc
        self.gene_list = set()
        self.snp_list = set()
        if gene_list:
            self.gene_list = set(gene_list)
        if snp_list:
            self.snp_list = set(snp_list)

    def __len__(self):
        return len(self.gene_list)

    def summary(self):
        print("\n".join([ x for x  in [self.id,self.name,self.type,self.desc] if x is not '']))
        print("Num SNPs: {}".format(len(self.snp_list)))
        print("Num Genes: {}".format(len(self.gene_list)))

    def snps(self,effective=True,window=50000):
        pass

    def add_gene(self,gene):
        self.gene_list.add(gene)

    def add_snp(self,snp):
        self.snp_list.add(snp) 
    
    def flanking_genes(self,refgen,gene_limit=4,chain=True,window_size=None):
        ''' returns flanking genes based on some set of arbitrary rules from a refgen '''
        if chain:
            return set(itertools.chain.from_iterable(
                [refgen.flanking_genes(x,gene_limit=gene_limit,window_size=window_size if window_size is not None else x.window) for x in self.snp_list]
            ))
        else:
            return [refgen.flanking_genes(x,gene_limit=gene_limit,window_size=window_size) for x in self.snp_list]

    def flanking_snps(self,gene,window_size=100000):
        ''' returns any nearby Term SNPs to a gene '''
        return [snp for snp in self.snp_list if abs(gene-snp) <= window_size]

    def bootstrap_flanking_genes(self,refgen,window_size=100000,gene_limit=4,chain=True):
        ''' returns random flanking genes, with similar properties, based on an arbitrary set of rules'''
        if chain:
            return set(itertools.chain(*[refgen.bootstrap_flanking_genes(x,gene_limit=gene_limit,window_size=window_size) for x in self.snp_list]))
        else:   
            return [refgen.bootstrap_flanking_genes(x,gene_limit=gene_limit,window_size=window_size) for x in self.snp_list]

    def effective_snps(self):
        ''' Collapse down snps that have overlapping windows'''
        snp_list = sorted(self.snp_list)
        collapsed = [snp_list.pop(0)]
        for snp in snp_list:
            # if they have overlapping windows, collapse
            if snp in collapsed[-1]:
                collapsed[-1] = collapsed[-1].collapse(snp)
            else:
                collapsed.append(snp)
        return collapsed

    def bootstrap_effective_flanking_genes(self,refgen,chain=True,window_size=10000,gene_limit=4):
        ''' Sometimes, especially with GWAS, SNPs will cluster and map to the same
        flanking gene set. This method accounts for this by collapsing down SNPs
        into meta-SNPs which have varying window sizes.'''
        # bootstrap flanking genes with specific window sizes
        if chain:
            return set(itertools.chain(*[
                refgen.bootstrap_flanking_genes(snp,gene_limit=gene_limit,window_size=snp.window) for snp in collapsed
            ]))
        else:   
            return [refgen.bootstrap_flanking_genes(snp,gene_limit=gene_limit,window_size=snp.window) for snp in collapsed]

    def print_stats(self,cob_list,file=sys.stdout, window_limit=100000, gene_limit=4,num_bootstrap=50,bootstrap_density=2):
        for cob in cob_list:
            # MAKE SURE WE ARE DEALING WITH A SET!!!!!!!!!! Lists will have duplicates in it!
            flanks = self.flanking_genes(cob.refgen,window_limit=window_limit,gene_limit=gene_limit)
            log("On Term {} with {}. {}/{} genes from {} SNPs",self.id,cob.name,len(flanks),len(self.gene_list),len(self.snp_list))
            density  = cob.density(flanks,min_distance=window_size)
            locality = cob.locality(flanks,min_distance=window_size)
            len_LCC = len(cob.lcc(flanks,min_distance=window_size).vs)
            print("{}\t{}_NumSNPS\t{}".format(self.id,cob.name,len(self.snp_list)),file=file)
            print("{}\t{}_NumGenes\t{}".format(self.id,cob.name,len(flanks)),file=file)
            print("{}\t{}_TransDensity\t{}".format(self.id,cob.name,density),file=file)
            print("{}\t{}_Locality\t{}".format(self.id,cob.name,locality),file=file)
            print("{}\t{}_LCC\t{}".format(self.id,cob.name,len_LCC),file=file)
            if density > bootstrap_density:
                log("Density > 2; Boostrapping!")
                # Calculate BootStrap 
                bs_density = []
                bs_local = []
                bs_lcc = []
                for x in range(num_bootstrap):
                    bs_flanks = list(chain.from_iterable(
                        [cob.refgen.bootstrap_flanking_genes(x,gene_limit=gene_limit,window_size=window_size) for x in self.snp_list]
                    ))
                    bs_density.append(cob.density(bs_flanks,min_distance=window_size))           
                    bs_local.append(cob.locality(bs_flanks,min_distance=window_size))
                    bs_lcc.append(len(cob.lcc(bs_flanks,min_distance=window_size).vs))
                print("{}\t{}_BS_TransDensity\t{}".format(self.id,cob.name,sum([x >= density for x in bs_density])),file=file)
                print("{}\t{}_BS_Locality\t{}".format(self.id,cob.name,sum([x >= locality for x in bs_local])),file=file)
                print("{}\t{}_BS_LCC\t{}".format(self.id,cob.name,sum([x >= len_LCC for x in bs_lcc])),file=file)

    def __str__(self):
        return "Term: {}, {} genes, {} SNPs".format(self.id,len(self.gene_list),len(self.snp_list))

    def __repr__(self):
        return str(self.id)

class Ontology(Camoco):
    ''' An Ontology is just a collection of terms. Each term is just a collection of genes. 
        Sometimes terms are related or nested within each other, sometimes not. Simple enough.  
    '''
    def __init__(self,name,description=None,basedir="~/.camoco"):
        super().__init__(name,description,type='Ontology',basedir=basedir)
        if self.refgen:
            self.refgen = RefGen(self.refgen)

    @classmethod
    def create(cls,name,description,refgen,basedir="~/.camoco"):
        self = cls(name,description=description,basedir=basedir)
        self._global('refgen',refgen.name)
        self.refgen = RefGen(self.refgen)
        self._create_tables()
        self._build_indices()
        return self

    def __getitem__(self,item):
        return self.term(item) 

    def __len__(self):
        return self.db.cursor().execute("SELECT COUNT(*) FROM terms;").fetchone()[0]


    def term(self,term_id,window_size=100000):
        ''' retrieve a term by name '''
        try:
            id,name,type,desc = self.db.cursor().execute(
                'SELECT id,name,type,desc from terms WHERE id = ?',(term_id,)
            ).fetchone()
            term_genes = list(self.refgen.from_ids([ x[0] for x in self.db.cursor().execute(
                'SELECT gene from gene_terms WHERE term = ?',(id,)).fetchall()
            ]))
            term_snps = [SNP(*x,window=window_size) for x in self.db.cursor().execute(
                'SELECT chrom,pos FROM snp_terms WHERE term = ?',(id,)
            ).fetchall()]
            return Term(id,name,type,desc,gene_list=term_genes,snp_list=term_snps)
        except TypeError as e: # Not in database
            return None

    def term_ids(self,like="%"):
        return [self.term(x[0]) for x in self.db.cursor().execute('SELECT id FROM terms WHERE id LIKE ?',(like,)).fetchall()]

    def term_names(self,like="%"):
        return [self.term(x[0]) for x in self.db.cursor().execute('SELECT id FROM terms WHERE name LIKE ?',(like,)).fetchall()]

    def iter_terms(self):
        for id, in self.db.cursor().execute("SELECT id FROM terms"):
            yield self.term(id)

    def terms(self):
        return list(self.iter_terms())

    def compare_COBs(self,COB1,COB2,min_term_size=10,max_term_size=300,plot=True):
        # get pairwise densities for terms
        self.log("Extracting Densities")
        densities = pd.DataFrame([
            (COB1.density(term.gene_list),COB2.density(term.gene_list)) for term in self.iter_terms() if len(term) > min_term_size and len(term) < max_term_size
        ],columns=['x','y'])
        if plot:
            self.log('Plotting')
            plt.clf()
            plt.figure(figsize=(8,8),dpi=200)
            plt.scatter(densities.x,densities.y)
            plt.xlim([0,150])
            plt.ylim([0,150])
            plt.xlabel(COB1.name + " Densities")
            plt.ylabel(COB2.name + " Densities")
            plt.savefig("{}_vs_{}_{}.png".format(COB1.name,COB2.name,self.name))
        self.log("Done")
        return densities

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
        return "Ontology: name:{} - desc: {} - contains {} terms for {}".format(self.name,self.description,len(self),self.refgen)
   

    def del_term(self,term):
        cur = self.db.cursor()
        cur.execute('DELETE FROM terms WHERE id = ?',(term.id,))
        cur.execute('DELETE FROM gene_terms WHERE term = ?',(term.id,))
        cur.execute('DELETE FROM snp_terms WHERE term = ?',(term.id,))

    def add_term(self,term,overwrite=True):
        ''' This will add a single term to the ontology '''
        cur = self.db.cursor()
        if overwrite:
            self.del_term(term)
        cur.execute('BEGIN TRANSACTION')
        # Add info to the terms tabls
        cur.execute('INSERT OR REPLACE INTO terms (id,name,type,desc) VALUES (?,?,?,?)',(term.id,term.name,term.type,term.desc))
        for gene in term.gene_list:
            cur.execute('INSERT OR REPLACE INTO gene_terms (gene,term) VALUES (?,?)',(gene.id,term.id))
        for snp in term.snp_list:
            cur.execute('INSERT OR REPLACE INTO snp_terms (chrom,pos,term) VALUES (?,?,?)',(snp.chrom, snp.pos, term.id))
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
        cur.execute('''
            CREATE TABLE IF NOT EXISTS snp_terms (
                chrom TEXT,
                pos TEXT,
                id TEXT,
                term TEXT
            );
        ''')
        self._build_indices()
    
