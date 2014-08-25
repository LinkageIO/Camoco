#!/usr/bin/python3

from collections import defaultdict

from Camoco import Camoco
from Locus import Gene,non_overlapping
from Chrom import Chrom
from Genome import Genome
from Tools import memoize

import itertools
import random

class RefGen(Camoco):
    def __init__(self,name,basedir="~/.camoco"):
        # initialize camoco instance
        super().__init__(name,type="RefGen",basedir=basedir)
        self.genome = Genome(
            self.type+self.name,
            chroms = [Chrom(*x) for x in  self.db.cursor().execute(''' 
                SELECT id,length FROM chromosomes
            ''')] 
        )
    @memoize
    def num_genes(self):
        ''' returns the number of genes in the dataset '''
        return self.db.cursor().execute(''' SELECT COUNT(*) FROM genes''').fetchone()[0]

    def random_gene(self):
        return Gene(*self.db.cursor().execute(''' 
            SELECT chromosome,start,end,strand,id from genes WHERE rowid = ?
            ''',(random.randint(1,self.num_genes()),)).fetchone()
        )

    def iter_chromosomes(self):
        ''' returns chrom object iterator '''
        return ( Chrom(*x) for x in self.db.cursor().execute(''' 
            SELECT id,length FROM chromosomee
        '''))

    def iter_genes(self,gene_filter=None):
        ''' iterates over genes in refgen, only returns genes within gene filter '''
        if not gene_filter:
            gene_filter = self
        return ( Gene(*x,build=self.build,organism=self.organism) for x in self.db.cursor().execute('''
            SELECT chromosome,start,end,strand,id FROM genes
            ''') if Gene(*x) in gene_filter
        )

    def from_ids(self,gene_list,gene_filter=None):
        ''' returns gene object iterable from an iterable of id strings  '''
        if not gene_filter:
            gene_filter = self
        return ( Gene(*x,build=self.build,organism=self.organism) for x in self.db.cursor().execute(''' 
            SELECT chromosome,start,end,strand,id FROM genes WHERE id IN ('{}')
            '''.format("','".join(gene_list))) if Gene(*x) in gene_filter
        )

    @memoize
    def __getitem__(self,item):
        try:
            gene_data = self.db.cursor().execute('''
                SELECT chromosome,start,end,strand,id FROM genes WHERE id = ?
            ''',(item,)).fetchone()
            return Gene(*gene_data,build=self.build,organism=self.organism)
        except Exception as e:
            pass
        try:
            _ = (x for x in item)
            return list(self.from_ids(list(_)))
        except TypeError as e:
            self.log('not iterable: {}',e)
            pass
        return None
    
    def chromosome(self,id):
        ''' returns a chromosome object '''
        try:
            return Chrom(*self.db.cursor().execute(
                '''SELECT id,length FROM chromosomes WHERE id = ?''',(id,)).fetchone()
            ) 
        except Exception as e:
            self.log("No chromosome where id = {}. Error: {}",id,e)

    def within_gene(self,locus,gene_filter=None):
        ''' returns the gene the locus is within, or None '''
        if not gene_filter:
            gene_filter = self
        try:
            x = [Gene(*x,build=self.build,organism=self.organism) for x in self.db.cursor().execute(''' 
                SELECT chromosome,start,end,strand,id FROM genes 
                WHERE chromosome = ?
                AND start < ?
                AND end > ?
            ''',(locus.chrom,locus.start,locus.start)) if Gene(*x) in gene_filter][0]
            return x
        except Exception as e:
            return None

    def genes_within(self,locus,gene_filter=None):
        ''' Returns the genes within a locus, or None '''
        if not gene_filter:
            gene_filter = self
        try:
            return [Gene(*x,build=self.build,organism=self.organism) for x in self.db.cursor().execute(''' 
                SELECT chromosome,start,end,strand,id FROM genes
                WHERE chromosome = ?
                AND start > ?
                AND end < ?
            ''',(locus.chrom,locus.start,locus.end)) if Gene(*x) in gene_filter]
        except Exception as e:
            return None
        
    def upstream_genes(self,locus,gene_limit=1000,pos_limit=10e6,gene_filter=None):
        ''' returns genes upstream of a locus. Genes are ordered so that the
            nearest genes are at the beginning of the list.'''
        if not gene_filter:
            gene_filter = self
        return [Gene(*x,build=self.build,organism=self.organism) for x in self.db.cursor().execute(''' 
            SELECT chromosome,start,end,strand,id FROM genes
            WHERE chromosome = ?
            AND start <= ?
            AND start >= ?
            ORDER BY start DESC
            LIMIT ?
            ''',(locus.chrom, locus.start, locus.start - int(pos_limit), 100*int(gene_limit))
        ) if Gene(*x) in gene_filter][0:int(gene_limit)]
    
    def downstream_genes(self,locus,gene_limit=1000,pos_limit=10e6,gene_filter=None):
        ''' returns genes downstream of a locus. Genes are ordered so that the 
            nearest genes are at the beginning of the list. '''
        if not gene_filter:
            gene_filter = self
        return [Gene(*x,build=self.build,organism=self.organism) for x in self.db.cursor().execute(''' 
            SELECT chromosome,start,end,strand,id FROM genes
            WHERE chromosome = ?
            AND start > ?
            AND start < ?
            ORDER BY start ASC
            LIMIT ?
            ''',(locus.chrom, locus.start, locus.start + int(pos_limit), 100*int(gene_limit))
        ) if Gene(*x) in gene_filter][0:int(gene_limit)]


    def flanking_genes(self, locus, gene_limit=2, pos_limit=10e6,chain=True,gene_filter=None):
        ''' Returns genes upstream and downstream from a locus
            ** including genes locus is within **
        '''
        up_genes = self.upstream_genes(locus, gene_limit=gene_limit/2, pos_limit=pos_limit, gene_filter=gene_filter )
        down_genes = self.downstream_genes(locus, gene_limit=gene_limit/2, pos_limit=pos_limit, gene_filter=gene_filter )
        if chain:
            return list(itertools.chain.from_iterable((up_genes,down_genes)))
        else:
            return (up_genes,down_genes)

    def candidate_genes(self, locus, gene_limit=2, pos_limit=10e6,chain=True,gene_filter=None):
        ''' if locus is within gene, returns that gene, 
            otherwise returns flanking genes '''
        if not gene_filter:
            gene_filter = self
        within = self.within_gene(locus,gene_filter=gene_filter)
        if not within:
            return self.flanking_genes(
                locus,gene_limit,pos_limit,chain=chain,gene_filter=gene_filter
            )
        else:
            return [within] # candidate genes as in plural, return iterable

    def nearest_gene(self,locus,gene_filter=None):
        ''' return the gene nearest the locus '''
        if not gene_filter:
            gene_filter = self
        candidates = self.flanking_genes(locus,gene_filter=gene_filter)
        val,idx  = min((val,idx) for (idx,val) in enumerate([abs(locus-candidate) for candidate in candidates]))
        return candidates[idx]

    def bootstrap_flanking_genes(self,locus,gene_limit=4,pos_limit=50000,gene_filter=None):
        ''' Returns a randoms set of non overlapping flanking genes calculated from 
            SNPs flanking genes which is the same number of genes as would be calculated
            from the actual flanking genes from SNP list'''
        flanking_genes_index = self.flanking_genes_index(gene_limit=gene_limit,pos_limit=pos_limit)
        num_genes = len(self.flanking_genes(locus,gene_limit,pos_limit,gene_filter=gene_filter))
        if num_genes == 0:
            return []
        else:
            return random.choice(flanking_genes_index[num_genes] )

    @memoize
    def flanking_genes_index(self,pos_limit=50000,gene_limit=4):
        ''' Generate an index of flanking genes useful for bootstrapping (i.e. we can get rid of while loop '''
        # iterate over genes keeping track of number of flanking genes 
        self.log("Generating flanking gene index...")
        index = defaultdict(list)
        for gene in self.iter_genes():
            flanks = self.flanking_genes(gene,gene_limit=gene_limit,pos_limit=pos_limit)
            index[len(flanks)].append(flanks)
        for num in index:
            self.log("Found {} genes with {} flanking genes",len(index[num]),num)
        return index
 
           
    def summary(self): 
        return ("\n".join([
            'Reference Genome: {} - {} - {}',
            '{} genes',
            'Genome:',
            '{}']).format(self.organism,self.build,self.name,self.num_genes(),self.genome))

    def __repr__(self):
        return 'Reference Genome: {} - {} - {}'.format(self.organism,self.build,self.name)

        

    def __contains__(self,obj):
        ''' flexible on what you pass into the 'in' function '''
        try:
            # you can pass in a gene object (this should ALWAYS be true if you 
            # created gene object from this RefGen)
            if self.db.cursor().execute('''
                SELECT COUNT(*) FROM genes WHERE id = ?''',(obj.id,)).fetchone()[0] == 1:
                return True
            else:
                return False
        except Exception as e:
            pass
        try:
            # Can be a string object
            if self.db.cursor().execute('''
                SELECT COUNT(*) FROM genes WHERE id = ?''',(str(obj),)).fetchone()[0] == 1:
                return True
            else:
                return False
        except Exception as e:
            pass
        self.log("object {} not correct type to test membership in {}",obj,self.name)


class RefGenBuilder(Camoco):
    def __init__(self,name,description,build,organism,basedir="~/.camoco"):
        # Initiate CAmoco instance
        super().__init__(name,description=description,type="RefGen",basedir=basedir)
        self._create_tables()
        # Set the table globals
        self._global('build',build)
        self._global('organism',organism)

    def _build_indices(self):
        cur = self.db.cursor()
        cur.execute('''
            CREATE INDEX IF NOT EXISTS genepos ON genes (chromosome,start);
            CREATE INDEX IF NOT EXISTS geneid ON genes (id);
        ''')
    def add_gene(self,gene):
        self.db.cursor().execute(''' 
            INSERT OR IGNORE INTO genes VALUES (?,?,?,?,?)
        ''',(gene.id,gene.chrom,gene.start,gene.end,gene.strand))

    def add_chromosome(self,chrom):
        ''' adds a chromosome object to the class '''
        self.db.cursor().execute('''
            INSERT OR REPLACE INTO chromosomes VALUES (?,?)
        ''',(chrom.id,chrom.length))

    def import_from_txt(self,filename,sep="\t"):
        ''' assumes columns are in the order: id, chromosome, start, end, strand ''' 
        genes = []
        with open(filename,'r') as TXT:
            for line in TXT:
                fields = line.strip().split("\t")
                genes.append(fields)
        cur = self.db.cursor()
        cur.execute('BEGIN TRANSACTION')
        cur.executemany('''INSERT INTO genes (id,chromosome,start,end,strand) 
            VALUES (?,?,?,?,?)''',genes)
        cur.execute('END TRANSACTION')

    def import_from_gtf(self,filename):
        pass

    @classmethod
    def import_from_RefGen(cls,name,description,refgen,gene_filter=None,basedir="~/.camoco"):
        ''' copies from a previous instance of refgen, making sure each gene is within gene filter '''
        self = cls(name,description,refgen.build,refgen.organism,basedir)
        if not gene_filter:
            gene_filter = refgen
        for chrom in refgen.iter_chromosomes():
            self.add_chromosome(chrom)
        for gene in refgen.iter_genes(gene_filter=gene_filter):
            self.add_gene(gene)
        self._build_indices()
        return self.name
       
    def _create_tables(self):
        cur = self.db.cursor()
        cur.execute("PRAGMA page_size = 1024;")
        cur.execute("PRAGMA cache_size = 100000;")
        cur.execute(''' 
            CREATE TABLE IF NOT EXISTS chromosomes (
                id TEXT NOT NULL UNIQUE,
                length INTEGER NOT NULL
            );
            CREATE TABLE IF NOT EXISTS genes (
                id TEXT NOT NULL,
                chromosome TEXT NOT NULL,
                start INTEGER,
                end INTEGER,
                strand INTEGER
            );
        ''');

