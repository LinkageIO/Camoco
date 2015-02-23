#!/usr/bin/python3
import pyximport; pyximport.install() 
import camoco.RefGenDist as RefGenDist

from collections import defaultdict

from camoco.Camoco import Camoco
from camoco.Locus import Gene,non_overlapping
from camoco.Chrom import Chrom
from camoco.Genome import Genome
from camoco.Tools import memoize

import itertools
import random
import pandas as pd

class RefGen(Camoco):
    def __init__(self,name,description=None,basedir="~/.camoco"):
        # initialize camoco instance
        super().__init__(name,type="RefGen",description=description,basedir=basedir)
        self._create_tables()
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
            SELECT id,length FROM chromosomes
        '''))

    def iter_genes(self,gene_filter=None):
        ''' iterates over genes in refgen, only returns genes within gene filter '''
        if not gene_filter:
            gene_filter = self
        return ( Gene(*x,build=self.build,organism=self.organism) for x in self.db.cursor().execute('''
            SELECT chromosome,start,end,strand,id FROM genes
            ''') if Gene(*x,build=self.build,organism=self.organism) in gene_filter
        )

    def from_ids(self,gene_list,gene_filter=None,check_shape=False,enumerated=False):
        ''' returns gene object list from an iterable of id strings '''
        if not gene_filter:
            gene_filter = self
        genes = [ Gene(*x,build=self.build,organism=self.organism) for x in self.db.cursor().execute(''' 
            SELECT chromosome,start,end,strand,id FROM genes WHERE id IN ('{}')
            '''.format("','".join(map(str.upper,gene_list)))) if Gene(*x) in gene_filter
        ]
        if check_shape and len(genes) != len(gene_list):
            raise ValueError('Some input ids do not have genes in reference')
        return genes

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


    def flanking_genes(self, locus, gene_limit=4, window_size=100000,chain=True,gene_filter=None):
        ''' Returns genes upstream and downstream from a locus
            ** including genes locus is within **
        '''
        import math
        upstream_gene_limit = math.ceil(gene_limit/2)
        downstream_gene_limit = math.floor(gene_limit/2)
        up_genes = self.upstream_genes(locus, gene_limit=upstream_gene_limit, pos_limit=window_size/2, gene_filter=gene_filter )
        down_genes = self.downstream_genes(locus, gene_limit=downstream_gene_limit, pos_limit=window_size/2, gene_filter=gene_filter )
        if chain:
            return list(itertools.chain.from_iterable((up_genes,down_genes)))
        else:
            return (up_genes,down_genes)

    def candidate_genes(self, locus, gene_limit=2, window_size=100000,chain=True,gene_filter=None):
        ''' if locus is within gene, returns that gene, 
            otherwise returns flanking genes '''
        if not gene_filter:
            gene_filter = self
        within = self.within_gene(locus,gene_filter=gene_filter)
        if not within:
            return self.flanking_genes(
                locus,gene_limit,window_size=window_size,chain=chain,gene_filter=gene_filter
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

    def bootstrap_flanking_genes(self,locus,gene_limit=4,window_size=100000,gene_filter=None):
        ''' Returns a randoms set of non overlapping flanking genes calculated from 
            SNPs flanking genes which is the same number of genes as would be calculated
            from the actual flanking genes from SNP list'''
        flanking_genes_index = self.flanking_genes_index(gene_limit=gene_limit,window_size=window_size)
        num_genes = len(self.flanking_genes(locus,gene_limit=gene_limit,window_size=window_size,gene_filter=gene_filter))
        if num_genes == 0:
            return []
        else:
            return random.choice(flanking_genes_index[num_genes] )

    def pairwise_distance(self, gene_list=None):
        ''' returns a vector containing the pairwise distances between genes 
            in gene_list in vector form. See np.squareform for matrix conversion.
        '''
        if gene_list is None:
            gene_list = list(self.iter_genes())
        query = ''' 
                SELECT genes.id, chrom.rowid, start FROM genes 
                LEFT JOIN chromosomes chrom ON genes.chromosome = chrom.id      
                WHERE genes.id in ("{}")
                ORDER BY genes.id
        '''.format('","'.join([g.id for g in gene_list]))
        # extract chromosome row ids and gene start positions for each gene
        positions = pd.DataFrame(
            # Grab the chromosomes rowid because its numeric
            self.db.cursor().execute(query).fetchall(),
            columns=['gene','chrom','pos']
        ).sort('gene')
        assert len(positions) == len(gene_list), 'Some genes in dataset not if RefGen'
        assert all(positions.gene == [g.id for g in gene_list]), 'Genes are not in the correct order!'
        distances = RefGenDist.gene_distances(positions.chrom.values,positions.pos.values) 
        return distances


    @memoize
    def flanking_genes_index(self,window_size=100000,gene_limit=4):
        ''' Generate an index of flanking genes useful for bootstrapping (i.e. we can get rid of while loop '''
        # iterate over genes keeping track of number of flanking genes 
        self.log("Generating flanking gene index window size: {}, gene_limit {}",window_size,gene_limit)
        index = defaultdict(list)
        for gene in self.iter_genes():
            flanks = self.flanking_genes(gene,gene_limit=gene_limit,window_size=window_size)
            index[len(flanks)].append(flanks)
        for num in index:
            self.log("Found {} genes with {} flanking genes",len(index[num]),num)
        return index
 
           
    def summary(self): 
        print ("\n".join([
            'Reference Genome: {} - {} - {}',
            '{} genes',
            'Genome:',
            '{}']).format(self.organism,self.build,self.name,self.num_genes(),self.genome))

    def __repr__(self):
        return 'Reference Genome: {} - {} - {}'.format(self.organism,self.build,self.name)


    def __len__(self):
        return self.num_genes()        

    def __contains__(self,obj):
        ''' flexible on what you pass into the 'in' function '''
        try:
            # you can pass in a gene object (this expression should ALWAYS be true if you 
            # created gene object from this RefGen)
            if self.db.cursor().execute(
                '''SELECT COUNT(*) FROM genes WHERE id = ?''',(obj.id.upper(),)).fetchone()[0] == 1:
                return True
            else:
                return False
        except Exception as e:
            pass
        try:
            # Can be a string object
            if self.db.cursor().execute('''
                SELECT COUNT(*) FROM genes WHERE id = ?''',(str(obj).upper(),)).fetchone()[0] == 1:
                return True
            else:
                return False
        except Exception as e:
            pass
        self.log("object {} not correct type to test membership in {}",obj,self.name)



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

    @classmethod
    def from_gff(cls,filename,name,description,build,organism,basedir='~/.camoco'):
        ''' imports refgen from gff file '''
        self = cls(name,description,basedir=basedir)
        self._global('build',build)
        self._global('organism',organism)
        with open(filename,'r') as IN:
            for line in IN:
                #skip comment lines
                if line.startswith('#'):
                    continue
                chrom,source,feature,start,end,score,strand,frame,attributes = line.strip().split()
                attributes = dict([(field.split('=')) for field in attributes.strip(';').split(';')])
                if feature == 'chromosome':
                    self.add_chromosome(Chrom(attributes['ID'],end))
                if feature == 'gene':
                    self.add_gene(Gene(chrom,start,end,strand,attributes['ID'],build,organism))
        return self

    @classmethod
    def from_RefGen(cls,name,description,refgen,gene_filter=None,basedir="~/.camoco"):
        ''' copies from a previous instance of refgen, making sure each gene is within gene filter '''
        self = cls(name,description,basedir)
        self._global('build',refgen.build)
        self._global('organism',refgen.organism)
        if not gene_filter:
            gene_filter = refgen
        for chrom in refgen.iter_chromosomes():
            self.add_chromosome(chrom)
        for gene in refgen.iter_genes():
            if gene in gene_filter:
                self.add_gene(gene)
        self._build_indices()
        return self
       
    def _create_tables(self):
        cur = self.db.cursor()
        cur.execute(''' 
            CREATE TABLE IF NOT EXISTS chromosomes (
                id TEXT NOT NULL UNIQUE,
                length INTEGER NOT NULL
            );
            CREATE TABLE IF NOT EXISTS genes (
                id TEXT NOT NULL UNIQUE,
                chromosome TEXT NOT NULL,
                start INTEGER,
                end INTEGER,
                strand INTEGER
            );
        ''');

