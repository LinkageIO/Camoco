#!/usr/bin/python4
import pyximport; pyximport.install() 
import camoco.RefGenDist as RefGenDist

from collections import defaultdict
import matplotlib.pylab as plt

from camoco.Camoco import Camoco
from camoco.Locus import Gene,Locus
from camoco.Chrom import Chrom
from camoco.Genome import Genome
from camoco.Tools import memoize

import itertools
import collections
import random
import pandas as pd
import numpy as np
import math

class RefGen(Camoco):
    def __init__(self,name):
        # initialize camoco instance
        super().__init__(name,type="RefGen")
        self._create_tables()

    @property
    def genome(self):
        return Genome(
            self.type + self.name,
            chroms = [Chrom(*x) for x in  self.db.cursor().execute(''' 
                SELECT id,length FROM chromosomes
            ''')] 
        )

    def Gene(self,chrom,start,end,name,window=0,sub_loci=None,**kwargs):
        attrs = dict(self.db.cursor().execute(''' 
            SELECT key,val FROM gene_attrs WHERE id = ?
        ''',(name,)).fetchall())
        return Gene(chrom,start,end,name,window,
                sub_loci,**kwargs).update(
                    attrs        
                )

    @memoize
    def num_genes(self):
        ''' 
            Returns the number of genes in the dataset 
        '''
        return self.db.cursor().execute(
           ''' SELECT COUNT(*) FROM genes'''
        ).fetchone()[0]

    def random_gene(self,**kwargs):
        '''
            Returns a random gene within the reference genome.
            Also allows passing of keyword arguments to Locus 
            constructor method allowing for flexible generation.
            See Locus.__init__ for more details.

            Parameters
            ----------
            **kwargs : key,value pairs
                Extra parameters passed onto the locus init method.

            Returns
            -------
            A locus object (Gene)

        '''
        return self.Gene(*self.db.cursor().execute(''' 
            SELECT chromosome,start,end,id from genes WHERE rowid = ?
            ''',(random.randint(1,self.num_genes()),)).fetchone(),
            **kwargs
        )

    def iter_chromosomes(self):
        ''' returns chrom object iterator '''
        return ( Chrom(*x) for x in self.db.cursor().execute(''' 
            SELECT id,length FROM chromosomes
        '''))

    def iter_genes(self):
        ''' iterates over genes in refgen, 
            only returns genes within gene filter '''
        return ( 
            self.Gene(*x,build=self.build,organism=self.organism) \
            for x in self.db.cursor().execute('''
                SELECT chromosome,start,end,id FROM genes
            ''')
        )

    def from_ids(self, gene_list, check_shape=False):
        ''' 
            Returns a list of gene object from an iterable of id strings 
            OR from a single gene id string.
            
            Parameters
            ----------
            gene_list : str OR iterable of str
                ID(s) of the genes you want to pull out
            check_shape : bool (default: False)
                Check if you get back the same number of ids you
                pass in. If false (default), just give back what
                you find, ignoring erronous ids.

            Returns
            -------
            A list of locus objects if you pass in an iterable,
            otherwise a single gene
                
        '''
        if isinstance(gene_list,str):
            # Handle when we pass in a single id
            gene_id = gene_list.upper()
            try:
                return self.Gene(
                    *self.db.cursor().execute('''
                        SELECT chromosome,start,end,id FROM genes WHERE id = ? 
                        ''',(gene_id,)
                    ).fetchone(),
                    build=self.build,
                    organism=self.organism
                ) 
            except TypeError as e:
                raise ValueError('{} not in {}'.format(gene_id,self.name))
        genes = [ 
            self.Gene(*x,build=self.build,organism=self.organism) \
            for x in self.db.cursor().execute(''' 
                SELECT chromosome,start,end,id FROM genes WHERE id IN ('{}')
            '''.format("','".join(map(str.upper,gene_list))))
        ]
        if check_shape and len(genes) != len(gene_list):
            raise ValueError('Some input ids do not have genes in reference')
        return genes

    def __getitem__(self,item):
        '''
            A convenience method to extract loci from the reference geneome. 
        '''
        return self.from_ids(item)
    
    def chromosome(self,id):
        ''' 
            returns a chromosome object 
        '''
        try:
            return Chrom(*self.db.cursor().execute(
                '''SELECT id,length FROM chromosomes WHERE id = ?''',
                (id,)).fetchone()
            ) 
        except Exception as e:
            self.log("No chromosome where id = {}. Error: {}",id,e)

    def genes_within(self,loci,chain=True):
        ''' 
            Returns the genes within a locus, or None 
        '''
        if isinstance(loci,Locus):
            return [
                self.Gene(*x,build=self.build,organism=self.organism) \
                for x in self.db.cursor().execute(''' 
                    SELECT chromosome,start,end,id FROM genes
                    WHERE chromosome = ?
                    AND start > ?
                    AND end < ?
                ''',(loci.chrom,loci.start,loci.end))]
        else:
            iterator = iter(loci)   
            genes = [self.genes_within(locus,chain=chain) for locus in iterator]
            if chain:
                genes = list(itertools.chain(*genes))
            return genes

    def upstream_genes(self,locus,gene_limit=1000):
        ''' 
            returns genes upstream of a locus. Genes are ordered so that the
            nearest genes are at the beginning of the list.
        '''
        return [
            self.Gene(*x,build=self.build,organism=self.organism) \
            for x in self.db.cursor().execute(''' 
                SELECT chromosome,start,end,id FROM genes
                WHERE chromosome = ?
                AND start <= ?
                AND start >= ?
                ORDER BY start DESC
                LIMIT ?
            ''',(locus.chrom, locus.start, locus.upstream, gene_limit)
        )]
    
    def downstream_genes(self,locus,gene_limit=1000):
        '''
            returns genes downstream of a locus. Genes are ordered so that the 
            nearest genes are at the beginning of the list. 
        '''
        return [
            self.Gene(*x,build=self.build,organism=self.organism) \
            for x in self.db.cursor().execute(''' 
                SELECT chromosome,start,end,id FROM genes
                WHERE chromosome = ?
                AND start > ?
                AND start < ?
                ORDER BY start ASC
                LIMIT ?
            ''',(locus.chrom, locus.end, locus.downstream, gene_limit)
        )]

    def flanking_genes(self, loci, gene_limit=4,chain=True):
        ''' 
            Returns genes upstream and downstream from a locus
            ** including genes locus is within **
        '''
        if isinstance(loci,Locus):
            # If we cant iterate, we have a single locus
            locus = loci
            upstream_gene_limit = math.ceil(gene_limit/2)
            downstream_gene_limit = math.floor(gene_limit/2)
            up_genes = self.upstream_genes(
                locus, gene_limit=upstream_gene_limit
            )
            down_genes = self.downstream_genes(
                locus, gene_limit=downstream_gene_limit
            )
            if chain:
                return list(itertools.chain(up_genes,down_genes))
            return (up_genes,down_genes)
        else:
            iterator = iter(loci)   
            genes = [self.flanking_genes(locus,gene_limit=gene_limit) \
                for locus in iterator
            ]
            if chain:
                genes = list(itertools.chain(*genes))
            return genes

    def bootstrap_candidate_genes(self,loci,gene_limit=4,chain=True):
        '''
            Returns candidate genes which are random, but conserves 
            total number of overall genes.

            Parameters
            ----------
            loci : camoco.Locus (also handles an iterable containing Loci)
                a camoco locus or iterable of loci 
            gene_limit : int (default : 4)
                The total number of flanking genes 
                considered a candidate surrounding a locus
            chain : bool (default : true)
                Calls itertools chain on results before returning  

            Returns
            -------
            a list of candidate genes (or list of lists if chain is False)

        '''
        if isinstance(loci,Locus):
            # We now have a single locus
            locus = loci
            # grab the actual candidate genes
            num_candidates = len(self.candidate_genes(
                locus,gene_limit=gene_limit,chain=True)
            )
            if num_candidates == 0:
                return []
            # Snps a random genes from the genome
            random_gene = self.random_gene()
            # Extend the window to something crazy
            random_gene.window = 10e10
            # Snag the same number of candidates
            random_candidates = self.upstream_genes(
                random_gene,gene_limit=num_candidates
            )
            if len(random_candidates) != num_candidates:
                # somehow we hit the end of a chromosome 
                # or something, just recurse
                return self.bootstrap_candidate_genes(
                    locus,gene_limit=gene_limit,chain=chain)
            assert len(random_candidates) == num_candidates
            return random_candidates
        else:
            # Sort the loci so we can collapse down 
            locus_list = sorted(loci) 
            seen = set()
            bootstraps = list()
            for locus in locus_list:
                # compare downstream of last locus to current locus
                target_len = len(self.candidate_genes(
                    locus,gene_limit=gene_limit)
                )
                genes = self.bootstrap_candidate_genes(
                    locus, gene_limit=gene_limit, chain=chain
                )
                # If genes randomly overlap, resample
                while np.any([x in seen for x in itertools.chain(genes,)]):
                    genes = self.bootstrap_candidate_genes(
                        locus, gene_limit=gene_limit, chain=chain
                    )
                # Add all new bootstrapped genes to the seen list 
                [seen.add(x) for x in itertools.chain(genes,)]
                assert target_len == len(genes)
                bootstraps.append(genes)
            if chain:
                bootstraps = list(set(itertools.chain(*bootstraps)))
            self.log("Found {} bootstraps",len(bootstraps))
            return bootstraps

    def candidate_genes(self, loci, gene_limit=4,chain=True):
        ''' 
            SNP to Gene mapping.
            Return Genes between locus start and stop, plus additional 
            flanking genes (up to gene_limit)

            Parameters
            ----------
            loci : camoco.Locus (also handles an iterable containing Loci)
                a camoco locus or iterable of loci 
            gene_limit : int (default : 4)
                The total number of flanking genes 
                considered a candidate surrounding a locus
            chain : bool (default : true)
                Calls itertools chain on results before returning  

            Returns
            -------
            a list of candidate genes (or list of lists if chain is False)

        '''
        if isinstance(loci,Locus):
            # If not an iterator, its a single locus
            locus = loci
            genes_within = self.genes_within(locus)
            up_genes,down_genes = self.flanking_genes(
                locus,gene_limit=gene_limit,chain=False
            )
            if chain:
                return list(itertools.chain(up_genes,genes_within,down_genes))
            return (up_genes,genes_within,down_genes)
        else:
            iterator = iter(sorted(loci))
            genes = [
                self.candidate_genes(
                    locus,gene_limit=gene_limit,chain=chain
                ) for locus in iterator
            ]
            if chain:
                genes = list(set(itertools.chain(*genes)))
            self.log("Found {} candidates",len(genes))
            return genes


    def pairwise_distance(self, gene_list=None):
        ''' 
            returns a vector containing the pairwise distances between genes 
            in gene_list in vector form. See np.squareform for matrix 
            conversion.
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
        # chromosome needs to be floats
        positions.chrom = positions.chrom.astype('float')
        assert len(positions) == len(gene_list), \
            'Some genes in dataset not if RefGen'
        assert all(positions.gene == [g.id for g in gene_list]), \
            'Genes are not in the correct order!'
        distances = RefGenDist.gene_distances(
            positions.chrom.values,positions.pos.values
        ) 
        return distances

    def summary(self): 
        print ("\n".join([
            'Reference Genome: {} - {} - {}',
            '{} genes',
            'Genome:',
            '{}']).format(
                self.organism,self.build,
                self.name,self.num_genes(),
                self.genome
            )
        )

    def plot_loci(self,loci,filename,gene_limit=4):
        ''' 
            Plots the loci, windows and candidate genes
            
            Parameters
            ----------
            loci : iterable of co.Loci
            filename : output filename
        '''
        plt.clf()
        # Each chromosome gets a plot
        chroms = set([x.chrom for x in loci])
        f, axes = plt.subplots(len(chroms),figsize=(10,4*len(chroms)))
        # Loci Locations
        chromloci = defaultdict(list)
        for locus in sorted(loci):
            chromloci[locus.chrom].append(locus)

        # iterate over Loci
        seen_chroms = set([loci[0].chrom])
        voffset = 1 # Vertical Offset
        hoffset = 0 # Horizonatal Offset
        current_chrom = 0
        for i,locus in enumerate(loci):
            # Reset the temp variables in necessary
            if locus.chrom not in seen_chroms:
                seen_chroms.add(locus.chrom)
                current_chrom += 1
                voffset = 1
                hoffset = 0  
            # Do the access things
            cax = axes[current_chrom]
            cax.set_ylabel('Chrom: '+ locus.chrom)
            cax.set_xlabel('Loci')
            cax.get_yaxis().set_ticks([])
            #cax.get_xaxis().set_ticks([])
            # shortcut for current axis
            cax.hold(True)
            # place marker for start window
            cax.scatter(hoffset,voffset,marker='>') 
            # place marker for start snp
            cax.scatter(hoffset+locus.window,voffset,marker='.',color='blue')
            # place marker for stop snp
            cax.scatter(hoffset+locus.window+len(locus),voffset,marker='.',color='blue')
            # place marker for stop snp
            cax.scatter(hoffset+locus.window+len(locus)+locus.window,voffset,marker='<')

            # place markers for sub snps
            for subsnp in locus.sub_loci:
                cax.scatter(hoffset+subsnp.start-locus.start+locus.window,voffset,marker='.',color='blue')

            # place a block for interlocal distance
            cax.barh(
                bottom=voffset,
                width=50,
                height=1,
                left=hoffset+locus.window+len(locus)+locus.window,
                color='red'
            )
            # grab the candidate genes
            for gene in self.candidate_genes(locus,gene_limit=gene_limit):
                cax.barh(
                    bottom=voffset,
                    width = len(gene),
                    height= 1,
                    left=gene.start-locus.start+locus.window,
                    color='red'
                )
            voffset += 5

        plt.savefig(filename)
        del f

    def __repr__(self):
        return 'Reference Genome: {} - {} - {}'.format(
            self.organism,self.build,self.name
        )

    def __len__(self):
        return self.num_genes()        

    def __contains__(self,obj):
        ''' flexible on what you pass into the 'in' function '''
        if isinstance(obj,Locus):
            # you can pass in a gene object (this expression 
            # should ALWAYS be true if you 
            # created gene object from this RefGen)
            if self.db.cursor().execute(
                '''SELECT COUNT(*) FROM genes WHERE id = ?''',
                (obj.id.upper(),)).fetchone()[0] == 1:
                return True
            else:
                return False
        elif isinstance(obj,str):
            # Can be a string object
            if self.db.cursor().execute('''
                SELECT COUNT(*) FROM genes WHERE id = ?''',
                (str(obj).upper(),)).fetchone()[0] == 1:
                return True
            else:
                return False
        else:
            raise TypeError('Cannot test for containment for {}'.format(obj))

    def _build_indices(self):
        cur = self.db.cursor()
        cur.execute('''
            CREATE INDEX IF NOT EXISTS genepos ON genes (chromosome,start);
            CREATE INDEX IF NOT EXISTS geneid ON genes (id);
            CREATE INDEX IF NOT EXISTS geneattr ON gene_attrs (id);
        ''')

    def add_gene(self,gene):
        if isinstance(gene,Locus):
            self.db.cursor().execute(''' 
            INSERT OR IGNORE INTO genes VALUES (?,?,?,?)
            ''',(gene.name,gene.chrom,gene.start,gene.end))
            self.db.cursor().executemany('''
            INSERT OR IGNORE INTO gene_attrs VALUES (?,?,?)
            ''',[(gene.id,key,val) for key,val in gene.attr.items()])
        else:
            # support adding lists of genes
            genes = list(gene)
            self.log('Adding Gene base info to database')
            self.db.cursor().executemany('''
            INSERT OR IGNORE INTO genes VALUES (?,?,?,?)
            ''',((gene.name,gene.chrom,gene.start,gene.end) for gene in genes))
            self.log('Adding Gene attr info to database')
            self.db.cursor().executemany('''
            INSERT OR IGNORE INTO gene_attrs VALUES (?,?,?)
            ''',[(gene.id,key,val) for gene in genes for key,val in gene.attr.items()])

    def add_chromosome(self,chrom):
        ''' adds a chromosome object to the class '''
        self.db.cursor().execute('''
            INSERT OR REPLACE INTO chromosomes VALUES (?,?)
        ''',(chrom.id,chrom.length))

    @classmethod
    def from_gff(cls,filename,name,description,build,organism):
        ''' 
            Imports RefGen object from a gff (General Feature Format) file.
            See more about the format here: 
            http://www.ensembl.org/info/website/upload/gff.html

            Parameters
            ----------

            filename : str
                The path to the GFF file.
            name : str
                The name if the RefGen object to be stored in the core
                camoco database.
            description : str
                A short description of the RefGen for future reference
            build : str
                A string designating the genome build, used for comparison
                operations, genes may share IDS but are different across build.
            organism : str
                A short string describing the organims this RefGen is coming 
                from. Again, is used in comparing equality among genes which 
                may have the same id or name.

        '''
        self = cls.create(name,description,type='RefGen')
        self._global('build',build)
        self._global('organism',organism)
        genes = list()
        with open(filename,'r') as IN:
            for line in IN:
                #skip comment lines
                if line.startswith('#'):
                    continue
                (chrom,source,feature,start,
                 end,score,strand,frame,attributes) = line.strip().split()
                attributes = dict([(field.split('=')) \
                    for field in attributes.strip(';').split(';')
                ])
                if feature == 'chromosome':
                    self.log('Found a chromosome: {}',attributes['ID'])
                    self.add_chromosome(Chrom(attributes['ID'],end))
                if feature == 'gene':
                    genes.append(
                        Gene(
                            chrom,int(start),int(end),
                            attributes['ID'].upper(),strand=strand,
                            build=build,organism=organism
                        ).update(attributes)
                    )
        self.add_gene(genes)
        return self

    @classmethod
    def filtered_refgen(cls,name,description,refgen,gene_list):
        ''' 
            Copies from a previous instance of refgen, making sure 
            each gene is within gene list 
        '''
        self = cls.create(name,description,'RefGen')
        self._global('build',refgen.build)
        self._global('organism',refgen.organism)
        # Should have the same chromosomes
        for chrom in refgen.iter_chromosomes():
            self.add_chromosome(chrom)
        # Add the genes from genelist
        self.add_gene(gene_list)
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
                end INTEGER
            );
            CREATE TABLE IF NOT EXISTS gene_attrs (
                id TEXT NOT NULL,
                key TEXT,
                val TEXT
            )
        ''');

    ''' ----------------------------------------------------------------------
            Unimplemented
    '''

    def within_gene(self,locus):
        ''' 
            Returns the gene the locus is within, or None 
        '''
        try:
            x = [self.Gene(*x,build=self.build,organism=self.organism) \
                for x in self.db.cursor().execute(''' 
                    SELECT chromosome,start,end,id FROM genes 
                    WHERE chromosome = ?
                    AND start < ?
                    AND end > ?
                ''',(locus.chrom,locus.start,locus.start))][0]
            return x
        except Exception as e:
            return None

    def nearest_gene(self,locus):
        ''' return the gene nearest the locus '''
        candidates = self.flanking_genes(locus)
        val,idx  = min((val,idx) for (idx,val) \
            in enumerate([abs(locus-candidate) for candidate in candidates]))
        return candidates[idx]

    @memoize
    def flanking_genes_index(self,gene_limit=4,window=50000):
        ''' 
            Generate an index of flanking genes useful for 
            bootstrapping (i.e. we can get rid of while loop 
        '''
        # iterate over genes keeping track of number of flanking genes 
        self.log("Generating flanking gene index gene_limit {}",gene_limit)
        index = defaultdict(list)
        for gene in self.iter_genes():
            gene.window = window
            flanks = self.flanking_genes(gene,gene_limit=gene_limit)
            index[len(flanks)].append(flanks)
        for num in index:
            self.log(
                "Found {} genes with {} flanking genes",
                len(index[num]),num
            )
        return index
