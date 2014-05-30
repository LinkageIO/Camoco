#!/usr/bin/python3

from Camoco import Camoco
from Locus import Gene
from Chrom import Chrom
from Genome import Genome

class RefGen(Camoco):
    def __init__(self,name,basedir="~/.camoco"):
        super().__init__(name,type="RefGen",basedir=basedir)

    def from_ids(self,id_list):
        ''' returns gene objects from an id  '''
        return [ Gene(*x,build=self.build,organism=self.organism) for x in self.db.cursor().execute(''' 
            SELECT chromosome,start,end,id,strand FROM genes WHERE id IN ('?')
            '''.format("'.'".join(id_list)))
        ]

    def genome(self):
        ''' returns a genome object based on loaded chromosome object'''
        return Genome(
            self.type+self.name,
            chroms = [Chrom(*x) for x in  self.db.cursor().execute(''' 
                SELECT id,length FROM chromosomes
            ''')] 
        )
    def chromosome(self,id):
        ''' returns a chromosome object '''
        try:
            return Chrom(*self.db.cursor().execute(
                '''SELECT id,length FROM chromosomes WHERE id = ?''',(id,)).fetchone()
            ) 
        except Exception as e:
            self.log("No chromosome where id = {}. Error: {}",id,e)
        

    def upstream_genes(self,locus,gene_limit=1000,pos_limit=10e6):
        ''' returns genes upstream of a locus '''
        return [Gene(*x,build=self.build,organism=self.organism) for x in self.db.cursor().execute(''' 
            SELECT chromosome,start,end,strand,id FROM genes
            WHERE chromosome = ?
            AND start < ?
            AND start > ?
            ORDER BY start ASC
            LIMIT ?
        ''',(locus.chrom, locus.start, locus.start - int(pos_limit), int(gene_limit)))]
    
    def downstream_genes(self,locus,gene_limit=1000,pos_limit=10e6):
        ''' returns genes downstream of a locus '''
        return [Gene(*x,build=self.build,organism=self.organism) for x in self.db.cursor().execute(''' 
            SELECT chromosome,start,end,strand,id FROM genes
            WHERE chromosome = ?
            AND start > ?
            AND start < ?
            ORDER BY start ASC
            LIMIT ?
        ''',(locus.chrom, locus.start, locus.start + int(pos_limit), int(gene_limit)))]

    def flanking_genes(self, locus, gene_limit=2, pos_limit=10e6):
        return (self.upstream_genes(locus,gene_limit=gene_limit/2,pos_limit=pos_limit),
                self.downstream_genes(locus,gene_limit=gene_limit/2,pos_limit=pos_limit))

    def nearest_gene(self,locus):
        ''' return the gene nearest the locus '''
        up,down = self.flanking_genes(locus,gene_limit=2,pos_limit=100e6)
        if locus in up[0]:
            return up[0]
        elif locus in down[0]:
            return down[0]
        elif abs(up[0] - locus) < abs(down[0] - locus):
            return up[0]
        else:
            return down[0] 

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
    def add_gene(self,id,chromosome,start,end,strand):
        self.db.cursor().execute(''' 
            INSERT OR IGNORE INTO genes VALUES (?,?,?,?,?)
        ''',(id,chromosome,start,end,strand))

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

    def import_from_fasta(self,filename):
        pass
       
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
