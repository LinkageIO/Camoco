#!/usr/bin/python3
import apsw as lite
import os as os
import time as time
import sys
import numpy as np
from scipy.stats import pearsonr

from camoco.Camoco import Camoco
from camoco.Chrom import Chrom

from matplotlib import pylab

class HapMap(Camoco):
    def __init__(self,name,basedir="~/.camoco"):
        # initialize super class
        super(self.__class__,self).__init__(name,type='HapMap',basedir=basedir)

    @property
    def num_accessions(self):
        return self.db.cursor().execute("SELECT COUNT(*) FROM accessions").fetchone()[0]

    @property
    def num_snps(self):
        return self.db.cursor().execute("SELECT COUNT(*) FROM snps").fetchone()[0]

    def accessions(self):
        ''' returns genotypes for HapMap object '''
        return [ x[0] for x in self.db.cursor().execute('SELECT * FROM accessions ORDER BY rowid').fetchall()]

    def snps(self,chromo,start,end):
        ''' Returns SNPs within range '''
        return [SNP(*x) for x in self.db.cursor().execute(
            '''SELECT chromosome,position,id FROM snps WHERE chromosome = ? AND position >= ? AND position <= ?''',
            (chromo,start,end)).fetchall()]

    def genotypes(self,snp,accessions=None):
        ''' returns the genotypes as a list based on position or snpid,
            returns None if snp not in database '''
        genos = [list(map(float,x[0].split("\t"))) for x in self.db.cursor().execute(
        ''' SELECT alleles FROM snps JOIN genotypes ON genotypes.snpRowID = snps.rowid
            WHERE snps.chromosome = ? AND snps.position = ?''', (snp.chrom,snp.pos)
        )]
        if len(genos) == 0:
            return None
        if accessions is not None:
            try:
                return np.array(genos[0])[list(map(self.accessions().index,accessions))]
            except ValueError as e:
                self.log(e)
        else:
            return np.array(genos[0])

    def genotype_missing_mask(self,snp):
        ''' returns a boolean mask on missing genotypes '''
        return self.genotypes(snp) == 9

    def ld(self,snp1,snp2,nearest=False):
        ''' returns the pearson r2 between two snps, if snps are not in the database,
            will use the closest snp if nearest=True (default False)'''
        if nearest:
            snp1 = self.nearest_snp(snp1)
            snp2 = self.nearest_snp(snp2)
        geno1 = self.genotypes(snp1)
        geno2 = self.genotypes(snp2)
        if geno1 is None or geno1 is None:
            return (np.nan,np.nan)
        non_missing = (geno1!=9)&(geno2!=9)
        return (pearsonr(geno1[non_missing],geno2[non_missing])[0]**2,sum(non_missing),snp2-snp1)

    def downstream_snps(self,snp,snp_limit=1000,pos_limit=1e8):
        ''' returns snps downstream of snp '''
        return [SNP(*x) for x in self.db.cursor().execute(''' 
            SELECT chromosome,position,id FROM snps 
            WHERE chromosome = ?
            AND position > ?
            AND position < ?
            ORDER BY position ASC
            LIMIT ?
            ''',(snp.chrom,snp.start,snp.start+pos_limit,int(snp_limit)))]
    def upstream_snps(self,snp,snp_limit=1000,pos_limit=1e8):
        ''' returns snps upstream of input '''
        return [SNP(*x) for x in self.db.cursor().execute(''' 
            SELECT chromosome,position,id FROM snps 
            WHERE chromosome = ?
            AND position < ?
            AND position > ?
            ORDER BY position ASC
            LIMIT ?
            ''',(snp.chrom,snp.start,snp.start-pos_limit,int(snp_limit)))]

    def flanking_snps(self,snp,snp_limit=100,pos_limit=1e8):
        ''' returns snps flanking '''
        return (self.upstream_snps(snp,snp_limit=snp_limit,pos_limit=pos_limit),
                self.downstream_snps(snp,snp_limit=snp_limit,pos_limit=pos_limit))

    def nearest_snp(self,snp,pos_limit=1e8):
        ''' retruns the nearest upstream or downstream SNP '''
        if snp in self:
            return snp
        nearest = (self.upstream_snps(snp,snp_limit=1,pos_limit=pos_limit) +
                   self.downstream_snps(snp,snp_limit=1,pos_limit=pos_limit))
        if len(nearest) == 0:
            self.log("No nearest SNPs within {}",pos_limit)
        elif len(nearest) == 1:
            return nearest[0]
        elif abs(snp - nearest[0]) < abs(snp - nearest[1]):
            return nearest[0]
        else:
            return nearest[1]
        return nearest
        
    def ld_interval(self,snp,edge_cutoff = 0.5,snp_limit=1000,pos_limit=10000):
        ''' Returns the interval around a locus/snp based on 
            a lowess regression of LD ~ distance-from-locus.
            This implementation calls the edge when the best fit
            line drops below the edge_cutoff parameter '''
        # grab upstream snp
        upstream,downstream = self.flanking_snps(snp,snp_limit=snp_limit,pos_limit=pos_limit) 
        # mask out snps with < 80 points with filter
        ld = [self.ld(snp,x) for x in upstream+downstream]
        # filter out nans
        n = np.array([x[1] for x in ld])
        mask = n>(max(n)/2)
        distances = np.array([snp.pos+x[2] for x in ld])[mask]
        r2 = np.array([x[0] for x in ld])[mask]
        flanking = np.array(upstream+downstream)[mask]
        r2[np.isnan(r2)] = 0
        # fit a lowess 
        a=(distances,n,r2,flanking)#,intercept,slope)
        pylab.clf()
        fig = pylab.figure()
        ax = fig.add_subplot(1,1,1)
        ax.scatter(distances,r2,alpha=0.5,label='pairwise ld')
        ax.ticklabel_format(axis='x',style='plain')
        try:
            yest = lowess(distances,r2,iter=10)
            ax.plot(distances,yest,label='lowess')
        except Exception as e:
            pass
        ax.legend()
        ax.set_ylim([0,1])
        ax.axvline(snp.pos,color='k')
        #ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=45)
        fig.tight_layout()
        fig.savefig('{}-{}:{}.png'.format(snp.id,snp.chrom,snp.pos))
        return a

    def sliding_window(self,interval=1000000):
        ''' Calculates average pairwise LD between SNPS within sliding window '''
        pass 


    def __contains__(self,obj):
        if isinstance(obj,SNP):
            (count,) = self.db.cursor().execute('''
                SELECT COUNT(*) FROM snps WHERE chromosome = ? AND position = ? ''',
                (obj.chrom,obj.pos)
            ).fetchone()
            if int(count) == 1:
                return True
            else:
                return False
        else:
            raise Exception("{} must be a SNP object".format(obj))

    def __str__(self):
        return '''
            HapMap: {}
            desc: {}
            #SNPs: {}
            #Accessions: {}
        '''.format(self.name,self.description,self.num_snps,self.num_accesions)

    def __repr__(self):
        return str(self)

class HapMapBuilder(Camoco):
    def __init__(self,name,description="",basedir="~/.camoco"):
        # add entry to camoco database
        # add yourself to the database
        super().__init__(name,description=description,type='HapMap',basedir=basedir)
        self._create_tables()
        self.snp_rows = []
        self.genotype_rows = []
        self.snpRowID = 1
        record_dump_size = 100000

    @property
    def num_accessions(self):
        return self.db.cursor().execute("SELECT COUNT(*) FROM accessions").fetchone()[0]
   
    def add_accession(self,accession_name):
        cur = self.db.cursor()
        try:    
            cur.execute('BEGIN TRANSACTION')
            cur.execute("INSERT INTO accessions VALUES(?)",(accession_name,)) 
            cur.execute('END TRANSACTION')
            return True
        except Exception as e:
            self.log('Something Bad Happened: {}',e)
            cur.execute("ROLLBACK")
            return False
             
    def add_snp(self,chrom,pos,snpid,alt,ref,qual,tabbed_genotypes):
        cur = self.db.cursor()
        if len(tabbed_genotypes.split("\t")) != self.num_accessions:
            self.log("The number of genotypes in does not match number of accessions")
            self.log("{} vs {}",len(tabbed_genotypes.split("\t")),self.num_accessions)
        try:    
            cur.execute('BEGIN TRANSACTION')
            cur.execute('''INSERT INTO snps
                VALUES (?,?,?,?,?,?)''',(int(chrom),int(pos),snpid,alt,ref,float(qual))) 
            snp_row_id = self.db.last_insert_rowid()
            cur.execute('''INSERT INTO genotypes VALUES (?,?)''',(snp_row_id,tabbed_genotypes))
            cur.execute('END TRANSACTION')
            return True
        except Exception as e:
            self.log('Something Bad Happened: {}',e)
            cur.execute("ROLLBACK")
            return False

    def import_vcf(self,filename):
        ''' imports hapmap information from a VCF '''
        self.log("Importing VCF: {}",filename)
        with open(filename,'r') as INVCF:
            genotype_rows = []
            snp_rows = []
            snpRowID = 1 # 
            record_dump_size = 100000
            cur = self.db.cursor()
            cur.execute('PRAGMA synchronous = off')
            cur.execute('PRAGMA journal_mode = memory')
            self._drop_indices()
            for i,line in enumerate(INVCF):
                line = line.strip()
                if line.startswith("##"):
                    next
                elif line.startswith("#"):
                    # we are at the header section, extract accessions here
                    (chrom,pos,id,ref,alt,qual,fltr,info,fmt,*accessions) = line.split()
                    for accession in accessions:
                        cur.execute('''
                            INSERT INTO accessions VALUES('{}');
                        '''.format(accession)
                        )
                else:
                    # we are at records, insert snp line by line but bulk insert genotypes
                    (chrom,pos,snpid,ref,alt,qual,fltr,info,fmt,*genotypes) = line.split()
                    if qual == '.':
                        qual = -1
                    snp_rows.append((chrom,int(pos),snpid,ref,alt,float(qual)))
                    GT_ind = fmt.split(":").index("GT") # we need to find the GT field
                    def zygosity(alleles):
                        if '.' in alleles:
                            return '9'
                        elif alleles == '00':
                            return '0'
                        elif alleles == '01' or alleles == '10':
                            return '1'
                        else:
                            return '2'
                    alleles = "\t".join(
                        [zygosity(x.split(':')[GT_ind].replace("|","/").replace("/",""))
                             for x in genotypes
                        ]
                    )
                    genotype_rows.append((snpRowID,alleles))
                    snpRowID += 1 # <- Do NOT fuck this up, genotypes need to reference correct SNPId
                if i % record_dump_size == 0 and i > 0:
                    self.log(" {} records reached. dumping data",i)
                    cur.execute("BEGIN TRANSACTION")
                    start_time = time.time()
                    # Bulk import the genotype values
                    self.log("\tInserting snps...")
                    cur.executemany(''' 
                        INSERT INTO snps (chromosome,position,id,alt,ref,qual)
                         VALUES(?,?,?,?,?,?)
                    ''',snp_rows)
                    self.log("\tInserting genotypes...")
                    cur.executemany('''
                        INSERT INTO genotypes VALUES(?,?)
                    ''',genotype_rows)
                    self.log("Total execution time for block: {} seconds : {} per second",
                        time.time()-start_time,
                        (len(genotype_rows)+len(snp_rows)) / (time.time()-start_time)
                    )
                    cur.execute("END TRANSACTION")
                    genotype_rows = []
                    snp_rows = []
        self.log("Creating Indices")
        self._build_indices()
        self.log("Import finished")
                    
    def _build_indices(self):
        cur = self.db.cursor()
        cur.execute('''
            CREATE INDEX IF NOT EXISTS snploc ON snps (chromosome,position);
        ''')
        self.log("Building genotype index")
        cur.execute('''
            CREATE INDEX IF NOT EXISTS genoSnpID ON genotypes (snpRowID);
        ''')
        self.log("Done")

    def _drop_indices(self):
        cur = self.db.cursor()
        cur.execute('''
            DROP INDEX IF EXISTS snpChrom;
            DROP INDEX IF EXISTS snpPos;
            DROP INDEX IF EXISTS genpSnpID;
        ''')

    def _create_tables(self):
        cur = self.db.cursor()
        cur.execute(''' 
            CREATE TABLE IF NOT EXISTS snps (
                chromosome TEXT NOT NULL,
                position INTEGER NOT NULL,
                id TEXT NOT NULL,
                alt TEXT,
                ref TEXT,
                qual REAL
            );
        ''')
        cur.execute(''' 
            CREATE TABLE IF NOT EXISTS snp_info (
                snpID INTEGER,
                key TEXT,
                val TEXT
            )
        ''')
        # sqlite has rowids by default, use that for numbering accessions
        cur.execute('''
            CREATE TABLE IF NOT EXISTS accessions (
                name TEXT NOT NULL
            );  
        ''')    
        cur.execute(''' 
            CREATE TABLE IF NOT EXISTS genotypes (
                snpRowID INTEGER,
                alleles TEXT
            );
        ''')
