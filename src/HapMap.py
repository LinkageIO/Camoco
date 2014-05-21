#!/usr/bin/python3
from __future__ import print_function,division
#import sqlite3 as lite
import apsw as lite
import os as os
import time as time
import sys
import numpy as np
from scipy.stats import pearsonr

from Camoco import Camoco
from Locus import SNP

from matplotlib import pylab

class HapMap(Camoco):
    def __init__(self,name,basedir="~/.camoco"):
        # initialize super class
        super(self.__class__,self).__init__(basedir)
        # get meta information
        (ID,name,description,type,added) = self.database('camoco').cursor().execute(
            "SELECT rowid,* FROM datasets WHERE name = ? AND type = ?",(name,"HapMap")).fetchone()
        self.name = name
        self.description = description
        self.type = "HapMap"
        self.basedir = basedir
        self.db = self.database(".".join([self.type,self.name]))
    
    def snps(self,chromo,start,end):
        ''' Returns SNPs within range '''
        return self.db.cursor().execute(
            '''SELECT rowid,* FROM snps WHERE chromosome = ? AND position >= ? AND position <= ?''',
            (chromo,start,end)
        )
    def genotypes(self,snp):
        ''' returns the genotypes as a list based on position or snpid,
            returns None if snp not in database '''
        genos = [list(map(int,list(x[0]))) for x in self.db.cursor().execute(
        ''' SELECT alleles FROM snps JOIN genotypes ON snps.rowid = genotypes.snpRowID
            WHERE snps.chromosome = ? AND snps.position = ?''', (snp.chrom,snp.pos)
        )]
        if len(genos) == 0:
            return None
        return np.array(genos[0])

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

    def downstream_snps(self,snp,snp_limit=1000,pos_limit=1000000):
        ''' returns snps downstream of snp '''
        return [SNP(*x) for x in self.db.cursor().execute(''' 
            SELECT chromosome,position,id FROM snps 
            WHERE chromosome = ?
            AND position >= ?
            AND position <= ?
            ORDER BY position ASC
            LIMIT ?
            ''',(snp.chrom,snp.pos,snp.pos+pos_limit,int(snp_limit)))]
    def upstream_snps(self,snp,snp_limit=1000,pos_limit=1000000):
        ''' returns snps upstream of input '''
        return [SNP(*x) for x in self.db.cursor().execute(''' 
            SELECT chromosome,position,id FROM snps 
            WHERE chromosome = ?
            AND position <= ?
            AND position >= ?
            ORDER BY position DESC
            LIMIT ?
            ''',(snp.chrom,snp.pos,snp.pos-pos_limit,int(snp_limit)))]

    def flanking_snps(self,snp,limit=100):
        ''' returns snps flanking '''
        return self.upstream_snps(snp,limit=limit/2) + self.downstream_snps(snp,limit=limit/2)

    def nearest_snp(self,snp):
        ''' retruns the nearest upstream or downstream SNP '''
        nearest = self.upstream_snps(snp,snp_limit=1) + self.downstream_snps(snp,snp_limit=1)
        if len(nearest) == 1:
            return nearest[0]
        elif abs(snp - nearest[0]) < abs(snp - nearest[1]):
            return nearest[0]
        else:
            return nearest[1]
        return nearest
        
    def ld_interval(self,snp,edge_cutoff = 0.5,snp_limit=10,pos_limit=10000):
        ''' Returns the interval around a locus/snp based on 
            a linear regression of LD ~ distance-from-locus.
            This implementation calls the edge when the best fit
            line drops below the edge_cutoff parameter '''
        # grab upstream snp
        downstream = self.downstream_snps(snp,snp_limit=snp_limit,pos_limit=pos_limit) 
        # mask out snps with < 80 points with filter
        ld = list(filter(lambda x: x[1] > 95, [self.ld(snp,x) for x in downstream]))
        # filter out nans
        distances = [x[2] for x in ld]
        r2 = [x[0] for x in ld]
        # fit a line  
        #best_fit = np.polyfit(distances,r2,deg=1,full=True)
        #slope = best_fit[0][0]
        #intercept = best_fit[0][1]
        a=(distances,r2)#,intercept,slope)
        pylab.clf()
        pylab.scatter(a[0],a[1])
        #pylab.plot(np.arange(max(a[0])),a[3]*np.arange(max(a[0]))+a[2])     
        pylab.savefig('tmp.png')


class HapMapBuilder(Camoco):
    def __init__(self,name,description="",basedir="~/.camoco"):
        # add entry to camoco database
        self.name = name
        self.description = description
        self.type = "HapMap"
        self.basedir = basedir
        # add yourself to the database
        super().__init__(basedir)
        self.add_dataset(name=name,description=description,type="HapMap")
        self.db = self.database(".".join([self.type,self.name]))
        self._create_tables()

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
                    snp_rows.append((chrom,pos,snpid,ref,alt,qual))
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
                    alleles = "".join(
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
        cur.execute("PRAGMA page_size = 262144;")
        cur.execute("PRAGMA cache_size = 100000;")
        cur.execute(''' 
            CREATE TABLE IF NOT EXISTS snps (
                chromosome INTEGER NOT NULL,
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
