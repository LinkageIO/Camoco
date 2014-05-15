#!/usr/bin/python3
from __future__ import print_function,division
#import sqlite3 as lite
import apsw as lite
import os as os
import time as time
import sys

from Camoco import Camoco

class HapMap(Camoco):
    def __init__(self,name,basedir="~/.camoco"):
        # initialize super class
        super(self.__class__,self).__init__(basedir)
        # get meta information
        (ID,name,description,type,added) = self.database('camoco').execute(
            "SELECT rowid,* FROM datasets WHERE name = '{}' AND type = 'HapMap'",).fetchone()
        self.db = self.database(".".join([self.type,self.name]))
        


class HapMapBuilder(Camoco):
    def __init__(self,name,description="",basedir="~/.camoco"):
        # add entry to camoco database
        self.name = name
        self.description = description
        self.type = "HapMap"
        self.basedir = basedir
        # add yourself to the database
        super(HapMapBuilder,self).__init__(basedir)
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
            cur.execute("BEGIN TRANSACTION")
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
                    for accessionRowID,genotype in enumerate(genotypes):
                        (allele1,allele2) = genotype.split(":")[GT_ind].replace("|","/").split("/")
                        genotype_rows.append((snpRowID,accessionRowID+1,allele1,allele2))
                    snpRowID += 1 # <- Do NOT fuck this up, genotypes need to reference correct SNPId
                if i % record_dump_size == 0 and i > 0:
                    self.log(" {} records reached. dumping data",i)
                    start_time = time.time()
                    # Bulk import the genotype values
                    self.log("\tInserting snps...")
                    cur.executemany(''' 
                        INSERT INTO snps (chromosome,position,id,alt,ref,qual)
                         VALUES(?,?,?,?,?,?)
                    ''',snp_rows)
                    self.log("\tInserting genotypes...")
                    cur.executemany('''
                        INSERT INTO genotypes VALUES(?,?,?,?)
                    ''',genotype_rows)
                    self.log("Total execution time for block: {} seconds : {} per second",
                        time.time()-start_time,
                        (len(genotype_rows)+len(snp_rows)) / (time.time()-start_time)
                    )
                    genotype_rows = []
                    snp_rows = []
        cur.execute("END TRANSACTION")
        self.log("Creating Indices")
        self._build_indices()
        self.log("Import finished")
                    
    def _build_indices(self):
        self._drop_indices()
        cur = self.db.cursor()
        self.log("Building snp chromo index")
        cur.execute('''
            CREATE INDEX IF NOT EXISTS snpChrom ON snps (chromosome);
        ''')
        self.log("Building snp position index")
        cur.execute('''
            CREATE INDEX IF NOT EXISTS snpPos ON snps (position);
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
            FROP INDEX IF EXISTS genpSnpID;
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
                accessionRowID INTEGER,
                allele1 INTEGER,
                allele2 INTEGER
            );
        ''')
