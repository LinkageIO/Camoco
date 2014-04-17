#!/usr/bin/python
from __future__ import print_function
import pandas as pd
import pandas.io.sql as psql
import numpy as np
import MySQLdb
import sys
import os as os
import time
import itertools
import tempfile

class COBDatabase(object):
    def sensitive(sens):
        ''' retrieves sensitive variables from the 'sensitive' directory in the 
            base dir. Useful for when you want to post code online and not have
            people trash your database... '''
        with open(os.path.join("/heap/cobDev/pyCOB/sensitive",sens),'r') as sens_file:
            return sens_file.readline().strip()
    # The Database Connection needs to be shared between all classes which inherit COBDatabase
    # Otherwise we open too many connections
    db = MySQLdb.connect(
        host   = sensitive("MySQLHost"),
        user   = sensitive("MySQLUser"),
        passwd = sensitive("MySQLPasswd"),
        db     = sensitive("MySQLDB")
    )
    def __init__(self, log_file=sys.stderr):
        if log_file != sys.stderr:
            self.log_file = open(log_file,'w')
        else:
            self.log_file = log_file
    def query(self, query, *variables):
        ''' Returns a pandas dataframe '''
        query = query.format(*variables)
        return psql.frame_query(query,self.db)
    def execute(self,query,*variables):
        ''' perform a database execution '''
        cur = self.db.cursor()
        cur.execute(query.format(*variables))
        cur.close()
    def fetchone(self,query,*variables):
        ''' Used for queries which only should have one row '''
        cur = self.db.cursor()
        cur.execute(query.format(*variables))
        return cur.fetchone()

    def log(self, *args):
        print("[COB LOG] ",time.ctime(),"-",*args,file=self.log_file)        

class COBDatabaseBuilder(COBDatabase): 
    def __init__(self):
        super(COBDatabaseBuilder,self).__init__()
        
    def add_gene(self,GrameneID,chromosome,chromo_start,chromo_end,strand,build):
        '''adds the gene name to the database '''
        self.execute('''
            INSERT IGNORE INTO genes (GrameneID,chromosome,chromo_start,chromo_end,strand,build) 
            VALUES ('{}',{},{},{},{},'{}')'''.format(GrameneID,chromosome,chromo_start,chromo_end,strand,build)
        )
    def add_COBGene(self,cg):
        ''' Adds an entry based on a COBGene '''
        self.add_gene(cg.GrameneID,cg.chromosome,cg.chromo_start,cg.chromo_end,cg.strand,cg.build)

    def add_accession(self,dataset_id, name, short_desc='', long_desc=''):
        self.execute('''INSERT IGNORE INTO accessions (DatasetID, name, class, description) 
            VALUES ({}, {}, {}, {})'''.format(dataset_id, name, short_desc, long_desc))

    def add_coex(self,dataset_id, gene_a, gene_b, score):
        self.execute('''INSERT IGNORE INTO coex (DatasetID, gene_a, gene_b, score) 
            VALUES ({}, {}, {}, {});'''.format(dataset_id, gene_a, gene_b, score)
        )


    def add_dataset(self,dataset,transform=np.arctanh,significance_thresh=3): 
        ''' Imports a COBDataset into the Database '''
        self.execute('''INSERT INTO datasets (name, description, FPKM, gene_build) 
        VALUES ('{}', '{}', {}, '{}')''',dataset.name, dataset.description, dataset.FPKM, dataset.gene_build)
        # This fetches the primary key for dataset
        dataset.id = self.db.insert_id()
        # Check that the genes are all in the database
        assert len(set([x.GrameneID for x in dataset.genes]) - set(self.query("SELECT GrameneID FROM genes").GrameneID))  == 0, "You must import all genes before you import the Dataset"
        # Create a Dataframe in memory
        self.log("Calculating Z-Scores for {}".format(dataset.name))
        scores = dataset.coex()
        # Calculate fisher transform
        scores = transform(scores)
        # Normalize to standard normal dist with mean == 0 and std == 1
        scores = (scores-scores.mean())/scores.std()
        # Create the dataframe
        tbl = pd.DataFrame(list(itertools.combinations([gene.ID for gene in dataset.genes],2)),columns=['gene_a','gene_b'])
        tbl['DatasetID'] = dataset.id
        tbl['score'] = scores
        tbl['significant'] = np.array(tbl['score'] >= significance_thresh,dtype=int)
        # Disable Keys on the coex table
        self.log('Adding Raw DataFrame to Database')
        self.execute("ALTER TABLE coex DISABLE KEYS")
        tmp_file = '/tmp/tmp{}.txt'.format(dataset.id)
        with open(tmp_file,'w') as FOUT:
            tbl[['DatasetID','gene_a','gene_b','score','significant']].to_csv(FOUT,sep="\t",index=False,header=False)
        self.execute("LOAD DATA INFILE '{}' INTO TABLE coex FIELDS TERMINATED BY '\t';".format(tmp_file))
        self.execute("ALTER TABLE coex ENABLE KEYS")
        
        self.log('Done Adding to Database')

    def clear_datasets(self):
        self.execute('''
            TRUNCATE TABLE coex; 
            TRUNCATE TABLE datasets;
            DELETE FROM coex;
            DELETE FROM datasets;
            ALTER TABLE datasets AUTO_INCREMENT = 1;
        ''')

    def __create_tables__(self): 
        self.execute(''' 
        CREATE TABLE  IF NOT EXISTS datasets(
            ID INT UNSIGNED AUTO_INCREMENT,
            name varchar(256),
            description varchar(4096),
            FPKM BOOL,
            gene_build ENUM('4a.53','5a','5b'),
            PRIMARY KEY(ID) 
        )
        ''')
        self.execute('''
        CREATE TABLE IF NOT EXISTS genes ( 
            ID INT UNSIGNED AUTO_INCREMENT, 
            GrameneID varchar(128), INDEX USING BTREE (GrameneID), 
            chromosome TINYINT UNSIGNED, 
            chromo_start INT UNSIGNED, 
            chromo_end INT UNSIGNED, 
            strand TINYINT, 
            build ENUM('4a.53','5a','5b'), 
            PRIMARY KEY(ID)
        );
        ''')
        self.execute(''' 
        CREATE TABLE IF NOT EXISTS accessions (
            ID INT UNSIGNED AUTO_INCREMENT,
            DatasetID INT UNSIGNED NOT NULL,
            name varchar(128) NOT NULL,
            class varchar(128),
            description TEXT NOT NULL,
            PRIMARY KEY(ID)
        );
        ''')
        self.execute(''' 
        CREATE TABLE IF NOT EXISTS expression (
            GeneID INT UNSIGNED,
            AccessionID INT UNSIGNED,
            DatasetID INT UNSIGNED,
            value FLOAT,
            PRIMARY KEY(GeneID,AccessionID,DatasetID)
        );
        ''')
        self.execute(''' 
        CREATE TABLE IF NOT EXISTS coex (
            DatasetID INT UNSIGNED, INDEX USING BTREE(DatasetID),
            gene_a INT UNSIGNED, INDEX USING BTREE(gene_a),
            gene_b INT UNSIGNED, INDEX USING BTREE(gene_b),
            score FLOAT, INDEX USING BTREE(score),
            significant BOOL, INDEX USING BTREE(significant),
            PRIMARY KEY(DatasetID,gene_a,gene_b)
        );
        ''')
        self.execute('''
        CREATE TABLE IF NOT EXISTS degree (
            DatasetID INT UNSIGNED, INDEX USING BTREE(DatasetID),
            GeneID INT UNSIGNED, INDEX USING BTREE(GeneID),
            Degree INT UNSIGNED
        )
        ''')
        self.execute('''
        CREATE TABLE IF NOT EXISTS gene_ontology(
            GeneID INT UNSIGNED,
            OntologyID INT UNSIGNED,
            PRIMARY KEY(GeneID,OntologyID)
        )
        ''')
        self.execute('''
        CREATE TABLE IF NOT EXISTS ontology (
            ID INT UNSIGNED AUTO_INCREMENT,
            name varchar(128),
            short_desc TEXT,
            long_desc MEDIUMTEXT,
            PRIMARY KEY(ID)
        );
        ''')

