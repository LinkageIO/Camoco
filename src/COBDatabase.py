#!/usr/bin/python

import pandas as pd
import pandas.io.sql as psql
import MySQLdb
import sys
import os.path as path

class COBDatabase(object):
    def sensitive(sens):
        ''' retrieves sensitive variables from the 'sensitive' directory in the 
            base dir. Useful for when you want to post code online and not have
            people trash your database... '''
        with open(path.join("/heap/cobDev/pyCOB/sensitive",sens),'r') as sens_file:
            return sens_file.readline().strip()
    # The Database Connection needs to be shared between all classes which inherit COBDatabase
    # Otherwise we open too many connections
    db = MySQLdb.connect(
        host   = sensitive("MySQLHost"),
        user   = sensitive("MySQLUser"),
        passwd = sensitive("MySQLPasswd"),
        db     = sensitive("MySQLDB")
    )
    def __init__(self):
        self.__create_tables__()
    def query(self, query, *variables):
        ''' Returns a pandas dataframe '''
        query = query.format(*variables)
        return psql.frame_query(query,self.db)
    def execute(self,query,*variables):
        ''' perform a database execution '''
        cur = self.db.cursor()
        cur.execute(query.format(*variables))
    def add_gene(self,GrameneID,chromosome,chromo_start,chromo_end,strand,build):
        '''adds the gene name to the database '''
        cur = self.db.cursor()
        cur.execute('''
            INSERT IGNORE INTO genes (GrameneID,chromosome,chromo_start,chromo_end,strand,build) 
            VALUES ('{}',{},{},{},{},'{}')'''.format(GrameneID,chromosome,chromo_start,chromo_end,strand,build)
        )
    def add_accession(self):
        pass

    def __create_tables__(self): 
        self.execute(''' 
        CREATE TABLE  IF NOT EXISTS datasets(
            ID INT UNSIGNED AUTO_INCREMENT,
            name varchar(256),
            description varchar(4096),
            FPKM BOOL,
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
            build ENUM('4.a53','5a','5b'), 
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
            DatasetID INT UNSIGNED,
            gene_a INT UNSIGNED, INDEX USING BTREE(gene_a),
            gene_b INT UNSIGNED, INDEX USING BTREE(gene_b),
            score FLOAT, INDEX USING BTREE(score),
            PRIMARY KEY(DatasetID,gene_a,gene_b)
        );
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
 
