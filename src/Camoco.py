#!/usr/bin/python
from __future__ import print_function,division
import sqlite3 as lite
import os as os
import time as time
import sys
import tempfile

class Camoco(object):
    def __init__(self,basedir="~/.camoco"):
        # Set up our base directory
        self.basedir = os.path.realpath(os.path.expanduser(basedir))
        if not os.path.exists(self.basedir):
            # If it doesn't exists, set up first time directories
            try:    
                os.mkdir(self.basedir)
                os.mkdir(os.path.join(self.basedir,"logs"))
                os.mkdir(os.path.join(self.basedir,"databases"))
                os.mkdir(os.path.join(self.basedir,"analyses"))
                os.mkdir(os.path.join(self.basedir,"tmp"))
            except Exception as e:
                print("[CAMOCO]",time.ctime(), '-', *args,file=sys.stderr)
        self.log_file = open(os.path.join(self.basedir,"logs","logfile.txt"),'w')


    def log(self,msg,*args):
        print("[CAMOCO]",time.ctime(), '-', msg.format(*args),file=sys.stderr)
        print("[CAMOCO]",time.ctime(), '-', msg.format(*args),file=self.log_file)
           
    def _resource(self,type,filename):
        return os.path.join(self.basedir,type,filename)

    def database(self,dbname):
        # return a connection if exists
        return lite.connect(self._resource("databases",str(dbname)+".db"))
    def tmp_file(self):
        # returns a handle to a tmp file
        return tempfile.NamedTemporaryFile(dir=os.path.join(self.basedir,"tmp"))
        

    def available_datasets(self):
        cur=self.database("camoco").cursor()
        return cur.execute("SELECT rowid,* FROM datasets;").fetchall()
    

    def add_dataset(self,name,description='',type=''):
        con = self.database("camoco")
        cur = con.cursor()
        try:
            # update information table
            cur.execute('''
                INSERT INTO datasets (name,description,type)
                VALUES ('{}','{}','{}')'''.format(name,description,type))
            con.commit()
            # create new sqlite file for dataset
            self.database(".".join([type,name]))
        except lite.IntegrityError as e:
            self.log("{} already added to database.",name)

    def del_dataset(self,name,type):
        con = self.database('camoco')
        con.execute(''' DELETE FROM datasets WHERE name = '{}' and type = '{}';'''.format(name,type))
        con.commit() 
        os.remove(self._resource("databases",".".join([type,name])+".db"))

    def _create_tables(self):
        camocodb = self.database("camoco")
        camocodb.execute(''' 
            CREATE TABLE IF NOT EXISTS datasets (
                name TEXT NOT NULL,
                description TEXT,
                type TEXT,
                added datetime DEFAULT CURRENT_TIMESTAMP,
                PRIMARY KEY(name,type)
            )
        ''')


