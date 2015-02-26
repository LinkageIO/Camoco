#!/usr/bin/python
from __future__ import print_function,division
import apsw as lite
import os as os
import time as time
import sys
import tempfile
import pandas as pd

from camoco.Tools import log

class Camoco(object):
    def __init__(self,name="Camoco",description=None,type='Camoco',basedir="~/.camoco"):
        # Set up our base directory
        self.basedir = os.path.realpath(os.path.expanduser(basedir))
        self._create_base_tables()
        self.log = log()
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
        self.log_file = open(os.path.join(self.basedir,"logs","logfile.txt"),'a')
        if description is not None:
            # create new sqlite file for dataset
            try:
                self.database('Camoco.Camoco').cursor().execute('''
                    INSERT OR FAIL INTO datasets (name,description,type)
                    VALUES (?,?,?)''',(name,description,type))
            except Exception as e:
                self.log.warn('CAUTION! {}.{} Database already exists. Dont do anything rash.',name,type)
        try:
            # A dataset already exists, return it
            self.db = self.database(".".join([type,name]))
            (self.ID,self.name,self.description,
                self.type,self.added) = self.database('Camoco.Camoco').cursor().execute(
                "SELECT rowid,* FROM datasets WHERE name = ? AND type = ?",
                (name,type)
            ).fetchone()
            cur = self.db.cursor()  
            cur.execute('PRAGMA main.page_size = 8192;')
            cur.execute('PRAGMA main.cache_size= 1000000;')
            cur.execute('PRAGMA main.temp_store = MEMORY;')
            #cur.execute('PRAGMA main.locking_mode=EXCLUSIVE;')
            cur.execute('PRAGMA main.synchronous=OFF;')
            cur.execute('PRAGMA count_changes=OFF')
            cur.execute('PRAGMA main.journal_mode=MEMORY;')
            cur.execute('''
                CREATE TABLE IF NOT EXISTS globals (
                    key TEXT,
                    val TEXT
                );
                CREATE UNIQUE INDEX IF NOT EXISTS uniqkey ON globals(key)
                ''')
        except Exception as e:
            raise NameError("Camoco dataset {} does not exist: {}".format(name,e))

    def _resource(self,type,filename):
        return os.path.expanduser(os.path.join(self.basedir,type,filename))

    def database(self,dbname):
        # return a connection if exists
        return lite.Connection(self._resource("databases",str(dbname)+".db"))

    def _tmp_file(self):
        # returns a handle to a tmp file
        return tempfile.NamedTemporaryFile(dir=os.path.join(self.basedir,"tmp"))

    def _global(self,key,val=None):
        # set the global for the dataset
        if val is not None:
            self.db.cursor().execute('''
                INSERT OR REPLACE INTO globals (key,val)VALUES (?,?)''',(key,val)
            )
        else:
            try:
                return self.db.cursor().execute(
                    '''SELECT val FROM globals WHERE key = ?''',(key,)
                ).fetchone()[0]
            except Exception as e:
                return None
            
    def _create_base_tables(self):
        camocodb = self.database("camoco")
        camocodb.cursor().execute(''' 
            CREATE TABLE IF NOT EXISTS datasets (
                name TEXT NOT NULL,
                description TEXT,
                type TEXT,
                added datetime DEFAULT CURRENT_TIMESTAMP,
                PRIMARY KEY(name,type)
            )
        ''')

    def __getattr__(self,name):
        return self._global(name)

