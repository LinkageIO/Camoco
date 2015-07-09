#!/usr/bin/env python3
import apsw as lite
import os as os
import time as time
import sys
import tempfile
import pandas as pd

from camoco.Tools import log
from camoco.Config import cf
from apsw import ConstraintError

class Camoco(object):

    '''
        Cache Money
    '''

    def __init__(self,name,type='Camoco',basedir="~/.camoco"):
        # Set up our base directory
        self.log = log()
        self.type = type
        # A dataset already exists, return it
        self.db = self._database(name)
        try:
            (self.ID,self.name,self.description,self.type,self.added) = \
            self._database('Camoco',type='Camoco') \
                .cursor().execute(
                "SELECT rowid,* FROM datasets WHERE name = ? AND type = ?",
                (name,type)
            ).fetchone()
            cur = self.db.cursor()  
            cur.execute('''
                CREATE TABLE IF NOT EXISTS globals (
                    key TEXT,
                    val TEXT
                );
                CREATE UNIQUE INDEX IF NOT EXISTS uniqkey ON globals(key)
                ''')
        except TypeError as e:
            raise TypeError('{}.{} does not exist'.format(type,name))


    def _database(self,dbname,type=None):
        # return a connection if exists
        if type is None: # This lets us grab databases for other types
            type = self.type
        return lite.Connection(
            os.path.expanduser(
                os.path.join(
                    cf.get('options','basedir'),
                    'databases',
                    "{}.{}.db".format(type,dbname)
                )
            )
        )

    def _hdf5(self,dbname,type=None):
        if type is None:
            type = self.type
        # return a connection if exists
        return pd.HDFStore(
            os.path.expanduser(
                os.path.join(
                    cf.get('options','basedir'),
                    'databases',
                    "{}.{}.hd5".format(type,dbname)
                )
            )
        )

    def _tmpfile(self):
        # returns a handle to a tmp file
        return tempfile.NamedTemporaryFile(
            dir=os.path.join(
                cf.get('options','basedir'),
                "tmp"
            )
        )

    def _global(self,key,val=None):
        # set the global for the dataset
        if val is not None:
            self.db.cursor().execute('''
                INSERT OR REPLACE INTO globals 
                (key,val)VALUES (?,?)''',(key,val)
            )
        else:
            try:
                return self.db.cursor().execute(
                    '''SELECT val FROM globals WHERE key = ?''',(key,)
                ).fetchone()[0]
            except Exception as e:
                return None

    def __getattr__(self,name):
        return self._global(name)


    @classmethod
    def create(cls,name,description,type='Camoco'):
        '''
            This is a class method to create a new camoco type object.
            It initializes base directory hierarchy 
        '''
        basedir = os.path.realpath(
            os.path.expanduser(cf.get('options','basedir'))
        )

        # Create the basedir if not exists
        try:    
            os.makedirs(basedir,exist_ok=True)
            os.makedirs(os.path.join(basedir,"logs"),exist_ok=True)
            os.makedirs(os.path.join(basedir,"databases"),exist_ok=True)
            os.makedirs(os.path.join(basedir,"analyses"),exist_ok=True)
            os.makedirs(os.path.join(basedir,"tmp"),exist_ok=True)
        except Exception as e:
            log(' Could not create files in {}',basedir)
            raise
        try:
        # Create the base camoco database
            lite.Connection(
                os.path.join(basedir,'databases','Camoco.Camoco.db')
            ).cursor().execute(''' 
                CREATE TABLE IF NOT EXISTS datasets (
                    name TEXT NOT NULL,
                    description TEXT,
                    type TEXT,
                    added datetime DEFAULT CURRENT_TIMESTAMP,
                    PRIMARY KEY(name,type)
                );
                INSERT OR IGNORE INTO datasets (name,description,type)
                VALUES ('Camoco','Camoco base','Camoco');
                INSERT OR FAIL INTO datasets (name,description,type)
                VALUES (?,?,?)''',(name,description,type)
            )
        except ConstraintError as e:
            log.warn('CAUTION! {}.{} Database already exists.',name,type)
        self = cls(name) 
        return self


