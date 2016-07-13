#!/usr/bin/env python3
import apsw as lite
import os as os
import tempfile
import pandas as pd
import feather as ft

from .Tools import log
from .Config import cf
from .Exceptions import CamocoExistsError
from apsw import ConstraintError


class Camoco(object):

    '''
        Cache Money
    '''

    def __init__(self, name, type='Camoco'):
        # Set up our base directory
        self.log = log()
        self.type = type
        # A dataset already exists, return it
        self.db = self._database(name)
        try:
            (self.ID, self.name, self.description, self.type, self.added) = \
                self._database('Camoco', type='Camoco') \
                .cursor().execute(
                "SELECT rowid, * FROM datasets WHERE name = ? AND type = ?",
                (name, type)
            ).fetchone()
            cur = self.db.cursor()
            cur.execute('''
                CREATE TABLE IF NOT EXISTS globals (
                    key TEXT,
                    val TEXT
                );
                CREATE UNIQUE INDEX IF NOT EXISTS uniqkey ON globals(key)
                ''')
        except TypeError:
            raise TypeError('{}.{} does not exist'.format(type, name))

    def _database(self, dbname, type=None):
        # This lets us grab databases for other types
        if type is None:
            type = self.type
        # return a connection if exists
        return lite.Connection(
            os.path.expanduser(
                os.path.join(
                    cf.options.basedir,
                    'databases',
                    "{}.{}.db".format(type, dbname)
                )
            )
        )

    def _ft(self, tblname, dbname=None, type=None, df=None):
        if type is None:
            type = self.type
        if dbname is None:
            dbname = self.name
        if df is None:
            # return the dataframe if it exists 
            df = ft.read_dataframe(
                os.path.expanduser(
                    os.path.join(
                        cf.options.basedir,
                        'databases',
                        "{}.{}.{}.ft".format(type, dbname, tblname)
                    )
                )
            )
            if 'idx' in df.columns.values:
                df.set_index('idx', drop=True, inplace=True)
                df.index.name = None
            return df
        
        else:
            if not(df.index.dtype_str == 'int64') and not(df.empty):
                df = df.copy()
                df['idx'] = df.index
            ft.write_dataframe(df, 
                os.path.expanduser(
                    os.path.join(
                        cf.options.basedir,
                        'databases',
                        "{}.{}.{}.ft".format(type, dbname, tblname)
                    )
                )
            )
            if 'idx' in df.columns.values:
                del df
            return 

    def _tmpfile(self):
        # returns a handle to a tmp file
        return tempfile.NamedTemporaryFile(
            dir=os.path.expanduser(
                os.path.join(
                    cf.options.basedir,
                    "tmp"
                )   
            )
        )

    def _global(self, key, val=None):
        # set the global for the dataset
        try:
            if val is not None:
                self.db.cursor().execute('''
                    INSERT OR REPLACE INTO globals
                    (key, val)VALUES (?, ?)''', (key, val)
                )
            else:
                return self.db.cursor().execute(
                    '''SELECT val FROM globals WHERE key = ?''', (key, )
                ).fetchone()[0]
        except TypeError:
            # It pains me to do, but but return none if key isn't in global
            # TODO: replace returning None with an exception
            return None

    def __getattr__(self, name):
        return self._global(name)

    @classmethod
    def create(cls, name, description, type='Camoco'):
        '''
            This is a class method to create a new camoco type object.
            It initializes base directory hierarchy
        '''
        basedir = os.path.realpath(
            os.path.expanduser(cf.options.basedir)
        )

        # Create the basedir if not exists
        try:
            os.makedirs(basedir, exist_ok=True)
            os.makedirs(os.path.join(basedir, "logs"), exist_ok=True)
            os.makedirs(os.path.join(basedir, "databases"), exist_ok=True)
            os.makedirs(os.path.join(basedir, "analyses"), exist_ok=True)
            os.makedirs(os.path.join(basedir, "tmp"), exist_ok=True)
        except Exception:
            log('Could not create files in {}', basedir)
            raise
        try:
            # Create the base camoco database
            lite.Connection(
                os.path.join(basedir, 'databases', 'Camoco.Camoco.db')
            ).cursor().execute('''
                CREATE TABLE IF NOT EXISTS datasets (
                    name TEXT NOT NULL,
                    description TEXT,
                    type TEXT,
                    added datetime DEFAULT CURRENT_TIMESTAMP,
                    PRIMARY KEY(name, type)
                );
                INSERT OR IGNORE INTO datasets (name, description, type)
                    VALUES ('Camoco', 'Camoco base', 'Camoco');
                INSERT OR IGNORE INTO datasets (name, description, type)
                    VALUES (?, ?, ?)''', (name, description, type)
            )
        except ConstraintError:
            raise CamocoExistsError(
                'CAUTION! {}.{} Database already exists.', name, type
            )
        self = cls(name)
        return self
