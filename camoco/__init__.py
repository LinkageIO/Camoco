""" Camoco Library
"""

__license__ = """ 
    Creative Commons Non-Commercial 4.0 Generic
    http://creativecommons.org/licenses/by-nc/4.0/
"""

from camoco.Camoco import Camoco
from camoco.Expr import Expr
from camoco.COB import COB
from camoco.RefGen import RefGen
from camoco.Ontology import Ontology,Term
from camoco.HapMap import HapMap
from camoco.Locus import *
from camoco.Tools import *
from camoco.GEO import *
from camoco.Annotation import *
import pandas as pd
import os


def available_datasets(type=None, basedir="~/.camoco"):
    cur = Camoco("Camoco").db.cursor()
    if type:
        datasets = cur.execute("SELECT type,name,description,added FROM datasets WHERE type = ? ORDER BY type;",(type,)).fetchall() 
    else:
        datasets = cur.execute("SELECT type,name,description,added FROM datasets ORDER BY type;").fetchall()
    if datasets:
        return pd.DataFrame(datasets,columns=["Type","Name","Description","Date Added"])
    else:
        return pd.DataFrame(columns=["Type","Name","Description","Date Added"])
def del_dataset(type,name,safe=True,basedir="~/.camoco"):
    c = Camoco("Camoco")
    if safe:
        c.log("Are you sure you want to delete {}",name)
        if input("[Y/n]:") != 'Y':
            c.log("Nothing Deleted")
            return
    try:
        c.log("Deleting {}",name)
        c.db.cursor().execute(''' DELETE FROM datasets WHERE name = '{}' and type = '{}';'''.format(name,type))
        os.remove(c._resource("databases",".".join([type,name])+".db"))
        if type == 'Expr':
            # also have to remove the COB specific refgen
            del_dataset('RefGen',name+'RefGen')
    except FileNotFoundError as e:
        c.log('Database Not Found')

def mv_dataset(type,name,new_name,basedir="~/.camoco"):
    c = Camoco("Camoco")
    c.db.cursor().execute("UPDATE datasets SET name = ? WHERE name = ? and type = ?",(new_name,name,type))
    os.rename(c._resource('databases','.'.join([type,name])+".db"),c._resource('databases',".".join([type,new_name])+".db"))

def redescibe_dataset(type,name,new_desc,basedir="~/.camoco"):
    c = Camoco("Camoco")
    c.db.cursor().execute("UPDATE datasets SET description = ? WHERE name = ? and type = ?",(new_desc,name,type))
