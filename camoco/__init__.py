""" 

Camoco Library - CoAnalysis of Molecular Components 

CacheMoneyCorn

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


def available_datasets(type=None,name=None):
    cur = Camoco("Camoco",type='Camoco').db.cursor()
    if type:
        datasets = cur.execute("SELECT type,name,description,added FROM datasets WHERE type = ? ORDER BY type;",(type,)).fetchall() 
    else:
        datasets = cur.execute("SELECT type,name,description,added FROM datasets ORDER BY type;").fetchall()
    if datasets:
        datasets = pd.DataFrame(datasets,columns=["Type","Name","Description","Date Added"])
    else:
        datasets = pd.DataFrame(columns=["Type","Name","Description","Date Added"])
    # Check to see if we are looking for a specific dataset
    if name is not None:
        return True if name in datasets['Name'].values else False
    else:
        return datasets

def del_dataset(type,name,safe=True):
    c = Camoco("Camoco")
    if safe:
        c.log("Are you sure you want to delete {}",name)
        if input("[Y/n]:") != 'Y':
            c.log("Nothing Deleted")
            return
    c.log("Deleting {}",name)
    c.db.cursor().execute(''' DELETE FROM datasets WHERE name = '{}' and type = '{}';'''.format(name,type))
    try:
        os.remove(c._resource("databases",".".join([type,name])+".db"))
    except FileNotFoundError as e:
        c.log('Database Not Found: {}'.format(e))
    try:
        os.remove(c._resource("databases",".".join([name])+".hd5"))
    except FileNotFoundError as e:
        c.log('Database Not Found: {}'.format(e))
    if type == 'Expr':
        # also have to remove the COB specific refgen
        del_dataset('RefGen','Filtered'+name,safe=safe)

def mv_dataset(type,name,new_name):
    c = Camoco("Camoco")
    c.db.cursor().execute("UPDATE datasets SET name = ? WHERE name = ? and type = ?",(new_name,name,type))
    os.rename(c._resource('databases','.'.join([type,name])+".db"),c._resource('databases',".".join([type,new_name])+".db"))

def redescibe_dataset(type,name,new_desc):
    c = Camoco("Camoco")
    c.db.cursor().execute("UPDATE datasets SET description = ? WHERE name = ? and type = ?",(new_desc,name,type))
