""" Camoco Library
"""

__license__ = """ 
The MIT License (MIT)
Copyright (c) 2013 Rob Schaefer
Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from camoco.Camoco import Camoco
from camoco.COB import COB
from camoco.RefGen import RefGen
from camoco.Ontology import Ontology
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
        if input("[Y/n]:") == 'Y':
            c.db.cursor().execute(''' DELETE FROM datasets WHERE name = '{}' and type = '{}';'''.format(name,type))
            os.remove(c._resource("databases",".".join([type,name])+".db"))
        else:
            c.log("Nothing Deleted")
    else:
        c.db.cursor().execute(''' DELETE FROM datasets WHERE name = '{}' and type = '{}';'''.format(name,type))
        os.remove(c._resource("databases",".".join([type,name])+".db"))
def mv_dataset(type,name,new_name,basedir="~/.camoco"):
    c = Camoco("Camoco")
    c.db.cursor().execute("UPDATE datasets SET name = ? WHERE name = ? and type = ?",(new_name,name,type))
    os.rename(c._resource('databases','.'.join([type,name])+".db"),c._resource('databases',".".join([type,new_name])+".db"))

def redescibe_dataset(type,name,new_desc,basedir="~/.camoco"):
    c = Camoco("Camoco")
    c.db.cursor().execute("UPDATE datasets SET description = ? WHERE name = ? and type = ?",(new_desc,name,type))
