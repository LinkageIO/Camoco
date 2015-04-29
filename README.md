Camoco
======
[![Build Status](https://travis-ci.org/schae234/Camoco.svg?branch=master)](https://travis-ci.org/schae234/Camoco)

Installation
------------
Newest tagged releases are available on PyPi: `pip install camoco`. Otherwise, you should be able to run 
`python setup.py install` once you have the following dependencies resolved:

Required:
+ pandas>=0.16
+ tables>=3.0.0 (hdf5)
+ apsw>=3.8
+ scipy>=0.15
+ numpy>=1.9.1
+ cython>=0.16    
+ matplotlib>=1.4.3

Optional:
+ igraph>=0.1.5
+ termcolor>=1.1.0


Co-Analysis of Molecular Components
-----------------------------------

Camoco is a python library modeled after the popular package pandas for building 
and analyzing co-expression networks. It exposes a nice API and is effecient for
working with such large data structures. It was created to be used interactively 
through something like iPython or normal scripting. All it takes to get started 
is an import:

```python
import camoco as co
```

This will expose many of the cool classes camoco implements. Most of the classes
in the library inherit directly from the Camoco class. This is because many of the
data structures the library implements are **expensive** to make and the library
itself relies heavily on caching. 

Camoco 
------
The base camoco class is almost entirely abstract. This means that you never really
interact with a camoco class directly, instead many other classes inherit properties
from the camoco class. This is nice because, as I said before, many of the data structures
are expensive to compute. What camoco does is to require datasets to have names and descriptions,
so after they are built, they can be cached on disk using sqlite and hdf5. This makes it fast
to come back! On disk objects are stored in a directory specified as `basedir` in the configuration file:
`~/.camoco.conf `. This is all done automatically, so if camoco is doing its job, you shouldn't
ever have to interact with it directly. It should just work! (not guaranteed ... ).

The other methods which are exposed when you import camoco are housekeeping methods. You can view, delete,
move, and rename camoco datasets using the following methods:

```python
# See docstrings for more information
co.available_datasets()
co.del_dataset()
co.mv_dataset()
co.redescribe_dataset() 
```

Other than that, Camoco has three main classes which each have some companion classes.
These classes/subclasses are as follows:

+ RefGen and Locus
+ COB and Expr
+ Ontology and Term


RefGen and Locus
----------------
Python has a fantastic range of built in data types. It lacks however in a data type which succinctly 
represents a locus. Loci cover a large range of biological data types; a gene is a locus, as is a SNP,
or a QTL, or an entire chromosome. It all depends on context. The 


CacheMoneyCorn
--------------
