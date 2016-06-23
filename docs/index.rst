.. Camoco documentation master file, created by
   sphinx-quickstart on Thu Feb 11 21:29:35 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Camoco
======

Contents:

.. toctree::
   :maxdepth: 2



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


Camoco 
------ 
The base camoco class is almost entirely abstract. Camoco require datasets 
to have names and descriptions, so after they are built,
they can be cached on disk using sqlite and hdf5. This makes it fast to come
back! On disk objects are stored in a directory specified as ``basedir`` in the
configuration file: ``~/.camoco.conf``. This is all done automatically, so if
camoco is doing its job, you shouldn't ever have to interact with it directly.
It should just work! (not guaranteed ... ).

The other methods which are exposed when you import camoco are housekeeping
methods. You can view, delete, move, and rename camoco datasets using the
following methods:

```
# See docstrings for more information
co.available_datasets()
co.del_dataset()
co.mv_dataset()
co.redescribe_dataset() 
```

Other than that, Camoco has three main classes which each have some companion
classes.  These classes/subclasses are as follows:

+ RefGen and Locus
+ COB and Expr
+ Ontology and Term


RefGen and Locus
----------------
Python has a fantastic range of built in data types. It lacks however in a data
type which succinctly represents a locus. Loci cover a large range of
biological data types; a gene is a locus, as is a SNP, or a QTL, or an entire
chromosome. It all depends on context. The Locus and RefGen classes implement
a basic locus data type in python. These objects are extensively used as parameters
to functions throughout Camoco.

