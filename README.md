Camoco
======
[![Build Status](https://travis-ci.org/schae234/Camoco.svg?branch=master)](https://travis-ci.org/schae234/Camoco)

Co-Analysis of Molecular Components
-----------------------------------

## Docker Install

Once you have docker installed and running on your computer (directions are
availible on docker's website), you simply run this command (possibly as sudo):

```
$ docker run -it --name='camoco' -v cobdata:/cobdata -v <folder with data to analyse>:/data monprin/camoco
```

If you would like to change where the databases created by analysis is stored, 
change 'cobdata' to a folder on your machine. Also replace the < > section to
the folder with you data. This will appear in the folder /data/ inside of the
container.

This will download the image, set where all the data should be, and then drops
into a shell, from here you can interact with camoco on the command line or type:

```
$ ipython3
```

And then you can interact with camoco in the interactive ipython shell. To end 
this session simply hit Ctrl+d until you are back at your home shell prompt.

To start the session again, just type:

```
$ docker start -ai camoco
```

And you will be right, where you left off!

There will be more documentaion later once I can work out some more kinks and run
this on a real computer.

## Introduction

Camoco is a python library modeled after the popular package pandas for
building and analyzing co-expression networks. It exposes a nice API and is
effecient for working with such large data structures. It was created to be
used interactively through something like iPython or normal scripting. All it
takes to get started is an import:

```python 
import camoco as co 
```

This will expose many of the cool classes camoco implements. Most of the
classes in the library inherit directly from the Camoco class. This is because
many of the data structures the library implements are **expensive** to make
and the library itself relies heavily on caching. sqlite3 and hdf5 files will
be created from you input datasets in a directory specified in `~/.camoco.conf`

Camoco 
------ 
The base camoco class is almost entirely abstract. This means that you
never really interact with a camoco class directly, instead many other classes
inherit properties from the camoco class. This is nice because, as I said
before, many of the data structures are expensive to compute. What camoco does
is to require datasets to have names and descriptions, so after they are built,
they can be cached on disk using sqlite and hdf5. This makes it fast to come
back! On disk objects are stored in a directory specified as `basedir` in the
configuration file: `~/.camoco.conf `. This is all done automatically, so if
camoco is doing its job, you shouldn't ever have to interact with it directly.
It should just work! (not guaranteed ... ).

The other methods which are exposed when you import camoco are housekeeping
methods. You can view, delete, move, and rename camoco datasets using the
following methods:

```python
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
chromosome. It all depends on context.

Installation
------------
Tagged releases are available on PyPi: `pip install camoco`. Otherwise,
you should be able to run `python setup.py install` once you have the following
dependencies resolved:

Required: python>=3.3

Hints on how we routinely build camoco can be found in our Travis-CI
configuration file, `.travis.yml`. Following these steps should result in a
successful build (at least on our test platform).

We decided to use APSW as an sqlite3 python wrapper, instructions on building
apsw can be found
[here](http://rogerbinns.github.io/apsw/build.html#recommended-build)

Required (install in this order to save pain):
+ cython>=0.16    
+ numpy>=1.9.1
+ scipy>=0.15
+ ipython>=3.1.0
+ pandas>=0.16
+ numexpr>=2.0.0
+ tables>=3.0.0 (hdf5)
+ apsw>=3.8
+ matplotlib>=1.4.3
+ statsmodels>=0.6.1
+ termcolor>=1.1.0
+ pyaml==15.6.3


Tests
-----
Unit tests within Camoco are implemented using the pytest framework. Unit tests
as well as Build tests can be found in the test directory.



CacheMoneyCorn
--------------
