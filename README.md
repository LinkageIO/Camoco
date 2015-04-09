Camoco
======

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
`~/.camoco.conf `. This is all done automatically, and if camoco is doing its job, you shouldn't
ever have to interact with it directly. It should just work!

The other methods which are exposed when you import camoco are housekeeping methods. You can view, delete,
move, and rename camoco datasets using the following methods:

```python
    # See docstrings for more information
    co.available_datasets()
    co.del_dataset()
    co.mv_dataset()
    co.redescribe_dataset() 
```




CacheMoneyCorn
--------------

