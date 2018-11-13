.. Camoco documentation master file, created by
   sphinx-quickstart on Tue Dec 20 15:10:35 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _overview:
.. currentmodule:: camoco

####################################################
Camoco - gene co-expression network analysis toolkit 
####################################################

**Camoco** is a python library for building and analyzing gene co-expression
networks.  Networks are built from tables of gene expression data typically derived 
from RNA-seq or micro-array experiments. Once networks are built, there are several 
tools available to validate or check the *health* of the co-expression network using 
annotated ontologies such as the `Gene Ontology <http://www.geneontology.org/>`__.
Co-expression can then be calculated among sets of genes using several different
metrics.


.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Camoco
======

.. automodule:: camoco

COB
===
.. autoclass:: COB
    :members:
