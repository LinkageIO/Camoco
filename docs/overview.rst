
.. _overview:
.. currentmodule:: camoco
.. ipython:: python
    :suppress:

    import camoco as co

########
Overview
########

Camoco handles the creation and the co-analysis/integration of several different
components. We follow a *build once, use many times* philosophy as many of the data
sets are computationally expensive to build. Camoco also was designed so that many 
of the standard analyses and procedures are available through the CLI. These procedures
include building networks and other datasets, listing available datasets, or calculating 
the co-expression overlap between a network and GWAS data. There are also more advanced 
work flows that can be created by importing Camoco as a python library or by analyzing 
data interactively using `ipython <https://ipython.readthedocs.io/en/stable/>`__.

Once Camoco datasets are built, they can be co-analyzed together to determine if there 
is significant *overlap* between the datasets. This overlap vocabulary is used throughout
this documentation to mean a statistical relationship between different types of data. For
example, one analysis that Camoco offers is to calculate the *overlap* between GWAS data
and a co-expression network. This analysis includes mapping GWAS SNPs to nearby genes, 
calculating the co-expression among the genes, then comparing the co-expression score to
randomized sets of genes to determine the statistical significance of the overlap. Likewise,
there are methods to calculate the overlap between a set of genes (potentially from a GWAS) 
and the clusters of a network. This *overlap* is calculated using an enrichment statistic: 
the hyper-geometric distribution.

In order to better understand the analysis being computed, its necessary to first define
some terms to represent the types of data Camoco can co-analyze. This section will mainly
describe the data-types at a high level while code and examples can be found in the tutorial
section.


Camoco Data Types
==================
There are several different types of data that Camoco build and in turn, compare. Many of these
data types are computationally expensive to generate and thus, Camoco stores them in using internal
databases. Several of the data types represent collections of other data types. For instance, there
are Camoco objects to represent a single gene, as well as a named set of genes. There is also a 
distinction between what objects are persistent (stored in the internal database) and objects that
are disposable. 

In general, if the object is *named* it is stored in a database. Also, persistent objects, in general,
generate disposable objects. This is better illustrated with an example. Camoco determines the 
significance of co-expression among a set of genes by comparing to the co-expression among random sets
of genes. To get these random sets of genes, Camoco generates them from a reference genome object. As
reference genomes include a large amount of data and are computationally expensive to generate,
RefGen objects are created once and are stored in a database. From these RefGen objects, random sets
of genes can be generated very quickly. 


Camoco Data Types include:

Locus
  The most fundamental element that Camoco uses is a :class:`Locus` object. A Locus is general purpose
  object that encodes a coordinate in the genome. It consists of a three mandatory elements:
  a chromosome name, a start position, and an end position. Locus objects are also the basis of 
  other genome coordinate objects such as SNPs and genes, which are just Loci with additional
  information.

RefGen*
  This object contains a reference set of Locus object, such as those that are defined in a `GFF
  File <https://uswest.ensembl.org/info/website/upload/gff.html>`__. Typically, RefGen objects are
  used by other Camoco objects to generate genes either by ID or by genome location.

Term
  A :class:`Term` object is collection of Loci objects with a common form or function. This is 
  commonly used as a container for loci in a :class:`Ontology` or :class:`GWAS`.

Ontology*
  This is a object containing a set of Terms and their hierarchical relationships. This is used 
  to store data related to things like `Gene Onotologies (GO) <http://www.geneontology.org/>`__ 
  or `MapMan <https://mapman.gabipd.org/>`__ pathways.
 
GWAS*
  Similar to an Ontology, the terms in a GWAS are SNPs. A GWAS object handles SNP to gene mapping
  by referencing a RefGen object.

COB*
  This is the co-expression browser object. A COB is the network object in Camoco and used to
  perform all the network analysis components.


\* Designates a persistent data type

Building data sets
==================

Datasets can be built in two main ways, the first is from raw data and the second is from
existing Camoco data sets. For example, a RefGen object is typically built from a GFF file.
From a RefGen object, Loci objects can be created. Similarly Ontologies and GWAS are built
from raw data files while Terms are created from Ontologies and GWAS objects. 

COB objects are built from raw gene expression datasets, however, they also require that 
a valid RefGen object is built. During this process, a filtered RefGen object is automatically
created containing only the genes in the expression data.

When a persistent Camoco object is built, a user specified name is required. This name is used
in analyses where specific datasets need to be referenced. For example, in maize, there are 
two gene references: a full gene set as well as a filtered gene set containing only high 
confidence genes. Perhaps you might be interested in what a co-expression network might look
like if it was built using genes in the full vs filtered gene set. First, you'd build the two
RefGen objects from the raw datasets. This requires naming the dataset as they need to be 
referenced when the networks are built later. If we named these datasets 'ZmFull' and 'ZmFiltered', 
we could build two networks, one from each RefGen dataset.

Another case where naming datasets can come in handy is when there are quality control steps in
the build process. Building a network requires several parameters that control quality control
including:filtering missing data in genes/accession, setting minimum threshold for gene expression,
and specifying a minimal level of variance in gene expression. While the parameters used to build
the datasets are stored internally in the Camoco dataset and are retrievable, it might be useful
to name your datasets to reflect how they were created. For instance: 'ZmStrict' and 'ZmLenient' 

