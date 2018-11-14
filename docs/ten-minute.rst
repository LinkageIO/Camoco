
.. _ten-minute: 

Tutorial
########
This is a whirlwind tutorial to get you started with Camoco. This tutorial assumes that 
you've successfully completed the :ref:`installation steps <installation>`. In this 
tutorial we will be recreating the datasets used in `this publication 
<https://doi.org/10.1105/tpc.18.00299>`_ which examined the co-expression among
genes related to the maize kernel ionome. 

First, we will need to download some data.

.. note:: 
  Data distributed here are subject to licensing terms specified by the original 
  publications. If used, please give proper attribution to both the source article
  (if raw data is used) as well as the `Camoco publication <https://doi.org/10.1105/tpc.18.00299>`_
  if Camoco was utilized.


Getting Data
============
Included with the source code on GitHub are test data sets as well as some raw input that 
can be used to create some Camoco objects. Here is a breakdown of some of the files:

RefGen Data
-----------

`ZmB73_5b_FGS.gff.gz <https://github.com/LinkageIO/Camoco/raw/master/tests/raw/RefGen/ZmB73_5b_FGS.gff.gz>`_
  This is the raw data needed to create the RefGen object. The format of this file is 
  a GFF, which you can read more about `here <https://uswest.ensembl.org/info/website/upload/gff.html>`_.


Network (COB) Data
------------------

`Hirsch2014_PANGenomeFPKM.txt.gz <https://github.com/LinkageIO/Camoco/raw/master/tests/raw/Expr/RNASEQ/Hirsch2014_PANGenomeFPKM.txt.gz>`_
  This is a FPKM table for 503 seedlings used to create the ZmPAN network from `Schaefer et al. <https://doi.org/10.1105/tpc.18.00299>`__ 
  using data originally generated from `Hirsch et al. <https://doi.org/10.1105/tpc.113.119982>`_


`Stelpflug2018_B73_Tissue_Atlas.txt.gz <https://github.com/LinkageIO/Camoco/raw/master/tests/raw/Expr/RNASEQ/Stelpflug2018_B73_Tissue_Atlas.txt.gz>`_
  This is a FPKM tables from 76 different tissues/time-points in the maize accession B73 described in
  `Stelpflug et al. <https://doi.org/10.3835/plantgenome2015.04.0025>`_. This is the basis of the ZmSAM
  network in `Schaefer et al. <https://doi.org/10.1105/tpc.18.00299>`_.

`Schaefer2018_ROOTFPKM.tsv.gz <https://github.com/LinkageIO/Camoco/raw/master/tests/raw/Expr/RNASEQ/Schaefer2018_ROOTFPKM.tsv.gz>`_
  This dataset contains FPKM values for 46 genotypically diverse maize root samples. It is the
  described as the ZmRoot network in `Schaefer et al. <https://doi.org/10.1105/tpc.18.00299>`_.


Ontology Data
-------------

`go.obo.gz <https://github.com/LinkageIO/Camoco/raw/master/tests/raw/GOnt/go.obo.gz>`_
  This contains a list of all Gene Ontology (GO) Terms. The .obo file contains the generic
  GO structure and the relationships between terms. A description of the data can
  be found `here <http://www.geneontology.org/page/download-ontology>`_.

`zm_go.tsv.gz <https://github.com/LinkageIO/Camoco/raw/master/tests/raw/GOnt/zm_go.tsv.gz>`_
  This file contains the mapping between maize genes and GO terms. 

GWAS Data
---------

`ZmIonome.allLocs.csv.gz <https://github.com/LinkageIO/Camoco/raw/master/tests/raw/GWAS/SchaeferPlantCell/ZmIonome.allLocs.csv.gz>`_
  This file contains GWAS for a study looking at nutritional quality (elemental composition) for
  17 traits in a maize mapping population (NAM) commonly referred to as the ionome. 

`Wallace_etal_2014_PlosGenet_GWAS.gz <https://github.com/LinkageIO/Camoco/raw/master/tests/raw/GWAS/WallacePLoSGenet/Wallace_etal_2014_PLoSGenet_GWAS_hits-150112.txt.gz>`_
  GWAS results from another maize mapping population (NAM) for 41 different traits. Details
  can be found in `Wallace et al. <https://doi.org/10.1371/journal.pgen.1004845>`_


To download all the data, you can use the following command: 

.. code:: bash

  wget \
   https://github.com/LinkageIO/Camoco/raw/master/tests/raw/RefGen/ZmB73_5b_FGS.gff.gz \
   https://github.com/LinkageIO/Camoco/raw/master/tests/raw/Expr/RNASEQ/Hirsch2014_PANGenomeFPKM.txt.gz \
   https://github.com/LinkageIO/Camoco/raw/master/tests/raw/Expr/RNASEQ/Stelpflug2018_B73_Tissue_Atlas.txt.gz \
   https://github.com/LinkageIO/Camoco/raw/master/tests/raw/Expr/RNASEQ/Schaefer2018_ROOTFPKM.tsv.gz \
   https://github.com/LinkageIO/Camoco/raw/master/tests/raw/GOnt/go.obo.gz \
   https://github.com/LinkageIO/Camoco/raw/master/tests/raw/GOnt/zm_go.tsv.gz \
   https://github.com/LinkageIO/Camoco/raw/master/tests/raw/GWAS/SchaeferPlantCell/ZmIonome.allLocs.csv.gz \
   https://github.com/LinkageIO/Camoco/raw/master/tests/raw/GWAS/WallacePLoSGenet/Wallace_etal_2014_PLoSGenet_GWAS_hits-150112.txt.gz

.. note::
  Please always be cautious when pasting shell commands from the internet! 

Now uncompress the files, so we can easily look inside them ::

  $ gunzip ./*.gz

Running the CLI
===============
The easiest way to get started with Camoco is to use the command line interface. This
can be accessed using the :code:`camoco` command from the shell::

  $ camoco

You should see the following output: ::

  $ camoco
  usage: camoco [-h] [--debug] [--interactive] [--force] Available Commands ...

        ___           ___           ___           ___           ___           ___      
       /  /\         /  /\         /__/\         /  /\         /  /\         /  /\     
      /  /:/        /  /::\       |  |::\       /  /::\       /  /:/        /  /::\    
     /  /:/        /  /:/\:\      |  |:|:\     /  /:/\:\     /  /:/        /  /:/\:\   
    /  /:/  ___   /  /:/~/::\   __|__|:|\:\   /  /:/  \:\   /  /:/  ___   /  /:/  \:\  
   /__/:/  /  /\ /__/:/ /:/\:\ /__/::::| \:\ /__/:/ \__\:\ /__/:/  /  /\ /__/:/ \__\:\ 
   \  \:\ /  /:/ \  \:\/:/__\/ \  \:\~~\__\/ \  \:\ /  /:/ \  \:\ /  /:/ \  \:\ /  /:/ 
    \  \:\  /:/   \  \::/       \  \:\        \  \:\  /:/   \  \:\  /:/   \  \:\  /:/  
     \  \:\/:/     \  \:\        \  \:\        \  \:\/:/     \  \:\/:/     \  \:\/:/   
      \  \::/       \  \:\        \  \:\        \  \::/       \  \::/       \  \::/    
       \__\/         \__\/         \__\/         \__\/         \__\/         \__\/ 

  Camoco (Co-analysis of Molecular Components) inter-relates and co-analyzes different 
  levels of genomic data. Namely it integrates genes present near and around GWAS loci
  using unbiased, functional information derived from co-expression networks.

  optional arguments:
    -h, --help          show this help message and exit
    --debug             Drop into ipdb when something bad happens.
    --interactive       Initiate an ipdb session right before exiting.
    --force             Overwrite output files from previous analyses.

  Camoco CLI program:
    Use --help with each command for more info

    Available Commands
      help              Prints this help message
      build-gwas        build a GWAS dataset
      build-go          Build a Gene Ontology (GO)
      build-refgen      Build a Reference Genome.
      build-cob         Build a Co-expression network.
      list (ls)         List camoco datasets.
      rm                Remove camoco dataset.
      overlap           Calculate network overlap among GWAS results. See
                        --method for details.
      health            Generate network health statistics
      snp2gene          Generate candidate genes and accompanying information
                        from GWAS SNPs
      neighbors         Generate significant gene neighbors from largest to
                        smallest Z-score

  version: 0.6.1
  src:/home/schae234/miniconda3/envs/camoco/lib/python3.6/site-packages/camoco/__init__.py
  Cache. Money. Corn. 


Building Camoco Objects
=======================

The first camoco object we are going to build is the RefGen object. This is needed because 
most of the other objects need a reference in order to properly interpret gene IDs. For example,
if you look at the first few lines of the file go gene mapping ::

  $ head zm_go.tsv.gz



