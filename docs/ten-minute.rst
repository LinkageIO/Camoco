
.. _ten-minute: 

Ten Minute Tutorial
###################
This is a whirlwind tutorial to get you started with Camoco. This tutorial assumes that 
you've successfully completed the :ref:`installation steps <installation>`.


Running the CLI
---------------
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

This is the CLI help screen.


