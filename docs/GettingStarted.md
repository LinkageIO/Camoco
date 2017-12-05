# Getting started


## Running the command line

Follow the installation instructions. Make sure that Camoco is installed by activating the virtual env
then running the camoco CLI command:

```
$ source activate camoco
$ camoco
```

You should see:
```shell
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
    rm                remove camoco dataset.
    plotGWAS          Plot a GWAS dataset
    overlap           Calculate network overlap among GWAS results. See
                      --method for details.
    health            Generate COB Health Stats
    crossref          cross reference a set of networks based on edge strength
                      and go term density
    simulateGWAS      Evaluate co-expression in simulated GWAS datasets
    cistrans          compare the density of interactions between cis
                      interactionsto trans interactions
    snp2gene          Generate candidate genes and accompanying information
                      from GWAS SNPs
    neighbors         Generate significant gene neighborsfrom largest to
                      smallest Z-score

version: 0.5.0
src:/home/rob/Codes/Camoco/camoco/__init__.py
Cache. Money. Corn.
```

Download some data:
```shell
# A Gene Reference Genome
wget https://s3.msi.umn.edu/schae234/LocusPocus/ZmB73_5b_FGS.gff.gz

# A Gene Expression Matrix
wget https://s3.msi.umn.edu/schae234/Camoco/MaizeRNASeqTissue.tsv.bz2
```

You could use the CLI to create a reference genome object, but lets do things interactively from Ipython:
```
$ ipython
```
```python
import camoco as co
import pandas as pd
```

Check out the function signature for importing a GFF file:
```python
co.RefGen.from_gff?

Signature: co.RefGen.from_gff(filename, name, description, build, organism, chrom_feature='chromosome', gene_feature='gene', ID_attr='ID', attr_split='=')
Docstring:
Imports RefGen object from a gff (General Feature Format) file.
See more about the format here:
http://www.ensembl.org/info/website/upload/gff.html

Parameters
----------

filename : str
    The path to the GFF file.
name : str
    The name if the RefGen object to be stored in the core
    camoco database.
description : str
    A short description of the RefGen for future reference
build : str
    A string designating the genome build, used for comparison
    operations, genes may share IDS but are different across build.
organism : str
    A short string describing the organims this RefGen is coming
    from. Again, is used in comparing equality among genes which
    may have the same id or name.
chrom_feature : str (default: chromosome)
    The name of the feature (in column 3) that designates a
    a chromosome.
gene_feature : str (default: gene)
    The name of the feature (in column 2) that designates a
    gene. These features will be the main object that the RefGen
    encompasses.
ID_attr : str (default: ID)
    The key in the attribute column which designates the ID or
    name of the feature.
attr_split : str (default: '=')
    The delimiter for keys and values in the attribute column
File:      ~/Codes/Camoco/camoco/RefGen.py
Type:      method

```

## Create a ReferenceGenome object

```python
co.RefGen.from_gff("./ZmB73_5b_FGS.gff.gz","CornGenes",'Test!','ZeaMays')
```

Bind the object to a variable:

```python
Zm5bFGS = co.RefGen("CornGenes")
# Get a random gene
Zm5bFGS.random_gene()
<GRMZM5G849243>9:41675156-41675491+0(0)


# get 10 random genes
Zm5bFGS.random_genes(n=10)

{<GRMZM2G030016>1:46126597-46129448+0(0),
 <GRMZM2G101140>1:220827363-220828136+0(0),
 <GRMZM2G140594>1:247621972-247626310+0(0),
 <GRMZM5G898755>10:4243325-4244555+0(0),
 <AC205405.3_FG001>2:158945280-158945702+0(0),
 <GRMZM2G703265>2:178811871-178813335+0(0),
 <GRMZM2G350662>4:233270559-233271155+0(0),
 <GRMZM2G178074>5:214716081-214717618+0(0),
 <GRMZM2G582044>6:28074763-28074969+0(0),
 <GRMZM2G076169>8:97992609-97995058+0(0)}

```

## Create a COB (co-expresssion browser) object
```python
co.COB.from_table("./MaizeRNASeqTissue.tsv.bz2",name='MaizeTissue',description='Test!',refgen=Zm5bFGS,rawtype='RNASEQ')
```

Bind the object to a variable
```pythonpython
co.COB.from_table("./MaizeRNASeqTissue.tsv.bz2",name='MaizeTissue',description='Test!',refgen=Zm5bFGS,rawtype='RNASEQ')

[LOG] Tue Dec  5 10:55:23 2017 - Loading Expr table
[LOG] Tue Dec  5 10:55:23 2017 - Building Expr Index
[LOG] Tue Dec  5 10:55:23 2017 - Loading RefGen
[LOG] Tue Dec  5 10:55:23 2017 - RefGen for MaizeTissue not set!
[LOG] Tue Dec  5 10:55:23 2017 - Loading Coex table
[LOG] Tue Dec  5 10:55:23 2017 - MaizeTissue is empty
[LOG] Tue Dec  5 10:55:23 2017 - Loading Global Degree
[LOG] Tue Dec  5 10:55:23 2017 - MaizeTissue is empty
[LOG] Tue Dec  5 10:55:23 2017 - Loading Clusters
[LOG] Tue Dec  5 10:55:23 2017 - Clusters not loaded for: MaizeTissue ()
[LOG] Tue Dec  5 10:55:23 2017 - Resetting raw expression data
[LOG] Tue Dec  5 10:55:23 2017 - Resetting expression data
[LOG] Tue Dec  5 10:55:23 2017 - Extracting raw expression values
[LOG] Tue Dec  5 10:55:23 2017 - Importing Raw Expression Values
[LOG] Tue Dec  5 10:55:23 2017 - Trans. Log: raw->RawRNASEQ
[LOG] Tue Dec  5 10:55:23 2017 - Resetting expression data
[LOG] Tue Dec  5 10:55:23 2017 - Extracting raw expression values
[LOG] Tue Dec  5 10:55:23 2017 - Performing Quality Control on genes
[LOG] Tue Dec  5 10:55:23 2017 - ------------Quality Control
[LOG] Tue Dec  5 10:55:23 2017 - Raw Starting set: 39456 genes 18 accessions
[LOG] Tue Dec  5 10:55:23 2017 - Found out 3050 genes not in Reference Genome: ZeaMays - 5bFGS - CornGenes
[LOG] Tue Dec  5 10:55:23 2017 - Filtering expression values lower than 0.01
[LOG] Tue Dec  5 10:55:27 2017 - Found 22752 genes with > 0.2 missing data
[LOG] Tue Dec  5 10:55:33 2017 - Found 11295 genes which do not have one sample above 5
[LOG] Tue Dec  5 10:55:34 2017 - Found 0 accessions with > 0.3 missing data
[LOG] Tue Dec  5 10:55:34 2017 - Genes passing QC:
has_id                 39456
pass_membership        36406
pass_missing_data      16704
pass_min_expression    28161
PASS_ALL               16074
dtype: int64
[LOG] Tue Dec  5 10:55:34 2017 - Accessions passing QC:
has_id               18
pass_missing_data    18
PASS_ALL             18
dtype: int64
[LOG] Tue Dec  5 10:55:35 2017 - Genes passing QC by chromosome:
         has_id  pass_membership  pass_missing_data  pass_min_expression  PASS_ALL
chrom                                                                             
1          5574             5574               2576                 4189      2550
10         2506             2506               1087                 1857      1076
2          4445             4445               1912                 3280      1895
3          3837             3837               1687                 2886      1673
4          3866             3866               1625                 2799      1609
5          4165             4165               2032                 3180      2010
6          3034             3034               1358                 2254      1340
7          2919             2919               1280                 2196      1263
8          3245             3245               1499                 2432      1476
9          2768             2768               1191                 2003      1176
UNKNOWN      47               47                  6                   25         6
[LOG] Tue Dec  5 10:55:35 2017 - Kept: 16074 genes 18 accessions
[LOG] Tue Dec  5 10:55:35 2017 - Trans. Log: raw->quality_control
[LOG] Tue Dec  5 10:55:35 2017 - Performing Raw Expression Normalization
[LOG] Tue Dec  5 10:55:35 2017 - ------------ Normalizing
[LOG] Tue Dec  5 10:55:35 2017 - Trans. Log: raw->quality_control->arcsinh
[LOG] Tue Dec  5 10:55:35 2017 - Filtering refgen: CornGenes
[LOG] Tue Dec  5 10:55:35 2017 - Building Indices
[LOG] Tue Dec  5 10:55:35 2017 - Building Indices
[LOG] Tue Dec  5 10:55:35 2017 - Adding 16074 Genes info to database
[LOG] Tue Dec  5 10:55:36 2017 - Adding Gene attr info to database
[LOG] Tue Dec  5 10:55:36 2017 - Building Indices
[LOG] Tue Dec  5 10:55:36 2017 - Calculating Coexpression
[LOG] Tue Dec  5 10:56:31 2017 - Applying Fisher Transform
[LOG] Tue Dec  5 10:56:35 2017 - Calculating Mean and STD
[LOG] Tue Dec  5 10:56:37 2017 - Finding adjusted scores
[LOG] Tue Dec  5 10:56:37 2017 - Build the dataframe and set the significance threshold
[LOG] Tue Dec  5 10:56:38 2017 - Calculating Gene Distance
Calculating for 16074 genes
[LOG] Tue Dec  5 10:56:52 2017 - Done
[LOG] Tue Dec  5 10:56:52 2017 - Building Degree
[LOG] Tue Dec  5 10:56:52 2017 - Calculating Gene degree
[LOG] Tue Dec  5 10:56:54 2017 - Calculating Leaves using single
[LOG] Tue Dec  5 10:57:00 2017 - Finding the leaves
[LOG] Tue Dec  5 10:57:00 2017 - Pulling the scores for the .dat
[LOG] Tue Dec  5 10:57:01 2017 - Finding the IDs
[LOG] Tue Dec  5 10:57:01 2017 - Writing the .dat
[LOG] Tue Dec  5 10:57:02 2017 - Done
[LOG] Tue Dec  5 10:57:02 2017 - running MCL: mcl /home/rob/.camoco/tmp/tmpct2m_ua2 --abc -scheme 7 -I 2.0 -o -
[LOG] Tue Dec  5 10:57:02 2017 - waiting for MCL to finish...
.........[mcl] new tab created
[mcl] pid 8729
 ite -------------------  chaos  time hom(avg,lo,hi) m-ie m-ex i-ex fmv
  1  ...................  85.10  0.53 1.09/0.02/7.52 8.93 6.87 6.87  40
  2  ................... 186.64  3.71 0.73/0.08/3.34 4.96 1.24 8.52  77
  3  ................... 166.84  4.15 0.56/0.02/3.88 4.56 0.61 5.17  86
  4  ...................  95.42  1.77 0.62/0.02/3.76 4.31 0.41 2.14  82
  5  ...................  37.74  0.42 0.77/0.04/3.20 2.82 0.25 0.54  57
  6  ...................   9.77  0.06 0.93/0.21/1.54 1.89 0.12 0.07  26
  7  ...................   1.93  0.00 0.97/0.45/1.42 1.22 0.48 0.03   0
  8  ...................   1.22  0.00 0.98/0.51/1.39 1.03 0.80 0.03   0
  9  ...................   0.63  0.00 0.99/0.51/1.32 1.01 0.90 0.02   0
 10  ...................   0.38  0.00 1.00/0.56/1.00 1.00 0.95 0.02   0
 11  ...................   0.25  0.00 1.00/0.76/1.00 1.00 0.97 0.02   0
 12  ...................   0.25  0.00 1.00/0.76/1.00 1.00 0.99 0.02   0
 13  ...................   0.25  0.00 1.00/0.76/1.00 1.00 1.00 0.02   0
 14  ...................   0.25  0.00 1.00/0.76/1.00 1.00 1.00 0.02   0
 15  ...................   0.25  0.00 1.00/0.76/1.00 1.00 1.00 0.02   0
 16  ...................   0.14  0.00 1.00/0.86/1.00 1.00 1.00 0.02   0
 17  ...................   0.12  0.00 1.00/0.95/1.00 1.00 1.00 0.02   0
 18  ...................   0.19  0.00 1.00/0.86/1.00 1.00 1.00 0.02   0
 19  ...................   0.25  0.00 1.00/0.76/1.00 1.00 1.00 0.02   0
 20  ...................   0.17  0.00 1.00/0.84/1.00 1.00 1.00 0.02   0
 21  ...................   0.03  0.00 1.00/0.97/1.00 1.00 1.00 0.02   0
 22  ...................   0.00  0.00 1.00/1.00/1.00 1.00 1.00 0.02   0
 23  ...................   0.00  0.00 1.00/1.00/1.00 1.00 1.00 0.02   0
[mcl] jury pruning marks: <97,96,97>, out of 100
[mcl] jury pruning synopsis: <96.8 or sensational> (cf -scheme, -do log)
[mcl] output is in -
[mcl] 1502 clusters found
[mcl] output is in -

Please cite:
    Stijn van Dongen, Graph Clustering by Flow Simulation.  PhD thesis,
    University of Utrecht, May 2000.
       (  http://www.library.uu.nl/digiarchief/dip/diss/1895620/full.pdf
       or  http://micans.org/mcl/lit/svdthesis.pdf.gz)
OR
    Stijn van Dongen, A cluster algorithm for graphs. Technical
    Report INS-R0010, National Research Institute for Mathematics
    and Computer Science in the Netherlands, Amsterdam, May 2000.
       (  http://www.cwi.nl/ftp/CWIreports/INS/INS-R0010.ps.Z
       or  http://micans.org/mcl/lit/INS-R0010.ps.Z)

[LOG] Tue Dec  5 10:57:13 2017 - MCL done, Reading results.
[LOG] Tue Dec  5 10:57:13 2017 - Building cluster dataframe
[LOG] Tue Dec  5 10:57:13 2017 - Creating Cluster Ontology
[LOG] Tue Dec  5 10:57:14 2017 - Adding 1501 terms to the database.
[LOG] Tue Dec  5 10:57:14 2017 - Building the indices.
[LOG] Tue Dec  5 10:57:14 2017 - Your gene ontology is built.
[LOG] Tue Dec  5 10:57:14 2017 - Finished finding clusters
Out[25]: <COB: MaizeTissue>

```

## Object work together!
Get the expression among 10 random genes

```python
MaizeTissue = co.COB('MaizeTissue')

MaizeTissue.expr(Zm5bFGS.random_genes(n=10))
```
