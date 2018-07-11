Camoco
======
[![Build Status](https://travis-ci.org/LinkageIO/Camoco.svg?branch=master)](https://travis-ci.org/LinkageIO/Camoco)
[![Coverage Status](https://coveralls.io/repos/github/schae234/Camoco/badge.svg?branch=master)](https://coveralls.io/github/schae234/Camoco?branch=master)
[![DOI](https://zenodo.org/badge/15055703.svg)](https://zenodo.org/badge/latestdoi/15055703)

Co-Analysis of Molecular Components
-----------------------------------
![Camoco Workflow](https://s3-us-west-2.amazonaws.com/camoco/CamocoWorkflow.png)
Camoco is a python library for building and analyzing co-expression networks.
Currently, networks are generally built from gene expression data, but given
that the input to Camoco is simply a gene-by-sample expression matrix, there is
not reason that the analysis couldn't include things like protein abundance or
metabolites. Hence the name: co-analysis of molecular components. Very briefly,
Camoco creates co-expression networks using table formatted expression data and
and a few common genome data files.

Installation
------------
#### Virtual Environment
The easiest way to install Camoco is to use the included installation script.
This script creates an anaconda based virtual environment which takes care of
all required python packages as well as several of the external binaries and 
programs used by Camoco. Installation should be as easy as:

```
# Download and run install script
git clone https://github.com/schae234/Camoco.git
cd Camoco
# Follow instructions on screen
./install.sh
# Activate the virtual environment
source activate camoco
# Use CLI or run script
camoco --help
# deactivate virtual environment
source deactivate
```

You will need to add a few lines to your .bashrc in order for the conda
environment to be available from your shell.

e.g.:
```
# Assuming your installed camoco to ~/.camoco
export LD_LIBRARY_PATH=~/.camoco/lib/:$LD_LIBRARY_PATH
export PATH=$BASE/bin:~/.camoco/conda/bin/:$PATH
```

The installation script accepts simple arguments: `-b (default: ~/.camoco)`:
the base directory to install camoco including storage for databases

Required Files:
+ FPKM (or equivalent) CSV/TSV file
+ GFF File

Optional Files:
+ Gene Ontology (.obo and gene mapping)

Once co-expression networks are built, you can interact with the data using
the included command line:

```
# If you are using the conda virtual environment
source activate camoco
# Run the command line interface to get an idea of how Camoco works
camoco --help
```
See more examples below in the CLI section.

Camoco is built almost entirely in Python and was designed to be modular and 
object based. This means that Camoco can also be used interactively from the
shell using an interactive python session such as iPython. Just import the
package!
```
import camoco as co
help(co.COB)
```
Object classes are well documented for inputs and outputs.

CLI
---
Once installed, the `camoco` command will be available through the terminal.
See `camoco --help` for options!


Tests
-----
Unit tests within Camoco are implemented using the pytest framework. Unit tests
as well as Build tests can be found in the test directory.

Tests can be run by executing py.test in the tests dir:

```
source activate camoco
cd tests
py.test
```

Code Of Conduct
---------------
We expect users to be excellent to each other as well as to provide supportive
and engaging conversations regardless of each others backgrounds or identities.
If you'd like to know specifics, Camoco shares the same code of conduct as the
[Mozilla Science Lab](https://science.mozilla.org/code-of-conduct). Follow the 
link to learn more.

**CacheMoneyCorn**
