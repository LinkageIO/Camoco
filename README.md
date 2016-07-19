Camoco
======
[![Build Status](https://travis-ci.org/monprin/Camoco.svg?branch=master)](https://travis-ci.org/monprin/Camoco)

Co-Analysis of Molecular Components
-----------------------------------
Camoco is a python library for building and analyzing co-expression networks.
Currently, networks are generally built from gene expression data, but given
that the input to Camoco is simply a gene-by-sample expression matrix, there is
not reason that the analysis couldn't include things like protein abundance or
metabolites. Hence the name: co-analysis of molecular components. Very breifly,
Camoco creates co-expression networks using table formatted expression data and
and a few common genome data files:

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

Installation
------------
#### Docker
The easiest way to install Camoco is to use the included DockerFile. Download and
install [Docker](https://www.docker.com) on Linux, windows or Mac OS X. Using
docker will ensure that Camoco is run in an environment closest to what it was
developed in. After you have installed and tested docker, you should be able to
build and run camoco:


```
# Download the repo
git clone https://github.com/schae234/Camoco.git
# Build the docker file (grab some coffee)
cd Camoco
docker build -t camoco .
# Run the docker image interactively
docker run -it camoco
# Inside docker image
souce activate camoco
camoco --help
source deactivate
# Exit the docker image
exit
``` 

#### Virtual Environment
Alternatively, you can install Camoco using the included installation script.
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

Tagged releases are available on PyPi: `pip install camoco`.

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

Documentation
-------------
We acknowledge a lack of formal manual. We are writing this. However, the
function definitions themselves are well commented throughout the code base.
Questions and concerns about usage issued through github will be addressed!
Contact us!

CacheMoneyCorn
--------------
