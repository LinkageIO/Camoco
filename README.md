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
Camoco was developed using python >=3.6 and requires numpy to be installed.
Otherwise install the latest tagged version with pip:
```
$ pip install numpy
$ pip install camoco
```
Or install the latest developmental version by cloning the github repository:
```
$ pip install numpy
$ git clone https://github.com/LinkageIO/Camoco.git
$ cd Camoco
$ python setup.py install
```

CLI
---
Once installed, the `camoco` command will be available through the terminal.
See `camoco --help` for options!


Code Of Conduct
---------------
We expect users to be excellent to each other as well as to provide supportive
and engaging conversations regardless of each others backgrounds or identities.
If you'd like to know specifics, Camoco shares the same code of conduct as the
[Mozilla Science Lab](https://science.mozilla.org/code-of-conduct). Follow the 
link to learn more.

**CacheMoneyCorn**
