"""

Camoco Library - CoAnalysis of Molecular Components

CacheMoneyCorn

"""

__license__ = """

The "MIT" License

Copyright (c) 2017-2019 Robert Schaefer

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.

"""

__version__ = '0.6.3'

import sys
import os
import numpy

import pyximport
pyximport.install(setup_args={
    "include_dirs":numpy.get_include() 
})

import matplotlib


from .Config import cf
from .Camoco import Camoco
from .Expr import Expr
from .COB import COB
from .RefGen import RefGen

from .Ontology import Ontology,Term
from .GWAS import GWAS
from .Locus import Locus
from .GOnt import GOnt
from .Overlap import Overlap

# Create yourself
Camoco.create('Camoco','Mother Database')
