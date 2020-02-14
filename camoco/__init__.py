"""

    Camoco Library - CoAnalysis of Molecular Components

"""

__license__ = """

The "MIT" License

Copyright (c) 2017-2020 Robert Schaefer

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

__version__ = "1.0.0"

import sys
import os
import numpy
import logging

log = logging.getLogger('camoco')
# One of DEBUG, INFO, WARNING, ERROR, CRITICAL
log.setLevel(logging.INFO)
# Set up the console handler
ch = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s - %(message)s',datefmt='%d-%b-%y %H:%M:%S')
ch.setFormatter(formatter)
ch.setLevel(logging.INFO)
log.addHandler(ch)

from .coex import Coex
