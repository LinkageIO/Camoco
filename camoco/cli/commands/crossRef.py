#!/usr/bin/env python3

import camoco as co
import matplotlib.pylab as plt

def crossref(args):
    cobs = [co.COB(x) for x in args.cobs]
