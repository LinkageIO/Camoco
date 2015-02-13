import os
import sys
import time
import re
import functools

from termcolor import colored,cprint
from itertools import chain
from camoco.Locus import *

import camoco as co
import matplotlib.pylab as pylab
import numpy as np
import pandas as pd

def memoize(obj):
    cache = obj.cache = {}
    @functools.wraps(obj)
    def memoizer(*args, **kwargs):
        # Give us a way to clear the cache
        if 'clear_cache' in kwargs:
            cache.clear()
        # This wraps the calling of the memoized object
        key = str(args) + str(kwargs)
        if key not in cache:
            cache[key] = obj(*args, **kwargs)
        return cache[key]
    return memoizer


class log(object):
    def __init__(self,msg=None,*args,color='green'):
        if msg is not None:
            print(colored(" ".join(["[LOG]",time.ctime(), '-', msg.format(*args)]),color=color),file=sys.stderr)
    @classmethod
    def warn(cls,msg,*args):
        cls(msg,*args,color='red')

    def __call__(self,msg,*args,color='green'):
        print(colored(" ".join(["[LOG]",time.ctime(), '-', msg.format(*args)]),color=color),file=sys.stderr)

def ext(filename):
    return os.path.join(os.path.expanduser("~/MetaboloCOB/"+filename))

def B73_eq_Mo17(snp,HM):
    genotypes = HM.genotypes(snp,accessions=['B73','MO17'])
    if genotypes[0] == genotypes[1]:
        return True
    else:
        return False

def plot_local_global_degree(term,filename=None,bootstraps=1):
    ROOT = co.COB("ROOT")
    RZM = ROOT.refgen # use root specific for bootstraps
    hood = ROOT.neighborhood(term.flanking_genes(RZM))
    bshood = pd.concat([ROOT.neighborhood(term.bootstrap_flanking_genes(RZM)) for _ in range(0,bootstraps)])
    pylab.clf()
    pylab.scatter(bshood['local'],bshood['global'],alpha=0.05)
    pylab.scatter(hood['local'],hood['global'],c='r')
    pylab.xlabel('Local Degree')
    pylab.ylabel('Global Degree')
    pylab.title('{} Locality'.format(term.id))
    if filename is None:
        filename = "{}_locality.png".format(term.id)
    pylab.savefig(filename)

def plot_local_vs_cc(term,filename=None,bootstraps=1):
    RZM = co.COB('ROOT').refgen # use root specific for bootstraps
    pylab.clf()
    for _ in range(0,bootstraps):
        graph = co.COB('ROOT').graph(term.bootstrap_flanking_genes(RZM))
        degree = np.array(graph.degree())
        cc = np.array(graph.transitivity_local_undirected(weights='weight'))
        nan_mask = np.isnan(cc)
        pylab.scatter(degree[~nan_mask],cc[~nan_mask],alpha=0.05)
    # plot empirical
    graph = co.COB('ROOT').graph(term.flanking_genes(RZM))
    degree = np.array(graph.degree())
    cc = np.array(graph.transitivity_local_undirected(weights='weight'))
    nan_mask = np.isnan(cc)
    pylab.scatter(degree[~nan_mask],cc[~nan_mask])
    pylab.xlabel('Local Degree')
    pylab.ylabel('Clustering Coefficient')
    if filename is None:
        filename = "{}_cc.png".format(term.id)
    pylab.savefig(filename)
   
