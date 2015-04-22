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
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std

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

def plot_flanking_vs_inter(cob):
    import numpy as np
    from scipy import stats
    import statsmodels.api as sm
    import matplotlib.pyplot as plt
    from statsmodels.distributions.mixture_rvs import mixture_rvs
    log('Getting genes')
    genes = sorted(list(cob.refgen.iter_genes()))
    flanking = np.array([cob.coexpression(genes[i],genes[i-1]).score for i in  range(1,len(genes))])
    inter = cob.coex[~np.isfinite(cob.coex.distance)].score.values
    log('Getting flanking KDE')
    # get the KDEs
    flanking_kde = sm.nonparametric.KDEUnivariate(flanking)
    flanking_kde.fit()
    log('Getting Inter KDE')
    inter_kde = sm.nonparametric.KDEUnivariate(inter)
    inter_kde.fit()
    log('Plotting')
    plt.clf()
    fig = plt.figure(figsize=(8,4))
    fig.hold(True)
    ax = fig.add_subplot(1,1,1)
    ax.set_xlim([-4,4])
    ax.set_ylim([0,0.5])
    ax.plot(flanking_kde.support,flanking_kde.density,lw=2,color='black',alpha=1)
    ax.fill(flanking_kde.support,flanking_kde.density,color='red',alpha=0.3,label='Cis Interactions')
    ax.scatter(np.median(flanking),-0.05,marker='D',color='red')
    ax.set_xlim([-4,4])
    ax.set_ylim([0,0.5])
    ax.plot(inter_kde.support,inter_kde.density,lw=2,color='black')
    ax.fill(inter_kde.support,inter_kde.density,color='blue',alpha=0.3,label='Trans Interactions')
    ax.scatter(np.median(inter),-0.05,marker='D',color='blue')
    ax.set_xlabel('CoExpression Interaction (Z-Score)')
    ax.set_ylabel('Distribution Density')
    fig.tight_layout()
    fig.savefig("{}_flank_inter.png".format(cob.name))


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

