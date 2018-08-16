import os
import sys
import time
import re
import functools
import glob
import shutil

from termcolor import colored, cprint
from itertools import chain
from collections import OrderedDict

from .Locus import Locus
from .Config import cf
from apsw import CantOpenError

import camoco as co

import matplotlib.pylab as pylab
import numpy as np
import pandas as pd
import statsmodels.api as sm

import gzip
import bz2

def mean_confidence_interval(data): # pragma no cover
    '''
        Convenience function to return both the mean as well as the 
        confidence interval on data. Good for chaining so inline data
        does not need to be evaluated twice (once for mean and another 
        for cf)
    '''
    return np.mean(data), confidence_interval(data)

def confidence_interval(data, confidence=0.95): # pragma no cover
    '''
        Returns the confidence (default 95%) for data.
    '''
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return 1.96*se

class NearestDict(OrderedDict): # pragma no cover
    '''
        This extension overrides the get item method 
        of dict where if a key does not exist, it returns
        the nearst key which does exist.
    '''
    def __getitem__(self,key): # pragma no cover
        'Returns the nearest key which exists'
        return dict.__getitem__(self,min(self.keys(),key=lambda x: abs(x-key)))



def available_datasets(type='%', name='%'): # pragma no cover
    try:
        cur = co.Camoco("Camoco", type='Camoco').db.cursor()
        datasets = cur.execute('''
            SELECT type, name, description, added
            FROM datasets 
            WHERE type LIKE ?
            AND name LIKE ?
            ORDER BY type;''', (type,name)).fetchall()
        if datasets:
            datasets = pd.DataFrame(
                datasets, 
                columns=["Type", "Name", "Description", "Date Added"],
            ).set_index(['Type'])
        else:
            datasets = pd.DataFrame(
                columns=["Type", "Name", "Description", "Date Added"]
            )
        # Check to see if we are looking for a specific dataset
        if '%' not in type and '%' not in name:
            return True if name in datasets['Name'].values else False
        else:
            return datasets
    except CantOpenError as e:
        raise e
        

def available(type=None,name=None): # pragma no cover
    # Laaaaaaaaazy
    return available_datasets(type=type,name=name)

def del_dataset(type, name, force=False): # pragma no cover
    try:
        c = co.Camoco("Camoco")
    except CantOpenError:
        return True
    if force == False:
        c.log("Are you sure you want to delete:\n {}.{}", type, name)
        if input("[Y/n]").upper() != 'Y':
            c.log("Nothing Deleted")
            return
    c.log("Deleting {}", name)
    try:
        c.db.cursor().execute('''
            DELETE FROM datasets 
            WHERE name LIKE '{}' 
            AND type LIKE '{}';'''.format(name, type)
        )
    except CantOpenError:
        pass
    try:
        dfiles = glob.glob(
            os.path.join(
                cf.options.basedir,
                'databases',
                '{}.{}.*'.format(type,name)
            )
        )
        for f in dfiles:
            c.log('Removing {}',f)
            try:
                os.remove(f)
            except IsADirectoryError:
                shutil.rmtree(f)
    except FileNotFoundError as e:
        pass
    if type == 'Expr':
        # also have to remove the COB specific refgen
        del_dataset('RefGen', 'Filtered'+name, force=force)
        del_dataset('Ontology', name+'MCL', force=force)
    return True

def mv_dataset(type,name,new_name): # pragma no cover
    c = co.Camoco("Camoco")
    c.db.cursor().execute('''
        UPDATE datasets SET name = ? 
        WHERE name = ? AND 
        type = ?''',(new_name,name,type)
    )
    os.rename(
        c._resource('databases','.'.join([type,name])+".db"),
        c._resource('databases',".".join([type,new_name])+".db")
    )

class rawFile(object): # pragma no cover
    def __init__(self,filename): # pragma no cover
        self.filename = filename
        if filename.endswith('.gz'): # pragma no cover
            self.handle = gzip.open(filename,'rt')
        elif filename.endswith('bz2'): # pragma no cover
            self.handle = bz2.open(filename,'rt')
        else:
            self.handle = open(filename,'r')
    def __enter__(self): # pragma no cover
        return self.handle
    def __exit__(self,type,value,traceback): # pragma no cover
        self.handle.close()


def redescribe_dataset(type,name,new_desc): # pragma no cover
    c = co.Camoco("Camoco")
    c.db.cursor().execute('''
        UPDATE datasets SET description = ? 
        WHERE name = ? AND type = ?''',
        (new_desc,name,type)
    )

def memoize(obj): # pragma no cover
    cache = obj.cache = {}
    @functools.wraps(obj)
    def memoizer(*args, **kwargs): # pragma no cover
        # Give us a way to clear the cache
        if 'clear_cache' in kwargs:
            cache.clear()
        # This wraps the calling of the memoized object
        key = str(args) + str(kwargs)
        if key not in cache:
            cache[key] = obj(*args, **kwargs)
        return cache[key]
    return memoizer


class log(object): # pragma no cover
    def __init__(self, msg=None, *args, color='green'): # pragma no cover
        if msg is not None and cf.logging.log_level == 'verbose':
            print(
                colored(
                    " ".join(["[LOG]", time.ctime(), '-', msg.format(*args)]), 
                    color=color
                ), file=sys.stderr
            )

    @classmethod
    def warn(cls, msg, *args): # pragma no cover
        cls(msg, *args, color='red')

    def __call__(self, msg, *args, color='green'): # pragma no cover
        if cf.logging.log_level == 'verbose':
            print(
                colored(
                    " ".join(["[LOG]", time.ctime(), '-', msg.format(*args)]), 
                    color=color
                ),
            file=sys.stderr
        )


def plot_flanking_vs_inter(cob): # pragma no cover
    import numpy as np
    from scipy import stats
    import statsmodels.api as sm
    import matplotlib.pyplot as plt
    from statsmodels.distributions.mixture_rvs import mixture_rvs
    log('Getting genes')
    genes = sorted(list(cob.refgen.iter_genes()))
    flanking = np.array([cob.coexpression(genes[i], genes[i-1]).score for i in  range(1, len(genes))])
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
    fig = plt.figure(figsize=(8, 4))
    fig.hold(True)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlim([-4, 4])
    ax.set_ylim([0, 0.5])
    ax.plot(flanking_kde.support, flanking_kde.density, lw=2, color='black', alpha=1)
    ax.fill(flanking_kde.support, flanking_kde.density, color='red', alpha=0.3, label='Cis Interactions')
    ax.scatter(np.median(flanking), -0.05, marker='D', color='red')
    ax.set_xlim([-4, 4])
    ax.set_ylim([0, 0.5])
    ax.plot(inter_kde.support, inter_kde.density, lw=2, color='black')
    ax.fill(inter_kde.support, inter_kde.density, color='blue', alpha=0.3, label='Trans Interactions')
    ax.scatter(np.median(inter), -0.05, marker='D', color='blue')
    ax.set_xlabel('CoExpression Interaction (Z-Score)')
    ax.set_ylabel('Distribution Density')
    fig.tight_layout()
    fig.savefig("{}_flank_inter.png".format(cob.name))


def plot_local_global_degree(term, filename=None, bootstraps=1): # pragma no cover
    ROOT = co.COB("ROOT")
    RZM = ROOT.refgen # use root specific for bootstraps
    hood = ROOT.neighborhood(term.flanking_genes(RZM))
    bshood = pd.concat([ROOT.neighborhood(term.bootstrap_flanking_genes(RZM)) for _ in range(0, bootstraps)])
    pylab.clf()
    pylab.scatter(bshood['local'], bshood['global'], alpha=0.05)
    pylab.scatter(hood['local'], hood['global'], c='r')
    pylab.xlabel('Local Degree')
    pylab.ylabel('Global Degree')
    pylab.title('{} Locality'.format(term.id))
    if filename is None:
        filename = "{}_locality.png".format(term.id)
    pylab.savefig(filename)

def plot_local_vs_cc(term, filename=None, bootstraps=1): # pragma no cover
    RZM = co.COB('ROOT').refgen # use root specific for bootstraps
    pylab.clf()
    for _ in range(0, bootstraps): # pragma no cover
        graph = co.COB('ROOT').graph(term.bootstrap_flanking_genes(RZM))
        degree = np.array(graph.degree())
        cc = np.array(graph.transitivity_local_undirected(weights='weight'))
        nan_mask = np.isnan(cc)
        pylab.scatter(degree[~nan_mask], cc[~nan_mask], alpha=0.05)
    # plot empirical
    graph = COB('ROOT').graph(term.flanking_genes(RZM))
    degree = np.array(graph.degree())
    cc = np.array(graph.transitivity_local_undirected(weights='weight'))
    nan_mask = np.isnan(cc)
    pylab.scatter(degree[~nan_mask], cc[~nan_mask])
    pylab.xlabel('Local Degree')
    pylab.ylabel('Clustering Coefficient')
    if filename is None:
        filename = "{}_cc.png".format(term.id)
    pylab.savefig(filename)

def read_density(path): # pragma no cover
    dfs = []
    for x in glob.glob(path): # pragma no cover
        df = pd.read_table(x,sep=',')
        dfs.append(df)
    df = pd.concat(dfs)
    df.insert(4,'TraitType','Element')
    df.loc[[x.startswith('Log') for x in df.Term],'TraitType'] = 'Log'
    df.loc[[x.startswith('PCA') for x in df.Term],'TraitType'] = 'PCA'
    df.loc[[x.startswith('Trans') for x in df.Term],'TraitType'] = 'Trans'
    return df.set_index(['Ontology','COB','Term','WindowSize','FlankLimit'])

def read_FDR(glob_path,sep=','): # pragma no cover
    dfs = []
    for x in glob.glob(glob_path): # pragma no cover
        df = pd.read_table(x,sep=sep)
        dfs.append(df)
    df = pd.concat(dfs)
    # I guess we forgot this before
    df.insert(4,'TraitType','Element')
    df.loc[[x.startswith('Log') for x in df.Term],'TraitType'] = 'Log'
    df.loc[[x.startswith('PCA') for x in df.Term],'TraitType'] = 'PCA'
    df.loc[[x.startswith('Trans') for x in df.Term],'TraitType'] = 'Trans'
    return df.set_index(['Ontology','COB','Term','WindowSize','FlankLimit'])

def zmax(a): # pragma no cover
    if len(a) == 0:
        return 0
    else:
        return np.max(a)

def zmin(a): # pragma no cover
    if len(a) == 0:
        return 0
    else:
        return np.min(a)

def groupedFDR(df): # pragma no cover
    def grouped_agg(x): # pragma no cover
        return pd.Series(
            {
                'Tot':  sum(x.numReal),
                'FDR10'   :zmax(x[x.FDR<=0.10].numReal),
                'FDR10_Z' :zmin(x[x.FDR<=0.10].zscore),
                'FDR35'   :zmax(x[x.FDR<=0.35].numReal),
                'FDR35_Z' :zmin(x[x.FDR<=0.35].zscore),
                'FDR50'   :zmax(x[x.FDR<=0.50].numReal),
                'FDR50_Z' :zmin(x[x.FDR<=0.50].zscore)
            }
        )
    groups = ['Ontology','COB','WindowSize','FlankLimit','TraitType','Term']
    return df.reset_index().groupby(groups).apply(grouped_agg)


class DummyRefGen(object): # pragma no cover
    '''
        This is a dummy refgen that will always return True
    '''
    def __init__(self): # pragma no cover
        self.name = 'DummyRefGen'
    def __contains__(self,x): # pragma no cover
        return True
