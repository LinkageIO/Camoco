import sys
import os

import pandas as pd
import numpy as np
import scipy as sp
import camoco as co

import matplotlib.pylab as plt

try:
    from sklearn.neighbors import KernelDensity
except ImportError as e:
    raise ImportError('This command requires sklearn')

def cistrans(args):
    cob = co.COB(args.cob) 
    if args.out == None:
        args.out = '{}_cistrans'.format(cob.name)
    # np.newaxis adds an empty axis in that position of the slice
    # the sklearn module requires the values to be in the rows:
    # http://scikit-learn.org/stable/auto_examples/neighbors/plot_kde_1d.html
    cis = cob.coex \
            .score[cob.coex.distance <= args.cis_distance]\
            .values[:,np.newaxis]
    trans = cob.coex\
            .score[np.isinf(cob.coex.distance)]\
            .values[:,np.newaxis]
    X_plot = np.linspace(-10,10,1000)[:,np.newaxis]
    print(
            'Found {:,} cis interactions and {:,} trans interactions'.format(
        cis.shape[0],
        trans.shape[0]
    ))
    # Fit the kernel
    kd=KernelDensity(bandwidth=0.2)
    kd.fit(cis)
    cis_kde = np.exp(kd.score_samples(X_plot))
    plt.fill(X_plot,cis_kde,alpha=0.5,label='Cis Interactions')
    # Fit the trans 
    kd.fit(trans[0:50000])
    trans_kde = np.exp(kd.score_samples(X_plot))
    plt.fill(X_plot,trans_kde,alpha=0.5,label='Trans Interactions')
    plt.legend()
    plt.title('Cis vs Trans Density: {}'.format(cob.name))
    # Calculate the mann whitney U test
    u,pval = sp.stats.mannwhitneyu(cis[:,0],trans[:,0]) 
    print('P-val: {}'.format(pval))
    print('Figure saved: {}'.format(args.out+'.png'))
    plt.savefig(args.out+'.png')
