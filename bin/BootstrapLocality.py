#!/usr/bin/env python3

import argparse
import sys
import os

import numpy as np
import scipy as sp
import scipy.stats

import camoco as co
import pandas as pd
import matplotlib.pylab as plt
import statsmodels.api as sm

lowess = sm.nonparametric.lowess

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return m, h


class NearestDict(dict):
    '''
        This extension overrides the get item method 
        of dict where if a key does not exist, it returns
        the nearst key which does exist.
    '''
    def __getitem__(self,key):
        'Returns the nearest key which exists'
        return dict.__getitem__(self,min(self.keys(),key=lambda x: abs(x-key)))

def main(args):
    # Initiate for args 
    ont = co.Ontology(args.ontology)
    term = ont[args.term]
    cob = co.COB(args.cob)
    # Generate outpue dirs
    os.makedirs(os.path.join(args.cob,args.ontology,args.term),exist_ok=True)
    cwd = os.getcwd()
    os.chdir(os.path.join(args.cob,args.ontology,args.term))

    '''---------------------------------------------------
        Empirical and SD Calculations
    '''

    # Find the empirical Locality
    loc = cob.locality(cob.refgen.candidate_genes(term.effective_snps(window_size=50000),gene_limit=10))
    # Find the windows 
    window_ticks = max(loc['global'])/args.window_size
    loc['window'] = [int(x/window_ticks) for x in loc['global']]
    
    # Find the Bootstrapped Locality
    bsloc = pd.concat(
            [cob.locality(
                cob.refgen.bootstrap_candidate_genes(term.effective_snps(window_size=50000),gene_limit=10),
                bootstrap_name=x
            ) for x in range(args.num_bootstraps)]
    ).sort('global')
    # Remove global degree larger than empirical
    bsloc = bsloc[bsloc['global'] <= max(loc['global'])]
    
    # Assign into windows
    bsloc['window'] = [int(x/window_ticks) for x in bsloc['global']]

    # Calculate the s.d. of the residuals in each window
    win_std = NearestDict(
        bsloc.groupby('window').apply(lambda df: df['resid'].std()).to_dict()
    )
    # divide empirical residuals by the s.d. in their respective window
    loc['bs_std'] = [win_std[x['window']] for i,x in loc.iterrows()]
    loc['zscore'] = [x['resid']/win_std[x['window']] for i,x in loc.iterrows()]
    # Do the same for bs
    bsloc['bs_std'] = [win_std[x['window']] for i,x in bsloc.iterrows()]
    bsloc['zscore'] = [x['resid']/win_std[x['window']] for i,x in bsloc.iterrows()]

    '''---------------------------------------------------
        FDR Calculations
    '''
    
    # Repeat bootstraps to assess global FDR
    fdr = pd.concat(
            [cob.locality(
                cob.refgen.bootstrap_candidate_genes(term.effective_snps(window_size=50000),gene_limit=10),
                bootstrap_name=x
            ) for x in range(args.num_bootstraps)]
    ).sort('global')
    # Remove global degree larger than empirical
    fdr = fdr[fdr['global'] <= max(loc['global'])]
    fdr['window'] = [int(x/window_ticks) for x in fdr['global']]

    # calculate z-scores for the global 
    fdr['bs_std'] = [win_std[x['window']] for i,x in fdr.iterrows()]
    fdr['zscore'] = [x['resid']/win_std[x['window']] for i,x in fdr.iterrows()]

    '''---------------------------------------------------
        Plotting
    '''
    # Make plots for later
    fig,ax = plt.subplots(figsize=(8,6)) 
    fig.hold(True)
    # Y axis is local degree (what we are TRYING to predict)
    ax.set_ylim(0,max(loc['local']))
    ax.set_xlim(0,max(loc['global']))
    ax.set_xlabel('Number Global Interactions')
    ax.set_ylabel('Number Local Interactions')

    # UGH! map lowess 
    fdrlowess = lowess(fdr['local'],fdr['global'],frac=0.2,it=5,delta=0.1*len(bsloc),is_sorted=False)
    # plot the bootstrap points
    plt.plot(fdr['global'],fdr['local'],'ro',alpha=0.05,label='Bootstraps')
    # plot the OLS lowess line
    plt.plot(fdrlowess[:,0],fdrlowess[:,1],'g--',label='Bootstrap OLS')

    #plot the empirical data
    plt.plot(loc['global'],loc['local'],'bo',alpha=1,label='Empirical')
    plt.plot(loc['global'],loc['fitted'],'k:',label='Empirical OLS')

    # plot the std envelope
    stdlowess = lowess(
        fdr['fitted']+(args.sd_envelope*fdr['bs_std']),
        fdr['global'], frac=0.2, it=5, delta=0.1*len(fdr),
        is_sorted=False
    )
    ax.plot(
        stdlowess[:,0],stdlowess[:,1],'r--'
        ,label='{} s.d. envelope'.format(args.sd_envelope)
    )

    # Output the datafiles
    loc.sort('zscore',ascending=False).to_csv('Empirical.csv')
    bsloc.to_csv('ZscoreBootstraps.csv')
    fdr.sort('zscore',ascending=False).to_csv('FDRBootstraps.csv')
    # finish plots
    legend = ax.legend(loc='best')
    plt.savefig("{}.png".format(args.term)) 


    '''---------------------------------------------------
        FDR Plotting
    '''
    
    fig,ax = plt.subplots(figsize=(8,6)) 
    fig.hold(True)

    # Plot the empirical Z-score distributions
    zscores = list(range(0,8))
    zcdf = [sum(loc['zscore'] >= x) for x in zscores]
    ax.plot(zscores,zcdf,'bo',label='Empirical Zscore > x')
    # plot the fdr scores spread
    zcdf = pd.DataFrame(
        [ mean_confidence_interval(
            fdr.groupby('bootstrap_name').apply(
                lambda df: sum(df['zscore'] >= x )
            )
        ) for x in zscores ],
        columns=['mean','ci']
    )
    ax.errorbar(zscores,zcdf['mean'],yerr=zcdf['ci'],label='Bootstrap Z-scores',color='red')  
    ax.set_xlabel('Z-Score')
    ax.set_ylabel('Number of Genes > Z')
    legend = ax.legend(loc='best')
    plt.savefig('FDR_ZScore.png')

    # Do the same thing, but with residuals

    fig,ax = plt.subplots(figsize=(8,6)) 
    fig.hold(True)

    resids = list(range(0,int(max(loc['resid'])+2)))
    rcdf = [sum(loc['resid'] >= x) for x in resids]
    ax.plot(resids,rcdf,'bo',label='Empirical Residuals')
    rcdf = pd.DataFrame(
         [ mean_confidence_interval(
            fdr.groupby('bootstrap_name').apply(
                lambda df: sum(df['resid'] >= x )
            )
        ) for x in resids ],
        columns=['mean','ci']
    )
    ax.errorbar(resids,rcdf['mean'],yerr=rcdf['ci'],label='Bootstrap Residuals',color='red')  
    ax.set_xlabel('Residual')
    ax.set_ylabel('Number of Genes > X')
    legend = ax.legend(loc='best')
    plt.savefig('FDR_Resid.png')



    # Go back to beginning
    os.chdir(cwd)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()    
    parser.add_argument('--cob')
    parser.add_argument('--ontology')
    parser.add_argument('--term')
    parser.add_argument('--window-size',default=100,type=int)
    parser.add_argument('--num-bootstraps',default=50,type=int)
    parser.add_argument('--sd-envelope',default=2,type=int)
    args = parser.parse_args()
    sys.exit(main(args))
    


