#!/usr/bin/env python3

import argparse
import sys
import os

import numpy as np
import scipy as sp
import scipy.stats

from collections import OrderedDict

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


class NearestDict(OrderedDict):
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
    # Generate output dirs
    cwd = os.getcwd()
    os.makedirs(args.out,exist_ok=True)
    os.chdir(args.out)
    os.makedirs('CSV',exist_ok=True)

    '''---------------------------------------------------
        Empirical and SD Calculations
    '''

    # Find the empirical Locality
    loc = cob.locality(
        cob.refgen.candidate_genes(
            term.effective_snps(window_size=args.candidate_window_size),
            gene_limit=args.candidate_gene_limit
        ),
        include_regression=True
    )
    
    # Find the Bootstrapped Locality
    bsloc = pd.concat(
            [cob.locality(
                cob.refgen.bootstrap_candidate_genes(
                    term.effective_snps(window_size=args.candidate_window_size),
                    gene_limit=args.candidate_gene_limit
                ),
                bootstrap_name=x
            ) for x in range(args.num_bootstraps)]
    )
    # We need to perform regression for the entire bootstrap dataset
    OLS = sm.OLS(bsloc['local'],bsloc['global']).fit()
    bsloc['fitted'] = OLS.fittedvalues
    bsloc['resid'] = OLS.resid
    # Remove global degree larger than empirical
    bslowess = lowess(
            bsloc['resid'],
            bsloc['fitted'],
            frac=0.15,
            it=5,
            delta=0.1*len(bsloc),
            is_sorted=False
        )

    # Windowing
    bsloc = bsloc.sort('fitted')
    # Find out how many tick there are with X items per window
    num_windows = len(bsloc)//args.window_size
    window_ticks = len(bsloc)//num_windows
    bsloc['window'] = [int(i/window_ticks) for i in range(len(bsloc))]
    # If there are not many in the last window, change it second to last
    max_window = max(bsloc['window'])
    if sum(bsloc['window'] == max_window) < args.window_size / 2:
        bsloc.loc[bsloc['window'] == max_window, 'window'] = max_window-1
    # create a dictionary so we can map the empirical data later
    win_map = NearestDict({
        # Good god this is a hack -- 
        # group the pandas df by window and calculate the mean fitted value, 
        # create dict from that and reverser keys and values
        fitted:window for window,fitted in bsloc.groupby('window').apply(
                lambda df: np.max(df['fitted'])
            ).to_dict().items() 
    })
    # create a dict of window to std mapping
    win_std = bsloc.groupby('window').apply(lambda df: df['resid'].std())
    # perform lowess on the std_windows
    win_std = NearestDict(
        win_std.to_dict()
    )
    fit_std = {f:win_std[w]for f,w in win_map.items()}
    
   
    # Create a dict where keys are fitted valus and values are that fitted values std
    fit_std = NearestDict(
        pd.DataFrame(
            lowess(np.array(list(fit_std.values())),np.array(list(fit_std.keys())),is_sorted=False),columns=['fitted','sd']
        ).sort('fitted').set_index('fitted').to_dict()['sd']
    )


    # Calculate the s.d. of the residuals in each window
    # divide empirical residuals by the s.d. in their respective window
    loc['bs_std'] = [fit_std[x] for x in loc['fitted']]
    loc['zscore'] = [x['resid']/x['bs_std'] for i,x in loc.iterrows()]

    '''---------------------------------------------------
        FDR Calculations
    '''
    
    # Repeat bootstraps to assess global FDR
    fdr = pd.concat(
            [cob.locality(
                cob.refgen.bootstrap_candidate_genes(
                    term.effective_snps(window_size=args.candidate_window_size),
                    gene_limit=args.candidate_gene_limit
                ),
                bootstrap_name=x
            ) for x in range(args.num_bootstraps)]
    ).sort('global')
    OLS = sm.OLS(fdr['local'],fdr['global']).fit()
    fdr['fitted'] = OLS.fittedvalues
    fdr['resid'] = OLS.resid

    # Remove global degree larger than empirical
    fdr = fdr[fdr['global'] <= max(loc['global'])]
    fdr['window'] = [int(x/window_ticks) for x in fdr['global']]

    # calculate z-scores for the global 
    fdr['bs_std'] = [fit_std[x] for x in fdr['fitted']]
    fdr['zscore'] = [x['resid']/x['bs_std'] for i,x in fdr.iterrows()]

    # Output the datafiles
    loc.sort('zscore',ascending=False).to_csv(
        os.path.join(
            'CSV',
            '{}.csv'.format('Empirical_'.join(map(str,[args.term,args.cob,args.window_size,args.term])))
        )
    )
    bsloc.to_csv(
        os.path.join(
            'CSV',
            '{}.csv'.format('Bootstrap_'.join(map(str,[args.term,args.cob,args.window_size,args.term])))
        )
    )
    fdr.sort(['bootstrap_name','zscore'],ascending=False).to_csv(
        os.path.join(
            'CSV',
            '{}.csv'.format('FDR_'.join(map(str,[args.term,args.cob,args.window_size,args.term])))
        )
    )


    ''' ---------------------------------------------------
        Plotting
    '''
    # Make plots for later
    gs = plt.GridSpec(3,2)
    fig = plt.figure(figsize=(16,16))
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[1,0])
    ax4 = fig.add_subplot(gs[1,1])
    ax5 = fig.add_subplot(gs[2,:1])
    ax1.hold(True)
    # Y axis is local degree (what we are TRYING to predict)
    ax1.set_ylim(0,max(loc['local']))
    ax1.set_xlim(0,max(loc['global']))
    ax1.set_xlabel('Number Global Interactions')
    ax1.set_ylabel('Number Local Interactions')
    ax1.set_title('OLS Fitting')

    # UGH! map lowess 
    fdrlowess = lowess(fdr['local'],fdr['global'],frac=0.15,it=5,delta=0.1*len(fdr),is_sorted=False)
    # plot the bootstrap points
    ax1.plot(fdr['global'],fdr['local'],'ro',alpha=0.05,label='Bootstraps')
    # plot the OLS lowess line
    ax1.plot(fdr['global'],fdr['fitted'],'g--',label='Bootstrap OLS')

    #plot the empirical data
    ax1.plot(loc['global'],loc['local'],'bo',alpha=1,label='Empirical')
    ax1.plot(loc['global'],loc['fitted'],'k-',alpha=1,label='Empirical OLS')
    # finish plots
    legend = ax1.legend(loc='best') 

    '''---------------------------------------------------
        FDR Plotting
    '''
    

    # Plot the empirical Z-score distributions
    zscores = list(range(1,8))
    zcdf = [sum(np.logical_and(loc['zscore'] >= x , loc['local'] >= args.min_fdr_degree)) for x in zscores]
    ax2.plot(zscores,zcdf,'bo',label='Empirical Zscore > x')
    # plot the fdr scores spread
    zcdf = pd.DataFrame(
        [ mean_confidence_interval(
            fdr.groupby('bootstrap_name').apply(
                lambda df: sum(np.logical_and(df['zscore'] >= x, df['local'] >= args.min_fdr_degree ))
            )
        ) for x in zscores ],
        columns=['mean','ci']
    )
    ax2.errorbar(zscores,zcdf['mean'],yerr=zcdf['ci'],label='Bootstrap Z-scores',color='red')  
    ax2.set_xlabel('Z-Score')
    ax2.set_ylabel('Number of Genes > Z')
    ax2.set_title('Z Score FDR')
    legend = ax2.legend(loc='best')

    # Do the same thing, but with residuals


    resids = list(range(0,int(max(loc['resid'])+2)))
    rcdf = [sum(loc['resid'] >= x) for x in resids]
    ax3.plot(resids,rcdf,'bo',label='Empirical Residuals')
    rcdf = pd.DataFrame(
         [ mean_confidence_interval(
            fdr.groupby('bootstrap_name').apply(
                lambda df: sum(df['resid'] >= x )
            )
        ) for x in resids ],
        columns=['mean','ci']
    )
    ax3.errorbar(resids,rcdf['mean'],yerr=rcdf['ci'],label='Bootstrap Residuals',color='red')  
    ax3.set_xlabel('Residual')
    ax3.set_ylabel('Number of Genes > X')
    ax3.set_title('Residual FDR')
    legend = ax3.legend(loc='best')

    # Plot the predicted vs residuals

    '''---------------------------------------------------
        Residual Plotting
    '''

    ax4.set_xlabel('Fitted')
    ax4.set_ylabel('Residual')
    ax4.set_title('Model Fitting')
    ax4.plot(loc['fitted'],loc['resid'],'bo')
    # plot the Z-scores for the std
    ax4.plot(list(fit_std.keys()),list(fit_std.values()),'g.',marker='.')
    # plot where windows are

    ax5.xaxis.set_visible(False)
    ax5.yaxis.set_visible(False)
    ax5.text(.1,.1,'''
        COB: {}
        Ontology: {}
        Term: {}
        Window Size: {}
        Num. Bootstraps: {}
        Min FDR Degree: {}
        Candidate Gene Window Size: {}
        Candidate Gene Limit: {}
        Num SNPs: {}
    '''.format(args.cob,args.ontology,args.term,args.window_size,args.num_bootstraps,args.min_fdr_degree,args.candidate_window_size,args.candidate_gene_limit,args.num_snps))

    plt.savefig(
        '{}_{}_{}_{}_{}_{}_{}_{}_{}.png'.format(args.term,args.cob,args.ontology,args.window_size,args.num_bootstraps,args.min_fdr_degree,args.candidate_window_size,args.candidate_gene_limit,args.num_snps)
    )
    os.chdir(cwd)

    if args.debug:
        import ipdb; ipdb.set_trace()



if __name__ == '__main__':
    parser = argparse.ArgumentParser()    
    parser.add_argument('--cob', help='The camoco network to use.')
    parser.add_argument('--ontology', help='The camoco ontology to use.')
    parser.add_argument('--term', help='The term within the ontology to use.')
    parser.add_argument('--window-size', default=15, type=int, help='The number of items within a window.')
    parser.add_argument('--num-bootstraps', default=50, type=int, help='The number of bootstraps to perform in order to estimate a null distribution.')
    parser.add_argument('--sd-envelope', default=2, type=int, help='The number or standard deviations to use for the regression window.')
    parser.add_argument('--min-fdr-degree', default=2, type=int, help='The miniumum degree to be included as true positive in FDR calculation.')
    parser.add_argument('--candidate-window-size',default=250000,type=int,help='The window size for mapping effective SNPs to genes.')
    parser.add_argument('--candidate-gene-limit',default=2,type=int,help='The number of flanking genes included in SNP to gene mapping')
    parser.add_argument('--num-snps',default=None,type=int,help='The number of SNPs to ')
    parser.add_argument('--out',default='Output',type=str,help='Name of output directory')
    parser.add_argument('--debug',default=False,action='store_true',help='Drop into an ipdb debugger before quitting.')
    args = parser.parse_args()
    sys.exit(main(args))
    


