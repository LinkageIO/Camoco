#!/usr/bin/env python3

import argparse
import sys
import os
import copy

import numpy as np
import scipy as sp
import scipy.stats

from collections import OrderedDict

import camoco as co
import pandas as pd
import matplotlib.pylab as plt
import statsmodels.api as sm

lowess = sm.nonparametric.lowess

def mean_confidence_interval(data):
    return np.mean(data), confidence_interval(data)

def confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return 1.96*se

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
    # Generate output dirs

    ont = co.Ontology(args.ontology)
    term = ont[args.term]
    cob = co.COB(args.cob)

    # Create a plot comparing spread of specified shit
    xaxes_key,yaxes_key = [x.replace('-','_') for x in args.axes]
    xaxes = getattr(args,xaxes_key)
    yaxes = getattr(args,yaxes_key)

    #Open a log file
    gs = plt.GridSpec(len(xaxes),len(yaxes))
    fig = plt.figure(figsize=(4*len(yaxes),4*len(xaxes)))

    # This gets a list of functions so we can map the string passed into
    # args.plot to a plotting function
    possibles = globals().copy()
    possibles.update(locals())

    # Add in text for axes
    for plot in args.plot:
        for i,xaxis in enumerate(xaxes):
            for j,yaxis in enumerate(yaxes):
                ax = fig.add_subplot(gs[i,j]) 
                print("Looking at {}:{} x {}:{}".format(
                    xaxes_key,xaxis,
                    yaxes_key,yaxis,
                ))
                perm_args = copy.deepcopy(args) 
                # Set the values we aren't permuting to their first value
                # I hate this
                for arg in perm_args.permutable:
                    if arg not in [xaxes_key,yaxes_key]:
                        setattr(perm_args,arg,getattr(args,arg)[0])
                # Change the relevent things in the permuted args 
                setattr(perm_args,xaxes_key,xaxis)
                setattr(perm_args,yaxes_key,yaxis)
                # Generate data using permuted arguments
                loc,bsloc,fdr = generate_data(cob,term,perm_args) 
                # Get the plot function based on args
                plot_function = possibles.get('plot_'+plot)
                # Plot the data
                plot_function(perm_args,loc,bsloc,fdr,ax)
                if i == 0:
                    ax.set_title(yaxis)
                if j == 0:
                    ax.set_ylabel(xaxis)
                if i == 0 and j == 0:
                    ax.set_title(str(yaxes_key)+' '+str(yaxis))
                    ax.set_ylabel(str(xaxes_key)+' '+str(xaxis))
        plt.tight_layout()
        plt.savefig("{}_{}.png".format(
            plot,
            args.out.replace('.png','')
        ))
    
    if args.debug:
        import ipdb; ipdb.set_trace()


def generate_data(cob,term,args):
    '''
        Generates the data according to parameters in args
    '''

    effective_snps = term.effective_snps(
            window_size=args.candidate_window_size
    ) 
    candidate_genes = cob.refgen.candidate_genes(
        effective_snps,
        gene_limit=args.candidate_gene_limit
    )
    # Find the empirical Locality
    loc = cob.locality(
        candidate_genes,
        include_regression=True
    )
    
    # Find the Bootstrapped Locality
    bsloc = pd.concat(
            [cob.locality(
                cob.refgen.bootstrap_candidate_genes(
                    term.effective_snps(
                        window_size=args.candidate_window_size
                    ),
                    gene_limit=args.candidate_gene_limit
                ),
                bootstrap_name=x
            ) for x in range(args.num_bootstraps)]
    )

    '''---------------------------------------------------
        Empirical and SD Calculations
    '''

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
    
    # Create a dict where keys are fitted valus and 
    # values are that fitted values std
    fit_std = NearestDict(
        pd.DataFrame(
            lowess(
                np.array(list(fit_std.values())),
                np.array(list(fit_std.keys())),
                is_sorted=False
            ),columns=['fitted','sd']
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
                    term.effective_snps(
                        window_size=args.candidate_window_size
                    ),
                    gene_limit=args.candidate_gene_limit
                ),
                bootstrap_name=x,
                include_regression=True
            ) for x in range(args.num_bootstraps)]
    ).sort('global')
    OLS = sm.OLS(fdr['local'],fdr['global']).fit()

    # Remove global degree larger than empirical
    fdr = fdr[fdr['global'] <= max(loc['global'])]
    fdr['window'] = [int(x/window_ticks) for x in fdr['global']]

    # calculate z-scores for the global 
    fdr['bs_std'] = [fit_std[x] for x in fdr['fitted']]
    fdr['zscore'] = [x['resid']/x['bs_std'] for i,x in fdr.iterrows()]

    # Give em the gold
    return loc,bsloc,fdr


def plot_scatter(args,loc,bsloc,fdr,ax):

    ''' ---------------------------------------------------
        Plotting
    '''
    # Y axis is local degree (what we are TRYING to predict)
    ax.set_ylim(0,max(loc['local']))
    ax.set_xlim(0,max(loc['global']))
    ax.set_xlabel('Number Global Interactions')
    ax.set_ylabel('Number Local Interactions')

    # UGH! map lowess 
    fdrlowess = lowess(
        fdr['local'],fdr['global'],
        frac=0.15,it=5,delta=0.1*len(fdr),
        is_sorted=False
    )
    # plot the bootstrap points
    ax.plot(fdr['global'],fdr['local'],'ro',alpha=0.05,label='Bootstraps')
    # plot the OLS lowess line
    ci = fdr.groupby('window')['fitted','global'].agg(
        [np.mean,confidence_interval]
    )
    
    #for win,df in fdr.groupby('bootstrap_name'):
    #    ax.plot(df['global'],df['fitted'],alpha=1)
        
    ax.errorbar(
        ci['global','mean'],ci['fitted','mean'],
        yerr=ci['fitted','confidence_interval'],
        color='g',label='Bootstrap OLS'
    )

    #plot the empirical data
    ax.plot(loc['global'],loc['local'],'bo',alpha=1,label='Empirical')
    ax.plot(loc['global'],loc['fitted'],'k-',alpha=1,label='Empirical OLS')
    # finish plots
    #legend = ax.legend(loc='best') 

    '''---------------------------------------------------
        FDR Plotting
    '''
    

def plot_fdr(args,loc,bsloc,fdr,ax):
    # Plot the empirical Z-score distributions
    zscores = list(range(1,8))
    zloc = [
        sum(np.logical_and(
            loc['zscore'] >= x ,
            loc['local'] >= args.min_fdr_degree
        )) 
        for x in zscores
    ]
    ax.plot(zscores,zloc,'bo',label='Empirical Zscore > x')
    # plot the fdr scores spread
    zcdf = pd.DataFrame(
        [ mean_confidence_interval(
            fdr.groupby('bootstrap_name').apply(
                lambda df: sum(np.logical_and(
                    df['zscore'] >= x,
                    df['local'] >= args.min_fdr_degree 
                ))
            )
        ) for x in zscores ],
        columns=['mean','ci']
    )
    #zcdf['emp'] = zloc
    #zcdf['el'] = args.term
    #zcdf['min_fdr_degree'] = args.min_fdr_degree
    #zcdf['can_gene_limit'] = args.candidate_gene_limit
    ax.errorbar(
        zscores,
        zcdf['mean'],
        yerr=zcdf['ci'],
        label='Bootstrap Z-scores',
        color='red'
    )  
    ax.set_xlabel('Z-Score')
    ax.set_ylabel('Number of Genes > Z')
    ax.set_title('Z Score FDR')


def plot_resid(args,loc,bsloc,fdr,ax):
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
    ax.errorbar(
        resids,
        rcdf['mean'],
        yerr=rcdf['ci'],
        label='Bootstrap Residuals',
        color='red'
    )  
    ax.set_xlabel('Residual')
    ax.set_ylabel('Number of Genes > X')
    ax.set_title('Residual FDR')
    legend = ax.legend(loc='best')

    # Plot the predicted vs residuals

    '''---------------------------------------------------
        Residual Plotting
    '''

    ax.set_xlabel('Fitted')
    ax.set_ylabel('Residual')
    ax.set_title('Model Fitting')
    ax.plot(loc['fitted'],loc['resid'],'bo')
    # plot the Z-scores for the std
    ax.plot(list(fit_std.keys()),list(fit_std.values()),'g.',marker='.')
    # plot where windows are

def plot_data(args,loc,bsloc,fdr,ax):
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.text(.1,.1,'''
        COB: {}
        Ontology: {}
        Term: {}
        Num Raw SNPS: {}
        Num Effective SNPS: {}
        Num Genes: {}
        CI Window Size: {}
        Num. Bootstraps: {}
        Min FDR Degree: {}
        Candidate Gene Window Size: {}
        Candidate Gene Limit: {}
        Num SNPs: {}
    '''.format(
        args.cob,
        args.ontology,
        args.term,
        len(term.locus_list),
        len(effective_snps),
        len(candidate_genes),
        args.window_size,
        args.num_bootstraps,
        args.min_fdr_degree,
        args.candidate_window_size,
        args.candidate_gene_limit,
        args.num_snps)
    )



    out_string = '_'.join(map(str,[
        args.term,
        args.cob,
        args.window_size,
        args.min_fdr_degree,
        args.candidate_gene_limit,
        args.candidate_window_size
    ]))

    #zcdf.to_csv(
    #    os.path.join(
    #        'CSV',
    #        '{}_FDR.csv'.format(out_string)
    #    )        
    #)


   ## Output the datafiles
   #loc.sort('zscore',ascending=False).to_csv(
   #    os.path.join(
   #        'CSV',
   #        'Empirical_{}.csv'.format(out_string)
   #    )
   #)
   #bsloc.to_csv(
   #    os.path.join(
   #        'CSV',
   #        'Bootstrap_{}.csv'.format(out_string)
   #    )
   #)
   #fdr.sort(['bootstrap_name','zscore'],ascending=False).to_csv(
   #    os.path.join(
   #        'CSV',
   #        'FDR_{}.csv'.format(out_string)
   #    )
   #)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()    
    # Data Set Arguments
    parser.add_argument(
        '--cob', help='The camoco network to use.'
    )
    parser.add_argument(
        '--ontology', help='The camoco ontology to use.'
    )
    parser.add_argument(
        '--term', help='The term within the ontology to use.'
    )

    # Data Analysis parameters  
    parser.add_argument(
       '--axes', nargs='*',
       help='This are the x and y axis for the plot'
    )
    parser.add_argument(
       '--num-bootstraps', default=50,type=int, 
       help='''The number of bootstraps to perform in order
             to estimate a null distribution.'''
    )
    parser.add_argument(
       '--sd-envelope', default=2, type=int,
       help='''The number or standard deviations to 
           use for the regression window.'''
    )
    parser.add_argument(
       '--window-size', default=15, 
       type=int, help='The number of items within a window.'
    )

    parser.add_argument(
       '--permutable', default=[
            'min_fdr_degree',
            'candidate_window_size',
            'candidate_gene_limit'
        ]        
    )
    # Permutation parameters
    parser.add_argument(
       '--min-fdr-degree', default=2, nargs='*',
       type=int, help='''The miniumum degree to be included 
       as true positive in FDR calculation.'''
    )
    parser.add_argument(
       '--candidate-window-size',default=50000,
       type=int, nargs='*',  
       help='The window size for mapping effective SNPs to genes.'
    )
    parser.add_argument(
       '--candidate-gene-limit',default=2,
       type=int, nargs='*',
       help='The number of flanking genes included in SNP to gene mapping'
    )
    # Dont worry about this right now
    # parser.add_argument(
    #    '--num-snps',default=None, nargs='*',
    #    type=int,help='The number of SNPs to '
    # )

    # Data Formatting Parameters
    parser.add_argument(
       '--out',default='Output',
       type=str,help='Name of output directory'
    )
    parser.add_argument(
       '--debug',default=False,
       action='store_true',
       help='Drop into an ipdb debugger before quitting.'
    )

    # Make this a list  
    parser.add_argument(
        '--plot',default='fdr',
        action='store', nargs='*',
        help=("Designates which item to plot. Must be in: ['fdr','scatter']")
    )

    with open('command_log.txt','a') as LOG:
        print('{}'.format(' '.join(sys.argv)),file=LOG)
    args = parser.parse_args()
    sys.exit(main(args))



