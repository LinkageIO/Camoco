#!/usr/bin/env python3

import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable

import matplotlib.pylab as plt
matplotlib.style.use('ggplot')
#matplotlib.use('Agg')


import argparse
import sys
import os
import copy
import re

import numpy as np
import scipy as sp
import scipy.stats

from collections import OrderedDict

import camoco as co
import pandas as pd
pd.set_option('display.width',300)
import glob as glob


class OverlapAnalysis(object):
    def __init__(self):
        self.results = pd.DataFrame()

    def generate_output_name(self):
        # Handle the different output schemes
        if self.args.out is None:
            self.args.out = '{}_{}_{}_{}_{}'.format(
                self.cob.name,
                self.ont.name,
                self.args.candidate_window_size,
                self.args.candidate_flank_limit,
                ':'.join(self.args.terms)
            )
        self.args.out = os.path.splitext(self.args.out)[0] + '.{}.overlap.tsv'.format(self.args.method)
        if os.path.dirname(self.args.out) != '':
            os.makedirs(os.path.dirname(self.args.out),exist_ok=True)
        if os.path.exists(self.args.out) and self.args.force != True:
            print(
                "Output for {} exists! Skipping!".format(
                    self.args.out
                ),file=sys.stderr
            )
            return

    def snp2gene(self,term):
        if 'effective' in self.args.snp2gene:
            # Map to effective
            return term.effective_loci(
                window_size=self.args.candidate_window_size
            )
        elif 'strongest' in self.args.snp2gene:
            return term.strongest_loci(
                window_size=self.args.candidate_window_size,
                attr=self.args.strongest_attr,
                lowest=self.args.strongest_higher
            )

    def generate_bootstraps(self,loci,overlap):
        '''
            bootstrapping procedure. Our target here is to provide enough bootstraps
            to identify loci that are significant at n==1000 bootstraps. The auto 
            procedure will continue untill we meet n==1000 OR we find 50 bootstraps
            that are have higher score such that we will never be significant at 1000
            boostraps (0.05 * 1000 = 50).
        '''
        target_score = overlap.score.mean()
        max_bs = 1000
        num_bs = 0
        bs = []
        if self.args.num_bootstraps == 'auto':
            # Create a bullshit generator... err bootstraps
            bs_generator = (self.overlap(loci,bootstrap=True,iter_name=x) \
                for x in range(max_bs)
            )
            while num_bs <= 50 and len(bs) < 1000: 
                # Add 50 bootstraps to current bootstraps
                bs = [next(bs_generator) for x in range(50)] + bs
                # Find the number of bs more extreme than the empirical score
                num_bs = sum([df.score.mean() >= target_score for df in bs])
                self.cob.log(
                    "Iteration: {} -- current pval: {} {}% complete",
                    len(bs),num_bs/len(bs),num_bs/50
                )
        else:
            # Be a lil broke back noodle and explicitly bootstrap
            bs = [self.overlap(loci,bootstrap=True,iter_name=x) \
                for x in range(int(self.args.num_bootstraps))
            ]
        return pd.concat(bs)

    @classmethod
    def from_csv(cls,dir='./',sep='\t'):
        self = cls()
        # Read in Gene specific data
        self.results = pd.concat(
            [pd.read_table(x,sep=sep) \
                for x in glob.glob(dir+"/*.overlap.tsv") ]
        )
        return self

    @classmethod
    def from_CLI(cls,args):
        '''
            Implements an interface for the CLI to perform overlap analysis
        '''
        self = cls()
        # Build base camoco objects
        self.args = args
        self.cob = co.COB(args.cob)
        self.ont = co.GWAS(args.gwas)
        self.generate_output_name()

        # Generate a terms iterable
        if 'all' in self.args.terms:
            terms = self.ont.iter_terms()
        else:
            terms = [self.ont[term] for term in self.args.terms]
        all_results = list()

        results = []
        # Iterate through terms and calculate
        for term in terms:
            self.cob.log(
                "Calculating Overlap for {} of {} in {}",
                term.id,
                self.ont.name,
                self.cob.name
            )
            loci = self.snp2gene(term)
            overlap = self.overlap(loci)
            bootstraps = self.generate_bootstraps(loci,overlap)
            bs_mean = bootstraps.groupby('iter').score.apply(np.mean).mean()
            bs_std  = bootstraps.groupby('iter').score.apply(np.std).mean()
            # Calculate z scores for density
            overlap['zscore'] = (overlap.score-bs_mean)/bs_std
            bootstraps['zscore'] = (bootstraps.score-bs_mean)/bs_std
            # Calculate FDR
            overlap['fdr'] = np.nan
            max_zscore = int(overlap.zscore.max()) + 1
            for zscore in np.arange(0, max_zscore,0.25):
                num_random = bootstraps\
                        .groupby('iter')\
                        .apply(lambda df: sum(df.zscore >= zscore))\
                        .mean()
                num_real = sum(overlap.zscore >= zscore)
                # Calculate FDR
                if num_real != 0 and num_random != 0:
                    fdr = num_random / num_real
                elif num_real != 0 and num_random == 0:
                    fdr = 0
                else:
                    fdr = 1
                overlap.loc[overlap.zscore >= zscore,'fdr'] = fdr
                overlap.loc[overlap.zscore >= zscore,'num_real'] = num_real
                overlap.loc[overlap.zscore >= zscore,'num_random'] = num_random
                overlap.loc[overlap.zscore >= zscore,'bs_mean'] = bs_mean
                overlap.loc[overlap.zscore >= zscore,'bs_std'] = bs_std
                overlap.sort_values(by=['zscore'],ascending=False,inplace=True)
            overlap_pval = (
                (sum(bootstraps.groupby('iter').apply(lambda x: x.score.mean()) >= overlap.score.mean()))\
                 / len(bootstraps.iter.unique())
            )
            # This gets collated into all_results below
            overlap['COB'] = self.cob.name
            overlap['Ontology'] = self.ont.name
            overlap['Term'] = term.id
            overlap['WindowSize'] = self.args.candidate_window_size
            overlap['FlankLimit'] = self.args.candidate_flank_limit
            overlap['TermLoci'] = len(term.loci)
            overlap['TermCollapsedLoci'] = len(loci)
            overlap['TermPValue'] = overlap_pval
            overlap['NumBootstraps'] = len(bootstraps.iter.unique())
            overlap['Method'] = self.args.method
            results.append(overlap.reset_index())
        self.results = pd.concat(results)
        self.results.to_csv(self.args.out,sep='\t',index=None)

    def overlap(self,loci,bootstrap=False,iter_name=None):
        '''
            Calculate Network-Term Overlap based on the method in CLI args

            This method must return a DataFrame 
        '''
        # generate emirical and bootstrap values 
        if self.args.method == 'density':
            return self.cob.trans_locus_density(
                loci,
                flank_limit=self.args.candidate_flank_limit,
                by_gene = True,
                bootstrap = bootstrap,
                iter_name = iter_name
            ) 
        elif self.args.method == 'locality':
            return self.cob.trans_locus_locality(
                loci,
                flank_limit=self.args.candidate_flank_limit,
                by_gene=True,
                bootstrap=bootstrap,
                iter_name=iter_name,
                include_regression=True,
            ).rename(columns={'resid':'score'})


    def num_below_fdr(self,fdr_cutoff=0.3):
        return pd.pivot_table(
            x.results,
            index=['Method','COB','WindowSize','FlankLimit'],
            columns='Term',
            values='fdr',
            aggfunc=lambda x: sum(x<=fdr_cutoff)
        )
            
    def plot_pval_heatmap(self,filename=None,pval_cutoff=0.05):
        '''
            Generate a heatmap based on TermPVal
        '''
        methods = self.results.Method.unique()
        cobs = self.results.COB.unique()
        fig,axes = plt.subplots(len(cobs),len(methods))
        fig.set_size_inches((8,10))
        def get_axis(i,j,axes):
            if len(axes.shape) == 1:
                return axes[i]
            else:
                return axes[i,j]
        for i,cob in enumerate(cobs):
            for j,method in enumerate(methods):
                axis = get_axis(i,j,axes)
                data = pd.pivot_table(
                    self.results.ix[(self.results.COB==cob) & (self.results.Method==method)],
                    index=['WindowSize','FlankLimit'],
                    columns=['Term'],
                    values='TermPValue',
                    aggfunc=np.mean
                )
                #data[data > pval_cutoff] = np.nan
                #data[data < pval_cutoff] = 0
                axis.set_frame_on(False)
                cmap = plt.cm.Greens_r
                cmap.set_bad('white',1.0)
                im = axis.pcolormesh(
                    np.ma.masked_invalid(data),
                    cmap=cmap,
                    #linewidth=0.1,
                    edgecolor='lightgrey',
                    vmin=0,
                    vmax=0.05
                )
                # Make the layout more natural
                if j == len(methods) - 1 :
                    axis.set_ylabel(cob,fontsize=10)
                    axis.yaxis.set_label_position("right")
                if j == 0:
                    axis.set_yticklabels(
                        axis.get_yticklabels(),size=7
                    )
                    axis.set(
                        yticks=np.arange(len(data.index))+0.5
                    )
                else:
                    axis.set_yticks([])
                    axis.set_yticklabels([])
                axis.yaxis.set_ticks_position('left')
                axis.invert_yaxis()
                if i == 0:
                    axis.set_title(method,y=1.15,size=15)
                    axis.xaxis.tick_top()
                    axis.set_xticks(np.arange(len(data.columns))+0.5)
                    axis.set_xticklabels(
                        [re.sub('\d','',x) for x in data.columns.values], 
                        rotation=45,
                        size=7
                    )
                else:
                    axis.set_xticks([])
                    axis.set_xticklabels([])
                axis.set_yticklabels(data.index.values)
        # Create a new axis to append the color bar to
        #fig.subplots_adjust(right=0.8)
        #cax = fig.add_axes([.95,0.5,0.01,0.4])
        #divider = make_axes_locatable(get_axes(0,0,))
        #cax = divider.append_axes("bottom", size="5%", pad=0.05)
        #cbar = fig.colorbar(
        #    im, cax=cax, orientation='vertical', ticks=[0,0.025,0.05]
        #)
        #cbar.set_ticklabels(['0','0.025','≥0.05'])
        #plt.tight_layout(pad=0.4,w_pad=0.5, h_pad=1.0)
        if filename is not None:
            plt.savefig(filename,dpi=300)
            plt.close()
        return (fig,axes)

    def num_candidate_genes(self,fdr_cutoff=0.3,min_snp2gene_obs=2):
        df = self.results[np.isfinite(self.results.fdr)]
        df = df[df.fdr <= fdr_cutoff]
        # Filter out genes that do not occur at 2+ SNP to gene mappings, 
        # they will never be included
        df = df.groupby(
            ['COB','Term','Method','gene']
        ).filter(lambda df: len(df)>=min_snp2gene_obs).copy()
        # Find genes that are significant in both density and locality 
        # in the same element and the same network
        both_same_net = df.groupby(
            ['COB','Term','gene']
        ).filter(lambda df: len(df.Method.unique()) == 2).copy()
        both_same_net.loc[:,'Method'] = 'both_same_net'
        # Find genes that are significant in both density and locality 
        # in the same element but any network
        both_any_net = df.groupby(
            ['Term','gene']
        ).filter(lambda df: len(df.Method.unique()) == 2).copy()
        both_any_net.loc[:,'Method'] = 'both_any_net'
        both_any_net['COB'] = 'Any'
        # Concatenate
        df = pd.concat([df,both_same_net,both_any_net])
        # Pivot
        return pd.pivot_table(
            df,
            columns=['Method','COB'],
            index=['Term'],
            values='gene',
            aggfunc=lambda x: len(set(x))
        ).fillna(0).astype(int)



