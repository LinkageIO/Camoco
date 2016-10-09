#!/usr/bin/env python3

import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable

import matplotlib.pylab as plt
matplotlib.style.use('ggplot')
matplotlib.use('Agg')


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
pd.set_option('display.width',300)
import glob as glob


tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]    
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)
 

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
        cur_num = 50
        num_bs = 0
        bs = []
        if self.args.num_bootstraps == 'auto':
            while num_bs <= 50 and len(bs) < 1000: 
                # Add some bootstraps
                bs = [self.overlap(loci,bootstrap=True,iter_name=x) \
                    for x in range(cur_num)
                ] + bs
                num_bs = sum([df.score.mean() >= target_score for df in bs])
                self.cob.log("Iteration: {} -- current pval: {} {}% complete",len(bs),num_bs/len(bs),num_bs/50)
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
            
    def plot_pval_heatmap(self,filename=None,pval_cutoff=0.05):
        '''
            Generate a heatmap based on TermPVal
        '''
        cobs = self.results.COB.unique()
        fig,axes = plt.subplots(len(cobs),1,sharex=True)
        fig.set_size_inches(8,11)
        for i,cob in enumerate(cobs):
            data = pd.pivot_table(
                self.results.ix[self.results.COB==cob,self.results.method=='density'],
                index=['WindowSize','FlankLimit'],
                columns=['Term'],
                values='TermPValue',
                aggfunc=np.mean
            )
            #data[data > pval_cutoff] = np.nan
            #data[data < pval_cutoff] = 0
            axes[i].set_frame_on(False)
            cmap = plt.cm.Greens_r
            cmap.set_bad('white',1.0)
            im = axes[i].pcolormesh(
                np.ma.masked_invalid(data),
                cmap=cmap,
                #linewidth=0.1,
                edgecolor='lightgrey',
                vmin=0,
                vmax=0.05
            )
            # Make the layout more natural
            axes[i].set_ylabel(cob,fontsize=20)
            axes[i].yaxis.set_label_position("right")
            axes[i].invert_yaxis()
            axes[i].set(
                yticks=np.arange(len(data.index))+0.5
            )
            axes[i].yaxis.set_ticks_position('left')
            axes[i].xaxis.set_ticks_position('none')
            if i == 0:
                axes[i].xaxis.tick_top()
                axes[i].set_xticks(np.arange(len(data.columns))+0.5)
                axes[i].set_xticklabels(data.columns.values, rotation=90)
            axes[i].set_yticklabels(data.index.values)
        # Create a new axis to append the color bar to
        divider = make_axes_locatable(axes[len(cobs)-1])
        cax = divider.append_axes("bottom", size="5%", pad=0.05)
        cbar = fig.colorbar(
            im, cax=cax, orientation='horizontal', ticks=[0,0.025,0.05]
        )
        cbar.set_ticklabels(['0','0.025','≥0.05'])
        plt.tight_layout(pad=0.4,w_pad=0.5, h_pad=1.0)
        if filename is not None:
            plt.savefig(filename,dpi=300)
            plt.close()
        return (fig,axes)
