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

class OverlapAnalysis(object):

    def __init__(self,args):
        # Build base camoco objects
        self.args = args
        self.cob = co.COB(args.cob)
        self.ont = co.GWAS(args.gwas)
        # Handle the different output schemes
        if self.args.out is None:
            self.args.out = '{}_{}_{}_{}_{}'.format(
                self.cob.name,
                self.ont.name,
                self.args.candidate_window_size,
                self.args.candidate_flank_limit,
                ':'.join(args.terms)
            )
        self.args.out = os.path.splitext(args.out)[0] + '.{}.overlap.tsv'.format(self.args.method)
        if os.path.dirname(args.out) != '':
            os.makedirs(os.path.dirname(args.out),exist_ok=True)
        if os.path.exists(args.out) and self.args.force != True:
            print(
                "Output for {} exists! Skipping!".format(
                    self.args.out
                ),file=sys.stderr
            )
            return
        self.analyze()

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

    def analyze(self):
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
            bootstraps = pd.concat([self.overlap(loci,bootstrap=True,iter_name=x) \
                    for x in range(self.args.num_bootstraps)
            ])
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
                 / self.args.num_bootstraps
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
            overlap['NumBootstraps'] = self.args.num_bootstraps
            results.append(overlap.reset_index())
        results = pd.concat(results)
        results.to_csv(self.args.out,sep='\t',index=None)

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
            
