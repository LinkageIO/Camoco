#!/usr/bin/env python3

import argparse
import sys
import os
import copy

import numpy as np
import scipy as sp
import scipy.stats


import camoco as co
import pandas as pd
import matplotlib.pylab as plt

from camoco.Config import cf

def main(args):
    
    cob = co.COB(args.cob)
    ont = co.Ontology(args.ontology)
    cf['options']['log_level'] = 'quiet'
    

    terms = []

    if args.term == 'all':
        terms = ont.iter_terms()
    else:
        terms = [ont[args.term]]

    for term in terms:

        print('-'*90)
        print("Boostrapping Density for {} of {} in {}".format(
            term.name,
            ont.name,
            cob.name
        ))
    
        emp = cob.trans_locus_density(
            term.effective_snps(window_size=args.candidate_window_size),
            gene_limit=args.candidate_gene_limit,
            bootstrap=False
        )
    
        print("Empirical density Z-score for {}: {}".format(term.name,emp))
    
        if emp < args.min_score:
            print("Empirical below Z-Score")
            continue
            
    
        bootstraps = np.array([
            cob.trans_locus_density(
                term.effective_snps(window_size=args.candidate_window_size),
                gene_limit=args.candidate_gene_limit,
                bootstrap=True
            ) for x in range(args.num_bootstraps)
        ])
    
        pval = sum([bs >= emp for bs in bootstraps])
        if pval > 0:
            pval = pval / args.num_bootstraps
    
    
        print("p-value for {}: {} (n={})".format(
            term.name,pval,args.num_bootstraps)
        )
        print("Bootstrap mean/std: {},{} (n={})".format(
            bootstraps.mean(),
            bootstraps.std(),
            args.num_bootstraps
        ))


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
    parser.add_argument(
        '--min-score', help='Term with emp. less than this will be ignored',
        type=int, default=2
    )
    parser.add_argument(
       '--num-bootstraps', default=50,type=int, 
       help='''The number of bootstraps to perform in order
             to estimate a null distribution.'''
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

    args = parser.parse_args()
    sys.exit(main(args))
