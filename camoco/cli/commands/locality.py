#!/usr/bin/env python3

import sys
import os
import copy

import numpy as np
import scipy as sp
import scipy.stats

import camoco as co
import pandas as pd
import matplotlib.pylab as plt
import statsmodels.api as sm

from camoco.Tools import log as coblog

from collections import namedtuple

lowess = sm.nonparametric.lowess

def locality(args):
    log = coblog()
    log('\n'
        '-----------------------\n'
        '   Network Locality    \n'
        '-----------------------\n'
    )
    # Generate output dirs
    if args.out != sys.stdout:
        args.out = "{}_Locality.tsv".format(args.out.replace('.tsv',''))        
    if os.path.dirname(args.out) != '':
        os.makedirs(os.path.dirname(args.out),exist_ok=True)
    if os.path.exists("{}_Locality.tsv".format(args.out.replace('.tsv',''))):
        log(
            "{}_Locality.csv exists! Skipping!".format(
                args.out.replace('.tsv','')
            )
        )
        return None
    # Grab the COB object
    cob = co.COB(args.cob)
    gwas = co.GWAS(args.gwas)
    # If there is a different score for 'significant', update the COB object
    if args.sig_edge_zscore is not None:
        cob.set_sig_edge_zscore(args.sig_edge_zscore)
    # If all, grab a generater
    if 'all' in args.terms:
        terms = gwas.iter_terms()
    else:
        # Otherwise get the term out of the GWAS
        terms = ( gwas[x] for x in args.terms ) 

    # Add in text for axes
    locality = pd.DataFrame([
        generate_data(cob,x,args) for x in terms
    ])
    locality.to_csv(args.out, sep='\t',index=None)
        
        
def generate_data(cob,term,args):
    '''
        Generates the data according to parameters in args
    '''
    if args.snp2gene == 'effective':
        loci = sorted(term.effective_loci(
            window_size=args.candidate_window_size
        ))
    elif args.snp2gene == 'strongest':
        loci = term.strongest_loci(
            window_size=args.candidate_window_size,
            attr=args.strongest_attr,
            lowest=args.strongest_higher
        )
    else:
        raise ValueError('{} not valid snp2gene mapping'.format(args.snp2gene))
    
    loc = cob.trans_locus_locality(
        loci, args.candidate_flank_limit, bootstrap=False, 
        by_gene=False, iter_name=None, 
        include_regression=True
    ).resid.mean()

    loc_bs = np.array([
        cob.trans_locus_locality(
            loci, args.candidate_flank_limit, bootstrap=True, 
            by_gene=False, iter_name=None, 
            include_regression=True
        ).resid.mean()  \
        for x in range(args.num_bootstraps)
    ])

    
    pval = sum(loc_bs >= loc)/args.num_bootstraps
    
    record = namedtuple('Record',['COB','Term','WindowSize','FlankLimit','Locality','PVal','Size'])
    # Give em the gold
    return record(
        cob.name, term.id, args.candidate_window_size, 
        args.candidate_flank_limit, loc, pval, len(loci)
    )
