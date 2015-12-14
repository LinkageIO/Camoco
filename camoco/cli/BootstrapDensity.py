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

co.cf.logging.log_level = 'quiet'

def density(args):
    
    cob = co.COB(args.cob)
    ont = co.GWAS(args.gwas)
    args.out = os.path.splitext(args.out)[0] + '_Density.csv'

    if os.path.dirname(args.out) != '':
        os.makedirs(os.path.dirname(args.out),exist_ok=True)
    if os.path.exists(args.out):
        print(
            "Output for {} exists! Skipping!".format(
                args.out
            ),file=sys.stderr
        )
        return

    if 'all' in args.terms:
        terms = ont.iter_terms()
    else:
        terms = [ont[term] for term in args.terms]
    
    all_results = list()

    for term in terms:
        results = OrderedDict()

        print("Boostrapping Density for {} of {} in {}".format(
            term.id,
            ont.name,
            cob.name
        ),file=sys.stderr)
   

        results['Ontology'] = ont.name
        results['COB'] = cob.name
        results['Term'] = term.id
        results['NumSNPs'] = len(term.loci)
        results['WindowSize'] = args.candidate_window_size
        results['FlankLimit'] = args.candidate_flank_limit
        results['NumBootstraps'] = args.num_bootstraps
    
        if 'effective' in args.snp2gene:
            # Map to effective
            loci = term.effective_loci(
                window_size=args.candidate_window_size
            )
        elif 'strongest' in args.snp2gene:
            loci = term.strongest_loci(
                window_size=args.candidate_window_size,
                attr=args.strongest_attr,
                lowest=args.strongest_higher
            )
     
        # Grab those candidates
        candidates = cob.refgen.candidate_genes(
            loci, flank_limit=args.candidate_flank_limit
        )
        emp = cob.trans_locus_density(
            loci,
            flank_limit=args.candidate_flank_limit,
            bootstrap=False
        )

        # ----------------------------
        # Effective Loci Bootstraps
        bootstraps = np.array([
            cob.trans_locus_density(
                loci,
                flank_limit=args.candidate_flank_limit,
                bootstrap=True
            ) for x in range(args.num_bootstraps)
        ])
        
        effective_pval = sum([bs >= effective_emp for bs in bootstraps])
        if effective_pval > 0:
            effective_pval = effective_pval / args.num_bootstraps
        
        results['NumCollapsedSNPs'] = len(loci)
        results['NumCandidates'] = len(candidates)
        results['Density'] = emp
        results['PValue'] = effective_pval
        results['BSMean'] = effective_bootstraps.mean()
        results['BSStd'] = effective_bootstraps.std()
        results['SNP2GeneMethod'] = args.snp2gene
            
        all_results.append(results)
    
    # Return the results 
    all_results = pd.DataFrame(all_results,columns=all_results[0].keys())
    # Make sure the output directory exists
    all_results.to_csv(args.out,index=None,sep='\t')
