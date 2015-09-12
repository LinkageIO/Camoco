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

from camoco.Config import cf
cf.logging.log_level = 'quiet'

def density(args):
    
    cob = co.COB(args.cob)
    ont = co.GWAS(args.gwas)

    if 'all' in args.terms:
        terms = ont.iter_terms()
    else:
        terms = [ont[term] for term in args.terms]
    
    all_results = pd.DataFrame(
        columns=[
            'Term','Ontology','COB','NumSNPs',
            'WindowSize','FlankLimit',
            'NumCollapsedSNPs','NumBootstraps',
            # All the Effective Things
            'NumEffectiveCandidates',
            'EffectiveDensity',
            'EffectivePValue','EffectiveBSMean','EffectiveBSStd',
            # All of the Strongest Things
            'NumStrongestCandidates',
            'StrongestDensity',
            'StrongestPValue','StrongestBSMean','StrongestBSStd'
        ]
    )

    for term in terms:
        results = {}

        print('-'*90,file=sys.stderr)
        print("Boostrapping Density for {} of {} in {}".format(
            term.id,
            ont.name,
            cob.name
        ),file=sys.stderr)
   
        print("Num SNPs:\t{}".format(len(term.locus_list)),file=sys.stderr)

        results['Term'] = term.id
        results['Ontology'] = ont.name
        results['COB'] = cob.name
        results['NumSNPs'] = len(term.locus_list)
        results['WindowSize'] = args.candidate_window_size
        results['FlankLimit'] = args.candidate_flank_limit
        results['NumBootstraps'] = args.num_bootstraps
    
        if 'effective' in args.snp2gene:
            # Map to effective
            effective_loci = term.effective_loci(
                window_size=args.candidate_window_size
            )
    
            # Grab those candidates
            effective_candidates = cob.refgen.candidate_genes(
                effective_loci, flank_limit=args.candidate_flank_limit
            )
            effective_emp = cob.trans_locus_density(
                effective_loci,
                flank_limit=args.candidate_flank_limit,
                bootstrap=False
            )

            # ----------------------------
            # Effective Loci Bootstraps
            effective_bootstraps = np.array([
                cob.trans_locus_density(
                    effective_loci,
                    flank_limit=args.candidate_flank_limit,
                    bootstrap=True
                ) for x in range(args.num_bootstraps)
            ])
            
            effective_pval = sum([bs >= effective_emp for bs in effective_bootstraps])
            if effective_pval > 0:
                effective_pval = effective_pval / args.num_bootstraps
        

            print("Num Collapsed SNPs:\t{}".format(len(effective_loci)),file=sys.stderr)
            print('Num Effective Genes:\t{}'.format(len(effective_candidates)),file=sys.stderr)
            print("Empirical effective density Z-score for {}:\t{}".format(
                term.id,effective_emp
            ),file=sys.stderr)
            print("Effective p-value for {}:\t{} (n={})".format(
                term.id,effective_pval,args.num_bootstraps),
                file=sys.stderr
            )
            print("Effective Bootstrap mean/std:\t{},{} (n={})".format(
                effective_bootstraps.mean(),
                effective_bootstraps.std(),
                args.num_bootstraps
            ),file=sys.stderr)
            print('{}'.format('-'*90))
       
            results['NumCollapsedSNPs'] = len(effective_loci)
            results['NumEffectiveCandidates'] = len(effective_candidates)
            results['EffectiveDensity'] = effective_emp
            results['EffectivePValue'] = effective_pval
            results['EffectiveBSMean'] = effective_bootstraps.mean()
            results['EffectiveBSStd'] = effective_bootstraps.std()
            
       
        if 'strongest' in args.snp2gene:
            strongest_loci = term.strongest_loci(
                window_size=args.candidate_window_size,
                attr=args.strongest_attr,
                lowest=args.strongest_higher
            )
            strongest_candidates = cob.refgen.candidate_genes(
                strongest_loci, flank_limit=args.candidate_flank_limit
            )
            strongest_emp = cob.trans_locus_density(
                strongest_loci,
                flank_limit=args.candidate_flank_limit,
                bootstrap=False
            )
            print("Empirical strongest density Z-score for {}:\t{}".format(
                term.id,strongest_emp),file=sys.stderr
            )
            print('Num Strongest Genes:\t{}'.format(len(strongest_candidates)),file=sys.stderr)
            results['NumStrongestCandidates'] = len(strongest_candidates)
            results['StrongestDensity'] = strongest_emp
            # ----------------------------
            # strongest Loci Bootstraps
            strongest_bootstraps = np.array([
                cob.trans_locus_density(
                    strongest_loci,
                    flank_limit=args.candidate_flank_limit,
                    bootstrap=True
                ) for x in range(args.num_bootstraps)
            ])
            
            strongest_pval = sum([bs >= strongest_emp for bs in strongest_bootstraps])
            if strongest_pval > 0:
                strongest_pval = strongest_pval / args.num_bootstraps
        
            results['StrongestPValue'] = strongest_pval
        
            print("Strongest p-value for {}:\t{} (n={})".format(
                term.id,strongest_pval,args.num_bootstraps),
                file=sys.stderr
            )
            print("Bootstrap mean/std:\t{},{} (n={})".format(
                strongest_bootstraps.mean(),
                strongest_bootstraps.std(),
                args.num_bootstraps
            ),file=sys.stderr)
            print('{}'.format('-'*90))
       
            results['NumCollapsedSNPs'] = len(strongest_loci)
            results['NumStrongestCandidates'] = len(strongest_candidates)
            results['StrongestDensity'] = strongest_emp
            results['StrongestPValue'] = strongest_pval
            results['StrongestBSMean'] = strongest_bootstraps.mean()
            results['StrongestBSStd'] = strongest_bootstraps.std()
    
        # Return the results 
        all_results = all_results.append(pd.Series(results),ignore_index=True)
    all_results.to_csv(args.out,index=None,sep='\t')



