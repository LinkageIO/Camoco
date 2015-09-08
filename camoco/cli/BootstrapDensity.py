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
            'Term','Ontology','COB','WindowSize','CandidateLimit',
            'NumSNPs','NumCandidates','EmpDensity','PValue','NumBootstraps',
            'BootstrapMean','BootstrapStd'
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
        results['Term'] = term.id
        results['Ontology'] = ont.name
        results['COB'] = cob.name

        effective_loci = term.effective_loci(
            window_size=args.candidate_window_size
        )
        candidates = cob.refgen.candidate_genes(
            effective_loci, flank_limit=args.candidate_flank_limit
        )
   
        print("Num SNPs:\t{}".format(len(term.locus_list)),file=sys.stderr)
        print("Num Effective SNPs:\t{}".format(len(effective_loci)),file=sys.stderr)
        print('Num Candidate Genes:\t{}'.format(len(candidates)),file=sys.stderr)

        results['NumSNPs'] = len(term.locus_list)
        results['NumEffective'] = len(effective_loci)
        results['NumCandidates'] = len(candidates)

        emp = cob.trans_locus_density(
            effective_loci,
            flank_limit=args.candidate_flank_limit,
            bootstrap=False
        )

        results['WindowSize'] = args.candidate_window_size
        results['CandidateLimit'] = args.candidate_flank_limit

        results['EmpDensity'] = emp
    
        print("Empirical density Z-score for {}:\t{}".format(term.id,emp),file=sys.stderr)
    
        bootstraps = np.array([
            cob.trans_locus_density(
                term.effective_loci(window_size=args.candidate_window_size),
                flank_limit=args.candidate_flank_limit,
                bootstrap=True
            ) for x in range(args.num_bootstraps)
        ])
        
        pval = sum([bs >= emp for bs in bootstraps])
        if pval > 0:
            pval = pval / args.num_bootstraps
    
        results['PValue'] = pval
        results['NumBootstraps'] = args.num_bootstraps
    
        print("p-value for {}:\t{} (n={})".format(
            term.id,pval,args.num_bootstraps),
            file=sys.stderr
        )
        print("Bootstrap mean/std:\t{},{} (n={})".format(
            bootstraps.mean(),
            bootstraps.std(),
            args.num_bootstraps
        ),file=sys.stderr)
        print('{}'.format('-'*90))
    
        results['BootstrapMean'] = bootstraps.mean()
        results['BootstrapStd'] = bootstraps.std()
        # Return the results 
        all_results = all_results.append(pd.Series(results),ignore_index=True)
    all_results.to_csv(args.out,index=None,sep='\t')



