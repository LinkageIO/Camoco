import sys
import os

import pandas as pd
import numpy as np
import scipy as sp
import camoco as co

from itertools import chain

def candidates(args):
    '''
        Perform SNP (locus) to candidate gene mapping
    '''

    if args.out != sys.stdout:
        # Create any non-existant directories
        if os.path.dirname(args.out) != '':
            os.makedirs(os.path.dirname(args.out),exist_ok=True)
        if os.path.exists(args.out):
            print(
                "Output for {} exists! Skipping!".format(
                    args.out
                ),file=sys.stderr
            )
            return None

    cob = co.COB(args.cob)
    ont = co.GWAS(args.gwas)

    if 'all' in args.terms:
        terms = ont.iter_terms()
    else:
        terms = [ont[term] for term in args.terms]

    data = pd.DataFrame()
    results = []
    for term in terms:
        effective_loci = term.effective_loci(
            window_size=args.candidate_window_size
        )
        genes = pd.DataFrame([ x.attr for x in 
            cob.refgen.candidate_genes(
                effective_loci,
                flank_limit=args.candidate_flank_limit,
                attrs={'term':term.id},
                include_parent_locus=True,
                include_num_siblings=True,
                include_num_intervening=True,
                include_rank_intervening=True
            )
        ])
        #densities = cob.trans_locus_density(
        #    effective_loci,
        #    flank_limit=args.candidate_flank_limit,
        #    by_gene=True
        #)
        #densities.columns = ['TransDensity']
        #genes = pd.merge(genes,densities,left_on='Name',right_index=True)
        data = pd.concat([data,genes])
    data['FlankLimit'] = args.candidate_flank_limit
    data['WindowSize'] = args.candidate_window_size
    data.to_csv(args.out,index=None,sep='\t')
        
