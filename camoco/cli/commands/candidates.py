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

    results = []
    for term in terms:
        results.append(
            cob.refgen.candidate_genes(
                term.effective_loci(window_size=args.candidate_window_size),
                flank_limit=args.candidate_flank_limit,
                attrs={'term':term.id},
                include_parent_locus=True,
                include_num_siblings=True,
                include_num_intervening=True,
                include_rank_intervening=True
            )
        )
    pd.DataFrame(
        [x.attr for x in chain(*results)]
    ).to_csv(args.out,index=None,sep='\t')
        
