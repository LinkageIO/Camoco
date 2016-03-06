import sys
import os
import progressbar

import pandas as pd
import numpy as np
import camoco as co
from collections import OrderedDict

def simulateGWAS(args):
    # Make sure that output files are handled correctly 
    if args.out != sys.stdout:
        args.out = "{}_GWASSim.csv".format(args.out.replace('.csv',''))        
    if os.path.dirname(args.out) != '':
        os.makedirs(os.path.dirname(args.out),exist_ok=True)
    if os.path.exists("{}_GWASSim.csv".format(args.out.replace('.csv',''))):
        print(
            "{}_GWASSim.csv exists! Skipping!".format(
                args.out.replace('.csv','')
            )
        )
        return None

    go = co.GOnt(args.GOnt)
    cob = co.COB(args.cob)
    terms = (x for x in go.iter_terms() \
            if len(x.loci) >= args.min_term_size \
            and len(x.loci) <= args.max_term_size 
    )
    results = []
    for i,term in enumerate(terms):
        cob.log('On term {}',i)   
        window_size = args.candidate_window_size
        flank_limit =  args.candidate_flank_limit
        # Generate a series of densities for parameters
        num_genes = len(term.loci)
        eloci = term.effective_loci(
            window_size=window_size
        )
        if args.percent_fdr != None and args.percent_fdr > 0:                                            
            # replace some loci with random genes if FDR specified              
            num_fdr = int(len(eloci) * args.percent_fdr)                        
            fdr_loci = cob.refgen.random_genes(num_fdr,window=window_size) 
            # permute and truncate the loci then add fdr loci                   
            eloci = np.concatenate([                                            
                np.random.permutation(eloci)[0:-1*len(fdr_loci)],                   
                np.array(list(fdr_loci))
            ])   
        candidates = cob.refgen.candidate_genes(
            eloci,
            flank_limit=flank_limit
        )
        if len(candidates) == 0:
            density = np.NaN
            trans_density = np.NaN
            pval = np.NaN
        else:
            density = cob.density(
                candidates
            )
            trans_density = cob.trans_locus_density(
                eloci,
                flank_limit=flank_limit
            )
            bootstraps = np.array([
                cob.trans_locus_density(
                    eloci,
                    flank_limit=flank_limit,
                    iter_name = x,
                    bootstrap=True
                ) for x in range(args.num_bootstraps)
            ])
            pval = sum([x>=trans_density for x in bootstraps])/args.num_bootstraps
        results.append(OrderedDict([
            ('COB',cob.name),
            ('GO',term.id),
            ('NumGenes',num_genes),
            ('WindowSize',window_size),
            ('FlankLimit',flank_limit),
            ('FDR',args.percent_fdr),
            ('NumEffective',len(eloci)),
            ('NumCandidates',len(candidates)),
            ('Density',density),
            ('TransDensity',trans_density),
            ('Pval',pval)
        ]))
    results = pd.DataFrame.from_dict(results)
    results.to_csv(args.out)
