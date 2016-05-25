import sys
import os

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
    if 'all' in args.terms:
        terms = go.iter_terms()
    else:
        # Otherwise get the term out of the GWAS
        terms = [ go[x] for x in args.terms ] 
        terms = filter(
            lambda x: len(x.loci) >= args.min_term_size or len(x.loci) <= args.max_term_size,
            terms
        )

    results = []
    for i,term in enumerate(terms):
        cob.log('On term {}',i)   
        window_size = args.candidate_window_size
        flank_limit =  args.candidate_flank_limit
        # Generate a series of densities for parameters
        num_genes = len(term.loci)
        if num_genes == 0:
            continue
        eloci = term.effective_loci(
            window_size=window_size
        )
        if args.percent_fcr != None and args.percent_fcr > 0:                                            
            # replace some loci with random genes if FDR specified              
            num_fcr = int(len(eloci) * args.percent_fcr)                        
            fcr_loci = cob.refgen.random_genes(num_fcr,window=window_size) 
            # permute and truncate the loci then add fcr loci                   
            eloci = np.concatenate([                                            
                np.random.permutation(eloci)[0:-1*len(fcr_loci)],                   
                np.array(list(fcr_loci))
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
            # Density
            density = cob.density(
                candidates
            )
            trans_density = cob.trans_locus_density(
                eloci,
                flank_limit=flank_limit
            )
            density_bs = np.array([
                cob.trans_locus_density(
                    eloci,
                    flank_limit=flank_limit,
                    iter_name = x,
                    bootstrap=True
                ) for x in range(args.num_bootstraps)
            ])
            density_pval = sum([x>=trans_density for x in density_bs])/args.num_bootstraps
            # Locality
            locality = cob.locality(
                candidates, include_regression=True
            ).resid.mean()
            trans_locality = cob.trans_locus_locality(
               eloci,
               flank_limit=flank_limit,
               include_regression=True
            ).resid.mean()
            locality_bs = np.array([
                cob.trans_locus_locality(
                    eloci,
                    flank_limit=flank_limit,
                    iter_name=x,
                    bootstrap=True,
                    include_regression=True
                ).resid.mean() for x in range(args.num_bootstraps)    
            ])
            locality_pval = sum([x>=trans_locality for x in locality_bs])/args.num_bootstraps
        results.append(OrderedDict([
           ('COB',cob.name),
           ('GO',term.id),
           ('NumGenes',num_genes),
           ('WindowSize',window_size),
           ('FlankLimit',flank_limit),
           ('FCR',args.percent_fcr),
           ('NumEffective',len(eloci)),
           ('NumCandidates',len(candidates)),
           ('Density',density),
           ('TransDensity',trans_density),
           ('Density_Pval',density_pval),
           ('Locality',locality),
           ('TransLocality',trans_locality),
           ('Locality_Pval',locality_pval)
        ]))
    results = pd.DataFrame.from_dict(results)
    results.to_csv(args.out,sep='\t',index=False)
