import sys
import os
import math

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
        if os.path.exists("{}".format(args.out)) and not args.force:
            print(
                "{} exists! Skipping!".format(
                    args.out.replace('.csv','')
                )
            )
            return None
    # Load camoco objects
    go = co.GOnt(args.GOnt)
    cob = co.COB(args.cob)
    if 'all' in args.terms:
        terms = list(go.iter_terms(
            min_term_size=args.min_term_size,
            max_term_size=args.max_term_size
        ))
    elif os.path.exists(args.terms[0]):
        terms = [go[x.strip()] for x in open(args.terms[0]).readlines()]
    else:
        # Otherwise get the term out of the GWAS
        terms = [ go[x] for x in args.terms ] 
    cob.log("Simulating GWAS for {} GO Terms",len(terms))
    cob.log("All terms are between {} and {} 'SNPs'", args.min_term_size, args.max_term_size)

    results = []
    for i,term in enumerate(terms):
        cob.log('-'*75)
        window_size = args.candidate_window_size
        flank_limit =  args.candidate_flank_limit
        # Generate a series of densities for parameters
        num_genes = len([x for x in term.loci if x in cob])
        eloci = [ x for x  in term.effective_loci(
            window_size=window_size
        ) if x in cob ]
        cob.log(
            'Simulation {}: {} ({}/{} genes in {})',
            i,term.id,len(eloci),num_genes,cob.name
        )   
        if num_genes > args.max_term_size:
            cob.log("Too many genes... skipping")
            continue
        elif num_genes < args.min_term_size:
            cob.log("Too few genes... skipping")
            continue
        elif num_genes == 0:
            continue
        # Remove a percentage of SNPs to simulate false negatives
        if args.percent_SNPs_missed != None and args.percent_SNPs_missed > 0:
            # Calulate the index needed to hit percent missing
            missing_index = math.ceil(len(eloci) * (1-(args.percent_SNPs_missed/100)))
            if missing_index < 2:
                missing_index = 2 
            new_eloci = np.random.permutation(eloci)[0:missing_index]
            cob.log('Simulating {}% of SNPs missed by GWAS ({} SNPs -> {})',args.percent_SNPs_missed,len(eloci),len(new_eloci))
            eloci = new_eloci
        # Replace a percentage of SNPs with false positives
        if args.percent_fcr != None and args.percent_fcr > 0:                                            
            # replace some loci with random genes if FDR specified              
            num_fcr = math.ceil(len(eloci) * (args.percent_fcr/101))               
            fcr_loci = cob.refgen.random_genes(num_fcr,window=window_size) 
            cob.log('Simulating {}% of SNPs as false positive -> adding {} SNPs',args.percent_fcr,len(fcr_loci))
            # permute and truncate the loci then add fcr loci                   
            eloci = np.concatenate([                                            
                #np.random.permutation(eloci)[0:-1*len(fcr_loci)],                   
                eloci,
                np.array(list(fcr_loci))
            ])   
        candidates = cob.refgen.candidate_genes(
            eloci,
            flank_limit=flank_limit
        )
        cob.log(
            "SNP to gene mapping finds {} genes at window:{} bp, "
            "flanking:{} genes", len(candidates),
            args.candidate_window_size,
            args.candidate_flank_limit
        )
        if len(candidates) == 0:
            density = np.NaN
            trans_density = np.NaN
            density_pval = np.NaN
            locality = np.NaN
            trans_locality = np.NaN
            locality_pval = np.NaN
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
           ('TermSize',len(term)),
           ('NumRealGenes',num_genes),
           ('WindowSize',window_size),
           ('FlankLimit',flank_limit),
           ('FCR',args.percent_fcr),
           ('MCR',args.percent_SNPs_missed),
           ('NumEffective',len(eloci)),
           ('NumCandidates',len(candidates)),
           ('Density',density),
           ('TransDensity',trans_density),
           ('Density_pval',density_pval),
           ('Locality',locality),
           ('TransLocality',trans_locality),
           ('Locality_pval',locality_pval)
        ]))
    results = pd.DataFrame.from_dict(results)
    results.to_csv(args.out,sep='\t',index=False)
