import sys
import os

import pandas as pd
import numpy as np
import scipy as sp
import camoco as co

from itertools import chain

def snp2gene(args):
    '''
        Perform SNP (locus) to candidate gene mapping
    '''

    if args.out != sys.stdout:
        # Create any non-existant directories
        if os.path.dirname(args.out) != '':
            os.makedirs(os.path.dirname(args.out),exist_ok=True)
        if os.path.exists(args.out) and args.overlook == True:
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
        if 'effective' in args.snp2gene:
            # Map to effective
            effective_loci = term.effective_loci(
                window_size=args.candidate_window_size
            )
        elif 'strongest' in args.snp2gene:
            effective_loci = term.strongest_loci(
                window_size=args.candidate_window_size,
                attr=args.strongest_attr,
                lowest=args.strongest_higher
            )
 
        genes = pd.DataFrame([ x.as_dict() for x in 
            cob.refgen.candidate_genes(
                effective_loci,
                flank_limit=args.candidate_flank_limit,
                attrs={'Term':term.id},
                include_parent_locus=True,
                include_num_siblings=True,
                include_num_intervening=True,
                include_rank_intervening=True,
                include_parent_attrs=args.include_parent_attrs
            )
        ])
       
        data = pd.concat([data,genes])
    data['COB'] = cob.name
    data['FlankLimit'] = args.candidate_flank_limit
    data['WindowSize'] = args.candidate_window_size

    # Add data from gene info files
    original_number_genes = len(data)
    for info_file in args.gene_info:
        cob.log('Adding info for {}',info_file)
        # Assume the file is a table
        info = pd.read_table(info_file,sep='\t')
        # try to match as many columns as possible
        data = pd.merge(data,info,how='left')
        if len(data) != original_number_genes:
            cob.log.warn(
                'Info file did not contain unique gene mappings, '
                'beware of duplicate candidate gene entries!'
            )
    
    # Generate the output file
    data.to_csv(args.out,index=None,sep='\t')
        
