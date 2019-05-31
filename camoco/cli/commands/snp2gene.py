import sys
import os

import pandas as pd
import numpy as np
import scipy as sp
import camoco as co

from itertools import chain

from camoco.Tools import log
# Initialize a new log object
log = log()

def snp2gene(args):
    '''
        Perform SNP (locus) to candidate gene mapping
    '''

    if args.out != sys.stdout:
        # Create any non-existant directories
        if os.path.dirname(args.out) != '':
            os.makedirs(os.path.dirname(args.out),exist_ok=True)
        if os.path.exists(args.out) and not args.force:
            print(
                "Output for {} exists! Skipping!".format(
                    args.out
                ),file=sys.stderr
            )
            return None

    # Set a flag saying this is from a COB refgen
    from_cob = False
    # Create the refgen (option to create it from a COB)
    if co.Tools.available_datasets('Expr',args.refgen):
        refgen = co.COB(args.refgen).refgen
        from_cob = args.refgen 
    elif co.Tools.available_datasets('RefGen',args.refgen):
        refgen = co.RefGen(args.refgen)
    # Create the GWAS object
    ont = co.GWAS(args.gwas)

    if 'all' in args.terms:
        terms = ont.iter_terms()
    else:
        terms = [ont[term] for term in args.terms]

    data = pd.DataFrame()
    results = []
    for term in terms:
        for window_size in args.candidate_window_size:
            for flank_limit in args.candidate_flank_limit:
                if 'effective' in args.snp2gene:
                    # Map to effective
                    effective_loci = term.effective_loci(
                        window_size=window_size
                    )
                elif 'strongest' in args.snp2gene:
                    effective_loci = term.strongest_loci(
                        window_size=window_size,
                        attr=args.strongest_attr,
                        lowest=args.strongest_higher
                    )
                genes = pd.DataFrame([ x.as_dict() for x in 
                    refgen.candidate_genes(
                        effective_loci,
                        flank_limit=flank_limit,
                        include_parent_locus=True,
                        include_num_siblings=True,
                        include_num_intervening=True,
                        include_rank_intervening=True,
                        include_SNP_distance=True,
                        include_parent_attrs=args.include_parent_attrs,
                        attrs={'Term':term.id},
                    )
                ])
                genes['FlankLimit'] = flank_limit
                genes['WindowSize'] = window_size
                genes['RefGen'] = refgen.name
                if from_cob != False:
                    genes['COB'] = from_cob
                data = pd.concat([data,genes])

    # Add data from gene info files
    original_number_genes = len(data)
    for info_file in args.gene_info:
        log('Adding info for {}',info_file)
        # Assume the file is a table
        info = pd.read_table(info_file,sep='\t')
        if len(info.columns) == 1:
            info = pd.read_table(info_file,sep=',')
        # try to match as many columns as possible
        matching_columns = set(data.columns).intersection(info.columns)
        log("Joining SNP2Gene mappings with info file on: {}",','.join(matching_columns))
        data = pd.merge(data,info,how='left')
        if len(data) != original_number_genes:
            log.warn(
                'There were multiple info rows for some genes. '
                'Beware of potential duplicate candidate gene entries! '
            )
    
    # Generate the output file
    data.to_csv(args.out,index=None,sep='\t')

    log("Summary stats")
    print('-'*100)
    #print('With {}kb windows and up to {} flanking genes'.format(int(args.candidate_window_size/1000),args.candidate_flank_limit))
    print("Mapped {} SNPs to {} genes".format(len(data.parent_locus.unique()),len(data.ID.unique())))
    print("Number of candidate genes per term:")
    print(data.groupby('Term').apply(lambda df: len(df.ID)))
