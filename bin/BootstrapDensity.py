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

def main(args):
    
    cob = co.COB(args.cob)
    ont = co.GWAS(args.GWAS)


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
            effective_loci, gene_limit=args.candidate_gene_limit
        )
   
        print("Num SNPs:\t{}".format(len(term.locus_list)),file=sys.stderr)
        print("Num Effective SNPs:\t{}".format(len(effective_loci)),file=sys.stderr)
        print('Num Candidate Genes:\t{}'.format(len(candidates)),file=sys.stderr)

        results['NumSNPs'] = len(term.locus_list)
        results['NumEffective'] = len(effective_loci)
        results['NumCandidates'] = len(candidates)

        emp = cob.trans_locus_density(
            effective_loci,
            gene_limit=args.candidate_gene_limit,
            bootstrap=False
        )

        results['WindowSize'] = args.candidate_window_size
        results['CandidateLimit'] = args.candidate_gene_limit

        results['EmpDensity'] = emp
    
        print("Empirical density Z-score for {}:\t{}".format(term.id,emp),file=sys.stderr)
    
        if args.min_score is not None and emp < args.min_score:
            print("Empirical below Z-Score",file=sys.stderr)
            results['NumBootstraps'] = np.nan
            results['PValue'] = np.nan
            results['BootstrapMean'] = np.nan
            results['BootstrapStd'] = np.nan
        else:

            bootstraps = np.array([
                cob.trans_locus_density(
                    term.effective_loci(window_size=args.candidate_window_size),
                    gene_limit=args.candidate_gene_limit,
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


if __name__ == '__main__':
    parser = argparse.ArgumentParser()    
    # Data Set Arguments
    parser.add_argument(
        '--cob',
        help='The camoco network to use.'
    )
    parser.add_argument(
        '--GWAS', 
        help='The camoco GWAS to use.'
    )
    parser.add_argument(
        '--terms',
        nargs='*', 
        help='The term within the ontology to use. default: all',
        default=['all']
    )
    parser.add_argument(
        '--min-score', 
        help=('Term with emp. less than this will be ignored '
              'If None, bootstraps will be calculated regardless. '
              'Default:None'),
        type=float, default=None
    )
    parser.add_argument(
       '--num-bootstraps', default=50,type=int, 
       help=('The number of bootstraps to perform in order '
           'to estimate a null distribution. default: 50')
    )
    parser.add_argument(
       '--candidate-window-size',default=50000,
       type=int,
       help=('The window size (in bp) for mapping effective SNPs to genes. '
             'default: 50000')
    )
    parser.add_argument(
       '--candidate-gene-limit',default=2,
       type=int,
       help=('The number of flanking genes included in SNP to gene mapping. '
           'default: 2' )
    )
    parser.add_argument(
        '--out', default=sys.stdout,
        help='OutputFile Name (default: Standard Out)'
    )

    import sys
    from IPython.core import ultratb
    sys.excepthook = ultratb.FormattedTB(mode='Verbose',
                 color_scheme='Linux', call_pdb=1)

    args = parser.parse_args()
    if args.out is not sys.stdout:
        args.out = open(args.out,'w')
    sys.exit(main(args))
