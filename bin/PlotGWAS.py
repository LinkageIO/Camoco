#!/usr/bin/env python3

import argparse
import sys

from collections import defaultdict

import camoco as co
import matplotlib.pylab as plt

from camoco.Config import cf
cf.logging.log_level = 'quiet'

def main(args):
    # snag the appropriate COB
    cob = co.COB(args.cob)
    # snag the GWAS object
    gwas = co.GWAS(args.GWAS)
    parent_refgen = co.RefGen(args.RefGen)

    if 'all' in args.terms:
        terms = gwas.iter_terms()
    else:
        terms = [gwas[term] for term in args.terms]

    # Make a plot for each Term
    for term in terms:
        loci = list(term.loci) 
        # create a dictionary of Loci which we can refer to using ids
        locus_lookup = {x.id:x for x in loci}
        effective_loci = sorted(term.effective_loci(
            window_size=args.candidate_window_size
        ))

        # Each chromosome gets a plot
        chroms = set([x.chrom for x in loci])

        # Create a figure with a subplot for each chromosome 
        f, axes = plt.subplots(len(chroms),figsize=(10,4*len(chroms)))

        # iterate over Loci
        seen_chroms = set()
        voffset = 1 # Vertical Offset
        current_axis = 0
        for i,locus in enumerate(effective_loci):
            hoffset = -1 * locus.window
            # Reset the temp variables in necessary
            if locus.chrom not in seen_chroms:
                # Plot the chromosome
                seen_chroms.add(locus.chrom)
                current_axis = len(seen_chroms)-1
                voffset = 1
            # Get current axis
            cax = axes[current_axis]
            # Set up labels if first time one axis
            if voffset == 1:
                cax.set_ylabel('Chrom: '+ locus.chrom)
                cax.set_xlabel('Loci')
                cax.get_yaxis().set_ticks([])
                #cax.get_xaxis().set_ticks([])
            # shortcut for current axis
            cax.hold(True)

            # Plot the Sub Loci
            for id in locus.sub_loci:
                if id in locus_lookup:
                    sub_locus = locus_lookup[id]
                    cax.scatter(
                        hoffset+locus.window+abs(sub_locus.start-locus.start),
                        voffset,
                        marker='.',
                        color='blue'
                    )

            # Plot ALL Genes
            for gene in parent_refgen.candidate_genes(
                    locus,
                    gene_limit = 10e10
                ):
                cax.barh(
                    bottom=voffset,
                    width = len(gene),
                    height= 5,
                    left = hoffset+gene.start-locus.start+locus.window,
                    color='black'
                )
            # Plot the candidate genes
            for gene in cob.refgen.candidate_genes(
                    locus,
                    gene_limit = 10e10
                ):
                cax.barh(
                    bottom=voffset,
                    width = len(gene),
                    height= 5,
                    left = hoffset+gene.start-locus.start+locus.window,
                    color='grey'
                )
    
            # Plot the candidate genes
            for gene in cob.refgen.candidate_genes(
                    locus,
                    gene_limit=args.candidate_gene_limit
                ):
                cax.barh(
                    bottom=voffset,
                    width = len(gene),
                    height= 5,
                    left = hoffset+gene.start-locus.start+locus.window,
                    color='red'
                )

            # Plot the Effective Locus
            cax.scatter(hoffset,voffset,marker='>') # Upstream
            cax.scatter(hoffset+locus.window,voffset,marker='.',color='blue') # Start
            cax.scatter(hoffset+locus.window+len(locus),voffset,marker='.',color='blue') # Stop
            cax.scatter(hoffset+locus.window+len(locus)+locus.window,voffset,marker='<') # Downstream


            # place a block for interlocal distance
            voffset += 10
        plt.savefig(args.out)



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
        '--RefGen',
        help='The parent refgen shared by COB and GWAS data.'
    )
    parser.add_argument(
        '--terms',
        nargs='*', 
        help='The term within the GWAS ontology to use. default: all',
        default=['all']
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




    from IPython.core import ultratb
    sys.excepthook = ultratb.FormattedTB(mode='Verbose',
                 color_scheme='Linux', call_pdb=1)

    args = parser.parse_args()
    sys.exit(main(args))
