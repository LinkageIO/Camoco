#!/usr/bin/env python3

import argparse
import sys

from collections import defaultdict

import camoco as co
import matplotlib.pylab as plt

from camoco.Config import cf

def commify(value):
    ' Return the value with commas delimiting the thousands '
    return '{:,}'.format(value)

def plot_gwas(args):
    # snag the appropriate COB
    cob = co.COB(args.cob)
    # snag the GWAS object
    gwas = co.GWAS(args.gwas)

    if 'all' in args.terms:
        terms = gwas.iter_terms()
    else:
        terms = [gwas[term] for term in args.terms]

    # Make a plot for each Term
    for term in terms:
        loci = list(term.loci) 
        # create a dictionary of Loci which we can refer to using ids
        locus_lookup = {x.id:x for x in loci}
        # Each chromosome gets a plot
        chroms = set([x.chrom for x in loci])

        # Create a figure with a subplot for each chromosome 
        f, axes = plt.subplots(
            len(chroms),
            figsize=(15,4*len(chroms))
        )
        plt.title('{} Term'.format(term.id))
        # Pull out the snp to gene mappings
        if args.snp2gene == 'effective':
            loci = sorted(term.effective_loci(
                window_size=args.candidate_window_size
            ))
        elif args.snp2gene == 'strongest':
            loci = term.strongest_loci(
                window_size=args.candidate_window_size,
                attr=args.strongest_attr,
                lowest=args.strongest_higher
            )
        else:
            raise ValueError('{} not valid snp2gene mapping'.format(args.snp2gene))

        # iterate over Loci
        seen_chroms = set()
        voffset = 1 # Vertical Offset
        current_axis = 0
        y_labels = []
        y_ticks = []
        for i,locus in enumerate(loci):
            hoffset = -1 * locus.window
            # Reset the temp variables in necessary
            if locus.chrom not in seen_chroms:
                seen_chroms.add(locus.chrom)
                current_axis = len(seen_chroms)-1
                voffset = 1
                if len(y_labels) > 0 and current_axis > 0:
                    # Set the old labels in the current 
                    axes[current_axis-1].set_yticks(y_ticks)
                    axes[current_axis-1].set_yticklabels(y_labels)
                y_labels = []
                y_ticks = []
            # Get current axis
            cax = axes[current_axis]
            # Set up labels if first time one axis
            if voffset == 1:
                cax.set_ylabel('Chrom: '+ locus.chrom)
                cax.set_xlabel('Loci')
            # shortcut for current axis
            cax.hold(True)

            # Plot ALL Genes
            for gene in gwas.refgen.candidate_genes(
                    locus,
                    flank_limit = 10e10
                ):
                cax.barh(
                    bottom=voffset,
                    width = len(gene),
                    height= 5,
                    zorder=1,
                    left = hoffset+gene.start-locus.start+locus.window,
                    label='RefGen Genes',
                    color='grey'
                )
            # Plot the candidate genes
            for gene in cob.refgen.candidate_genes(
                    locus,
                    flank_limit = 10e10
                ):
                cax.barh(
                    bottom=voffset,
                    width = len(gene),
                    height= 5,
                    zorder=1,
                    left = hoffset+gene.start-locus.start+locus.window,
                    label = 'Gene Passed QC',
                    color='green'
                )
    
            # Plot the candidate genes
            for gene in cob.refgen.candidate_genes(
                    locus,
                    flank_limit=args.candidate_flank_limit
                ):
                cax.barh(
                    bottom=voffset,
                    width = len(gene),
                    height= 5,
                    zorder=1,
                    left = hoffset+gene.start-locus.start+locus.window,
                    label='Candidate Gene',
                    color='red'
                )

            # Plot the Effective Locus
            cax.scatter( # Upstream
                hoffset,voffset,
                marker='>',zorder=2
            )
            cax.scatter( # Start
                hoffset+locus.window,voffset,
                marker='.',color='blue',zorder=2
            )
            cax.scatter( # Stop
                hoffset+locus.window+len(locus),
                voffset,marker='.',color='blue',zorder=2
            )
            cax.scatter( # Downstream
                hoffset+locus.window+len(locus)+locus.window,
                voffset,marker='<',zorder=2 
            )
            # Plot the Sub Loci
            for id in locus.sub_loci:
                if id in locus_lookup:
                    sub_locus = locus_lookup[id]
                    cax.scatter(
                        hoffset+locus.window+abs(sub_locus.start-locus.start),
                        voffset,
                        zorder=2,
                        marker='.',
                        label='SNP',
                        color='blue'
                    )

            # place a block for interlocal distance
            y_labels.append(commify(locus.start))
            y_ticks.append(voffset)
            voffset += 10
        # Have to finish off the ticks on the last chromosome
        axes[current_axis].set_yticks(y_ticks)
        axes[current_axis].set_yticklabels(y_labels)
        # Save Plot
        plt.savefig(args.out.replace('.png','_{}.png'.format(term.id)))
        plt.close()
