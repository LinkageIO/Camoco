#!/usr/bin/python

from COB import *

# Only in extreme cases
cb = COBDatabase()



# Import 5b Genes 
genes = pd.read_table("genes/ZmB73_5b_FGS_info.txt")
genes.chromosome.replace({'chrMt':'chr11'},inplace=True)
genes.chromosome.replace({'chrPt':'chr12'},inplace=True)
genes.chromosome = genes.chromosome.apply(lambda x : int(x.replace('chr','')))
# Filter out genes to only contain the longest transcript
genes = genes.groupby('gene_id').apply(lambda a : a.ix[(a.transcript_end - a.transcript_start).idxmax()])

genes[['chromosome','gene_id','transcript_start','transcript_end','transcript_strand']].drop_duplicates().apply(
    lambda x: cb.add_gene(x.gene_id, x.chromosome, x.transcript_start, x.transcript_end, x.transcript_strand, '5b'),
    axis = 1
)


DevelRNASEQ = pd.read_table("")

