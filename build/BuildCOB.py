#!/usr/bin/python
import pandas as pd
import itertools
import COB as cob
# Only in extreme cases
cb = cob.COBDatabaseBuilder()
#cb.__create_tables__()
#cb.log('Clearing Database')
#cb.clear_datasets()
#cb.log('Done')

# -----------------------
# Genotype Network PNAS
#Load the Genotype Datasets
#cb.log("Reading in Genotype Dataset")
#GenRaw = pd.read_csv("./datasets/CghFilteredCombined.txt",sep="\t")
#Geno = cob.COBDataset("Genotype","Rewiring of the Maize Transcriptome",FPKM=False,gene_build='4a.53')
#Geno.from_DataFrame(GenRaw)
#cb.log("Adding Dataset: ", Geno)
#cb.add_dataset(Geno)
#cb.log("Done")

#------------------------
# Devel Network RNASEQ
# Load the transcript expression data from csv
# Create a new DataSet
Devel = cob.COBDataset("DevelRNASEQ","Developmental Atlas RNASEQ dataset",FPKM=True,gene_build='5b')
# Import 5b Genes 
Devel.log("Reading in Zm5B Gene Table")
Zm5bGenes = pd.read_table("./genes/ZmB73_4a.53_FGS_info.txt")
Zm5bGenes.chromosome.replace({'chrMt':'chr11'},inplace=True)
Zm5bGenes.chromosome.replace({'chrPt':'chr12'},inplace=True)
Zm5bGenes.chromosome = Zm5bGenes.chromosome.apply(lambda x : int(x.replace('chr','')))
# Filter out Zm5bGenes to only contain the longest transcript
Zm5bGenes = Zm5bGenes.groupby('gene_id').apply(lambda a : a.ix[(a.transcript_end - a.transcript_start).idxmax()])
# Set the index to be transcript ID
Zm5bGenes.index = Zm5bGenes.transcript_id
# Import the Devel Transcript FPKM
Devel.log('Reading in Devel FPKM')
DvlRaw = pd.read_csv("./datasets/DevelAtlasRNASEQ.csv",sep="\t")
# Filter out transcript to only contain the longest one
DvlRaw = DvlRaw.ix[Zm5bGenes['transcript_id']]
# Convert transcript IDs to gene IDs
DvlRaw.index = Zm5bGenes.ix[DvlRaw.index]['gene_id']
# Filter out rows which are do not have FPKM > 0 in at least one tissue
DvlRaw = DvlRaw[DvlRaw.apply(lambda  x: any(x > 0),axis=1)]
# Fill on data from DataFrame
Devel.from_DataFrame(DvlRaw)
Devel.log("Adding Dataset:",Devel)
cb.add_dataset(Devel)

exit()

#------------------------
# Genotype Network ProbeBlast PNAS
RawProbeBlast = pd.read_table("./datasets/afterBlast_ProbeDesigns.csv")
# Filter out probes with more than 1 Blast hit
ProbeBlast = RawProbeBlast[RawProbeBlast.NUMBER_OF_BLAST_HITS == 1]
# Figure  out why there are only 54 inds and if it matters
PBRaw = pd.read_table("./datasets/GSE30036_series_matrix.txt",mangle_dupe_cols=False).groupby(level=0,axis=1).mean()


# Compare Old database with new database

# Import 4a.53 genes
Zm4aGenes = pd.read_table("./genes/ZmB73_4a.53_FGS_info.txt")
Zm4aGenes.chromosome.replace({'chrMt':'chr11'},inplace=True)
Zm4aGenes.chromosome.replace({'chrPt':'chr12'},inplace=True)
Zm4aGenes.chromosome = Zm4aGenes.chromosome.apply(lambda x : int(x.replace('chr','')))
# Filter out Zm4aGenes to only contain the longest transcript
Zm4aGenes = Zm4aGenes.groupby('gene_id').apply(lambda a : a.ix[(a.transcript_end - a.transcript_start).idxmax()])
# Set the index to be transcript ID
Zm4aGenes.index = Zm4aGenes.transcript_id

#Zm4aGenes[['chromosome','gene_id','transcript_start','transcript_end','transcript_strand']].drop_duplicates().apply(
#    lambda x: cb.add_gene(x.gene_id, x.chromosome, x.transcript_start, x.transcript_end, x.transcript_strand, '4a.53'),
#    axis = 1
#)


# Import 5b Genes 
Zm5bGenes = pd.read_table("./genes/ZmB73_4a.53_FGS_info.txt")
Zm5bGenes.chromosome.replace({'chrMt':'chr11'},inplace=True)
Zm5bGenes.chromosome.replace({'chrPt':'chr12'},inplace=True)
Zm5bGenes.chromosome = Zm5bGenes.chromosome.apply(lambda x : int(x.replace('chr','')))
# Filter out Zm5bGenes to only contain the longest transcript
Zm5bGenes = Zm5bGenes.groupby('gene_id').apply(lambda a : a.ix[(a.transcript_end - a.transcript_start).idxmax()])
# Set the index to be transcript ID
Zm5bGenes.index = Zm5bGenes.transcript_id

#Zm5bGenes[['chromosome','gene_id','transcript_start','transcript_end','transcript_strand']].drop_duplicates().apply(
#    lambda x: cb.add_gene(x.gene_id, x.chromosome, x.transcript_start, x.transcript_end, x.transcript_strand, '5b'),
#    axis = 1
#)


