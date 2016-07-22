import camoco as co
import pandas as pd
from camoco.Tools import DummyRefGen

def build_cob(args):
    # Build the refgen
    if not args.ignore_refgen_membership:
        refgen = co.RefGen(args.refgen)
    else:
        refgen = DummyRefGen() 
    # Basically just pass all the CLI arguments to the COB class method  
    cob = co.COB.from_table(
        args.filename,
        args.name,
        args.description,
        refgen,
        # Optional arguments
        sep=args.sep,
        rawtype=args.rawtype,
        quantile=False,
        max_gene_missing_data=args.max_gene_missing_data,
        max_accession_missing_data=args.max_accession_missing_data,
        min_single_sample_expr=args.min_single_sample_expr,
        min_expr=args.min_expr,
        max_val=args.max_val,
        dry_run=args.dry_run,
        quality_control=args.skip_quality_control
    )
    print("Build successful!")
    print(cob.summary())

def build_refgen(args):
    co.RefGen.from_gff(
        args.filename, 
        args.name,
        args.description,
        args.build,
        args.organism,
        # Optional Arguments
        chrom_feature=args.chrom_feature,
        gene_feature=args.gene_feature,
        ID_attr=args.ID_attr,
        attr_split=args.attr_split
    )
    print("Build successful!")

def build_gont(args):
    refgen = co.RefGen(args.refgen)
    co.GOnt.from_obo(
        args.obo_filename,
        args.filename,
        args.name,
        args.description,
        refgen
    )
    print('Build Successful')

def build_GWAS(args):
    df = pd.DataFrame.from_csv(args.filename,sep=args.sep).reset_index()
    refgen = co.RefGen(args.refgen)
    gwas = co.GWAS.from_DataFrame(
        df,
        args.name,
        args.description,
        refgen,
        term_col=args.trait_col,
        chr_col=args.chrom_col,
        pos_col=args.pos_col
    )
    print("Build Successful:")
    print(gwas.summary())
