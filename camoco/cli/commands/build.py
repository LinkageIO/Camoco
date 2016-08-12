import camoco as co
import pandas as pd
from camoco.Tools import DummyRefGen
from camoco.Locus import Gene

def build_cob(args):
    # Build the refgen
    refgen = co.RefGen(args.refgen)
    # Check that the sep is likely right.
    if len(pd.read_table(args.filename,sep=args.sep).columns) == 1:
        print(
            ("Detected only 1 column in {}, are you sure "
            "colunms are separated by '{}'?").format(args.filename,args.sep)
        )
        return None
    if args.allow_non_membership:
        refgen = refgen.copy(
            '{}_tmp'.format(refgen.name), 
            'temp refgen'.format(refgen.name)
        )
        # Add non membership genes
        for gid in pd.read_table(args.filename,sep=args.sep).index:
            refgen.add_gene(Gene(None,None,id=gid))

    quality_control = False if args.skip_quality_control else True
    normalize = False if args.skip_normalization else True
        
    # Basically just pass all the CLI arguments to the COB class method  
    cob = co.COB.from_table(
        args.filename,
        args.name,
        args.description,
        refgen,
        # Optional arguments
        sep=args.sep,
        rawtype=args.rawtype,
        # Data Processing
        quality_control=quality_control,
        normalization=normalize,
        quantile=args.quantile,
        # Data processing parameters
        max_gene_missing_data=args.max_gene_missing_data,
        max_accession_missing_data=args.max_accession_missing_data,
        min_single_sample_expr=args.min_single_sample_expr,
        min_expr=args.min_expr,
        max_val=args.max_val,
        dry_run=args.dry_run
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
    if len(df.columns) == 1:
        raise ValueError("Only 1 column found, check --sep, see --help")
    print('Loading {}'.format(args.refgen))
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
