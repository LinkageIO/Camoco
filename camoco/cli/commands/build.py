import camoco as co
import pandas as pd
import sys
from camoco.Tools import DummyRefGen,log,available_datasets
from camoco.Locus import Gene

def build_cob(args):
    try:
        # Build the refgen
        refgen = co.RefGen(args.refgen)
        # Check that the sep is likely right.
        if len(pd.read_table(args.filename,sep=args.sep).columns) == 1:
            print(
                ("Detected only 1 column in {}, are you sure "
                "colunms are separated by '{}'?").format(args.filename,args.sep)
            )
            return None
        elif len(pd.read_table(args.filename,sep=args.sep).columns) < 20 and args.non_interactive != True:
            print(
                ("Detected fewer than 20 accessions in the expression matrix. "
                 "Calculating co-expression with this many datapoints is not advised")
            )
            if input('are you sure you want to continue? [y/n]: ').upper() == 'Y':
                pass
            else:
                sys.exit(1)
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
        quantile = False if args.skip_quantile else True

        # Check to see if this dataset is already built
        if available_datasets('Expr',args.name):
            print('Warning! This dataset has already been built.')
            co.Tools.del_dataset('Expr',args.name,force=args.force)
            
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
            quantile=quantile,
            # Data processing parameters
            max_gene_missing_data=args.max_gene_missing_data,
            max_accession_missing_data=args.max_accession_missing_data,
            min_single_sample_expr=args.min_single_sample_expr,
            min_expr=args.min_expr,
            max_val=args.max_val,
            dry_run=args.dry_run,
            zscore_cutoff=args.zscore_cutoff,
            index_col=args.index_col
        )
        print(cob.summary())
    except Exception as e:
        print(f"Build failed for {args.name}. Rolling back: removing corrupted files...")
        co.Tools.del_dataset('Expr',args.name,force=True)
        raise e

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
    # Check to see if this dataset is already built
    if available_datasets('GOnt',args.name):
        print('Warning! This dataset has already been built.')
        co.Tools.del_dataset('GOnt',args.name,force=args.force)
 
    go = co.GOnt.from_obo(
        args.obo_filename,
        args.filename,
        args.name,
        args.description,
        refgen,
        go_col=args.go_col,
        id_col=args.id_col,
        headers=args.gene_term_header
    )
    print("Done: {}".format(go.summary()))
    print('Build Successful')

def build_GWAS(args):
    df = pd.DataFrame.from_csv(args.filename,sep=args.sep).reset_index()
    if len(df.columns) == 1:
        raise ValueError("Only 1 column found, check --sep, see --help")
    # Add the trait if its missing
    if args.add_trait is not None:
        df[args.trait_col] = args.add_trait
    if args.query is not None:
        df = df.query(args.query)
    print('Loading {}'.format(args.refgen))
    refgen = co.RefGen(args.refgen)
    # Filter out traits that are in args.skip_trait
    df = df[[x not in args.skip_traits for x in df[args.trait_col]]]
    # Build
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
