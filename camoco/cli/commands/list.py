import camoco as co

def list_command(args):
    if args.type and args.name:
        if args.type == 'GWAS':
            gwas = co.GWAS(args.name)
            print('\n'.join([x.id for x in gwas.iter_terms()]))
        elif args.type =='GOnt':
            gont = co.GOnt(args.name)
            print('\n'.join([x.id for x in gont.iter_terms()]))
    else:
        args.type = '%'
        args.name = '%'
        print(co.available_datasets(args.type,args.name).to_string())
