import camoco as co

def list_command(args):
    print(co.available_datasets(args.type,args.name))
