import camoco as co

def remove(args):
    co.del_dataset(args.type,args.name,safe=args.force)
    print('Done')
