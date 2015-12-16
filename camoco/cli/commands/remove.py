import camoco as co

def remove(args):
    print(co.del_dataset(args.type,args.name,safe=args.force))
