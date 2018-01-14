import camoco as co
from camoco.Tools import del_dataset

def remove(args):
    del_dataset(args.type,args.name,force=args.force)
    print('Done')
