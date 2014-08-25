import os
import sys
import time
import re
import functools

from itertools import chain
from Locus import *

def memoize(obj):
    cache = obj.cache = {}

    @functools.wraps(obj)
    def memoizer(*args, **kwargs):
        # This wraps the calling of the memoized object
        key = str(args) + str(kwargs)
        if key not in cache:
            cache[key] = obj(*args, **kwargs)
        return cache[key]
    return memoizer

def log(msg,*args):
    print("[LOG]",time.ctime(), '-', msg.format(*args),file=sys.stderr)

def ext(filename):
    return os.path.join(os.path.expanduser("~/MetaboloCOB/"+filename))
def snp_from_str(strs):
    snps = []
    for s in strs:
        m = re.match("S(\d\d?)_(\d+)_.*",s)
        snps.append(SNP(m.group(1),m.group(2)))
    return snps
def B73_eq_Mo17(snp,HM):
    genotypes = HM.genotypes(snp,accessions=['B73','MO17'])
    if genotypes[0] == genotypes[1]:
        return True
    else:
        return False


def TermStats(term,cob_list,OUT):
    for cob in cob_list:
        # Create an output dir
        log("On Term {} with {}",term.id,cob.name)
        flanks = [cob.refgen.flanking_genes(x,gene_limit=4,pos_limit=50000) for x in term.snp_list]
        print("{}\t{}_NumSNPS\t{}".format(term.id,cob.name,len(term.snp_list)),file=OUT)
        print("{}\t{}_NumGenes\t{}".format(term.id,cob.name,len(list(chain(*flanks)))),file=OUT)
        density  = cob.density(chain(*flanks),min_distance=100000)
        locality = cob.locality(chain(*flanks),min_distance=100000)
        len_LCC = cob.lcc(list(chain(*flanks)),min_distance=100000)
        print("{}\t{}_TransDensity\t{}".format(term.id,cob.name,density),file=OUT)
        print("{}\t{}_Locality\t{}".format(term.id,cob.name,cob.locality(chain(*flanks))),file=OUT)
        print("{}\t{}_LCC\t{}".format(term.id,cob.name,len_LCC),file=OUT)
        if density > 2:
            log("Density > 2; Boostrapping!")
            # Calculate BootStrap 
            bs_density = []
            bs_local = []
            bs_lcc = []
        
            for x in range(50):
                bs_flanks = list(chain.from_iterable([cob.refgen.bootstrap_flanking_genes(x,gene_limit=4,pos_limit=50000) for x in term.snp_list]))
                bs_density.append(cob.density(bs_flanks,min_distance=100000))           
                bs_local.append(cob.locality(bs_flanks,min_distance=100000))
                bs_lcc.append(cob.lcc(bs_flanks,min_distance=100000))
            print("{}\t{}_BS_TransDensity\t{}".format(term.id,cob.name,sum([x >= density for x in bs_density])),file=OUT)
            print("{}\t{}_BS_Locality\t{}".format(term.id,cob.name,sum([x >= locality for x in bs_local])),file=OUT)
            print("{}\t{}_BS_LCC\t{}".format(term.id,cob.name,sum([x >= len_LCC for x in bs_lcc])),file=OUT)
            cob.heatmap(cob.gene_expression_vals(chain(*flanks)),filename="{}_{}_Heatmap.png".format(term.id,cob.name)) 
            cob.plot(list(chain(*flanks)),filename="{}_{}_Network.png".format(term.id,cob.name),layout='kk',height=800,width=800)
            
        
    
        
    
