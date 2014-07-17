import os
from itertools import chain
import re
from Locus import *

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
def Bootrap(snps,networks,reference,ontologies):
    pass 



