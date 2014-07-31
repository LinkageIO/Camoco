import os
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

def log(msg,*args,basedir="~/.camoco"):
    if os.path.exists(os.path.expanduser(os.path.join(basedir,'logs','logfile.txt'))):
        print("[CAMOCO]",time.ctime(), '-', self.name, '-', msg.format(*args),file=self.log_file)
    print("[CAMOCO]",time.ctime(), '-', self.name, '-', msg.format(*args),file=sys.stderr)

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
