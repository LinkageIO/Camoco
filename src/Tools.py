import os
from itertools import chain

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
def flanking_genes(snps,ZM):
    return list(chain(*[list(chain(*ZM.flanking_genes(s,gene_limit=4,pos_limit=50000))) for s in snps]))
def Bootrap(snps,networks,reference,ontologies):
    pass 



TEST_GENES = ZM.from_ids('GRMZM2G024993','GRMZM2G138060','GRMZM2G015534','GRMZM2G348551')
