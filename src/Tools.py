import os

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
def flanking_genes(snps):
    return list(chain(*[list(chain(*ZM.flanking_genes(s,gene_limit=4,pos_limit=50000))) for s in snps]))
def heatmap(snps,network,ZM,filename):
    genes = ZM.from_ids(list(chain(
            *map(flanking_genes,map(snp_from_str,snps))
    )))
    network.heatmap(network.gene_expr_vals(genes),filename)


