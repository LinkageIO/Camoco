import pytest
import camoco as co
from camoco import cf

from camoco import Locus

'''
    Unit tests
'''

def test_from_ids(testRefGen):
    random_genes = sorted(testRefGen.random_genes(n=10))
    from_ids = sorted(testRefGen.from_ids([x.id for x in random_genes]))
    assert set(random_genes) == set(from_ids)

def test_get_item(testRefGen):
    random_gene = testRefGen.random_gene()
    assert random_gene == testRefGen[random_gene.id]

def test_get_items_from_list(testRefGen):
    random_genes = sorted(testRefGen.random_genes(n=10))
    from_ids = sorted(testRefGen[[x.id for x in random_genes]])
    assert set(random_genes) == set(from_ids)

def test_lowercase_get_item(testRefGen):
    random_gene = testRefGen.random_gene()
    random_id = random_gene.id
    # Stupid mutability
    random_id.lower()
    assert random_gene == testRefGen[random_id]

def test_genes_within(testRefGen):
    random_gene = testRefGen.random_gene()
    bigger_locus = Locus(
        random_gene.chrom,
        start=random_gene.start-100,
        end=random_gene.end+100
    )
    genes = testRefGen.genes_within(bigger_locus)
    assert random_gene in genes

def test_locus_not_in_upstream_downstream(testRefGen):
    '''
        Upstream and downstream should not include the gene of interest.
    '''
    random_gene = testRefGen.random_gene()
    upstream = testRefGen.upstream_genes(
        random_gene,window_size=50e5,gene_limit=5
    )
    downstream = testRefGen.downstream_genes(
        random_gene,gene_limit=5,window_size=50e5
    )
    assert random_gene not in upstream
    assert random_gene not in downstream

def test_upstream_downstream_genes(testRefGen):
    '''
        Take downstream of genes, then upstream genes of the 
        last gene in downstream. Tests that returns the same interval 
    '''
    # Grab downstream genes of random genes
    random_gene = testRefGen.random_gene()
    # Grab 10 downstream genes
    downstream_genes = testRefGen.downstream_genes(
        random_gene,gene_limit=11,window_size=50e10
    )
    assert len(downstream_genes) == 11
    # grab last gene
    last_gene = downstream_genes.pop(-1)
    # Grab upstream genes
    upstream_genes = testRefGen.upstream_genes(
        last_gene,gene_limit=10,window_size=50e10
    )
    assert sorted(downstream_genes) == sorted(upstream_genes)

def test_flanking_genes(testRefGen):
    random_gene = testRefGen.random_gene()
    downstream = testRefGen.downstream_genes(
        random_gene, window_size=50e6, gene_limit=5
    )
    upstream = testRefGen.upstream_genes(
        random_gene, window_size=50e6, gene_limit=5
    )
    flanking = testRefGen.flanking_genes(
        random_gene, window_size=50e6, flank_limit=5
    )
    assert sorted(flanking) == sorted(upstream + downstream)

def test_flanking_genes_includes_within_genes_for_SNPS(testRefGen):
    random_gene = testRefGen.random_gene()
    # test snp
    test_snp = Locus(random_gene.chrom,random_gene.start,window=50e5)
    flanking = testRefGen.flanking_genes(test_snp)
    assert random_gene not in flanking

def test_candidate_genes_from_SNP(testRefGen):
    random_gene = testRefGen.random_gene()
    # grab a bunch of downstream genes
    down1,down2 = testRefGen.downstream_genes(
        random_gene,gene_limit=2,window_size=50e6
    )
    # Create a Locus that is on gene 5
    test_snp = Locus(
        down1.chrom,
        down1.start-50,
        end=down2.end+50,
        window=50e6
    )
    candidates = testRefGen.candidate_genes(
        test_snp,flank_limit=5,chain=False
    )
    assert len(candidates) == 12 

def test_candidate_genes_from_gene_includes_gene(testRefGen):
    random_gene = testRefGen.random_gene()
    # grab a bunch of downstream genes
    downstream = testRefGen.downstream_genes(
        random_gene,gene_limit=10,window_size=50e6
    )
    # Create a Locus that is on gene 5
    candidates = testRefGen.candidate_genes(
        downstream[5],flank_limit=10,window_size=50e6
    )
    assert downstream[4] in candidates

def test_non_chained_candidates(testRefGen):
    random_genes = testRefGen.random_genes(n=10)
    # Create a Locus that is on gene 5
    candidates = testRefGen.candidate_genes(
        random_genes,flank_limit=10,window_size=50e6,chain=False
    )
    # test that we got candidates for each random locus
    assert len(candidates) == len(random_genes)
   

def test_flank_limit_for_candidate_genes(testRefGen):
    random_gene = testRefGen.random_gene()
    # Create a Locus that is on gene 5
    candidates = testRefGen.candidate_genes(
        random_gene,flank_limit=5,window_size=50e6,chain=True
    )
    assert len(candidates) == 11

def test_flank_limit_for_candidate_genes_from_SNP(testRefGen):
    random_gene = testRefGen.random_gene()
    downstream = testRefGen.downstream_genes(
        random_gene,gene_limit=10,window_size=50e6
    )
    test_snp = Locus(downstream[5].chrom,downstream[5].start,window=50e6)
    # Create a Locus that is on gene 5
    candidates = testRefGen.candidate_genes(
        test_snp,flank_limit=5,window_size=50e6
    )
    assert len(candidates) == 11

def test_bootstrap_candidate_length_equal_from_SNP(testRefGen):
    random_gene = testRefGen.random_gene()
    test_snp = Locus(random_gene.chrom,random_gene.start,window=50e6)
    candidates = testRefGen.candidate_genes(test_snp)
    bootstraps = testRefGen.bootstrap_candidate_genes(test_snp)
    assert len(candidates) == len(bootstraps)

def test_bootstrap_candidate_length_equal_from_gene(testRefGen):
    random_gene = testRefGen.random_gene()
    candidates = testRefGen.candidate_genes(random_gene,window_size=5e10)
    bootstraps = testRefGen.bootstrap_candidate_genes(random_gene,window_size=5e10)
    assert len(candidates) == len(bootstraps)

def test_refgen_length(testRefGen):
    # grab length from sqlite 
    from_sql = testRefGen.db.cursor().execute('''
        SELECT COUNT(*) FROM genes;
    ''').fetchone()[0]
    assert from_sql == len(testRefGen)

def test_filtered_refgen(testRefGen):
    co.del_dataset('RefGen','test_filtered_refgen',safe=False) 
    random_genes = set(testRefGen.random_genes(n=500))
    test_filtered_refgen = testRefGen.filtered_refgen(
        'test_filtered_refgen',
        'test. please ignore',
        testRefGen,
        random_genes
    )
    assert len(test_filtered_refgen) == len(random_genes)
    for x in random_genes:
        assert x in test_filtered_refgen
    co.del_dataset('RefGen','test_filtered_refgen',safe=False) 

def test_rowid_equals_1_after_refgen_rebuild(Zm5bFGS_duplicate):
    '''
        This was a regression bug where when a refgen was rebuilt
        the rowid was not reset resulting in weird random_gene 
        method which relies on rowid 
    '''
    assert Zm5bFGS_duplicate\
        .db.cursor().execute(
            "SELECT MIN(rowid) from genes"
        ).fetchone()[0] == 1


def test_random_genes_returns_correct_n(testRefGen):
    assert len(testRefGen.random_genes(n=cf.test.num)) == cf.test.num
