import pytest

def test_locus_initialization(simple_Locus):
    # numeric chromosomes
    self.assertIsInstance(a, co.Locus)


def test_candidate_vs_bootstrap_length():
   snps = self.Fe.effective_snps(window_size=50000)
   self.assertEqual(
       len(self.ZM.candidate_genes(snps)),
       len(self.ZM.bootstrap_candidate_genes(snps)),
       'Bootstrap vs emirical candidate gene sizes do not match'
   )


def test_get_attr_and_generating_from_ids():
   self.assertIsInstance(self.ZM[cf['test']['gene']], co.Locus)
   self.assertIsInstance(
       self.ZM.from_ids([cf['test']['gene']])[0],
       co.Locus
   )
