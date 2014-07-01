from HapMap import HapMap
from RefGen import RefGen,RefGenBuilder
from Ontology import Ontology
from Locus import *
from COB import COB,COBBuilder

TEST_LOCI_NON_OVERLAP = [Locus(1,1,100),Locus(1,101,200),Locus(1,201,300),Locus(2,1,100)]
assert non_overlapping(TEST_LOCI_NON_OVERLAP) == True
TEST_LOCI_OVERLAP = [Locus(1,1,100),Locus(1,50,150)]
assert non_overlapping(TEST_LOCI_OVERLAP) == False
