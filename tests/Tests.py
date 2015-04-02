#!/usr/bin/python3

import unittest
import os

import camoco as co
import pandas as pd
from camoco.Config import cf

# Set the basedir to the testdir
cf.set('options','basedir', cf.get('options','testdir'))

# write test case to import refgen from GFF

class LocusBaseTestCase(unittest.TestCase):
    def test_locus_initialization(self):
        # numeric chromosomes
        a = co.Locus(1,500)
        self.assertIsInstance(a,co.Locus)

class RefGenBaseTestCase(unittest.TestCase):
    def test_gff_creation(self):
        gff = os.path.join(cf.get('options','testdir'),'raw','TAIR10_GFF3_genes.gff')
        co.del_dataset('RefGen','T10',safe=False)
        T10 = co.RefGen.from_gff(gff,'T10','Tair 10','10','Arabidopsis')
        self.assertIsInstance(T10,co.RefGen)

    def test_Zm_Gff_creation(self):
        gff = os.path.join(cf.get('options','testdir'),'raw','ZmB73_5b_FGS.gff')
        co.del_dataset('RefGen','Zm5bFGS',safe=False)
        ZM = co.RefGen.from_gff(gff,'Zm5bFGS','Maize 5b Filtered Gene Set','5b','Zea Mays')
        self.assertIsInstance(ZM,co.RefGen)

class OntologyCreationTestCase(unittest.TestCase):
    def test_Ionome_import(self):
        csv = os.path.join(cf.get('options','testdir'),'raw','sigGWASsnpsCombinedIterations.longhorn.allLoc.csv')
        ZM = co.RefGen('Zm5bFGS')
        df = pd.DataFrame.from_csv(csv,index_col=None)
        IONS  = co.Ontology.from_DataFrame(df,'ZmIonome','Maize Ionome',ZM,term_col='el',chr_col='chr',pos_col='pos');


class COBCreationTestCase(unittest.TestCase):
    def setUp(self):
        pass

class COBBaseTestCase(unittest.TestCase):
    def setUp(self):
        pass

#   def test_coex_table_length(self):
#       pass
#       self.assertEqual(len(self.cob.coex), comb(self.cob.num_genes(),2), 
#               "coex dimentions mismatch")

#   def test_coexpress_concordance(self):
#       pass
#       for a,b in itertools.combinations([self.refgen.random_gene() for x in range(50)],2):
#           self.assertTrue(abs(self.coexpression(a,b).score - self._coex_concordance(a,b)) < 0.001)
#           dis_dif = abs(self.coexpression(a,b).distance - abs(a-b)) 
#           self.assertTrue( np.isnan(dis_dif) or dis_dif < 0.001)


if __name__ == '__main__':
    unittest.main()
