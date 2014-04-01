#!/usr/bin/python

from COBLocus import Locus

class Gene(Locus):
    def __init__(self,id,gene_build='4a.53'):
        info = self.query('''SELECT chromosome, chromo_start, chromo_end FROM mzn_gene_loci 
            WHERE gene_id = {} and gene_build = '{}' ''',id,gene_build)
        super(Gene,self).__init__(info.iloc[0]['chromosome'],info.iloc[0]['chromo_start'],info.iloc[0]['chromo_end'])
        self.id = int(id)

    @property
    def gramene_id(self):
        return self.query("SELECT gene_name FROM mzn_gene_name WHERE gene_id = {}",self.id).iloc[0]['gene_name']

    @property
    def common_name(self):
        info = self.query("SELECT common_name FROM mzn_gene_common WHERE common_id = {}",self.id)
        return None if info.shape[0] == 0 else info.iloc[0]['common_name']

    @property
    def arab_ortho(self):
        info = self.query('''SELECT * FROM mzn_arab_orthologs orth
            LEFT JOIN mzn_arab_gene info ON info.arab_id = orth.arab_id 
            LEFT JOIN mzn_arab_gene_types type ON info.type_id = type.type_id
            LEFT JOIN mzn_arab_short_desc short ON info.short_id = short.short_id
            LEFT JOIN mzn_arab_curator_desc cur ON info.curated_id = cur.curator_id
            LEFT JOIN mzn_arab_comp_desc comp ON info.comp_id = comp.comp_id
            WHERE gene_id = {} ''',self.id)
        return info

    @property
    def go_terms(self):
        info = self.query('''SELECT term_name, term_short_desc, term_long_desc, space_desc
            FROM mzn_gene_go_terms
            LEFT JOIN mzn_go_terms ON go_id = term_id
            LEFT JOIN mzn_go_space ON term_space = space_id
            WHERE gene_id = {}
        ''',self.id)
        return info


