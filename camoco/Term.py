
from .Tools import log

class Term(object):
    '''
        A Term is a just group of loci that are related.
    '''
    def __init__(self, id, desc='', locus_list=None, **kwargs):
        self.id = id
        self.desc = desc
        self.attr = {}
        self.locus_list = set()
        if locus_list:
            self.locus_list = set(locus_list)
        for key, val in kwargs.items():
            self.attrs[key] = val

    def __len__(self):
        '''
            Returns the number of loci in the term.
        '''
        return len(self.locus_list)

    def add_locus(self, locus):
        '''
            Adds a locus to the Term.
        '''
        self.locus_list.add(locus)

    def flanking_loci(self, gene, window_size=100000):
        '''
            returns any nearby Term SNPs to a gene
        '''
        return [locus for locus in self.locus_list if abs(gene-locus) <= window_size]

    def effective_loci(self, window_size=None, max_genes_between=1):
        '''
            Collapse down loci that have overlapping windows.
            Also collapses down snps that have
        '''
        locus_list = sorted(self.locus_list)
        if window_size is not None:
            for locus in locus_list:
                locus.window = window_size
        collapsed = [locus_list.pop(0)]
        for locus in locus_list:
            # if they have overlapping windows, collapse
            if locus in collapsed[-1]:
                # Collapse if the windows overlap
                collapsed[-1] = collapsed[-1] + locus
            else:
                collapsed.append(locus)
        log('{}: Found {} SNPs -> {} effective SNPs', self.id, len(self.locus_list), len(collapsed))
        return collapsed

    def __str__(self):
        return "Term: {}, Desc: {}, {} Loci".format(self.id, self.desc, len(self))

    def __repr__(self):
        return str(self.id)


