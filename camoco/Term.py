#!/usr/bin/python3

import numpy as np

from .Tools import log

class Term(object):
    '''
        A Term is a just named group of loci that are related.

        Parameters
        ----------
        id : unique identifier 
        desc: short description
        loci : iterable of loci that are related
        ** kwargs : dictionary of other term attributes

        Returns
        -------
        A Term Object

    '''
    def __init__(self, id, desc='', loci=None, **kwargs):
        self.id = id
        self.name = id
        self.desc = desc
        self.attrs = {}
        self.loci = set()
        if loci:
            self.loci = set(loci)
        for key, val in kwargs.items():
            self.attrs[key] = val

    #@property
    #def name(self):
    #    return self.id

    @property
    def locus_list(self): #pragma: no cover
        raise Exception('This is deprecated')

    def __len__(self):
        '''
            Returns the number of loci in the term.
        '''
        return len(self.loci)

    def __getitem__(self,key):
        return self.attrs[key]

    def add_locus(self, locus):
        '''
            Adds a locus to the Term.
        '''
        self.loci.add(locus)

    def flanking_loci(self, locus, window_size=100000):
        '''
            returns any nearby Term SNPs to a locus
        '''
        return [flank for flank in self.loci if abs(locus-flank) <= window_size]

    def copy(self,id=None,desc='',loci=None,**kwargs):
        '''
            Creates a copy of a term with the option to 
            expand loci and attrs. 

            Parameters
            ----------
            name : str
                A required name for the new term.
            desc : str
                An optional short description for the term.
            loci : iterable of co.Loci objects
                These loci will be added to the Term object
                in addition to the loci objects that were
                in the original Term.
            **kwargs : key value pairs
                Additional key value pairs will be added 
                as attributes to the term object.

            Returns
            -------
            A Term object.
        '''
        if id == None:
            id = self.id
        if loci == None:
            loci = set()
        loci = self.loci.union(loci)
        new_attrs = self.attrs.copy()
        new_attrs.update(**kwargs)
        copy = Term(
            id,
            desc=desc,
            loci=loci,
            **new_attrs
        )

        return copy
    
    def effective_loci(self, window_size=None):
        '''
            Collapse down loci that have overlapping windows into
            'effective' loci. Looks like:

            Locus1:    |--------o-------|
            Locus2:        |--------o--------|
            Locus3:                         |--------o--------|
            Effective: |--------o---+----------------o--------|

            Legend: '|' : Window edge, used to collapse
                    'o' : 'Locus' edge (SNPs in this case)
                    '+' : Sub loci, kept for downstream analysis

            Parameters
            ----------
            window_size : int (default: None)
                If not None, maps a new window size to each locus.      
        '''
        loci = sorted(self.loci)
        if window_size is not None:
            for locus in loci:
                locus.window = window_size
        collapsed = [loci.pop(0)]
        for locus in loci:
            # if they have overlapping windows, collapse
            if locus in collapsed[-1]:
                # Collapse if the windows overlap
                collapsed[-1] = collapsed[-1] + locus
            else:
                collapsed.append(locus)
        log('{}: Found {} SNPs -> {} effective SNPs with window size {} bp', 
            self.id, len(self.loci), len(collapsed), window_size
        )
        return collapsed

    def strongest_loci(self, attr, window_size=None,lowest=True):
        '''
            Collapses down loci that have overlapping windows,
            then returns the locus with the strongest 'attr'
            value. Looks like:

            Locus1:    |--------o-------| (attr: 7)
            Locus2:        |--------o--------| (attr: 2)
            Locus3:                             |--------o--------| (attr: 8)
            Strongest: |--------o-------|       |--------o--------|

            Legend: '|' : Window edge, used to collapse
                    'o' : 'Locus' edge (SNPs in this case)

            Parameters
            ----------
            attr : str
                The locus attribute to use to determine the 'strongest'
            window_size : int (default: None)
                If not None, maps a new window size to each locus.      
            lowest: bool (default: True)
                When sorting by attr, lowest is strongest (i.e. p-vals) 
        '''
        is_reverse = not lowest
        return [
            # sort by attr and take first item
            sorted(
                locus.sub_loci,
                key=lambda x: float(x.default_getitem(attr,np.inf)),
                reverse=is_reverse
            )[0] for locus in self.effective_loci(window_size=window_size)
        ]

    def __str__(self):
        return "Term: {}, Desc: {}, {} Loci".format(self.id, self.desc, len(self))

    def __repr__(self):
        return str(self.id)


