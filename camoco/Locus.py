#!/usr/bin/python3
from collections import defaultdict
from itertools import chain

import hashlib
import re

class Locus(object):
    def __init__(self, chrom, start, end=None, id=None, window=0, sub_loci=None, **kwargs):
        if id is None or id.startswith('<None>'):
            self._id = None
        else:
            self._id = id
        self.chrom = chrom
        self._start = int(start)
        self._end = int(end) if end is not None else int(start)
        self.window = int(window)
        self.attr = kwargs
        self.sub_loci = set(sub_loci) if sub_loci is not None else set()
        if len(self.sub_loci) == 0:
            self.sub_loci.add(self)
        #  Warn us if something seems off
        if self._start > self._end:
            raise ValueError("Wonky start and stop positions for: {}".format(self))

    def as_dict(self):
        a_dict = {
            'name'  : self.name,
            'chrom' : self.chrom,
            'start' : self.start,
            'end'   : self.end
        }
        a_dict.update(self.attr)
        return a_dict

    @property
    def id(self):
        if self._id is None:
            return '''<{}>{}:{}-{}'''.format(
                self._id, self.chrom,
                self.start, self.end
            )
        else:
            return self._id

    def update(self,dict):
        '''
            updates the attr attribute with values from the dictionary
        '''
        self.attr.update(dict)
        return self


    def as_record(self):
        return (self.chrom,self.start,self.end,self.name,self.window,self.id)

    @classmethod
    def from_record(cls,tpl):
        return cls(*tpl)

    def __setitem__(self,key,val):
        self.attr[key] = val

    def __getitem__(self,key):
        return self.attr[key]

    def default_getitem(self,key,default=None):
        ''' 
            Return a default value if the attr[key] value is None
        '''
        if self[key] is None:
            return default
        else:
            return self[key]

    @property
    def start(self):
        return max(0,int(self._start))

    @property
    def end(self):
        return int(self._end)

    @property
    def coor(self):
        ''' returns a tuple with start and stop '''
        return (self.start,self.end)

    @property
    def upstream(self):
        return self.start - self.window

    @property
    def downstream(self):
        return self.end + self.window

    @property
    def name(self):
        return self.id
        # Make sure this doesn't break anything
        # return "{}:{}:{}:{}".format(self.name,self.chrom,self.start,self.end)

    def __add__(self,locus):
        ''' 
            collapse two loci into a new 'meta' locus. The start is the min
            between the two loci and the window extends within as well as 1/2
            upstream and downstream of the original window sizes
        '''
        # must be on the same chromosome to collapse
        if self-locus == float('Inf'):
            raise TypeError('Loci must be on same chromosome to collapse.')
        new_start = int(max(0,min(self.start,locus.start)))
        new_end = int(max(self.end,locus.end))
        new_window = self.window
        # This can be a list, since set gets called anyways
        new_sub_loci = self.sub_loci | locus.sub_loci

        return Locus(
            self.chrom, new_start, new_end,
            window=new_window, sub_loci=new_sub_loci
        )

    def __eq__(self,locus):
        if (self.chrom == locus.chrom and
            self.start == locus.start and
            self.end == locus.end):
            return True
        else:
            return False

    def __contains__(self,locus):
        if (locus.chrom == self.chrom and
               # The locus has as 'start' position within the Locus window
               (( locus.upstream >= self.upstream and locus.upstream <= self.downstream)
               # The locus has an 'end' position within the Locus window
               or(locus.downstream <= self.downstream and locus.downstream >= self.upstream)
            )):
            return True
        else:
            return False


    def __len__(self):
        ''' inclusive length of locus '''
        if self.start == self.end:
            return 1
        else:
            return self.end - self.start

    def __cmp__(self,locus):
        if self.chrom == locus.chrom:
            return self.start - locus.start
        elif self.chrom > locus.chrom:
            return 1
        else:
            return -1

    def __lt__(self,locus):
        if self.chrom == locus.chrom:
            return self.start < locus.start
        else:
            return self.chrom < locus.chrom

    def __gt__(self,locus):
        if self.chrom == locus.chrom:
            return self.start > locus.start
        else:
            return self.chrom > locus.chrom

    def __sub__(self,other):
        if self.chrom != other.chrom:
            return float('Inf')
        if self == other:
            return 0
        else:
            # sort them
            a,b = sorted([self,other])
            return b.start - a.end

    def __str__(self):
        return '''<{}>{}:{}-{}+{}({})'''.format(
            self._id, self.chrom,
            self.start, self.end,
            self.window, len(self.sub_loci)-1
        )

    def summary(self):
        return '\n'.join([
            'ID: {}',
            'Chromosome: {}',
            'Start Position: {}',
            'End Position: {}',
            'Window Size: {}',
            'Sub Loci: {}'
        ]).format(
            self._id, self.chrom,
            self.start, self.end,
            self.window, len(self.sub_loci)-1
        )

    def __repr__(self):
        return str(self)

    def __hash__(self):
        digest = hashlib.md5(
            str.encode(str(self))
        ).hexdigest()
        return int(digest,base=16)

class Gene(Locus):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
