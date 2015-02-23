#!/usr/bin/python3
from collections import defaultdict
import re

class Locus(object):
    def __init__(self, chrom, start, end, id=None ,gene_build=None, organism=None, window=100000):
        self.id = id
        self.gene_build = gene_build
        self.organism = organism
        self.chrom = chrom
        try:
            # Sometimes we have underflows, so choose 0 in that case
            self.start = max(0,int(start))
        except TypeError as e:
            self.start = None
        try:
            self.end = int(end)
        except TypeError as e:
            self.end = None
        # Warn us if something seems off
        if self.start > self.end:
            raise ValueError("Wonky start and stop positions for: {}".format(self))


    @property
    def stop(self):
        ''' because im an idiot ''' 
        return self.end

    @property
    def coor(self):
        ''' returns a tuple with start and stop '''
        return (self.start,self.stop)

    def __eq__(self,locus):
        if (self.chrom == locus.chrom and
            self.start == locus.start and
            self.end == locus.end and
            self.gene_build == locus.gene_build and
            self.organism == locus.organism):
            return True
        else:
            return False

    def __contains__(self,locus):
        if (locus.chrom == self.chrom and
               (( locus.start >= self.start and locus.start <= self.end)
               or(locus.end   <= self.end   and locus.end   >= self.start)
            )):
            return True
        else:
            return False
    def __len__(self):
        ''' inclusive length of locus '''
        return self.stop - self.start + 1

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
        else:
            return self.start - other.start
    def __str__(self):
        return '''
            organism:{} 
            type: {}
            id: {}
            chromosome: {}
            start: {}
            end: {}'''.format(self.organism,self.__class__,self.id,self.chrom,self.start,self.stop)
    def __repr__(self):
        return str(self)
    def __hash__(self):
        try:
            return int("{}{}".format(self.chrom,self.start))
        except ValueError as e:
            return int("-1{}".format(abs(self.start)))

class SNP(Locus):
    def __init__(self, chrom, pos, id=None ,gene_build='5b', organism='Zea',window=0):
        super().__init__(
            chrom=str(chrom),
            start=int(pos)-int(window/2),
            end=int(pos)+int(window/2),
            id=id,
            gene_build=gene_build,
            organism=organism
        )
        self.pos = int(pos)
        self.window = window
        # This is the number of genes a snp tags
        self.tag_genes = 4

    def collapse(self,snp):
        ''' collapse two snps into a new 'meta' SNP. The pos is half way
        between the two snps and the window extends within as well as 1/2
        upstream and downstream of the original window sizes '''
        # must be on the same chromosome to collapse
        assert(self-snp < float('Inf'))
        new_pos = int((self.pos + snp.pos)/2)
        new_window = max(self.end,snp.end) - min(self.start,snp.start)
        new_id = str(self.id)+'-'+str(snp.id)
        new_tag_genes = self.tag_genes + snp.tag_genes
        return SNP(self.chrom,new_pos,window=new_window)

    def __str__(self):
        return '''
            organism:{} 
            type: {}
            id: {}
            chromosome: {}
            pos: {}
            window: {}'''.format(self.organism,self.__class__,self.id,self.chrom,self.pos,self.window)


    def summary(self):
        return "S{}:{}".format(self.chrom,self.start)

    @classmethod
    def from_str(cls,string,regex='.*(\d+)_(\d+).*'):
        m = re.search(regex,string)
        if m is not None:
            return cls(m.group(1),m.group(2))
       

class Gene(Locus):
    def __init__(self,chrom=None,start=None,end=None,strand=None,id=None,build='5b',organism='Zea'):
        super().__init__(chrom,start,end,id,build,organism) 
        self.strand = strand
        self.attributes = dict()
    def __str__(self):
        return self.id
    def __repr__(self):
        return str(self)
    def summary(self):
        return "Gene-{}-{}:{}-{}".format(self.id,self.chrom,self.start,self.end)

    def to_gtf(self):
        return "\t".join(
            [self.chrom,'Camoco','Gene',self.start,self.end,".",self.strand,'0',
            "ID={};Organism={};Build={}".format(self.id, self.organism, self.build)]
        )

    @classmethod
    def from_gtf(cls,line,attribute_id="ID"):
        # GTF fields MUST be tab seperated as per Ensemble rules
        seqname,source,feature,start,end,score,strand,frame,attributes = line.strip().split()
        attributes = attributes.replace('=',';').split(";")
        attributes = dict(zip(attributes[::2],attributes[1::2]))
        if 'build 'not in attributes:
            attributes['build'] = None
        if 'organism' not in attributes:
            attributes['organism'] = None
        self = cls(chrom=seqname,start=start,end=end,strand=strand,id=attributes[attribute_id],build=attributes['build'],organism=attributes['organism'])
        self.attributes = attributes
        return self

class QTL(Locus):
    def __init__(self,chrom,start,end,id=None,build='5b',organism='Zea'):
        if id == None:
            self.id = "QTL-chr{}:{}-{}".format(chrom,start,end)
        else:
            self.id = id
        super().__init__(str(chrom),start,end,self.id)
    def __str__(self):
        return "{}".format(self.id)


# -----------------------End of Class Definitions -------------------------------------------

# Useful abstractions for analyses
def bootstrap(loci,genome):
    ''' Returns a non-overlapping, randoms set of equal 
        size loci from the genome '''
    random_loci = [genome.rLocus(len(locus)) for locus in loci]
    if non_overlapping(random_loci):
        return random_loci
    else:
        return bootstrap(loci,genome)

def num_overlap(loci1,loci2,unique=False):
    ''' Returns the number of loci in list 1 which overlap list 2 
        If unique is True, function return the number of unique loci 
        overlapping between lists '''
    total = defaultdict(int)
    for x in loci1:
        for y in loci2:
            if x in y:
                total[str(y)] += 1
    if unique:
        return len(total)
    else:
        return sum(total.values())

def non_overlapping(loci_list):
    ''' Return True if loci do not overlap with one another
        otherwise returns False if any two loci overlap '''
    if num_overlap(loci_list,loci_list) > len(loci_list):
        return False
    else:
        return True    

def significant_overlap(loci1,loci2,genome,iterations=1000,unique=False):
    ''''Returns the number of bootstrapped overlaps larger than or equal to loci1 in loci2
        for i iterations '''
    emp_overlap = num_overlap(loci1,loci2,unique=unique)
    bs = [ num_overlap(loci1,bootstrap(loci2,genome),unique=unique) for x in range(iterations)]
    return sum([x >= emp_overlap for x in bs])/iterations
