#cython: boundscheck=False, wraparound=False, cdivision=True
import numpy as np
from scpiy.misc import comb

def gene_distances(double[:] chr, double[:] position):
    # loop
    cdef double[:] distances = np.empty(comb(chr.shape[0],2,exact=True))
    for i in range(chr.shape[0]):
        for j in range(i,chr.shape[0]):

