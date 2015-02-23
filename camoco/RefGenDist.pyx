#cython: boundscheck=False, wraparound=False, cdivision=True
import numpy as np
from scipy.misc import comb

def gene_distances(double[:] chr, long[:] position):
    # Create an array to put the results in
    cdef double[:] distances = np.empty(comb(chr.shape[0],2,exact=True))
    # to remember which permutation we are one
    cdef long i, j, index
    cdef int debug

    print("Calculating for {} genes".format(len(position)))

    # Loop through genes and calcualate distances
    index = 0
    for i in range(chr.shape[0]):
        for j in range(i+1,chr.shape[0]):
            if chr[i] == chr[j]:
                distances[index] = np.nan 
            else:
                distances[index] = abs(position[i] - position[j])
            index += 1
    return distances.base

