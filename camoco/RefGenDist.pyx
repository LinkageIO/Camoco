import numpy as np
from scipy.misc import comb

def gene_distances(double[:] chr, long[:] start, long[:] end):
    # Create an array to put the results in
    cdef double[:] distances = np.empty(comb(chr.shape[0],2,exact=True))
    # to remember which permutation we are one
    cdef long i, j, index

    print("Calculating for {} genes".format(len(chr)))

    # Loop through genes and calcualate distances
    index = 0
    for i in range(chr.shape[0]):
        for j in range(i+1,chr.shape[0]):
            # We Cant compare genes on different chromosomes
            if chr[i] != chr[j]:
                distances[index] = np.inf 
            elif start[i] < start[j]: # i is upstream of j
                distances[index] = float(start[j] - end[i])
            else: # i is upstream of j
                distances[index] = float(start[i] - end[j])
            index += 1
    assert index == len(distances)
    return distances.base

