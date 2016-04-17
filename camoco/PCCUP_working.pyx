import numpy as np
from scipy.misc import comb 

cdef extern from "math.h":
    bint isnan(double x)
    double sqrt(double x)

# input is a typed numpy memoryview (::1 means c contiguous array)
def pair_correlation(double[:, ::1] x):
    # Define a new memoryview on an empty gene X gene matrix
    cdef double[:] pccs = np.empty(comb(x.shape[0],2,exact=True))
    cdef double u, v
    cdef int i, j, k, count
    cdef double du, dv, d, n, r
    cdef double sum_u, sum_v, sum_u2, sum_v2, sum_uv
    cdef long index

    index = 0
    for i in range(x.shape[0]):
        for j in range(i+1, x.shape[0]):
            sum_u = sum_v = sum_u2 = sum_v2 = sum_uv = 0.0
            count = 0            
            for k in range(x.shape[1]):
                u = x[i, k]
                v = x[j, k]
                # skips if u or v are nans
                if u == u and v == v:
                    sum_u += u
                    sum_v += v
                    sum_u2 += u*u
                    sum_v2 += v*v
                    sum_uv += u*v
                    count += 1
            if count < 10:
                pccs[index] = np.nan
            else:
                um = sum_u / count
                vm = sum_v / count
                n = sum_uv - sum_u * vm - sum_v * um + um * vm * count
                du = sqrt(sum_u2 - 2 * sum_u * um + um * um * count) 
                dv = sqrt(sum_v2 - 2 * sum_v * vm + vm * vm * count)
                if (du * dv) == 0:
                    pccs[index] = np.nan
                else:
                    r = 1 - n / (du * dv)
                    pccs[index] = r
            index += 1
    # Return the base of the memory view
    return pccs.base

def coex_index(long[:] ids, int mi):
    '''
        Camoco stores the coexpression matrix in long form. This is 
        space efficient, but accessing elements in an [i,j] format 
        requires additional overhead. This function takes in the original
        indicies for the [i,j] matrix and returns the long form indices for each 
        pairwise combinations of ids.

        Parameters
        ----------
        ids : array of indices 
            gene indices from the Expr matrix
        mi : int
            The total number of genes in the Expr matrix

        Returns
        -------
        An array of indices you can extract from the coex table

    '''
    cdef long[::] indices = np.empty(comb(ids.shape[0],2,exact=True),dtype=np.long)
    cdef long count = 0
   
    for ix in range(ids.shape[0]):
        for jx in range(ix+1,ids.shape[0]):
            i = min(ids[ix],ids[jx])
            j = max(ids[ix],ids[jx])
            # Calculate what the index would be if it were a square matrix
            #k = ((i * mi) + j) 
            # Calculate the number of cells in the lower diagonal
            #ld = (((i+1)**2) - (i+1))/2
            # Calculate the number of items on diagonal
            #d = i + 1
            #indices[count] = k-ld-d
            indices[count] = square_to_vector(i,j,mi)
            count += 1
    return indices.base 

def square_to_vector(long i, long j, mi):
    '''
        Convert an index from its square form
        to its vector form
    '''
    k = ((i * mi) + j) 
    # Calculate the number of cells in the lower diagonal
    ld = (((i+1)**2) - (i+1))/2
    # Calculate the number of items on diagonal
    d = i + 1
    return k-ld-d

def coex_neighbors(long id, int mi):
    '''
        Calculate the indices for the neighbors of id.
        This is better by example. If i == 4, what are 
        the indices for the neighbors (designated as '^')
        in the picture. Diagonal is '-'.
        Looks like:

         0123456789   
         
      0  -000^00000
      1  0-00^00000
      2  00-0^00000
      3  000-^00000
      4  0000-^^^^^
      5  00000-0000
      6  000000-000
      7  0000000-00
      8  00000000-0
      9  000000000-

    '''

    cdef long[::] indices = np.empty(mi-1,dtype=np.long)
    cdef long count = 0
    cdef long pivot

    for i in range(id):
        indices[count] = square_to_vector(i,id,mi)
        count += 1
    pivot = square_to_vector(id,id+1,mi)
    for j in range(id+1,mi):
        indices[count] = pivot
        pivot += 1
        count += 1
    return indices.base
