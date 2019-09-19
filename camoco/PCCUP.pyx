import numpy as np
cimport numpy as np
from scipy.special import comb

cimport cython

from libc.math cimport sqrt
#from libc.math cimport isnan
cdef extern from "numpy/npy_math.h" nogil:
    long double NAN "NPY_NAN"
    bint isnan "npy_isnan"(long double)


#cdef extern from "math.h":
#    bint isnan(double x)
#    double sqrt(double x)

def pcc2(double[:] x, double[:] y):
    '''
    This implements the PCC as defined by Wikipedia
    '''
    cdef long k
    cdef long double xi, yi, n
    cdef long double xsum, ysum
    cdef long double xbar, ybar
    cdef long double xsquare, ysquare
    cdef long double xy

    assert x.shape[0] == y.shape[0]
    n = xsum = ysum = xbar = ybar = xy = xsquare = ysquare = 0.0
    for k in range(x.shape[0]):
        xi = x[k]
        yi = y[k]
        # Only iterate over values that are not nan in 
        # both x AND y (nan == nan returns False for a nan value)
        if xi == xi and yi == yi:
            xy += (xi*yi)
            xsum += xi
            ysum += yi
            xsquare += xi**2
            ysquare += yi**2
            n += 1
    # Short circuit if there are not enough values
    if n < 10:
        r = np.nan
    else:
        xbar = xsum / n
        ybar = ysum / n

        numerator   = xy - n*(xbar*ybar)
        denominator = sqrt(xsquare - (n*(xbar**2))) * sqrt(ysquare - (n*(ybar**2)))
        if denominator == 0:
            r = np.nan
        else:
            r = numerator / denominator
            if abs(r) > 1.1:
                raise ValueError(f"PCC out of range")
    return r

def pcc(double[:] x, double[:] y):
    # ref: https://stackoverflow.com/questions/25077080/calculate-special-correlation-distance-matrix-faster?rq=1
    # NOTE: This returns the pcc + 1 because its an implementation of the pdist function
    cdef float u, v
    cdef int k, count
    cdef long num_rows, num_cols
    cdef double du, dv, d, n
    cdef long double sum_u, sum_v, sum_u2, sum_v2, sum_uv
    cdef long index
    cdef float r

    num_cols = x.shape[0]
    sum_u = sum_v = sum_u2 = sum_v2 = sum_uv = 0.0
    count = 0            
    for k in range(num_cols):
        u = x[k]
        v = y[k]
        # skips if u or v are nans
        if u == u and v == v:
            sum_u += u
            sum_v += v
            sum_u2 += u*u
            sum_v2 += v*v
            sum_uv += u*v
            count += 1
    if count < 10:
        r = np.nan
    else:
        um = sum_u / count
        vm = sum_v / count
        n = sum_uv - sum_u * vm - sum_v * um + um * vm * count
        du = sqrt(sum_u2 - 2 * sum_u * um + um * um * count) 
        dv = sqrt(sum_v2 - 2 * sum_v * vm + vm * vm * count)
        if (du * dv) == 0:
            r = np.nan
        else:
            r = 1 - n / (du * dv)
    return r


# input is a typed numpy memoryview (::1 means c contiguous array)
def pair_correlation(double[:, ::1] x):
    # Define a new memoryview on an empty gene X gene matrix
    cdef float[::1] pccs = np.empty(comb(x.shape[0],2,exact=True)).astype('float32')
    cdef long i, j
    cdef long num_rows
    cdef long index
    cdef float r

    index = 0
    num_rows = x.shape[0]
    for i in range(num_rows):
        for j in range(i+1, num_rows):
            r = pcc2(x[i,:],x[j,:])
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
    cdef long ix, jx, i, j
    cdef long num_rows

    num_rows = ids.shape[0]
    for ix in range(num_rows):
        for jx in range(ix+1,num_rows):
            i = min(ids[ix],ids[jx])
            j = max(ids[ix],ids[jx])
            # Calculate what the index would be if it were a square matrix
            indices[count] = square_to_vector(i,j,mi)
            count += 1
    return indices.base 

cdef square_to_vector(long i, long j, mi):
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

def coex_expr_index(long[:] ids, int num_genes):
    '''
        Convert a list of coex indexes to a tuple of expr indexes
    '''
    cdef int num_rows = ids.shape[0]
    coors = np.zeros([num_rows,2], dtype=np.int32)
    if num_rows == 0:
        return coors
    cdef long idx, pos, i, j
    idx = 0
    pos = 0
    
    for i in range(num_genes):
        if (ids[idx] < (pos + (num_genes - (i+1)))):
            for j in range(i+1, num_genes):
                if ids[idx] == pos:
                    coors[idx, 0] = i
                    coors[idx, 1] = j
                    idx += 1
                if idx >= num_rows:
                    break
                pos += 1
        else:
            pos += (num_genes - (i+1))
        
        if idx >= num_rows:
            break
    
    return coors

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
