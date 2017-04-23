import numpy as np

def squareform_diagfill(arr1D):
    n = int(np.sqrt(arr1D.size*2))
    if (n*(n+1))//2!=arr1D.size:
        print("Size of 1D array not suitable for creating a symmetric 2D array!")
        return None
    else:
        R,C = np.triu_indices(n)
        out = np.zeros((n,n),dtype=arr1D.dtype)
        out[R,C] = arr1D
        out[C,R] = arr1D
    return out
