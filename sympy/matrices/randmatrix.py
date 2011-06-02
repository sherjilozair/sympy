from matrices import Matrix
from dokmatrix import DOKMatrix
from sympy.polys.specialpolys import random_poly
import random

"""
Code
01: Dense Matrix,   unsymmetric,    default sparsity 0.9,       default domain ZZ,  using Matrix class,     not guranteed to be invertible, but most probably is.
02: Dense Matrix,   symmetric,      default sparsity 0.9,                           using Matrix class,     not guranteed to be invertible, but most probably is.

03: Sparse Matrix,  unsymmetric,    default sparsity 0.1,                           using DOKMatrix class,  not guaranteed to be invertible,
04: Sparse Matrix,  symmetric,      default sparsity 0.1,                           using DOKMatrix class,  not guaranteed to be invertible,
05: Sparse Matrix,  unsymmetric,    default sparsity 1/sqrt(n),                     using DOKMatrix class,  guaranteed to be invertible,
06: Sparse Matrix,  symmetric,      default sparsity 1/sqrt(n),                     using DOKMatrix class,  guaranteed to be invertible,

"""

def randMatrix(m, n, min=0, max=99, seed=None, sparsity=1, symmetric=False, invertible=True, **kwargs):
    """Create random matrix m x n"""
    if seed == []:
        prng = random.Random()  # use system time
    else:
        prng = random.Random(seed)
    if sparsity > 0.8:
        rM = Matrix(m, n, lambda i,j: prng.randint(min,max) if prng.random() < sparsity else 0)
    else:
        rM = DOKMatrix(m, n, lambda i,j: prng.randint(min,max) if prng.random() < sparsity else 0)

    return rM
