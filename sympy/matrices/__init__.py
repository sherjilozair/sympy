"""A module that handles matrices.

Includes functions for fast creating matrices like zero, one/eye, random matrix etc.
"""
from matrices import (Matrix, zeros, ones, eye, diag,
     hessian, GramSchmidt, wronskian, casoratian,
     list2numpy, matrix2numpy, DeferredVector, block_diag, symarray, ShapeError,
     NonSquareMatrixError)

from randmatrix import randMatrix

from dokmatrix import *

from lilmatrix import *

from matrix_ import Matrix_
