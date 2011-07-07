from sympy.matrices.matrices import Matrix
from sympy.matrices.dokmatrix import DOKMatrix
from sympy.matrices.lilmatrix import LILMatrix

def _from_callable(*args):
    rows = args[0]
    cols = args[1]
    op = args[2]
    mat = {}
    for i in xrange(rows):
        for j in xrange(cols):
            val = op(i, j)
            if val != 0:
                mat[i, j] = val
    return rows, cols, mat

def _from_list(*args):
    rows = args[0]
    cols = args[1]
    inp = args[2]
    mat = {}
    for i in xrange(rows):
        for j in xrange(cols):
            val = mat[i * rows + j]
            if val != 0:
                mat[i, j] = val
    return rows, cols, mat

def _dict_to_densematrix(rows, cols, mat):
    return Matrix._from_dict(rows, cols, mat)

def _dict_to_dokmatrix(rows, cols, mat):
    return DOKMatrix._from_dict(rows, cols, mat)

def _dict_to_lilmatrix(rows, cols, mat):
    return LILMatrix._from_dict(rows, cols, mat)

def slice2bounds(self, key, defmax):
        """
        Takes slice or number and returns (min,max) for iteration
        Takes a default maxval to deal with the slice ':' which is (none, none)
        """
        if isinstance(key, slice):
            lo, hi = 0, defmax
            if key.start is not None:
                if key.start >= 0:
                    lo = key.start
                else:
                    lo = defmax+key.start
            if key.stop is not None:
                if key.stop >= 0:
                    hi = key.stop
                else:
                    hi = defmax+key.stop
            return lo, hi
        elif isinstance(key, int):
            if key >= 0:
                return key, key+1
            else:
                return defmax+key, defmax+key+1
        else:
            raise IndexError("Improper index type")

def vecs2matrix(vecs, repr='dense'):
    """Join column vectors to a matrix."""
    m = len(vecs[0])
    n = len(vecs)
    if repr == 'dok':
        A = DOKMatrix.zeros((m, n))
    elif repr == 'lil':
        A = LILMatrix.zeros((m, n))
    elif repr == 'dense':
        A = Matrix.zeros((m, n))
    else:
        raise Exception('repr not recognized')
    for i in xrange(m):
        for j in xrange(n):
            A[i,j] = vecs[j][i,0]
    return A

