from __future__ import division
import random
from sympy.matrices import Matrix
from sympy.printing import sstr, pretty
from sympy.simplify.simplify import simplify as sympy_simplify
from sympy import S

def _iszero(x):
    return x == 0

class LILMatrix(object):
    def __init__(self, *args, **kwargs):
        self.type = lambda i: i
        self._cache = {}
        if len(args) == 3 and callable(args[2]):
            op = args[2]
            if not isinstance(args[0], int) or not isinstance(args[1], int):
                raise TypeError("`args[0]` and `args[1]` must both be integers.")
            self.rows = args[0]
            self.cols = args[1]
            self.mat = [[] for i in xrange(self.rows)]
            for i in xrange(self.rows):
                for j in xrange(self.cols):
                    value = self.type(op(i,j))
                    if value != self.type(0):
                        self.mat[i].append((j, value))
        elif len(args)==3 and isinstance(args[0],int) and \
                isinstance(args[1],int) and isinstance(args[2], (list, tuple)):
            self.rows = args[0]
            self.cols = args[1]
            mat = args[2]
            self.mat = [[] for i in xrange(self.rows)]
            for i in range(self.rows):
                for j in range(self.cols):
                    value = self.type(mat[i*self.cols+j])
                    if value != self.type(0):
                        self.mat[i].append((j, value))
        elif len(args)==3 and isinstance(args[0],int) and \
                isinstance(args[1],int) and isinstance(args[2], dict):
            self.rows = args[0]
            self.cols = args[1]
            self.mat = [[] for i in xrange(self.rows)]
            # manual copy, copy.deepcopy() doesn't work
            for key in args[2].keys():
                val = args[2][key]
                if val != self.type(0):
                    self.mat[i].append((j, value))
        else:
            if len(args) == 1:
                mat = args[0]
            else:
                mat = args
            if not isinstance(mat[0], (list, tuple)):
                mat = [ [element] for element in mat ]
            self.rows = len(mat)
            self.cols = len(mat[0])
            self.mat = [[] for i in xrange(self.rows)]
            for i in range(self.rows):
                if len(mat[i]) != self.cols:
                    raise ValueError("All arguments must have the same length.")
                for j in range(self.cols):
                    value = self.type(mat[i][j])
                    if value != self.type(0):
                        self.mat[i].append((j, value))

    def __str__(self):
        return sstr(self.toMatrix())

    def __repr__(self):
        return sstr(self.toMatrix())

    def toMatrix(self):
        mat = [[self.type(0)] * self.cols for i in xrange(self.rows)]
        for i in xrange(self.rows):
            for j, val in self.mat[i]:
                mat[i][j] = val
        return Matrix(mat)

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

    def __getitem__(self, key):
        # if not (0, 0) <= key < (self.rows, self.cols):
        #     raise Exception
        i, j = key
        if type(i) is slice or type(j) is slice:
                return self.submatrix2(key)
        for j2, val in self.mat[i]:
            if j2 >= j:
                if j2 == j:
                    return val
                else:
                    return self.type(0)
        return self.type(0)

    def set_type(self, type):
        self.type = type
        self.applyfunc(type)

    def __setitem__(self, key, value):
        if not (0, 0) <= key < (self.rows, self.cols):
            raise Exception
        i, j = key
        for ind, (j2, val) in enumerate(self.mat[i]):
            if j2 >= j:
                if j2 == j:
                    if value == 0:
                        self.mat[i].pop(ind)
                    else:
                        self.mat[i][ind] = (j, value)
                    return
                else:
                    if value != 0:
                        self.mat[i].insert(ind, (j, value))
                    return
        if value != self.type(0):
            self.mat[i].append((j, value))

    def submatrix2(self, keys):
        if not isinstance(keys[0], slice) and not isinstance(keys[1], slice):
            raise TypeError("At least one element of `keys` must be a slice object.")

        rlo, rhi = self.slice2bounds(keys[0], self.rows)
        clo, chi = self.slice2bounds(keys[1], self.cols)
        if not ( 0<=rlo<=rhi<=self.rows and 0<=clo<=chi<=self.cols ):
            raise IndexError("Slice indices out of range: a[%s]"%repr(keys))
        outRows, outCols = rhi-rlo, chi-clo
        outMat = []
        for i in xrange(rlo, rhi):
            startset = False
            start = end = len(self.mat[i])
            for ind, (j, val) in enumerate(self.mat[i]):
                if j >= clo and not startset:
                    start = ind
                    startset = True
                if j >= chi:
                    end = ind
                    break
            outMat.append(self.mat[i][start:end])
        for i in xrange(len(outMat)):
            for ind, (j, val) in enumerate(outMat[i]):
                 outMat[i][ind] = (j - clo, val)
        A = LILMatrix.zeros((outRows, outCols))
        A.mat = outMat
        A.type = self.type
        return A

    def submatrix(self, keys):
        if not isinstance(keys[0], slice) and not isinstance(keys[1], slice):
            raise TypeError("At least one element of `keys` must be a slice object.")

        rlo, rhi = self.slice2bounds(keys[0], self.rows)
        clo, chi = self.slice2bounds(keys[1], self.cols)
        if not ( 0<=rlo<=rhi and 0<=clo<=chi ):
            raise IndexError("Slice indices out of range: a[%s]"%repr(keys))
        outRows, outCols = rhi-rlo, chi-clo
        outMat = []
        for i in xrange(rlo, rhi):
            startset = False
            start = 0
            end = len(self.mat[i])
            for ind, (j, val) in enumerate(self.mat[i]):
                if not startset and j >= clo:
                    start = ind
                    startset = True
                if j >= chi:
                    end = ind
                    break
            outMat.append(self.mat[i][start:end])
        A = LILMatrix.zeros((outRows, outCols))
        A.mat = outMat
        return A

    def copyin_matrix(self, key, value):
        rlo, rhi = self.slice2bounds(key[0], self.rows)
        clo, chi = self.slice2bounds(key[1], self.cols)
        if value.rows != rhi - rlo or value.cols != chi - clo:
            raise ShapeError("The Matrix `value` doesn't have the same dimensions " +
                "as the in sub-Matrix given by `key`.")
        raise NotImplemented

    def row_add(self, r1, r2, alpha):
        if r1 == r2:
            return
        row1 = self.mat[r1]
        row2 = self.mat[r2]
        self.mat[r1] = _row_add(row1, row2, alpha)

    def row_scale(self, r, alpha):
        for ind, (j, val) in enumerate(self.mat[r]):
            self.mat[r][ind] = (j, alpha * val)

    def row_add_bad(self, r1, r2, alpha):
        row2 = self.mat[r2]
        for elem in row2:
            self[r1,elem[0]] += alpha * elem[1]

    def row_functor(self, r, f):
        for i in xrange(len(self.mat[r])):
            self.mat[r][i] = self.mat[r][i][0],f(self.mat[r][i][1],self.mat[r][i][0])

    row = row_functor

    def row_swap(self, r1, r2):
        self.mat[r1], self.mat[r2] = self.mat[r2], self.mat[r1]


    @classmethod
    def zeros(cls, shape):
        if isinstance(shape, tuple):
            return cls(shape[0], shape[1], lambda i,j:0)
        else:
            return cls(shape, shape, lambda i, j: 0)

    def LUdecomposition_Simple(self, iszerofunc=_iszero):
        """
        Returns A comprised of L,U (L's diag entries are 1) and
        p which is the list of the row swaps (in order).
        """
        n = self.rows
        A = self[:,:]
        p = []
        # factorization
        for j in range(n):
            for i in range(j):
                for k in range(i):
                    A[i,j] = A[i,j] - A[i,k]*A[k,j]
            pivot = -1
            for i in range(j,n):
                for k in range(j):
                    A[i,j] = A[i,j] - A[i,k]*A[k,j]
                # find the first non-zero pivot, includes any expression
                if pivot == -1 and not iszerofunc(A[i,j]):
                    pivot = i
            if pivot < 0:
                # this result is based on iszerofunc's analysis of the possible pivots, so even though
                # the element may not be strictly zero, the supplied iszerofunc's evaluation gave True
                raise ValueError("No nonzero pivot found; inversion failed.")
            if pivot != j: # row must be swapped
                A.row_swap(pivot,j)
                p.append([pivot,j])
            scale = 1 / A[j,j]
            for i in range(j+1,n):
                A[i,j] = A[i,j] * scale
        return A, p

    def LUdecomposition(self, iszerofunc=_iszero):
        """
        Returns the decomposition LU and the row swaps p.

        Example:
        >>> from sympy import Matrix
        >>> a = Matrix([[4, 3], [6, 3]])
        >>> L, U, _ = a.LUdecomposition()
        >>> L
        [  1, 0]
        [3/2, 1]
        >>> U
        [4,    3]
        [0, -3/2]

        """
        combined, p = self.LUdecomposition_Simple(iszerofunc=_iszero)
        L = self.zeros(self.rows)
        U = self.zeros(self.rows)
        for i in range(self.rows):
            for j in range(self.rows):
                if i > j:
                    L[i,j] = combined[i,j]
                else:
                    if i == j:
                        L[i,i] = 1
                    U[i,j] = combined[i,j]
        return L, U, p

    def gauss_sparse(self):
        A = self[:, :]
        for i in xrange(A.rows):
            if A[i, i] == 0:
                print 'bad pivot, exchanging', i
                for k in xrange(i + 1, A.rows):
                    print '\tconsidering', k, i
                    if A[k, i] != 0:
                        print '\t\tfound', k, i
                        print pretty(A)
                        print

                        A.row_swap(k, i)

                        print pretty(A)
                        break
                if A[i, i] == 0:
                    print 'bad bad pivot', i
                    print pretty(A)
                    raise Exception
            ind = 0
            l = len(A.mat[i])
            while ind < l:
                j, v = A.mat[i][ind]
                if j < i:
                    print 'zeroing out', i, j, A.mat[i]
                    A.row_add(i, j, - v / A[j, j])
                    print 'zeroed out', i, j, A.mat[i]
                    l = len(A.mat[i])
                else:
                    break
        return A

    def gauss_col(self):
        "Gaussian elimnation, currently tested only on square matrices"
        A = self[:, :]
        for j in xrange(A.rows):
            rlist = A.nz_col_lower(j)
            if A[j, j] == self.type(0):
                if rlist:
                    A.row_swap(j, rlist[0])
                    rlist.pop(0)
                else:
                    continue
            for i in rlist:
                A.row_add(i, j, - A[i, j] / A[j, j])
        return A

    def _upper_triangular_solve(self, rhs):
        X = LILMatrix.zeros((rhs.rows, 1))
        for i in reversed(xrange(self.rows)):
            X[i, 0] = (rhs[i, 0] - sum(value * X[j, 0] for j, value in self.mat[i])) / self[i, i]
        return X

    def _lower_triangular_solve(self, rhs):
        X = LILMatrix.zeros((rhs.rows, 1))
        for i in xrange(self.rows):
            X[i, 0] = (rhs[i, 0] - sum(value * X[j, 0] for j, value in self.mat[i])) / self[i, i]
        return X

    def LUsolve(self, rhs):
        L, U, p = self.LU_sparse()
        rhs.permute(p)
        Y = L._lower_triangular_solve(rhs)
        return U._upper_triangular_solve(Y)

    def solve_gauss(self, rhs):
        big = self.join_rows(rhs)
        ref = big.gauss_col()
        U, b = ref[:, :self.cols], ref[:, self.cols:]
        return U._upper_triangular_solve(b)

    def solve_rref(self, rhs):
        big = self.join_rows(rhs)
        rref = big.rref()
        return rref[:, self.cols]

    def solve(self, rhs, method="GE"):
        if method == "GE":
            return self.solve_gauss(rhs)
        elif method == "RREF":
            return self.solve_rref(rhs)
        else:
            raise ValueError('Unrecognised method')

    def __mul__(self, other):
        if isinstance(other, LILMatrix):
            prod = self.toDOKMatrix() * other.toDOKMatrix()
            return LILMatrix(prod.rows, prod.cols, lambda i, j: prod[i, j]) 
        else:
            return scalar_product(self, other)

    def __rmul__(self, other):
        ## assume other is scalar
        return self.__mul__(other)

    def __sub__(self, other):
        return self + (-1 * other)

    def det(self):
        if not self.is_square():
            raise Exception
        ref = self.gauss_col()
        det = 1
        for i in xrange(ref.rows):
            det *= ref[i, i]
        return det
        
    def toDOKMatrix(self):
        from sympy import DOKMatrix
        Mat = DOKMatrix(self.rows, self.cols, {})
        for i in xrange(self.rows):
            for j, value in self.mat[i]:
                Mat[i, j] = value
        return Mat

    def rref2(self):
        pivot, r = 0, self[:,:]
        pivotlist = []
        for i in range(r.cols):
            if pivot == r.rows:
                break
            if r[pivot,i] == 0:
                for k in xrange(pivot + 1, r.rows):
                    if r[k, i] != 0:
                        r.row_swap(pivot, k)
                        break
                if r[pivot, i] == 0:
                    continue
            r.row_scale(pivot, self.type(1) / r[pivot, i])
            for j in r.nz_col(i):
                if j == pivot:
                    continue
                scale = r[j,i]
                r.row_add(j, pivot, - scale)
            pivotlist.append(i)
            pivot += 1
        return r, pivotlist

    def rref(self):
        "rref"
        A = self[:, :]
        for j in xrange(A.rows):
            rlist = A.nz_col_lower(j)
            if A[j, j] == self.type(0):
                if rlist:
                    A.row_swap(j, rlist[0])
                    rlist.pop(0)
                else:
                    continue
            A.row_scale(j, self.type(1) /  A[j, j])
            for i in A.nz_col(j):
                if i != j:
                    A.row_add(i, j, - A[i, j])
        return A
                    
    def nz_col(self, j):
        li = []
        for i in xrange(self.rows):
            if self[i, j] != self.type(0):
                li.append(i)
        return li

    def nz_col_lower(self, j):
        " Returns the row indices of non-zero elements in column j, below the diagonal"
        li = []
        for i in xrange(j + 1, self.rows):
            if self[i, j] != self.type(0):
                li.append(i)
        return li

    def is_upper(self):
        return all(j >= i for i in xrange(self.rows) for j, _ in self.mat[i])
            
    def applyfunc(self, f):
        for i in xrange(self.rows):
            for ind, (j, value) in enumerate(self.mat[i]):
                self.mat[i][ind] = (self.mat[i][ind][0], f(self.mat[i][ind][1]))

    def sparsity(self):
        return float(self.nnz()) / (self.rows * self.cols)

    def nnz(self):
        return sum(len(self.mat[i]) for i in xrange(self.rows))

    def join_rows(A, B):
        A = A[:, :]
        B = B[:, :]
        assert A.rows == B.rows
        for i in xrange(B.rows):
            for ind, (j, val) in enumerate(B.mat[i]):
                B.mat[i][ind] = (j + A.cols, val)
            A.mat[i].extend(B.mat[i])
        A.cols += B.cols
        return A

    @classmethod
    def eye(cls, n, one = 1, zero = 0):
        return cls(n, n, lambda i, j: one if i==j else zero)

    def inv_rref(self):
        aug = self.join_rows(LILMatrix.eye(self.rows, one = self.type(1), zero = self.type(0)))
        reduced = aug.rref()
        return reduced[:,self.rows:]

    def clone(self):
        return LILMatrix(self.rows, self.cols, lambda i, j: self[i, j])
        
    def __add__(self, other):
        if not isinstance(other, LILMatrix):
            if self.is_square():
                other = other * LILMatrix.eye(self.rows)
            else:
                raise Exception
        add = self.clone()
        for i in xrange(self.rows):
            add.mat[i] = _row_add(self.mat[i], other.mat[i], 1)
        return add

    def is_square(self):
        return self.rows == self.cols

    def transpose(self):
        T = LILMatrix.zeros(self.rows)
        for i in xrange(self.rows):
            for j, value in self.mat[i]:
                T.mat[j].append((i, value))
        return T

    @property
    def shape(self):
        return (self.rows, self.cols)
    
    def __eq__(self, other):
        if not self.shape == other.shape:
            return False
        return all(self.mat[i][ind] == other.mat[i][ind] for i in xrange(self.rows) for ind in xrange(len(self.mat[i]))) 

    def __ne__(self, other):
        if not self.shape == other.shape:
            return True
        return any(self.mat[i][ind] != other.mat[i][ind] for i in xrange(self.rows) for ind in xrange(len(self.mat[i]))) 

    def doolittle(self):
        A = self.clone()
        n = A.rows
        for i in xrange(n):
            for j in xrange(i):
                a = A[i, j]
                for p in xrange(j):
                    a -= A[i, p] * A[p, j]
                A[i, j] = a / A[j, j]
            for j in xrange(i, n):
                a = A[i, j]
                for p in xrange(i):
                    a -= A[i, p] * A[p, j]
                A[i, j] = a
        return A

    def LU_sparse(self):
        row_swaps = []
        A = self.clone()
        n = self.rows
        for k in xrange(n):

            rlist = A.nz_col_lower(k)

            # Pivoting
            if self[k, k] == 0:
                if not rlist:
                    print k
                    raise Exception('Singular')
                A.row_swap(k, rlist[0])
                row_swaps.append((k, rlist[0]))
                rlist.pop(0)
            assert A[k, k]

            # Algorithm
            for i in rlist:
                A[i, k] /= A[k, k]

            for j, val in A.mat[k]:
                if j <= k:
                    continue
                for i in rlist:
                    A[i, j] -= A[i, k] * val
    
        L = LILMatrix.eye(self.rows)
        for i in xrange(L.rows):
            for j in xrange(i):
                L[i, j] = A[i, j]
        U = LILMatrix.zeros(self.rows)
        for i in xrange(U.rows):
            for j in xrange(i, U.rows):
                U[i, j] = A[i, j]
        return L, U, row_swaps

    def permute(self, row_swaps):
        for r1, r2 in row_swaps:
            self.row_swap(r1, r2)
        
            
            
                
def _row_add(row1, row2, alpha):
    "li = row1 + alpha * row2 "
    li = []
    i1 = i2 = 0
    n1 = len(row1)
    n2 = len(row2)
    while i1 < n1 or i2 < n2:
        # print i1, i2, len(row1), len(row2)
        if i1 < n1 and (i2 >= n2 or row1[i1][0] < row2[i2][0]):
            li.append(row1[i1])
            i1 += 1
        elif i1 >= n1 or row1[i1][0] > row2[i2][0]:
            li.append((row2[i2][0], alpha * row2[i2][1]))
            i2 += 1
        else:
            val = row1[i1][1] + alpha * row2[i2][1]
            if val != 0:
                li.append((row1[i1][0], val))
            i1 += 1
            i2 += 1
    return li   

def scalar_product(self, scalar):
    prod = self.clone()
    for i in xrange(self.rows):
        for ind, (j, value) in enumerate(self.mat[i]):
            prod.mat[i][ind] = (j, scalar * value)
    return prod

def randInvLILMatrix(n, d, min=-5, max=10):
    A = LILMatrix(n, n, lambda i, j: random.randint(min, max) if abs(i - j) <= d-1 else 0)
    return A

def test(n, d):
    A = randInvLILMatrix(n, d)
    A.applyfunc(S)
    B = A.gauss_col()
    if not B.is_upper():
        return A, A.toMatrix().det()

def randLILMatrix(i, j, min=1, max=1, sparsity=0.5):
    return LILMatrix(i, j, lambda i, j: random.randint(min, max) if random.random() < sparsity else 0)
    
