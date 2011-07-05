from sympy.matrices.matrices import Matrix
from sympy.matrices.dokmatrix import DOKMatrix
from sympy.matrices.lilmatrix import LILMatrix

from sympy.matrices.matrixutils import _from_callable, _from_list, _dict_to_densematrix, _dict_to_dokmatrix, _dict_to_lilmatrix

# This is a high-level class.
# TODO: Make the internals low-level.

class Matrix_(object):

    # # # # # # # # # # # # # #
    # User-level constructor  #
    # # # # # # # # # # # # # #

    def __init__(self, *args, **kwargs):
        # Preprocessing data.
        # Making a dict of all non-zero key-value pairs.
        if len(args) == 3:
            if callable(args[2]):
                rows, cols, dok = _from_callable(*args)
            elif isinstance(args[2], (list, tuple)):
                rows, cols, dok = _from_list(*args)
            else: # Value passing only through lambda and list right now
                raise NotImplementedError("Data type not understood")
            # dok is now a dict of of all non-zero elements

            # Choosing a representation for the matrix.     
            if 'repr' in kwargs:
                # User-specified representation.
                if kwargs['repr'] == 'dok':
                    mat = _dict_to_dokmatrix(rows, cols, dok)
                elif kwargs['repr'] == 'lil':
                    mat = _dict_to_lilmatrix(rows, cols, dok)
                elif kwargs['repr'] == 'dense':
                    mat = _dict_to_densematrix(rows, cols, dok)
                else:
                    raise Exception('repr %s not recognized' % kwargs['repr'])
            else:
                # Representation chosen by sparsity check.
                sparsity = float(len(dok)) / (rows * cols)
                if sparsity > 0.5:
                    mat = _dict_to_densematrix(rows, cols, dok)
                else:
                    mat = _dict_to_dokmatrix(rows, cols, dok)

                # mat is now one of the low-level matrices, 
                # generated using the data given.

        elif len(args) == 1:
            # Matrix_(InternalMatrix) would work and efficiently.
            if isinstance(args[0], (DOKMatrix, LILMatrix, Matrix)):
                mat = args[0]
            else:
                raise TypeError('Data type not understood')
        # Put the internal matrix as self.mat
        self.mat = mat
        # We're done.

    # # # # # # # # # # #       
    # Element accessors #
    # # # # # # # # # # #
      
    
    def __setitem__(self, key, value):
        """

        """
        if type(key) is tuple:
            i, j = key
            if type(i) is slice or type(j) is slice:
                if isinstance(value, Matrix):
                    self.copyin_matrix(key, value)
                    return
                if isinstance(value, (list, tuple)):
                    self.copyin_list(key, value)
                    return
            else:
                # a2idx inlined
                if not type(i) is int:
                    try:
                        i = i.__index__()
                    except AttributeError:
                        raise IndexError("Invalid index a[%r]" % (key,))

                # a2idx inlined
                if not type(j) is int:
                    try:
                        j = j.__index__()
                    except AttributeError:
                        raise IndexError("Invalid index a[%r]" % (key,))


                if not (i>=0 and i<self.rows and j>=0 and j < self.cols):
                    raise IndexError("Index out of range: a[%s]" % (key,))
                else:
                    self.mat[i, j] = value
                    return

        else:
            # row-wise decomposition of matrix
            if type(key) is slice:
                raise IndexError("Vector slices not implemented yet.")
            else:
                if not type(key) is int:
                    try:
                        key = key.__index__()
                    except AttributeError:
                        raise IndexError("Invalid index a[%r]" % (key,))
                if key is not None:
                    self.mat[key] = value
                    return

        raise IndexError("Invalid index: a[%s]"%repr(key))

    def __getitem__(self,key):
        """
        """
        if type(key) is tuple:
            i, j = key
            if type(i) is slice or type(j) is slice:
                return self.mat.submatrix(key)

            else:
                if not isinstance(i, int):
                    try:
                        i = i.__index__()
                    except AttributeError:
                        raise IndexError("Invalid index a[%r]" % key)

                if not isinstance(i, int):
                    try:
                       j = j.__index__()
                    except AttributeError:
                        raise IndexError("Invalid index a[%r]" % key)

                if not (i>=0 and i<self.rows and j>=0 and j < self.cols):
                    raise IndexError("Index out of range: a[%s]" % key)
                else:
                    return self.mat[i, j]


        else:
            if type(key) is slice:
                return self.rowdecomp[key] # implement rowdecomp function in the internals
            else:
                if not type(key) is int:
                    try:
                        key = key.__index__()
                    except AttributeError:
                        raise IndexError("Invalid index a[%r]" % (key,))
                if key is not None:
                    return self.mat[key]

        raise IndexError("Invalid index: a[%s]" % repr(key))

    # # # # # # # # # # # # #
    # Relational operators  #
    # # # # # # # # # # # # #

    def __eq__(self, other):
        # Where to put __eq__ methods. 
        # Matrix_ or the internals or both ?
        # Code duplication is fine ?
        pass
        

    def __ne__(self, other):
        pass

    def __neg__(self):
        return self * (-1)
        

    def __add__(self, other):
        # if either is dense, the other will be converted to dense.
        # else if one is DOK, then other will be converted to DOK.
        if isinstance(other, Matrix_):
            # If other is a Matrix_ (not a scalar)
            if not self.mat.shape == other.mat.shape:
                raise ValueError('Matrices have different shape')
                # Make a ShapeError Exception
            if isinstance(self.mat, Matrix) or isinstance(other.mat, Matrix):
                # If either is dense
                return Matrix_(self.mat.to_dense() + other.mat.to_dense())
            # Convert internals to dok and add
            return Matrix_(self.mat.to_dokmatrix() + other.mat.to_dokmatrix())
        # `matrix + scalar` not supported
        return NotImplemented

    __radd__ = __add__
        
    def __sub__(self, other):
        '''
        Subtraction of Matrices.
        if either is dense, the other will be converted to dense.
        else if one is DOK, then other will be converted to DOK.
        '''
        if isinstance(other, Matrix_):
            # If other is a Matrix_ (not a scalar)
            if isinstance(self.mat, Matrix) or isinstance(other.mat, Matrix):
                # If either is dense
                return Matrix_(self.mat.to_dense() - other.mat.to_dense())
            # Convert internals to dok and subtract
            return Matrix_(self.mat.to_dokmatrix() - other.mat.to_dokmatrix())
        # `matrix - scalar` not supported
        return NotImplemented

    def __rsub__(self, other):
        # Do we want to support `scalar + matrix` ?
        # If not, then this is of no use.
        pass

    def __mul__(self, other):
        if isinstance(other, Matrix_):
            # If other is a Matrix_ (not a scalar)
            if not self.cols == other.rows:
                raise ValueError('Given matrices cannot be multiplied')
            if isinstance(self.mat, Matrix) or isinstance(other.mat, Matrix):
                # If either is dense
                return Matrix_(self.mat.to_dense() * other.mat.to_dense())
            # Convert internals to dok and multiply
            return Matrix_(self.mat.to_dokmatrix() * other.mat.to_dokmatrix())
        # Scalar multiplication
        return Matrix_(self.mat.scalar_multiply(other))

    def __rmul__(self, other):
        # Since left multiplicand couldn't handle the multiplication, 
        # other must be a scalar
        return Matrix_(self.mat.scalar_multiply(other))

    # # # # # # # # # # #
    # Rows and Columns  #
    # # # # # # # # # # #

    @property
    def rows(self):
        return self.mat.rows

    @property
    def cols(self):
        return self.mat.cols

    # # # # # # # # # # # # # #
    # Linear Algebra methods  #
    # # # # # # # # # # # # # #

    def det(self, method=None):
        # DOK/cholesky does not support det
        # LIL/GE is the only available method for det
        # default for each internal would be different.
        if not self.is_square():
            raise ValueError('Matrix should be square.')
        if isinstance(self.mat, Matrix) or method == 'bareis' or method == 'berkowitz':
            return self.mat.to_dense().det(method=method)
        return self.mat.to_lilmatrix().det(method=method)

    def solve(self, rhs, method=None):
        # All three matrices provide solves
        if not isinstance(rhs, Matrix_):
            raise TypeError('%s is not a vector' % repr(rhs))

        if self.cols != rhs.cols:
            raise TypeError('Wrong shape of rhs') #ShapeError

        if self.rows < self.cols:
                raise ValueError('Under-determined system')
        
        if method == 'CH':
            if isinsntance(self.mat, Matrix):
                return Matrix_(self.mat.cholesky_solve(rhs.mat.to_dense()))
            return Matrix_(self.mat.to_dokmatrix().cholesky_solve(rhs.mat.to_dokmatrix()))

        if method == 'QR':
            return Matrix_(self.mat.to_dense().QRsolve(rhs.mat.to_dense()))

        if method == 'LU':
            if isinsntance(self.mat, Matrix):
                return Matrix_(self.mat.LUsolve(rhs.mat.to_dense()))
            return Matrix_(self.mat.to_lilmatrix().LUsolve(rhs.mat.to_lilmatrix()))

        if not method:
            return Matrix_(self.mat.solve(rhs.mat))

        raise Exception('method not recognized')
            

    def inv(self, method=None):
        # All three matrices provide inverse
        # FIXME: Redo this !
        # TODO: Write Matrix.inverse_cholesky and DOKMatrix.inverse_cholesky

        if not self.is_square():
            raise ValueError('ShapeError') # ShapeError

        if isinstance(self.mat, Matrix):
            return Matrix_(self.mat.inv(method=method))

        if method == 'CH':
            if isinsntance(self.mat, Matrix):
                return Matrix_(self.mat.inv_cholesky())
            return Matrix_(self.mat.to_dokmatrix().inv_cholesky())
        
        if method == 'GE':
            if isinsntance(self.mat, Matrix):
                return Matrix_(self.mat.LUsolve(rhs.mat.to_dense()))
            return Matrix_(self.mat.to_lilmatrix().LUsolve(rhs.mat.to_lilmatrix()))

        if not method:
            return Matrix_(self.mat.solve(rhs.mat))

        raise Exception('method not recognized')
        

    def rref(self):
        # DOK would be converted to LIL for this.
        if isinstance(self.mat, Matrix):
            return Matrix_(self.mat.rref())
        return Matrix_(self.mat.to_lilmatrix().rref())
    
    def nullspace(self):
        # TODO write nullspace method for LILMatrix

        if not self.is_square():
            raise ValueError('ShapeError') # ShapeError

        if isinstance(self.mat, Matrix):
            return Matrix_(self.mat.nullspace())
        return Matrix_(self.mat.to_lilmatrix().nullspace())

    def factorize(self, method=None):
        # We have many.
        if method == 'LU':
            try:
                if isinstance(self.mat, Matrix):
                    L, U, p = self.mat.LUdecomposition()
                else:
                    L, U, p = self.mat.to_lilmatrix().LUdecomposition()
            except SingularMatrixException:
                raise Exception('Matrix cannot be LU factorized')
            else:
                return Matrix_(L), Matrix_(U), p

        if method == 'CH':
            if not self.is_symmetric():
                raise Exception('Matrix cannot be factorized by cholesky\'s method')
            try:
                if isinstance(self.mat, Matrix):
                    return Matrix_(self.mat.cholesky())
                return Matrix_(self.mat.to_dokmatrix().cholseky())
            except SingularMatrixException:
                raise Exception('Matrix cannot be factorized by cholesky\'s method')

        if method == 'LDL':
            if not self.is_symmetric():
                raise Exception('Matrix cannot be LDL factorized')
            try:
                if isinstance(self.mat, Matrix):
                    L, D = self.mat.LDLdecomposition()
                else:
                    L, D = self.mat.to_dokmatrix().LDLdecomposition()
            except SingularMatrixException:
                raise Exception('Matrix cannot be LDL factorized.')
            else:
                return Matrix_(L), Matrix_(D)

        if method == 'QR':
            if self.rows < self.cols:
                raise Exception('Matrix cannot be QR factorized')
            try:
                Q, R = self.mat.to_dense().QRdecomposition()
            except SingularMatrixException:
                raise Exception('Matrix cannot be QR factorized')
            else:
                return Matrix_(Q), Matrix_(R)
            
            
    # # # # # # # # # # # # 
    # Utility functions   #
    # # # # # # # # # # # #

    def __str__(self):
        return self.mat.__str__()

    def __repr__(self):
        return self.mat.__str__()

    def clone(self):
        return Matrix_(self.mat.clone())

    # # # # # # # # # # # # # 
    # .is_something methods #
    # # # # # # # # # # # # #

    def is_symmetric(self):
        return self.mat.is_symmetric()

    def is_lower(self):
        return self.mat.is_lower()

    def is_upper(self):
        return self.mat.is_upper()

    def is_square(self):
        return self.mat.is_square()

