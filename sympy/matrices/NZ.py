from sympy.matrices import *
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.core import Add, S
from sympy.functions.elementary.integers import floor
import random

def nfsch(self):
    Lstruc = cholesky_structure(self)
    Lstruc2 = []
    for i in Lstruc:
        Lstruc2.extend(i)
    Lstruc2.sort()
    L = DOKMatrix(self.rows, self.cols, lambda i,j:0)
    for i, j in Lstruc2:
        if i == j:
            L[j, j] = (self[j, j] - sum(L[j, k] ** 2 for k in xrange(j))) ** (S(1)/2)
        else:
            L[i, j] = (S(1) / L[j, j]) * (self[i, j] - sum(L[i, k] * L[j, k] for k in xrange(j)))
    return L

def _lower_nonzero_structure(self):
    n = self.rows
    NZlist = [[]] * n
    keys = [(i, j) for (j, i) in sorted((j,i) for (i, j) in self.mat.keys())]
    k = 0
    startset = False
    l = len(keys)
    for i in xrange(l):
        if startset and keys[i][1] == k + 1:
            NZlist[k] = keys[start:i]
            startset = False
        k = keys[i][1]
        if not startset and keys[i][1] == k and keys[i][1] <= keys[i][0]:
            start = i
            startset = True
    if (n-1, n-1) == keys[-1]:
        NZlist[-1].append((n-1,n-1))
    return NZlist

def NZ2_test():
    self = randDOKMatrix(15,15)
    n = self.rows
    NZlist = [[]] * n
    keys = [(i, j) for (j, i) in sorted((j,i) for (i, j) in self.mat.keys())] # Make this DS independent
    k = 0
    startset = False
    for i in range(len(keys)):
        if startset and keys[i][1] == k + 1:
            NZlist[k] = keys[start:i]
            startset = False
        k = keys[i][1]
        if not startset and keys[i][1] == k and keys[i][1] <= keys[i][0]:
            start = i
            startset = True
    if (n-1, n-1) == keys[-1]:
        NZlist[-1].append((n-1,n-1))
    return NZlist
    

def non_zero_pattern(self):
    n = self.rows
    A = [[]] * n
    keys = [(i, j) for (j, i) in sorted((j,i) for (i, j) in self.mat.keys())] # Make this DS independent
    oldindex = oldcol = k = 0
    startset = True
    for i in range(1, len(keys)):
        if keys[i][1] > oldcol and startset:
            A[k] = keys[oldindex:i]
            k += 1
            oldcol = keys[i][1]
            startset = False
        if keys[i][1] <= keys[i][0] and (keys[i-1][1] > keys[i-1][0] or keys[i][1]!=keys[i-1][1]):
            oldindex = i
            startset = True
    return A

NZ = non_zero_pattern


def sch(self):
    n = self.rows; A = NZ(self); L = [set([])] * n; pi = [0] * n;
    for i in range(n):
        L[i] = set(A[i])
        for j in range(n):
            if pi[j] == i:
                L[i] = L[i] | L[j] - set([j])
            pi[i] = setmin(L[i] - set([i]))
            print i, pi[i]
    return [sorted(list(i)) for i in L]

def setmin(a):
    try:
        return min(a)
    except:
        return 0

def cholesky_structure(self):
    n = self.rows
    parent = [0] * n
    E = NZ2(self)
    for k in range(n):
        parent[k] = setmin(i[0] for i in E[k])
        for i, j in E[k]:
            add(E, (i, parent[i]))
    return E

def add(E, ij):
    if ij not in E[ij[1]]:
        E[ij[1]].append(ij)
        E[ij[1]].sort()

def add2(E, ij):
    li = E[ij[1]]
    for k in xrange(len(li)):
        if li[k] > ij:
            if ij == li[k-1]:
                break
            li.insert(k,ij) 
            break

def testF(n):
    A = randDOKMatrix(n, n)
    B = fsch(A)
    C = DOKMatrix(n, n, lambda i,j: 1 if (i,j) in B[j] else 0)
    return A, C

def testNZ(n):
    A = randDOKMatrix(n, n)
    B = NZ(A)
    C = DOKMatrix(n, n, lambda i,j: 1 if (i,j) in B[j] else 0)
    for i in range(A.rows):
        for j in range(A.rows):
            if i >= j:
                if not A[i,j] == C[i, j]:
                    print "error", i, j
            else:
                if not C[i, j] == 0:
                    print "error", i, j
    return A, C

def testNZ2(n):
    A = randDOKMatrix(n, n)
    B = NZ2(A)
    C = DOKMatrix(n, n, lambda i,j: 1 if (i,j) in B[j] else 0)
    for i in range(A.rows):
        for j in range(A.rows):
            if i >= j:
                assert A[i,j] == C[i, j]
            else:
                assert C[i, j] == 0
    return A, C

def testmany(n, times):
    for i in xrange(times):
        testNZ2(n)

def randNZ2(n):
    NZ2(randDOKMatrix(n, n))

#randNZ2(15)

def randSymDOKMatrix(n):
    self = randMatrix(n,n); self = self.T*self; self = self.applyfunc(lambda i: floor(i/(n*1000)-2));B = DOKMatrix(self.rows, self.cols, self.mat); return B

def randSymMDOKMatrix(n):
    self = randMatrix(n,n); self1 = self.T*self; self1 = self1.applyfunc(lambda i: floor(i/(n*1000)-2));B = DOKMatrix(self1.rows, self1.cols, self1.mat); return self1, B

def randDiagMatrix(m, n, d, _min=0, _max=1):
    return DOKMatrix(m, n, lambda i, j: randint(_min, _max) if abs(i - j) <= d-1 else 0)

# def randInvDOKMatrix(n, d, _min=0, _max=1):
#     A = DOKMatrix(n, n, lambda i, j: 1 if i==j else random.randint(_min, _max) if abs(i - j) <= d-1 else 0)
#    if A.toMatrix().det() != 0:
#        return A
#    else:
#        return randInvDOKMatrix(n, d, _min, _max)

        

def code():
    # from line_profiler import LineProfiler; from NZ import *; profile = LineProfiler(nfsch)
    # from NZ import *; self = randMatrix(5,5);self = self.T * self; self = DOKMatrix(self.rows, self.cols, self.mat)
    # results = profile.runcall(nfsch, self)
    # profile.print_stats()
    pass
