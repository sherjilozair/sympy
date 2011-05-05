from sympy import *
from random import randint
from time import *

def setupM(num, n):
    mats=[]
    vectors=[]
    for i in range(num):
        mats.append(Matrix(n,n,lambda i,j : r()))
        vectors.append(Matrix(n,1,lambda i,j:r()))
    return mats, vectors

def testCholesky(mats, vectors):
    t1 = time()
    for A in mats:
        for B in vectors:
            assert A * A.ch(B) == B
    t2=time()
    print "Chloesky takes",
    print t2-t1,
    print "seconds"
    return (t2-t1)/len(mats)*len(vectors)

def testLDL(mats, vectors):
    t1 = time()
    for A in mats:
        for B in vectors:
            assert A * A.LDLSolve(B) == B
    t2=time()
    print "LDL takes",
    print t2-t1,
    print "seconds"
    return (t2-t1)/len(mats)*len(vectors)


def testLU(mats, vectors):
    t1 = time()
    for A in mats:
        for B in vectors:
            assert A * A.LUsolve(B) == B
    t2=time()
    print "LU takes",
    print t2-t1,
    print "seconds"
    return (t2-t1)/len(mats)*len(vectors)
    

def r():
    return randint(-100,100)

def start(num, n):
    M, V = setupM(num, n)
    print testCholesky(M, V)
    print testLDL(M,V)
    print testLU(M,V)
