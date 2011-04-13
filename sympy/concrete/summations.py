from sympy.core import (Expr, S, C, Symbol, Equality, Interval, sympify, Wild,
                        Tuple, Dummy)
from sympy.solvers import solve
from sympy.utilities import flatten

class Sum(Expr):
    """Represents unevaluated summation."""

    def __new__(cls, f, *symbols, **assumptions):
        f = sympify(f)

        if not symbols:
            raise ValueError("No symbols given.")

        else:
            limits = []

            for V in symbols:
                if isinstance(V, Symbol):
                    limits.append(Tuple(V))
                    continue
                elif isinstance(V, Equality):
                    if isinstance(V.lhs, Symbol):
                        if isinstance(V.rhs, Interval):
                            limits.append(Tuple(V.lhs, V.rhs.start, V.rhs.end))
                        else:
                            limits.append(Tuple(V.lhs, V.rhs))

                        continue
                elif isinstance(V, (tuple, list, Tuple)):
                    V = flatten(V)
                    if len(V) == 1:
                        if isinstance(V[0], Symbol):
                            limits.append(Tuple(V[0]))
                            continue
                    elif len(V) in (2, 3):
                        if isinstance(V[0], Symbol):
                            limits.append(Tuple(*map(sympify, V)))
                            continue

                raise ValueError("Invalid summation variable or limits")

        obj = Expr.__new__(cls, **assumptions)
        arglist = [f]
        arglist.extend(limits)
        obj._args = tuple(arglist)

        return obj

    @property
    def function(self):
        return self._args[0]

    @property
    def limits(self):
        return self._args[1:]

    @property
    def variables(self):
        """Return a list of the summation variables

        >>> from sympy import Sum
        >>> from sympy.abc import x, i
        >>> Sum(x**i, (i, 1, 3)).variables
        [i]
        """
        return [l[0] for l in self.limits]

    @property
    def free_symbols(self):
        """
        This method returns the symbols that will exist when the
        summation is evaluated. This is useful if one is trying to
        determine whether a sum is dependent on a certain
        symbol or not.

        >>> from sympy import Sum
        >>> from sympy.abc import x, y
        >>> Sum(x, (x, y, 1)).free_symbols
        set([y])
        """
        # analyze the summation
        # >>> Sum(x*y,(x,1,2),(y,1,3)).args
        # (x*y, Tuple(x, 1, 2), Tuple(y, 1, 3))
        # >>> Sum(x, x, y).args
        # (x, Tuple(x), Tuple(y))
        intgrl = self
        args = intgrl.args
        integrand, limits = args[0], args[1:]
        if integrand.is_zero:
            return set()
        isyms = integrand.free_symbols
        for ilim in limits:
            if len(ilim) == 1:
                isyms.add(ilim[0])
                continue
            # take out the target symbol
            if ilim[0] in isyms:
                isyms.remove(ilim[0])
            if len(ilim) == 3 and ilim[1] == ilim[2]:
                # if two limits are the same the sum is 0
                # and there are no symbols
                return set()
            # add in the new symbols
            for i in ilim[1:]:
                isyms.update(i.free_symbols)
        return isyms

    def doit(self, **hints):
        #if not hints.get('sums', True):
        #    return self
        f = self.function
        for i, a, b in self.limits:
            f = eval_sum(f, (i, a, b))
            if f is None:
                return self

        if hints.get('deep', True):
            return f.doit(**hints)
        else:
            return f

    def _eval_summation(self, f, x):
        return

    def euler_maclaurin(self, m=0, n=0, eps=0, eval_integral=True):
        """
        Return an Euler-Maclaurin approximation of self, where m is the
        number of leading terms to sum directly and n is the number of
        terms in the tail.

        With m = n = 0, this is simply the corresponding integral
        plus a first-order endpoint correction.

        Returns (s, e) where s is the Euler-Maclaurin approximation
        and e is the estimated error (taken to be the magnitude of
        the first omitted term in the tail):

            >>> from sympy.abc import k, a, b
            >>> from sympy import Sum
            >>> Sum(1/k, (k, 2, 5)).doit().evalf()
            1.28333333333333
            >>> s, e = Sum(1/k, (k, 2, 5)).euler_maclaurin()
            >>> s
            7/20 - log(2) + log(5)
            >>> from sympy import sstr
            >>> print sstr((s.evalf(), e.evalf()), full_prec=True)
            (1.26629073187416, 0.0175000000000000)

        The endpoints may be symbolic:

            >>> s, e = Sum(1/k, (k, a, b)).euler_maclaurin()
            >>> s
            -log(a) + log(b) + 1/(2*a) + 1/(2*b)
            >>> e
            Abs(-1/(12*b**2) + 1/(12*a**2))

        If the function is a polynomial of degree at most 2n+1, the
        Euler-Maclaurin formula becomes exact (and e = 0 is returned):

            >>> Sum(k, (k, 2, b)).euler_maclaurin()
            (-1 + b/2 + b**2/2, 0)
            >>> Sum(k, (k, 2, b)).doit()
            -1 + b/2 + b**2/2

        With a nonzero eps specified, the summation is ended
        as soon as the remainder term is less than the epsilon.
        """
        m = int(m)
        n = int(n)
        f = self.function
        if not len(self.limits) == 1:
            raise NotImplementedError("Euler-Maclaurin series for more than one variables not implemented.")
        i, a, b = self.limits[0]
        s = S.Zero
        if m:
            for k in range(m):
                term = f.subs(i, a+k)
                if (eps and term and abs(term.evalf(3)) < eps):
                    return s, abs(term)
                s += term
            a += m
        x = Dummy('x')
        I = C.Integral(f.subs(i, x), (x, a, b))
        if eval_integral:
            I = I.doit()
        s += I
        def fpoint(expr):
            if b is S.Infinity:
                return expr.subs(i, a), 0
            return expr.subs(i, a), expr.subs(i, b)
        fa, fb = fpoint(f)
        iterm = (fa + fb)/2
        g = f.diff(i)
        for k in xrange(1, n+2):
            ga, gb = fpoint(g)
            term = C.bernoulli(2*k)/C.Factorial(2*k)*(gb-ga)
            if (eps and term and abs(term.evalf(3)) < eps) or (k > n):
                break
            s += term
            g = g.diff(i, 2)
        return s + iterm, abs(term)

    def _eval_subs(self, old, new):
        if self == old:
            return new
        newlimits = []
        for lim in self.limits:
            if lim[0] == old:
                return self
            newlimits.append( (lim[0],lim[1].subs(old,new),lim[2].subs(old,new)) )

        return Sum(self.args[0].subs(old, new), *newlimits)


def summation(f, *symbols, **kwargs):
    """
    Compute the summation of f with respect to symbols.

    The notation for symbols is similar to the notation used in Integral.
    summation(f, (i, a, b)) computes the sum of f with respect to i from a to b,
    i.e.,

                                b
                              ____
                              \   `
    summation(f, (i, a, b)) =  )    f
                              /___,
                              i = a


    If it cannot compute the sum, it returns an unevaluated Sum object.
    Repeated sums can be computed by introducing additional symbols tuples::

    >>> from sympy import summation, oo, symbols, log
    >>> i, n, m = symbols('i n m', integer=True)

    >>> summation(2*i - 1, (i, 1, n))
    n**2
    >>> summation(1/2**i, (i, 0, oo))
    2
    >>> summation(1/log(n)**n, (n, 2, oo))
    Sum(log(n)**(-n), (n, 2, oo))
    >>> summation(i, (i, 0, n), (n, 0, m))
    m/3 + m**2/2 + m**3/6

    """
    return Sum(f, *symbols, **kwargs).doit(deep=False)

def telescopic_direct(L, R, n, (i, a, b)):
    """Returns the direct summation of the terms of a telescopic sum

    L is the term with lower index
    R is the term with higher index
    n difference between the indexes of L and R

    For example:

    >>> from sympy.concrete.summations import telescopic_direct
    >>> from sympy.abc import k, a, b
    >>> telescopic_direct(1/k, -1/(k+2), 2, (k, a, b))
    1/a + 1/(1 + a) - 1/(1 + b) - 1/(2 + b)

    """
    s = 0
    for m in xrange(n):
        s += L.subs(i,a+m) + R.subs(i,b-m)
    return s

def telescopic(L, R, (i, a, b)):
    '''Tries to perform the summation using the telescopic property

    return None if not possible
    '''
    if L.is_Add or R.is_Add:
        return None
    s = None
    #First we try to solve using match
    #Maybe this should go inside solve
    k = Wild("k")
    sol = (-R).match(L.subs(i, i + k))
    if sol and k in sol:
        if L.subs(i,i + sol[k]) == -R:
            #sometimes match fail(f(x+2).match(-f(x+k))->{k: -2 - 2x}))
            s = sol[k]
    #Then we try to solve using solve
    if not s or not s.is_Integer:
        m = Symbol("m")
        try:
            s = solve(L.subs(i, i + m) + R, m)[0]
        except IndexError:#(ValueError, IndexError):
            pass
    if s and s.is_Integer:
        if s < 0:
            return telescopic_direct(R, L, abs(s), (i, a, b))
        elif s > 0:
            return telescopic_direct(L, R, s, (i, a, b))
    return None

def eval_sum(f, (i, a, b)):
    if f.is_Number:
        if f is S.NaN:
            return S.NaN
        elif f is S.Zero:
            return S.Zero

    if not f.has(i):
        return f*(b-a+1)
    definite = a.is_Integer and b.is_Integer
    # Doing it directly may be faster if there are very few terms.
    if definite and (b-a < 100):
        return eval_sum_direct(f, (i, a, b))
    # Try to do it symbolically. Even when the number of terms is known,
    # this can save time when b-a is big.
    # We should try to transform to partial fractions
    value = eval_sum_symbolic(f.expand(), (i, a, b))
    if value is not None:
        return value
    # Do it directly
    if definite:
        return eval_sum_direct(f, (i, a, b))

def eval_sum_symbolic(f, (i, a, b)):
    if not f.has(i):
        return f*(b-a+1)
    # Linearity
    if f.is_Mul:
        L, R = f.as_two_terms()
        if not L.has(i):
            sR = eval_sum_symbolic(R, (i, a, b))
            if sR: return L*sR
        if not R.has(i):
            sL = eval_sum_symbolic(L, (i, a, b))
            if sL: return R*sL
    if f.is_Add:
        L, R = f.as_two_terms()
        lrsum = telescopic(L, R, (i, a, b))
        if lrsum: return lrsum
        lsum = eval_sum_symbolic(L, (i, a, b))
        rsum = eval_sum_symbolic(R, (i, a, b))
        if None not in (lsum, rsum):
            return lsum + rsum
    # Polynomial terms with Faulhaber's formula
    p = C.Wild('p')
    e = f.match(i**p)
    if e is not None:
        c = p.subs(e)
        B = C.bernoulli
        if c.is_integer and c >= 0:
            s = (B(c+1, b+1) - B(c+1, a))/(c+1)
            return s.expand()
    # Geometric terms
    c1 = C.Wild('c1', exclude=[i])
    c2 = C.Wild('c2', exclude=[i])
    c3 = C.Wild('c3', exclude=[i])
    e = f.match(c1**(c2*i+c3))
    if e is not None:
        c1 = c1.subs(e)
        c2 = c2.subs(e)
        c3 = c3.subs(e)
        # TODO: more general limit handling
        return c1**c3 * (c1**(a*c2) - c1**(c2+b*c2)) / (1 - c1**c2)
    return None

def eval_sum_direct(expr, (i, a, b)):
    s = S.Zero
    if expr.has(i):
        for j in xrange(a, b+1):
            s += expr.subs(i, j)
    else:
        for j in xrange(a, b+1):
            s += expr
    return s
