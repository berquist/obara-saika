from math import sqrt, pi, exp
from mpmath import gamma, gammainc
import sys

#-------------------------------------------------------------------------------

class X2:

    """The parent class for two-center integrals.

    Keyword arguments:
    scale -- The multiplicative scaling factor for the non-zero term contribution.
    prefactors -- A list of the multiplicative prefactors that appear in recursion expressions before each integral.
    q -- A list of 6 integers representing the angular momentum of each center [la, ma, na, lb, mb, nb].
    kind -- A one-letter string representing the type of integral (one of S, T, V, M, L, E, J).
    operator -- For integrals with variable operators (M, E), the 3 powers of the operator components [nx, ny, nz].
    d -- For integrals with multiple components (L, J), one of {0, 1, 2} specifying the Cartesian component.
    order -- For auxiliary integrals, the superscript present in the formulas.
    """

    def __init__(self,
                 scale=1,
                 prefactors=[],
                 q=[0,0,0,0,0,0],
                 kind='S',
                 operator=[0,0,0],
                 d=0,
                 order=0):
        self.scale = scale
        self.prefactors = prefactors
        self.q = q
        self.kind = kind
        self.operator = operator
        self.d = d
        self.order = order

    def __str__(self):
        if self.kind in ('S'):
            return "{}({})".format(self.kind, self.q)
        elif self.kind in ('T', 'V'):
            return "{}(q={}, order={})".format(self.kind,
                                               self.q,
                                               self.order)

    def __repr__(self):
        return "X2({}, {}, {}, {}, {}, {}, {})".format(self.scale,
                                                       self.prefactors,
                                                       self.q,
                                                       self.kind,
                                                       self.operator,
                                                       self.d,
                                                       self.order)

#-------------------------------------------------------------------------------

class X4:

    """The parent class for four-center integrals."""

    def __init__(self,
                 scale=1,
                 prefactors=[],
                 q=[0,0,0,0,0,0,0,0,0,0,0,0],
                 order=0):
        self.scale = scale
        self.prefactors = prefactors
        self.q = q
        self.order = order

#-------------------------------------------------------------------------------

def list_is_flat(l):
    """Determine if the given list is flat (only check for other lists,
    not other iterables).
    """
    for x in l:
        if isinstance(x, list):
            return False
    return True

#-------------------------------------------------------------------------------

def flatten_sub(l):
    """Flatten a list by one level."""
    return [item for sublist in l for item in sublist]

#-------------------------------------------------------------------------------

def flatten(l):
    """Flatten a list as much as possible."""
    l_out = l[:]
    for _ in range(50):
        if list_is_flat(l_out):
            return l_out
        else:
            l_out = flatten_sub(l_out)

    sys.stderr.write('ERROR: max depth reached in flatten\n')
    sys.exit(1)

#-------------------------------------------------------------------------------

def find_fun_to_lower(q, n):
    """Determine which basis function to lower in an n-center integral.

    Arguments:
    q -- The list of exponents, length == 3 * n.
    n -- The number of integral centers.

    >>> find_fun_to_lower([1, 0, 0, 0, 0, 0], 2)
    0
    >>> find_fun_to_lower([0, 1, 0, 0, 0, 0], 2)
    0
    >>> find_fun_to_lower([0, 0, 1, 0, 0, 0], 2)
    0
    >>> find_fun_to_lower([0, 0, 0, 1, 0, 0], 2)
    1
    >>> find_fun_to_lower([0, 0, 0, 0, 1, 0], 2)
    1
    >>> find_fun_to_lower([0, 0, 0, 0, 0, 1], 2)
    1
    >>> find_fun_to_lower([0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], 4)
    1
    >>> find_fun_to_lower([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], 4)
    3
    >>> find_fun_to_lower([1, 0, 0, 0, 0, 0, 0, 0, 1], 3)
    0
    >>> find_fun_to_lower([0, 0, 0, 0, 0, 0, 0, 0, 1], 3)
    2
    >>> find_fun_to_lower([0, 0, 0, 0, 1, 0, 0, 0, 1], 3)
    1
    >>> find_fun_to_lower([0, 0, 0, 0, 2, 0, 0, 0, 1], 3)
    2
    """

    # Determine the total angular momentum on each center.
    l = []
    for i in range(n):
        l.append(q[i*3] + q[i*3 + 1] + q[i*3 + 2])

    # find function to lower
    # start with lowest angular momentum above s
    fun = -1
    kmax = max(l) + 1
    for i in range(n):
        k = l[i]
        # If we're larger than an s-function...
        if k > 0:
            if k < kmax:
                kmax = k
                fun = i

    return fun

#-------------------------------------------------------------------------------

def find_component_to_lower(fun):
    """Determine which exponent/component of a basis function to lower.

    The first component with a non-zero value will be the result.

    Arguments:
    fun -- The coefficients for a single basis function (list of 3 integers).

    >>> find_component_to_lower([0, 0, 1])
    2
    >>> find_component_to_lower([0, 1, 1])
    1
    >>> find_component_to_lower([1, 0, 1])
    0
    """

    for i, c in enumerate(fun):
        if c > 0:
            return i

    return -1

#-------------------------------------------------------------------------------

def apply_os4(x4):

    if sum(x4.q) == 0:
        return [x4]

    fun = find_fun_to_lower(x4.q, 4)
    if fun == -1:
        return [x4]
    component = find_component_to_lower([x4.q[fun*3], x4.q[fun*3 + 1], x4.q[fun*3 + 2]])

    if component == 0:
        i1 = [  0,  1,  2,  3][fun]
        i2 = [  4,  4,  5,  5][fun]
    if component == 1:
        i1 = [  6,  7,  8,  9][fun]
        i2 = [ 10, 10, 11, 11][fun]
    if component == 2:
        i1 = [ 12, 13, 14, 15][fun]
        i2 = [ 16, 16, 17, 17][fun]

    bra = [0, 1]
    ket = [2, 3]

    pre = []
    pre.append(i1)
    pre.append(i2)
    if fun in bra:
        pre.append(18)
        pre.append(20)
        pre.append(18)
        pre.append(20)
    else:
        pre.append(19)
        pre.append(21)
        pre.append(19)
        pre.append(21)
    pre.append(22)
    pre.append(22)

    a = fun
    if a in bra:
       bra.remove(a)
       b = bra.pop()
       c = ket.pop()
       d = ket.pop()
    else:
       ket.remove(a)
       b = ket.pop()
       c = bra.pop()
       d = bra.pop()

    x4_copy = []
    scale = x4.scale
    order = x4.order
    for term in range(8):
        x4_new = X4(scale=scale,
                    prefactors=x4.prefactors[:],
                    q=x4.q[:],
                    order=order)
        x4_copy.append(x4_new)
        x4_copy[term].q[fun*3 + component] -= 1

    for term in [1, 3, 5, 6, 7]:
        x4_copy[term].order += 1

    x4_copy[2].q[a*3 + component] -= 1
    x4_copy[3].q[a*3 + component] -= 1
    x4_copy[4].q[b*3 + component] -= 1
    x4_copy[5].q[b*3 + component] -= 1
    x4_copy[6].q[c*3 + component] -= 1
    x4_copy[7].q[d*3 + component] -= 1

    n = []
    n.append(1)
    n.append(1)
    n.append(x4.q[a*3 + component] - 1)
    n.append(x4.q[a*3 + component] - 1)
    n.append(x4.q[b*3 + component])
    n.append(x4.q[b*3 + component])
    n.append(x4.q[c*3 + component])
    n.append(x4.q[d*3 + component])

    x4_list = []
    for term in range(8):
        if n[term] > 0:
            if all(i >= 0 for i in x4_copy[term].q):
                if n[term] > 1:
                    x4_copy[term].scale *= n[term]
                x4_copy[term].prefactors.append(pre[term])
                x4_list.append(x4_copy[term])

    x4_final = []
    for x4 in x4_list:
        if all(i == 0 for i in x4.q):
            x4_final.append([x4])
        else:
            x4_final.append(apply_os4(x4))

    return flatten(x4_final)

#-------------------------------------------------------------------------------

def apply_os2(x, kind):
    """Apply the Obara-Saika recursion scheme to a two-center integral.

    Arguments:
    x -- A two-center integral.
    kind -- The type of integral (string S, T, V, M, L, E, J).
    """

    orders = x.q + x.operator

    # base case
    if sum(orders) == 0:
        x.kind = kind
        return [x]

    # Determine which basis function and component to lower.
    # The component is one of (x, y, z).
    fun = find_fun_to_lower(orders, 3)
    # Make sure to not choose the operator vrr until q is exhausted.
    if fun == 2 and sum(x.q) > 0:
        fun = find_fun_to_lower(x.q, 2)
    if fun == -1:
        x.kind = kind
        return [x]
    component = find_component_to_lower([orders[fun*3], orders[fun*3 + 1], orders[fun*3 + 2]])

    # Determine the index of q to descend on.
    # where q = [xa, xb, ya, yb, za, zb].
    # Note: I thought the order of q == [xa, ya, za, xb, yb, zb]?
    if component == 0:
        i1 = [0, 1, 6][fun]
    if component == 1:
        i1 = [2, 3, 7][fun]
    if component == 2:
        i1 = [4, 5, 8][fun]

    if kind == 'S':
        # The vrr for overlap integrals consists of three 'terms'.
        pre = []
        pre.append(i1)
        pre.append(6)
        pre.append(7)

    if kind == 'T':
        pre = []
        pre.append(i1)
        pre.append(6)
        pre.append(7)
        pre.append(8)
        i2 = [9, 10][fun]
        pre.append(i2)

    if kind == 'V':
        if component == 0:
            i2 = 6
        if component == 1:
            i2 = 7
        if component == 2:
            i2 = 8
        pre = []
        pre.append(i1)
        pre.append(i2)
        pre.append(9)
        pre.append(10)
        pre.append(9)
        pre.append(10)

    if kind == 'M':
        pre = []
        pre.append(i1)
        pre.append(9)
        pre.append(10)
        pre.append(11)

    if kind == 'L':
        pre = []
        pre.append(i1)
        pre.append(6)
        pre.append(7)
        pre.append(8 + x.d)
        pre.append(11)
        pre.append(12)
        pre.append(13)

    # Determine which of the basis functions is ('a', 'b').
    l = [0, 1]
    if fun == 2:
        a, b = 2, 2
    else:
        a = fun
        l.remove(a)
        b = l[0]

    # These are the number of integrals that appear in the main
    # recursion equations for each kind.
    if kind == 'S':
        num_terms = 3
    elif kind == 'T':
        num_terms = 5
    elif kind == 'V':
        num_terms = 6
    elif kind == 'M':
        num_terms = 4
        if fun == 2:
            num_terms = 2
    elif kind == 'L':
        num_terms = 7
    # Not implemented yet...
    elif kind == 'E':
        num_terms = -1
    elif kind == 'J':
        num_terms = -1
    else:
        sys.stderr.write('ERROR: unexpected kind\n')
        sys.exit(1)

    # Make copies of the current integral to manipulate later,
    # one for each term in the recursion expression.
    x_copy = []
    scale = x.scale
    order = x.order # currently only used for V
    for term in range(num_terms):
        x_new = X2(scale=scale,
                   prefactors=x.prefactors[:],
                   q=x.q[:],
                   kind=kind,
                   operator=x.operator[:],
                   d=x.d,
                   order=order)
        x_copy.append(x_new)

    # These are the terms in [A19] with (m+1).
    if kind == 'V':
        for term in [1, 3, 5]:
            x_copy[term].order += 1

    # Look at the last line of [A12].
    if kind == 'T':
        x_copy[3].kind = 'S'
        x_copy[4].kind = 'S'

    # Look at the last two lines of [A31].
    if kind == 'L':
        x_copy[3].kind = 'S'
        x_copy[4].kind = 'S'
        x_copy[5].kind = 'S'
        x_copy[6].kind = 'S'

    # 1. Lower the target component for all three terms.
    # 2. Lower again on center 'a'.
    # 3. Lower again on center 'b'.
    if kind == 'S':
        x_copy[0].q[fun*3 + component] -= 1
        x_copy[1].q[fun*3 + component] -= 1
        x_copy[2].q[fun*3 + component] -= 1
        x_copy[1].q[a*3 + component] -= 1
        x_copy[2].q[b*3 + component] -= 1

    if kind == 'T':
        x_copy[0].q[fun*3 + component] -= 1
        x_copy[1].q[fun*3 + component] -= 1
        x_copy[2].q[fun*3 + component] -= 1
        # term 4 (x_copy[3]) has the same components but becomes an
        # overlap integral
        x_copy[4].q[fun*3 + component] -= 2
        x_copy[1].q[a*3 + component] -= 1
        x_copy[2].q[b*3 + component] -= 1

    if kind == 'V':
        x_copy[0].q[fun*3 + component] -= 1
        x_copy[2].q[fun*3 + component] -= 1
        x_copy[4].q[fun*3 + component] -= 1
        x_copy[2].q[a*3 + component] -= 1
        x_copy[4].q[b*3 + component] -= 1
        x_copy[1].q = x_copy[0].q[:]
        x_copy[3].q = x_copy[2].q[:]
        x_copy[5].q = x_copy[4].q[:]

    if kind == 'M':
        # For the case of the moment operator over s functions.
        if fun == 2:
            x_copy[0].operator[component] -= 1
            x_copy[1].operator[component] -= 2
        else:
            x_copy[0].q[fun*3 + component] -= 1
            x_copy[1].q[fun*3 + component] -= 1
            x_copy[2].q[fun*3 + component] -= 1
            x_copy[3].q[fun*3 + component] -= 1
            x_copy[1].q[a*3 + component] -= 1
            x_copy[2].q[b*3 + component] -= 1
            x_copy[3].operator[component] -= 1

    if kind == 'L':
        x_copy[0].q[fun*3 + component] -= 1
        x_copy[1].q[fun*3 + component] -= 1
        x_copy[2].q[fun*3 + component] -= 1
        x_copy[3].q[fun*3 + component] -= 1
        x_copy[4].q[fun*3 + component] -= 1
        x_copy[5].q[fun*3 + component] -= 1
        x_copy[6].q[fun*3 + component] -= 1
        x_copy[1].q[a*3 + component] -= 1
        x_copy[2].q[b*3 + component] -= 1
        x_copy[4].q[b*3 + 0] -= 1
        x_copy[5].q[b*3 + 1] -= 1
        x_copy[6].q[b*3 + 2] -= 1

    # Now that the descending part of the vrr has been performed, keep
    # track of which new terms are going to be zero/non-zero.
    n = []
    # The first term is always going to be non-zero, otherwise we
    # wouldn't even be in the vrr routine.
    n.append(1)

    if kind == 'S':
        n.append(x.q[a*3 + component] - 1)
        n.append(x.q[b*3 + component])

    if kind == 'T':
        n.append(x.q[a*3 + component] - 1)
        n.append(x.q[b*3 + component])
        n.append(1)
        n.append(x.q[a*3 + component] - 1)

    if kind == 'V':
        n.append(1)
        n.append(x.q[a*3 + component] - 1)
        n.append(x.q[a*3 + component] - 1)
        n.append(x.q[b*3 + component])
        n.append(x.q[b*3 + component])

    if kind == 'M':
        if fun == 2:
            n.append(x.operator[component] - 1)
        else:
            n.append(x.q[a*3 + component] - 1)
            n.append(x.q[b*3 + component])
            n.append(x.operator[component])

    if kind == 'L':
        n.append(x.q[a*3 + component] - 1)
        n.append(x.q[b*3 + component])
        n.append([1 if d == component else 0 for d in range(3)][x.d])
        n.append([1 if d == component and d == 0 else 0 for d in range(3)][x.d])
        n.append([1 if d == component and d == 1 else 0 for d in range(3)][x.d])
        n.append([1 if d == component and d == 2 else 0 for d in range(3)][x.d])

    # Generate a list of all non-zero terms for an expression.
    x_list = []
    for term in range(num_terms):
        if n[term] > 0:
            orders = x_copy[term].q + x_copy[term].operator
            if all(i >= 0 for i in orders):
                if n[term] > 1:
                    x_copy[term].scale *= n[term]
                x_copy[term].prefactors.append(pre[term])
                x_list.append(x_copy[term])

    # If we've hit the base case where all the components are zero,
    # terminate, otherwise recurse through each term.
    x_final = []
    for y in x_list:
        orders = y.q + y.operator
        if all(i == 0 for i in orders):
            x_final.append([y])
        else:
            x_final.append(apply_os2(y, kind=y.kind))

    return flatten(x_final)

#-------------------------------------------------------------------------------

def get_r12_squared(r1, r2):
    """Get the distance between two centers in Cartesian space."""
    return (r1[0] - r2[0])**2.0 + (r1[1] - r2[1])**2.0 + (r1[2] - r2[2])**2.0

#-------------------------------------------------------------------------------

def get_k(z1, z2, r1, r2):

    r12 = get_r12_squared(r1, r2)
    f0 = z1 + z2
    if r12 > 0.0:
        f1 = -z1*z2*r12/f0
        f2 = exp(f1)
    else:
        f2 = 1.0
    return sqrt(2.0)*f2*pi**(5.0/4.0)/f0

#-------------------------------------------------------------------------------

def get_rho(za, zb, zc, zd):

    z = za + zb
    n = zc + zd
    return z*n/(z + n)

#-------------------------------------------------------------------------------

def get_bi_center(z1, z2, r1, r2):

    z = z1 + z2
    rx = (z1*r1[0] + z2*r2[0])/z
    ry = (z1*r1[1] + z2*r2[1])/z
    rz = (z1*r1[2] + z2*r2[2])/z
    return [rx, ry, rz]

#-------------------------------------------------------------------------------

def boys(n, x):

    if x > 0.0:
        f = 2.0*x**(n + 0.5)
        g = gamma(n + 0.5)
        gi = 1.0 - gammainc(n + 0.5, x, regularized=True)
        return g*gi/f
    else:
        return 1.0/(n*2 + 1)

#-------------------------------------------------------------------------------

def get_aux(za, zb, zc, zd, ra, rb, rc, rd):

    k1 = get_k(za, zb, ra, rb)
    k2 = get_k(zc, zd, rc, rd)
    return k1*k2/sqrt(za + zb + zc + zd)

#-------------------------------------------------------------------------------

def get_coulomb(za, zb, zc, zd, ra, rb, rc, rd, c):

    rp = get_bi_center(za, zb, ra, rb)
    rq = get_bi_center(zc, zd, rc, rd)
    rw = get_bi_center((za + zb), (zc + zd), rp, rq)
    t  = get_rho(za, zb, zc, zd)*get_r12_squared(rp, rq)
    s  = get_aux(za, zb, zc, zd, ra, rb, rc, rd)

    z = za + zb
    n = zc + zd
    rho = z*n/(z + n)

    prefac = []
    prefac.append(rp[0] - ra[0])
    prefac.append(rp[0] - rb[0])
    prefac.append(rq[0] - rc[0])
    prefac.append(rq[0] - rd[0])
    prefac.append(rw[0] - rp[0])
    prefac.append(rw[0] - rq[0])
    prefac.append(rp[1] - ra[1])
    prefac.append(rp[1] - rb[1])
    prefac.append(rq[1] - rc[1])
    prefac.append(rq[1] - rd[1])
    prefac.append(rw[1] - rp[1])
    prefac.append(rw[1] - rq[1])
    prefac.append(rp[2] - ra[2])
    prefac.append(rp[2] - rb[2])
    prefac.append(rq[2] - rc[2])
    prefac.append(rq[2] - rd[2])
    prefac.append(rw[2] - rp[2])
    prefac.append(rw[2] - rq[2])
    prefac.append(0.5/z)
    prefac.append(0.5/n)
    prefac.append(-0.5*rho/(z*z))
    prefac.append(-0.5*rho/(n*n))
    prefac.append(0.5/(z + n))

    fun = X4(q=c)
    expansion = apply_os4(fun)
    integral = 0.0
    for i in range(sum(c) + 1):
        b = boys(i, t)*s
        for f in expansion:
            if f.order == i:
                g = 1.0
                for k in f.prefactors:
                    g *= prefac[k]
                integral += float(f.scale)*b*g
    return integral

#-------------------------------------------------------------------------------

def get_overlap(za, zb, ra, rb, c):
    """Compute the overlap integral
     ...
    where
     ...

    Recursion expression:
     (a + 1_{i} || b) = (P_{i} - A_{i})*(a || b) +
                        (0.5/z)*(a - 1_{i} || b) +
                        (0.5/z)*(a || b - 1_{i})

    Arguments:
    za, zb -- Exponent of each basis function center
    ra, rb -- Position of each basis function center in Cartesian space
    c -- A list of 6 integers for the angular momentum of each component of each bf [xa, ya, za, xb, yb, zb]
    """

    # Compute the base case:
    # equations [14], [13], product Gaussian center, distance between basis function centers, and (s||s)/[22] [Table VI, p3972]
    z = za + zb
    e = za*zb/(za + zb)
    rp = get_bi_center(za, zb, ra, rb)
    ab = get_r12_squared(ra, rb)
    aux = exp(-e*ab)*(pi/z)**1.5

    # Determine all possible prefactors for each term in the recursion
    # expression.
    # equation A2, p3971
    # The components of (P - B) appear because (a||b) == (b||a).
    prefac = []
    prefac.append(rp[0] - ra[0])
    prefac.append(rp[0] - rb[0])
    prefac.append(rp[1] - ra[1])
    prefac.append(rp[1] - rb[1])
    prefac.append(rp[2] - ra[2])
    prefac.append(rp[2] - rb[2])
    prefac.append(0.5/z)
    prefac.append(0.5/z)

    fun = X2(q=c)
    expansion = apply_os2(fun, kind='S')
    integral = 0.0
    for f in expansion:
        g = 1.0
        for k in f.prefactors:
            g *= prefac[k]
        integral += float(f.scale)*aux*g

    return integral

#-------------------------------------------------------------------------------

def get_kinetic(za, zb, ra, rb, c):

    # Compute the base case:
    # equations [14], [13], product Gaussian center, distance between basis function centers, and (s||s) [22]
    z = za + zb
    e = za*zb/(za + zb)
    rp = get_bi_center(za, zb, ra, rb)
    ab = get_r12_squared(ra, rb)
    aux = get_overlap(za, zb, ra, rb, [0, 0, 0, 0, 0, 0])

    prefac = []
    prefac.append(rp[0] - ra[0])
    prefac.append(rp[0] - rb[0])
    prefac.append(rp[1] - ra[1])
    prefac.append(rp[1] - rb[1])
    prefac.append(rp[2] - ra[2])
    prefac.append(rp[2] - rb[2])
    prefac.append(0.5/z)
    prefac.append(0.5/z)
    prefac.append(2.0*e)
    prefac.append(-e/za)
    prefac.append(-e/zb)

    fun = X2(q=c)
    expansion = apply_os2(fun, kind='T')
    integral = 0.0
    for f in expansion:
        # equation A13
        if f.kind == 'T':
            g = e*(3.0 - 2.0*e*ab)
            for k in f.prefactors:
                g *= prefac[k]
            integral += float(f.scale)*aux*g
        if f.kind == 'S':
            g = 1.0
            for k in f.prefactors:
                g *= prefac[k]
            integral += float(f.scale)*aux*g

    return integral

#-------------------------------------------------------------------------------

def get_nuclear(za, zb, ra, rb, rc, c):

    z = za + zb
    rp = get_bi_center(za, zb, ra, rb)
    pc = get_r12_squared(rp, rc)
    u = z*pc
    aux = -2.0*(z/pi)**0.5*get_overlap(za, zb, ra, rb, [0, 0, 0, 0, 0, 0])

    prefac = []
    prefac.append(rp[0] - ra[0])
    prefac.append(rp[0] - rb[0])
    prefac.append(rp[1] - ra[1])
    prefac.append(rp[1] - rb[1])
    prefac.append(rp[2] - ra[2])
    prefac.append(rp[2] - rb[2])
    prefac.append(-rp[0] + rc[0])
    prefac.append(-rp[1] + rc[1])
    prefac.append(-rp[2] + rc[2])
    prefac.append(0.5/z)
    prefac.append(-0.5/z)

    fun = X2(q=c)
    expansion = apply_os2(fun, kind='V')
    integral = 0.0
    for f in expansion:
        g = 1.0
        for k in f.prefactors:
            g *= prefac[k]
        integral += float(f.scale)*aux*g*boys(f.order, u)

    return integral

#-------------------------------------------------------------------------------

def get_moment(za, zb, ra, rb, rc, c, order):
    """Compute the Cartesian moment integral
     ...
    where
     ...

    Recursion expression:
     (a + 1_{i} | M(mu) | b) = (P_{i} - A_{i})*(a | M(mu) | b) +
                               (0.5/z)*(a - 1_{i} | M(mu) | b) +
                               (0.5/z)*(a | M(mu) | b - 1_{i}) +
                               (0.5/z)*(a | M(mu - 1_{i}) | b)

    Arguments:
    za, zb -- Exponent of each basis function center
    ra, rb -- Position of each basis function center in Cartesian space
    rc -- Postion of the moment operator in Cartesian space
    c -- A list of 6 integers for the angular momentum of each component of each bf [xa, ya, za, xb, yb, zb]
    order -- The powers of the three components of the moment operator [nx, ny, nz]
    """

    z = za + zb
    e = za*zb/(za + zb)
    rp = get_bi_center(za, zb, ra, rb)
    ab = get_r12_squared(ra, rb)
    aux = get_overlap(za, zb, ra, rb, [0, 0, 0, 0, 0, 0])

    prefac = []
    prefac.append(rp[0] - ra[0])
    prefac.append(rp[0] - rb[0])
    prefac.append(rp[1] - ra[1])
    prefac.append(rp[1] - rb[1])
    prefac.append(rp[2] - ra[2])
    prefac.append(rp[2] - rb[2])
    prefac.append(rp[0] - rc[0])
    prefac.append(rp[1] - rc[1])
    prefac.append(rp[2] - rc[2])
    prefac.append(0.5/z)
    prefac.append(0.5/z)
    prefac.append(0.5/z)

    fun = X2(q=c, operator=order)
    expansion = apply_os2(fun, kind='M')
    integral = 0.0
    for f in expansion:
        g = 1.0
        for k in f.prefactors:
            g *= prefac[k]
        integral += float(f.scale)*aux*g

    return integral

#-------------------------------------------------------------------------------

def get_angmom(za, zb, ra, rb, rc, c, d):
    """Compute the orbital angular momentum integral
     ...
    where
     ...

    Recursion expression:

    Arguments:
    za, zb -- Exponent of each basis function center
    ra, rb -- Position of each basis function center in Cartesian space
    rc -- Postion of the angular momentum operator in Cartesian space
    c -- A list of 6 integers for the angular momentum of each component of each bf [xa, ya, za, xb, yb, zb]
    d -- ... {0, 1, 2}
    """

    z = za + zb
    e = za*zb/(za + zb)
    rp = get_bi_center(za, zb, ra, rb)
    # Equation [A31]
    aux = 2*e*(ra[d]-rc[d])*(rb[d]-rc[d])*get_overlap(za, zb, ra, rb, [0, 0, 0, 0, 0, 0])

    prefac = []
    prefac.append(rp[0] - ra[0])
    prefac.append(rp[0] - rb[0])
    prefac.append(rp[1] - ra[1])
    prefac.append(rp[1] - rb[1])
    prefac.append(rp[2] - ra[2])
    prefac.append(rp[2] - rb[2])
    prefac.append(0.5/z)
    prefac.append(0.5/z)
    prefac.append((zb/z)*(rb[0] - rc[0]))
    prefac.append((zb/z)*(rb[1] - rc[1]))
    prefac.append((zb/z)*(rb[2] - rc[2]))
    prefac.append(0.5/z)
    prefac.append(0.5/z)
    prefac.append(0.5/z)

    fun = X2(q=c, d=d)
    expansion = apply_os2(fun, kind='L')
    integral = 0.0
    for f in expansion:
        g = 1.0
        for k in f.prefactors:
            g *= prefac[k]
        integral += float(f.scale)*aux*g

    return integral

#-------------------------------------------------------------------------------

def get_efield():
    """Compute the electric field integral
     ...
    where
     ...

    Recursion expression:

    Arguments:
    """

    z = za + zb
    e = za*zb/(za + zb)
    rp = get_bi_center(za, zb, ra, rb)
    ab = get_r12_squared(ra, rb)

    prefac = []

    fun = X2(q=c)
    expansion = apply_os2(fun, kind='E')
    integral = 0.0
    for f in expansion:
        g = 1.0
        for k in f.prefactors:
            g *= prefac[k]
        integral += float(f.scale)*aux*g

    return integral

#-------------------------------------------------------------------------------

def get_spinorb():
    """Compute the spin-orbit interaction integral
     ...
    where
     ...

    Recursion expression:

    Arguments:
    """

    z = za + zb
    e = za*zb/(za + zb)
    rp = get_bi_center(za, zb, ra, rb)
    ab = get_r12_squared(ra, rb)

    prefac = []

    fun = X2(q=c, d=d)
    expansion = apply_os2(fun, kind='J')
    integral = 0.0
    for f in expansion:
        g = 1.0
        for k in f.prefactors:
            g *= prefac[k]
        integral += float(f.scale)*aux*g

    return integral

#-------------------------------------------------------------------------------

def get_fermi(za, zb, ra, rb, rc, c):
    """Compute the Fermi contact integral
     ...
    """

    z = za + zb
    e = za*zb/(za + zb)
    rp = get_bi_center(za, zb, ra, rb)
    ab = get_r12_squared(ra, rb)
    aux = exp(-e*ab)*(pi/z)**1.5

    # Determine all possible prefactors for each term in the recursion
    # expression.
    # equation A2, p3971
    # The components of (P - B) appear because (a||b) == (b||a).
    prefac = []
    prefac.append(rp[0] - ra[0])
    prefac.append(rp[0] - rb[0])
    prefac.append(rp[1] - ra[1])
    prefac.append(rp[1] - rb[1])
    prefac.append(rp[2] - ra[2])
    prefac.append(rp[2] - rb[2])
    prefac.append(0.5/z)
    prefac.append(0.5/z)

    fun = X2(q=c)
    expansion = apply_os2(fun, kind='S')
    integral = 0.0
    for f in expansion:
        g = 1.0
        for k in f.prefactors:
            g *= prefac[k]
        integral += float(f.scale)*aux*g

    return integral
