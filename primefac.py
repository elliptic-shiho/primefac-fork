#! /usr/bin/env python

# Formatting note: this file's maximum line length is 128 characters.

# TODO: Python 2.x/3.x compatibility

from multiprocessing import Process, Queue as mpQueue
from itertools import count, takewhile
from random import randrange
from math import log

try: from gmpy2 import mpz; mpzv, inttypes =      2, (int, long, type(mpz(1)))
except ImportError:
 try: from gmpy import mpz; mpzv, inttypes =      1, (int, long, type(mpz(1)))
 except ImportError:   mpz, mpzv, inttypes = int, 0, (int, long)



def gcd(a, b):
    while b: a, b = b, a%b
    return abs(a)
if mpzv == 1: from gmpy import gcd
if mpzv == 2: from gmpy2 import gcd

def isqrt(n):
    if n == 0: return 0
    x, y = n, (n + 1) / 2
    while y < x: x, y = y, (y + n/y) / 2
    return x
if mpzv == 1: from gmpy import sqrt as isqrt
if mpzv == 2: from gmpy2 import isqrt

def introot(n, r=2):
    if n < 0: return None if r%2 == 0 else -introot(-n, r)
    if n < 2: return n
    if r == 2: return isqrt(n)
    lower, upper = 0, n
    while lower != upper - 1:
        mid = (lower + upper) / 2
        m = mid**r
        if   m == n: return  mid
        elif m <  n: lower = mid
        elif m >  n: upper = mid
    return lower
if mpzv == 1:
    from gmpy import root
    def introot(n, r=2):
        if n < 0: return None if r%2 == 0 else -introot(-n, r)
        return root(n, r)[0]
if mpzv == 2:
    from gmpy2 import iroot
    def introot(n, r=2):
        if n < 0: return None if r%2 == 0 else -introot(-n, r)
        return iroot(n, r)[0]

# Recursive sieve of Eratosthenes
def primegen():
    yield 2; yield 3; yield 5; yield 7; yield 11; yield 13;
    ps = primegen() # yay recursion
    p = ps.next() and ps.next()
    q, sieve, n = p**2, {}, 13
    while True:
        if n not in sieve:
            if n < q: yield n
            else:
                next, step = q + 2*p, 2*p
                while next in sieve: next += step
                sieve[next] = step
                p = ps.next()
                q = p**2
        else:
            step = sieve.pop(n)
            next = n + step
            while next in sieve: next += step
            sieve[next] = step
        n += 2

def primes(n): return list(takewhile(lambda p: p < n, primegen()))      # The primes STRICTLY LESS than n

def listprod(a): return reduce(lambda x,y:x*y, a, 1)

def nextprime(n):
    if n < 2: return 2
    if n == 2: return 3
    n = (n + 1) | 1    # first odd larger than n
    m = n % 6
    if m == 3:
        if isprime(n+2): return n+2
        n += 4
    elif m == 5:
        if isprime(n): return n
        n += 2
    for m in count(n, 6):
        if isprime(m  ): return m
        if isprime(m+4): return m+4

def pfactor(n):
    s, d, q = 0, n-1, 2
    while not d & q - 1: s, q = s+1, q*2
    return s, d / (q / 2)

def sprp(n, a, s=None, d=None):
    if n%2 == 0: return False
    if (s is None) or (d is None): s, d = pfactor(n)
    x = pow(a, d, n)
    if x == 1: return True
    for i in xrange(s):
        if x == n - 1: return True
        x = pow(x, 2, n)
    return False
if mpzv == 2:
    from gmpy2 import is_strong_prp
    def sprp(n, a, s=None, d=None): return is_strong_prp(n, a)

def jacobi(a, p):
    if (p%2 == 0) or (p < 0): return None # p must be a positive odd number
    if (a == 0) or (a == 1): return a
    a, t = a%p, 1
    while a != 0:
        while not a & 1:
            a /= 2
            if p & 7 in (3, 5): t *= -1
        a, p = p, a
        if (a & 3 == 3) and (p & 3) == 3: t *= -1
        a %= p
    return t if p == 1 else 0
if mpzv == 1: from gmpy import jacobi
if mpzv == 2: from gmpy2 import jacobi

def chain(n, u1, v1, u2, v2, d, q, m): # Used in SLPRP.  TODO: figure out what this does.
    k = q
    while m > 0:
        u2, v2, q = (u2*v2)%n, (v2*v2-2*q)%n, (q*q)%n
        if m%2 == 1:
            u1, v1 = u2*v1+u1*v2, v2*v1+u2*u1*d
            if u1%2 == 1: u1 = u1 + n
            if v1%2 == 1: v1 = v1 + n
            u1, v1, k = (u1/2)%n, (v1/2)%n, (q*k)%n
        m /= 2
    return u1, v1, k

def isprime(n, tb=(3,5,7,11), eb=(2,), mrb=()):      # TODO: more streamlining
    # tb: trial division basis
    # eb: Euler's test basis
    # mrb: Miller-Rabin basis
    
    # This test suite's first false positve is unknown but has been shown to be greater than 2**64.
    # Infinitely many are thought to exist.
    
    if n%2 == 0 or n < 13 or n == isqrt(n)**2: return n in (2, 3, 5, 7, 11) # Remove evens, squares, and numbers less than 13
    if any(n%p == 0 for p in tb): return n in tb                            # Trial division
    
    for b in eb:                                                            # Euler's test
        if b >= n: continue
        if not pow(b, n-1, n) == 1: return False
        r = n - 1
        while r%2 == 0: r /= 2
        c = pow(b, r, n)
        if c == 1: continue
        while c != 1 and c != n-1: c = pow(c, 2, n)
        if c == 1: return False
    
    s, d = pfactor(n)
    if not sprp(n, 2, s, d): return False
    if n < 2047: return True
    if n >= 3825123056546413051: # BPSW has two phases: SPRP with base 2 and SLPRP.  We just did the SPRP; now we do the SLPRP:
        d = 5
        while True:
            if gcd(d, n) > 1:
                p, q = 0, 0
                break
            if jacobi(d, n) == -1:
                p, q = 1, (1 - d) / 4
                break
            d = -d - 2*d/abs(d)
        if p == 0: return n == d
        s, t = pfactor(n + 2)
        u, v, u2, v2, m = 1, p, 1, p, t/2
        k = q
        while m > 0:
            u2, v2, q = (u2*v2)%n, (v2*v2-2*q)%n, (q*q)%n
            if m%2 == 1:
                u, v = u2*v+u*v2, v2*v+u2*u*d
                if u%2 == 1: u += n
                if v%2 == 1: v += n
                u, v, k = (u/2)%n, (v/2)%n, (q*k)%n
            m /= 2
        if (u == 0) or (v == 0): return True
        for i in xrange(1, s):
            v, k = (v*v-2*k)%n, (k*k)%n
            if v == 0: return True
        return False
    
    if not mrb:
        if   n <             1373653: mrb = [3]
        elif n <            25326001: mrb = [3,5]
        elif n <          3215031751: mrb = [3,5,7]
        elif n <       2152302898747: mrb = [3,5,7,11]
        elif n <       3474749660383: mrb = [3,5,6,11,13]
        elif n <     341550071728321: mrb = [3,5,7,11,13,17]   # This number is also a false positive for primes(19+1).
        elif n < 3825123056546413051: mrb = [3,5,7,11,13,17,19,23]   # Also a false positive for primes(31+1).
    return all(sprp(n, b, s, d) for b in mrb)                               # Miller-Rabin
if mpzv == 2: from gmpy2 import is_bpsw_prp as isprime

def ilog(x, b): # greatest integer l such that b**l <= x.
    l = 0
    while x >= b:
        x /= b
        l += 1
    return l

# Returns the largest integer that, when squared/cubed/etc, yields n, or 0 if no such integer exists.
# Note that the power to which this number is raised will be prime.
def ispower(n):
    for p in primegen():
        r = introot(n, p)
        if r is None: continue
        if r ** p == n: return r
        if r == 1: return 0

def pollardRho_brent(n):
    if isprime(n): return n
    g = n
    while g == n:
        y, c, m, g, r, q = randrange(1, n), randrange(1, n), randrange(1, n), 1, 1, 1
        while g==1:
            x, k = y, 0
            for i in xrange(r): y = (y**2 + c) % n
            while k < r and g == 1:
                ys = y
                for i in xrange(min(m, r-k)):
                    y = (y**2 + c) % n
                    q = q * abs(x-y) % n
                g, k = gcd(q, n), k+m
            r *= 2
        if g==n:
            while True:
                ys = (ys**2+c)%n
                g = gcd(abs(x-ys), n)
                if g > 1: break
    return g

# http://programmingpraxis.com/2010/04/27/modern-elliptic-curve-factorization-part-2/
def pollard_pm1(n, B1=100, B2=1000):       # TODO: What are the best default bounds and way to increment them?
    if isprime(n): return n
    m = ispower(n)
    if m: return m
    while True:
        pg = primegen()
        q = 2           # TODO: what about other initial values of q?
        p = pg.next()
        while p <= B1: q, p = pow(q, p**ilog(B1, p), n), pg.next()
        g = gcd(q-1, n)
        if 1 < g < n: return g
        while p <= B2: q, p = pow(q, p, n), pg.next()
        g = gcd(q-1, n)
        if 1 < g < n: return g
        # These bounds failed.  Increase and try again.
        B1 *= 10
        B2 *= 10

def mlucas(v, a, n):
    """ Helper function for williams_pp1().  Multiplies along a Lucas sequence modulo n. """
    v1, v2 = v, (v**2 - 2) % n
    for bit in bin(a)[3:]: v1, v2 = ((v1**2 - 2) % n, (v1*v2 - v) % n) if bit == "0" else ((v1*v2 - v) % n, (v2**2 - 2) % n)
    return v1

def williams_pp1(n):
    if isprime(n): return n
    m = ispower(n)
    if m: return m
    for v in count(1):
        for p in primegen():
            e = ilog(isqrt(n), p)
            if e == 0: break
            for _ in xrange(e): v = mlucas(v, p, n)
            g = gcd(v - 2, n)
            if 1 < g < n: return g
            if g == n: break

# http://programmingpraxis.com/2010/04/23/modern-elliptic-curve-factorization-part-1/
# http://programmingpraxis.com/2010/04/27/modern-elliptic-curve-factorization-part-2/
def ecadd(p1, p2, p0, n): # Add two points p1 and p2 given point P0 = P1-P2 modulo n
    x1,z1 = p1; x2,z2 = p2; x0,z0 = p0
    t1, t2 = (x1-z1)*(x2+z2), (x1+z1)*(x2-z2)
    return (z0*pow(t1+t2,2,n) % n, x0*pow(t1-t2,2,n) % n)
def ecdub(p, A, n): # double point p on A modulo n
    x, z = p; An, Ad = A
    t1, t2 = pow(x+z,2,n), pow(x-z,2,n)
    t = t1 - t2
    return (t1*t2*4*Ad % n, (4*Ad*t2 + t*An)*t % n)
def ecmul(m, p, A, n): # multiply point p by m on curve A modulo n
    if m == 0: return (0, 0)
    elif m == 1: return p
    else:
        q = ecdub(p, A, n)
        if m == 2: return q
        b = 1
        while b < m: b *= 2
        b /= 4
        r = p
        while b:
            if m&b: q, r = ecdub(q, A, n), ecadd(q, r, p, n)
            else:   q, r = ecadd(r, q, p, n), ecdub(r, A, n)
            b /= 2
        return r
def ecm(n, B1=10, B2=20):       # TODO: Determine the best defaults for B1 and B2 and the best way to increment them and iters
    # "Modern" ECM using Montgomery curves and an algorithm analogous to the two-phase variant of Pollard's p-1 method
    # TODO: We currently compute the prime lists from the sieve as we need them, but this means that we recompute them at every
    #       iteration.  While it would not be particularly efficient memory-wise, we might be able to increase time-efficiency
    #       by computing the primes we need ahead of time (say once at the beginning and then once each time we increase the
    #       bounds) and saving them in lists, and then iterate the inner while loops over those lists.
    if isprime(n): return n
    m = ispower(n)
    if m: return m
    iters = 1
    while True:
        for _ in xrange(iters):     # TODO: multiprocessing?
            # TODO: Do we really want to call the randomizer?  Why not have seed be a function of B1, B2, and iters?
            # TODO: Are some seeds better than others?
            seed = randrange(6, n)
            u, v = (seed**2 - 5) % n, 4*seed % n
            p = pow(u, 3, n)
            Q, C = (pow(v-u,3,n)*(3*u+v) % n, 4*p*v % n), (p, pow(v,3,n))
            pg = primegen()
            p = pg.next()
            while p <= B1: Q, p = ecmul(p**ilog(B1, p), Q, C, n), pg.next()
            g = gcd(Q[1], n)
            if 1 < g < n: return g
            while p <= B2:
                # "There is a simple coding trick that can speed up the second stage. Instead of multiplying each prime times Q,
                # we iterate over i from B1 + 1 to B2, adding 2Q at each step; when i is prime, the current Q can be accumulated
                # into the running solution. Again, we defer the calculation of the greatest common divisor until the end of the
                # iteration."                                                TODO: Implement this trick and compare performance.
                Q = ecmul(p, Q, C, n)
                g *= Q[1]
                g %= n
                p = pg.next()
            g = gcd(g, n)
            if 1 < g < n: return g
            # This seed failed.  Try again with a new one.
        # These bounds failed.  Increase and try again.
        B1 *= 3
        B2 *= 3
        iters *= 2


# legendre symbol (a|m)
# TODO: which is faster?
def legendre1(a, p): return ((pow(a, (p-1) >> 1, p) + 1) % p) - 1
def legendre2(a, p):                                                 # TODO: pretty sure this computes the Jacobi symbol
    if a == 0: return 0
    x, y, L = a, p, 1
    while 1:
        if x > (y >> 1):
            x = y - x
            if y & 3 == 3: L = -L
        while x & 3 == 0: x >>= 2
        if x & 1 == 0:
            x >>= 1
            if y & 7 == 3 or y & 7 == 5: L = -L
        if x == 1: return ((L+1) % p) - 1
        if x & 3 == 3 and y & 3 == 3: L = -L
        x, y = y % x, x
if mpzv == 0: legendre = legendre1
else:
    if mpzv == 1: from gmpy  import legendre as legendre0
    if mpzv == 2: from gmpy2 import legendre as legendre0
    def legendre(n, p): return legendre0(n, p) if (n > 0) and (p % 2 == 1) else legendre1(n, p)
    def legendre(n, p): return legendre0(n, p) if (n > 0) and (p % 2 == 1) else legendre1(n, p)

# modular sqrt(n) mod p
# p must be prime
def mod_sqrt(n, p):
    a = n%p
    if p%4 == 3: return pow(a, (p+1) >> 2, p)
    elif p%8 == 5:
        v = pow(a << 1, (p-5) >> 3, p)
        i = ((a*v*v << 1) % p) - 1
        return (a*v*i)%p
    elif p%8 == 1: # Shank's method
        q, e = p-1, 0
        while q&1 == 0:
            e += 1
            q >>= 1
        n = 2
        while legendre(n, p) != -1: n += 1
        w, x, y, r = pow(a, q, p), pow(a, (q+1) >> 1, p), pow(n, q, p), e
        while True:
            if w == 1: return x
            v, k = w, 0
            while v != 1 and k+1 < r:
                v = (v*v)%p
                k += 1
            if k == 0: return x
            d = pow(y, 1 << (r-k-1), p)
            x, y = (x*d)%p, (d*d)%p
            w, r = (w*y)%p, k
    else: return a # p == 2

# modular inverse of a mod m
def modinv(a, m):
    a, x, u = a%m, 0, 1
    while a: x, u, m, a = u, x - (m/a)*u, a, m%a
    return x

# Multiple Polynomial Quadratic Sieve
# Most of this function is copied verbatim from https://codegolf.stackexchange.com/questions/8629/9088#9088
def mpqs(n):
    # When the bound proves insufficiently large, we throw out all our work and start over.
    # TODO: When this happens, get more data, but don't trash what we already have.
    # TODO: Rewrite to get a few more relations before proceeding to the linear algebra.
    # TODO: When we need to increase the bound, what is the optimal increment?
    
    # Special cases: this function poorly handles primes and perfect powers:
    m = ispower(n)
    if m: return m
    if isprime(n): return n
    
    root_n, root_2n = isqrt(n), isqrt(2*n)
    bound = ilog(n**6, 10)**2  # formula chosen by experiment
    
    while True:
        try:
            prime, mod_root, log_p, num_prime = [], [], [], 0
            
            # find a number of small primes for which n is a quadratic residue
            p = 2
            while p < bound or num_prime < 3:
                leg = legendre(n%p, p)
                if leg == 1:
                    prime += [p]
                    mod_root += [mod_sqrt(n, p)]    # the rhs was [int(mod_sqrt(n, p))].  If we get errors, put it back.
                    log_p += [log(p, 10)]
                    num_prime += 1
                elif leg == 0: return p
                p = nextprime(p)
            
            x_max = len(prime)*60    # size of the sieve
            
            m_val = (x_max * root_2n) >> 1    # maximum value on the sieved range
            
            # fudging the threshold down a bit makes it easier to find powers of primes as factors
            # as well as partial-partial relationships, but it also makes the smoothness check slower.
            # there's a happy medium somewhere, depending on how efficient the smoothness check is
            thresh = log(m_val, 10) * 0.735
            
            # skip small primes. they contribute very little to the log sum
            # and add a lot of unnecessary entries to the table
            # instead, fudge the threshold down a bit, assuming ~1/4 of them pass
            min_prime = mpz(thresh*3)
            fudge = sum(log_p[i] for i,p in enumerate(prime) if p < min_prime)/4
            thresh -= fudge
            
            smooth, used_prime, partial = [], set(), {}
            num_smooth, num_used_prime, num_partial, num_poly, root_A = 0, 0, 0, 0, isqrt(root_2n / x_max)
            
            while num_smooth <= num_used_prime:
                # find an integer value A such that:
                # A is =~ sqrt(2*n) / x_max
                # A is a perfect square
                # sqrt(A) is prime, and n is a quadratic residue mod sqrt(A)
                while True:
                    root_A = nextprime(root_A)
                    leg = legendre(n, root_A)
                    if leg == 1: break
                    elif leg == 0: return root_A
                
                A = root_A**2
                
                # solve for an adequate B
                # B*B is a quadratic residue mod n, such that B*B-A*C = n
                # this is unsolvable if n is not a quadratic residue mod sqrt(A)
                b = mod_sqrt(n, root_A)
                B = (b + (n - b*b) * modinv(b + b, root_A))%A
                C = (B*B - n) / A        # B*B-A*C = n <=> C = (B*B-n)/A
                
                num_poly += 1
                
                # sieve for prime factors
                sums, i = [0.0]*(2*x_max), 0
                for p in prime:
                    if p < min_prime:
                        i += 1
                        continue
                    logp = log_p[i]
                    inv_A = modinv(A, p)
                    # modular root of the quadratic
                    a, b, k = mpz(((mod_root[i] - B) * inv_A)%p), mpz(((p - mod_root[i] - B) * inv_A)%p), 0
                    while k < x_max:
                        if k+a < x_max: sums[k+a] += logp
                        if k+b < x_max: sums[k+b] += logp
                        if k:
                            sums[k-a+x_max] += logp
                            sums[k-b+x_max] += logp
                        k += p
                    i += 1
                
                # check for smooths
                i = 0
                for v in sums:
                    if v > thresh:
                        x, vec, sqr = x_max-i if i > x_max else i, set(), []
                        # because B*B-n = A*C
                        # (A*x+B)^2 - n = A*A*x*x+2*A*B*x + B*B - n
                        #               = A*(A*x*x+2*B*x+C)
                        # gives the congruency
                        # (A*x+B)^2 = A*(A*x*x+2*B*x+C) (mod n)
                        # because A is chosen to be square, it doesn't need to be sieved
                        val = sieve_val = (A*x + 2*B)*x + C
                        if sieve_val < 0: vec, sieve_val = {-1}, -sieve_val
                        
                        for p in prime:
                            while sieve_val%p == 0:
                                if p in vec: sqr += [p] # track perfect sqr facs to avoid sqrting something huge at the end
                                vec ^= {p}
                                sieve_val = mpz(sieve_val / p)
                        if sieve_val == 1: # smooth
                            smooth += [(vec, (sqr, (A*x+B), root_A))]
                            used_prime |= vec
                        elif sieve_val in partial:
                            # combine two partials to make a (xor) smooth
                            # that is, every prime factor with an odd power is in our factor base
                            pair_vec, pair_vals = partial[sieve_val]
                            sqr += list(vec & pair_vec) + [sieve_val]
                            vec ^= pair_vec
                            smooth += [(vec, (sqr + pair_vals[0], (A*x+B)*pair_vals[1], root_A*pair_vals[2]))]
                            used_prime |= vec
                            num_partial += 1
                        else: partial[sieve_val] = (vec, (sqr, A*x+B, root_A))      # save partial for later pairing
                    i += 1
                
                num_smooth, num_used_prime = len(smooth), len(used_prime)
            
            used_prime = sorted(list(used_prime))
            
            # set up bit fields for gaussian elimination
            masks, mask, bitfields = [], 1, [0]*num_used_prime
            for vec, vals in smooth:
                masks += [mask]
                i = 0
                for p in used_prime:
                    if p in vec: bitfields[i] |= mask
                    i += 1
                mask <<= 1
            
            # row echelon form
            offset = 0
            null_cols = []
            for col in xrange(num_smooth):
                pivot = bitfields[col-offset] & masks[col] == 0 # This occasionally throws IndexErrors.
                # TODO: figure out why it throws errors and fix it.
                for row in xrange(col+1-offset, num_used_prime):
                    if bitfields[row] & masks[col]:
                        if pivot: bitfields[col-offset], bitfields[row], pivot = bitfields[row], bitfields[col-offset], False
                        else: bitfields[row] ^= bitfields[col-offset]
                if pivot:
                    null_cols += [col]
                    offset += 1
            
            # reduced row echelon form
            for row in xrange(num_used_prime):
                mask = bitfields[row] & -bitfields[row]        # lowest set bit
                for up_row in xrange(row):
                    if bitfields[up_row] & mask: bitfields[up_row] ^= bitfields[row]
            
            # check for non-trivial congruencies
            # TODO: if none exist, check combinations of null space columns...
            # if _still_ none exist, sieve more values
            for col in null_cols:
                all_vec, (lh, rh, rA) = smooth[col]
                lhs = lh   # sieved values (left hand side)
                rhs = [rh] # sieved values - n (right hand side)
                rAs = [rA] # root_As (cofactor of lhs)
                i = 0
                for field in bitfields:
                    if field & masks[col]:
                        vec, (lh, rh, rA) = smooth[i]
                        lhs += list(all_vec & vec) + lh
                        all_vec ^= vec
                        rhs += [rh]
                        rAs += [rA]
                    i += 1
                factor = gcd(listprod(rAs)*listprod(lhs) - listprod(rhs), n)
                if 1 < factor < n: return factor
        
        except IndexError: pass
        
        bound *= 1.2

def multifactor(n, methods=(pollardRho_brent, pollard_pm1, williams_pp1, ecm, mpqs), verbose=False):
    # Note that the multiprocing incurs relatively significant overhead.  Only call this if n is proving difficult to factor.
    def factory(method, n, output): output.put((method(n), str(method).split()[1]))
    factors = mpQueue()
    procs = [Process(target=factory, args=(m, n, factors)) for m in methods]
    for p in procs: p.start()
    (f, g) = factors.get()
    for p in procs: p.terminate()
    if verbose:
        names = {"pollardRho_brent":"prb", "pollard_pm1":"p-1", "williams_pp1":"p+1"}
        print "\033[1;31m" + (names[g] if g in names else g) + "\033[;m",
        stdout.flush()
    return f

def primefac(n, trial_limit=1000, rho_rounds=42000, verbose=False,
             methods=(pollardRho_brent, pollard_pm1, williams_pp1, ecm, mpqs)):
    # Obtains a complete factorization of n, yielding the prime factors as they are obtained.
    # If the user explicitly specifies a splitting method, use that method.  Otherwise,
    # 1.  Pull out small factors with trial division.
    # TODO: a few rounds of Fermat's method?
    # 2.  Do a few rounds of Pollard's Rho algorithm.
    # TODO: a few rounds of ECM by itself?
    # TODO: a certain amount of P-1?
    # 3.  Launch multifactor on the remainder.  Multifactor has enough overhead that we want to be fairly sure that rho isn't
    #     likely to yield new factors soon.  The default value of rho_rounds=42000 seems good for that but is probably overkill.
    
    if n < 2: return
    if isprime(n): yield n; return
    
    factors, nroot = [], isqrt(n)
    for p in primegen(): # Note that we remove factors of 2 whether the user wants to or not.
        if n%p == 0:
            while n%p == 0:
                yield p
                n /= p
            nroot = isqrt(n)
            if isprime(n):
                yield n
                return
        if p > nroot:
            if n != 1: yield n
            return
        if p >= trial_limit: break
    if isprime(n): yield n; return
    
    if rho_rounds == "inf":
        factors = [n]
        while len(factors) != 0:
            n = min(factors)
            factors.remove(n)
            f = pollardRho_brent(n)
            if isprime(f): yield f
            else: factors.append(f)
            n /= f
            if isprime(n): yield n
            else: factors.append(n)
        return
    
    factors, difficult = [n], []
    while len(factors) != 0:
        rhocount = 0
        n = factors.pop()
        try:
            g = n
            while g == n:
                x, c, g = randrange(1, n), randrange(1, n), 1
                y = x
                while g==1:
                    if rhocount >= rho_rounds: raise Exception
                    rhocount += 1
                    x = (x**2 + c) % n
                    y = (y**2 + c) % n
                    y = (y**2 + c) % n
                    g = gcd(x-y, n)
            # We now have a nontrivial factor g of n.  If we took too long to get here, we're actually at the except statement.
            if isprime(g): yield g
            else: factors.append(g)
            n /= g
            if isprime(n): yield n
            else: factors.append(n)
        except Exception: difficult.append(n) # Factoring n took too long.  We'll have multifactor chug on it.
    
    factors = difficult
    while len(factors) != 0:
        n = min(factors)
        factors.remove(n)
        f = multifactor(n, methods=methods, verbose=verbose)
        if isprime(f): yield f
        else: factors.append(f)
        n /= f
        if isprime(n): yield n
        else: factors.append(n)

def factorint(n, trial_limit=1000, rho_rounds=42000, methods=(pollardRho_brent, pollard_pm1, williams_pp1, ecm, mpqs)):
    out = {}
    for p in primefac(n, trial_limit=trial_limit, rho_rounds=rho_rounds, methods=methods): out[p] = out.get(p, 0) + 1
    return out

usage = """
This is primefac version 1.1.

USAGE:
    primefac [-vs|-sv] [-v|--verbose] [-s|--summary] [-t=NUM] [-r=NUM]
          [-m=[prb][,p-1][,p+1][,ecm][,mpqs]] rpn

    "rpn" is evaluated using integer arithmetic.  Each number that remains on
    the stack after evaluation is then factored.

    "-t" is the trial division limit.  Default == 1000.  Use "-t=inf" to use
    trial division exclusively.

    "-r" is the number of rounds of Pollard's rho algorithm to try before
    calling a factor "difficult".  Default == 42,000.  Use "-r=inf" to use
    Pollard's rho exclusively once the trial division is completed.

    If verbosity is invoked, we indicate in the output which algorithm produced
    which factors during the multifactor phase.

    If the summary flag is absent, then output is identical to the output of the
    GNU factor command, except possibly for the order of the factors and, if
    verbosity has been turned on, the annotations indicating which algorithm
    produced which factors.

    If the summary flag is present, then output is modified by adding a single
    newline between each item's output, before the first, and after the last.
    Each item's output is also modified by printing a second line of data
    summarizing the results by describing the number of decimal digits in the
    input, the number of decimal digits in each prime factor, and the factors'
    multiplicities.  For example:

    >>> user@computer:~$ primefac  -s   24 ! 1 -   7 !
    >>> 
    >>> 620448401733239439359999: 991459181683 625793187653
    >>> Z24  =  P12 x P12  =  625793187653 x 991459181683
    >>> 
    >>> 5040: 2 2 2 2 3 3 5 7
    >>> Z4  =  P1^4 x P1^2 x P1 x P1  =  2^4 x 3^2 x 5 x 7
    >>> 
    >>> user@computer:~$

    Note that the primes in the summary lines are listed in strictly-increasing
    order, regardless of the order in which they were found.
    
    The single-character versions of the verbosity and summary flags may be
    combined into a single flag, "-vs" or "-sv".

    The "-m" flag controls what methods are run during the multifactor phase.
    prb and ecm can be listed repeatedly to run multiple instances of these
    methods; running multiple instances of p-1, p+1, or mpqs confers no benefit,
    so repeated listings of those methods are ignored.

    This program can also be imported into your Python scripts as a module.


DETAILS:
    Factoring: 1.  Trial divide using the primes <= the specified limit.
               2.  Run Pollard's rho algorithm on the remainder.  Declare a
                   cofactor "difficult" if it survives more than the specified
                   number of rounds of rho.
               3.  Subject each remaining cofactor to five splitting methods in
                   parallel: Pollard's rho algorithm with Brent's improvement,
                             Pollard's p-1 method,
                             Williams' p+1 method,
                             the elliptic curve method,
                         and the multiple-polynomial quadratic sieve.
               Using the "verbose" option will cause primefac to report which of
               the various splitting methods separated which factors in stage 3.

    RPN:       The acceptable binary operators are + - * / % **.
               They all have the same meaning as they do in Python source code
               --- i.e., they are addition, subtraction, multiplication, integer
               division, remainder, and exponentiation.
               The acceptable unary operators are ! #.  They are the factorial
               and primorial, respectively.
               There are three aliases: x for *, xx for **, and p! for #.
               You may enclose the RPN expression in quotes if you so desire.


PERFORMANCE:


CREDITS:
    Not much of this code was mine from the start.
     * The MPQS code was copied mostly verbatim from
       https://codegolf.stackexchange.com/questions/8629/9088#9088
     * The functions to manipulate points in the elliptic curve method were
       copied from a reply to the Programming Praxis post at
       http://programmingpraxis.com/2010/04/23/

""" # TODO performance, credits

def rpn(instr):
    stack = []
    for token in instr.split():
        if set(token).issubset("1234567890"): stack.append(int(token))
        elif len(token) > 1 and token[0] == '-' and set(token[1:]).issubset("1234567890"): stack.append(int(token))
        elif token in ('+', '-', '*', '/', '%', '**', 'x', 'xx'):   # binary operators
            b = stack.pop()
            a = stack.pop()
            if   token == '+' : res = a  + b
            elif token == '-' : res = a  - b
            elif token == '*' : res = a  * b
            elif token == 'x' : res = a  * b
            elif token == '/' : res = a  / b
            elif token == '%' : res = a  % b
            elif token == '**': res = a ** b
            elif token == 'xx': res = a ** b
            stack.append(res)
        elif token in ('!', '#', 'p!'):                             # unary operators
            a = stack.pop()
            if   token == '!' : res = listprod(xrange(1, a+1))
            elif token == '#' : res = listprod(primes(a+1))
            elif token == 'p!': res = listprod(primes(a+1))
            stack.append(res)
        else: raise Exception, "Failed to evaluate RPN expression: not sure what to do with '{t}'.".format(t=token)
    return map(mpz, stack)

# TODO timeout?
if __name__ == "__main__":
    from sys import stdout, exit, argv
    if len(argv) == 1: exit(usage)
    start, rpx, tr, rr, veb, su = 1, [], 1000, 42000, False, False
    ms = {"prb":pollardRho_brent, "p-1":pollard_pm1, "p+1":williams_pp1, "ecm":ecm, "mpqs":mpqs}
    methods = (pollardRho_brent, pollard_pm1, williams_pp1, ecm, mpqs)
    try:
        for arg in argv[1:]:
            if arg in ("-v", "--verbose"): veb = True
            elif arg in ("-s", "--summary"): su = True
            elif arg in ("-vs", "-sv"): veb, su = True, True
            elif arg[:3] == "-t=": tr = "inf" if arg[3:] == "inf" else int(arg[3:])    # Maximum number for trial division
            elif arg[:3] == "-r=": rr = "inf" if arg[3:] == "inf" else int(arg[3:])    # Number of rho rounds before multifactor
            elif arg[:3] == "-m=": #methods = tuple(ms[x] for x in arg[3:].split(',') if x in ms)
                methods = []
                for x in arg[3:].split(','):
                    if x in ms:
                        if x in ("p-1", "p+1", "mpqs") and ms[x] in methods: continue
                        methods.append(ms[x])
            else: rpx.append(arg)
        nums = rpn(' '.join(rpx))
        for x in nums: assert isinstance(x, inttypes)
    except: exit("Error while parsing arguments")
    if su: print
    for n in nums:
        print "%d:" % n,
        f = {}
        for p in primefac(n, trial_limit=(n if tr == "inf" else tr), rho_rounds=rr, verbose=veb, methods=methods):
            f[p] = f.get(p, 0) + 1
            print p,
            stdout.flush()
            assert isprime(p) and n%p == 0, (n, p)
        print
        if su:
            print "Z%d  = " % len(str(n)),
            outstr = ""
            for p in sorted(f):
                if f[p] == 1: outstr += "P%d x " % len(str(p))
                else: outstr += "P%d^%d x " % (len(str(p)), f[p])
            outstr = outstr[:-2] + " = "
            for p in sorted(f):
                if f[p] == 1: outstr += " %d x" % p
                else: outstr += " %d^%d x" % (p, f[p])
            print outstr[:-2]
            print

# Fun examples:
# primefac -v 1489576198567193874913874619387459183543154617315437135656
#   On my system, the factor race is a bit unpredicatble on this number.  prb, ecm, p-1, and mpqs all show up reasonably often.
# primefac -v 12956921851257164598146425167654345673426523793463
#   Z50 = P14 x P17 x P20 = 24007127617807 x 28050585032291527 x 19240648901716863967.  p-1 gets the P14 and p+1 gets the rest.
# primefac -v 38 ! 1 +  -->  Z45 = P23 x P23 = 14029308060317546154181 x 37280713718589679646221
#   The MPQS (almost always) gets this one.  Depending on the system running things, this can take from 10 seconds to 3 minutes.

