#! /usr/bin/env python

from __future__ import print_function, division
from _primefac import *

# Formatting note: this file's maximum line length is 128 characters.

# TODO: Python 2.x/3.x compatibility
from six.moves import xrange, reduce
import six

from multiprocessing import Process, Queue as mpQueue
from itertools import count, takewhile
from random import randrange
from math import log

def pollardRho_brent(n):
    if isprime(n):
        return n
    g = n
    while g == n:
        y, c, m, g, r, q = randrange(1, n), randrange(1, n), randrange(1, n), 1, 1, 1
        while g == 1:
            x, k = y, 0
            for _ in xrange(r):
                y = (y**2 + c) % n
            while k < r and g == 1:
                ys = y
                for _ in xrange(min(m, r-k)):
                    y = (y**2 + c) % n
                    q = q * abs(x-y) % n
                g, k = gcd(q, n), k+m
            r *= 2
        if g == n:
            while True:
                ys = (ys**2+c) % n
                g = gcd(abs(x-ys), n)
                if g > 1:
                    break
    return g

# http://programmingpraxis.com/2010/04/27/modern-elliptic-curve-factorization-part-2/
def pollard_pm1(n, B1=100, B2=1000):       # TODO: What are the best default bounds and way to increment them?
    if isprime(n):
        return n
    m = ispower(n)
    if m:
        return m
    while True:
        pg = primegen()
        q = 2           # TODO: what about other initial values of q?
        p = six.next(pg)
        while p <= B1:
            q, p = pow(q, p**ilog(B1, p), n), six.next(pg)
        g = gcd(q-1, n)
        if 1 < g < n:
            return g
        while p <= B2:
            q, p = pow(q, p, n), six.next(pg)
        g = gcd(q-1, n)
        if 1 < g < n:
            return g
        # These bounds failed.  Increase and try again.
        B1 *= 10
        B2 *= 10

def mlucas(v, a, n):
    """ Helper function for williams_pp1().  Multiplies along a Lucas sequence modulo n. """
    v1, v2 = v, (v**2 - 2) % n
    for bit in bin(a)[3:]:
        v1, v2 = ((v1**2 - 2) % n, (v1*v2 - v) % n) if bit == "0" else ((v1*v2 - v) % n, (v2**2 - 2) % n)
    return v1

def williams_pp1(n):
    if isprime(n):
        return n
    m = ispower(n)
    if m:
        return m
    for v in count(1):
        for p in primegen():
            e = ilog(isqrt(n), p)
            if e == 0:
                break
            for _ in xrange(e):
                v = mlucas(v, p, n)
            g = gcd(v - 2, n)
            if 1 < g < n:
                return g
            if g == n:
                break

# http://programmingpraxis.com/2010/04/23/modern-elliptic-curve-factorization-part-1/
# http://programmingpraxis.com/2010/04/27/modern-elliptic-curve-factorization-part-2/
def ecadd(p1, p2, p0, n):  # Add two points p1 and p2 given point P0 = P1-P2 modulo n
    x1, z1 = p1
    x2, z2 = p2
    x0, z0 = p0
    t1, t2 = (x1-z1)*(x2+z2), (x1+z1)*(x2-z2)
    return (z0*pow(t1+t2, 2, n) % n, x0*pow(t1-t2, 2, n) % n)
def ecdub(p, A, n):  # double point p on A modulo n
    x, z = p
    An, Ad = A
    t1, t2 = pow(x+z, 2, n), pow(x-z, 2, n)
    t = t1 - t2
    return (t1*t2*4*Ad % n, (4*Ad*t2 + t*An)*t % n)
def ecmul(m, p, A, n):  # multiply point p by m on curve A modulo n
    if m == 0:
        return (0, 0)
    elif m == 1:
        return p
    else:
        q = ecdub(p, A, n)
        if m == 2:
            return q
        b = 1
        while b < m:
            b *= 2
        b //= 4
        r = p
        while b:
            if m & b:
                q, r = ecdub(q, A, n), ecadd(q, r, p, n)
            else:
                q, r = ecadd(r, q, p, n), ecdub(r, A, n)
            b //= 2
        return r
def ecm(n, B1=10, B2=20):
    # TODO: Determine the best defaults for B1 and B2 and the best way to increment them and iters
    # "Modern" ECM using Montgomery curves and an algorithm analogous to the two-phase variant of Pollard's p-1 method
    # TODO: We currently compute the prime lists from the sieve as we need them, but this means that we recompute them at every
    #       iteration.  While it would not be particularly efficient memory-wise, we might be able to increase time-efficiency
    #       by computing the primes we need ahead of time (say once at the beginning and then once each time we increase the
    #       bounds) and saving them in lists, and then iterate the inner while loops over those lists.
    if isprime(n):
        return n
    m = ispower(n)
    if m:
        return m
    iters = 1
    while True:
        for _ in xrange(iters):     # TODO: multiprocessing?
            # TODO: Do we really want to call the randomizer?  Why not have seed be a function of B1, B2, and iters?
            # TODO: Are some seeds better than others?
            seed = randrange(6, n)
            u, v = (seed**2 - 5) % n, 4*seed % n
            p = pow(u, 3, n)
            Q, C = (pow(v-u, 3, n)*(3*u+v) % n, 4*p*v % n), (p, pow(v, 3, n))
            pg = primegen()
            p = six.next(pg)
            while p <= B1:
                Q, p = ecmul(p**ilog(B1, p), Q, C, n), six.next(pg)
            g = gcd(Q[1], n)
            if 1 < g < n:
                return g
            while p <= B2:
                # "There is a simple coding trick that can speed up the second stage. Instead of multiplying each prime times Q,
                # we iterate over i from B1 + 1 to B2, adding 2Q at each step; when i is prime, the current Q can be accumulated
                # into the running solution. Again, we defer the calculation of the greatest common divisor until the end of the
                # iteration."                                                TODO: Implement this trick and compare performance.
                Q = ecmul(p, Q, C, n)
                g *= Q[1]
                g %= n
                p = six.next(pg)
            g = gcd(g, n)
            if 1 < g < n:
                return g
            # This seed failed.  Try again with a new one.
        # These bounds failed.  Increase and try again.
        B1 *= 3
        B2 *= 3
        iters *= 2


def fermat(n):
    x = isqrt(n) + 1
    y = isqrt(x**2 - n)
    while True:
        w = x**2 - n - y**2
        if w == 0:
            break
        if w > 0:
            y += 1
        else:
            x += 1
    return x+y


# Multiple Polynomial Quadratic Sieve
# Most of this function is copied verbatim from https://codegolf.stackexchange.com/questions/8629/9088#9088
def mpqs(n):
    # When the bound proves insufficiently large, we throw out all our work and start over.
    # TODO: When this happens, get more data, but don't trash what we already have.
    # TODO: Rewrite to get a few more relations before proceeding to the linear algebra.
    # TODO: When we need to increase the bound, what is the optimal increment?

    # Special cases: this function poorly handles primes and perfect powers:
    m = ispower(n)
    if m:
        return m
    if isprime(n):
        return n

    root_2n = isqrt(2*n)
    bound = ilog(n**6, 10)**2  # formula chosen by experiment

    while True:
        try:
            prime, mod_root, log_p, num_prime = [], [], [], 0

            # find a number of small primes for which n is a quadratic residue
            p = 2
            while p < bound or num_prime < 3:
                leg = legendre(n % p, p)
                if leg == 1:
                    prime += [p]
                    mod_root += [mod_sqrt(n, p)]    # the rhs was [int(mod_sqrt(n, p))].  If we get errors, put it back.
                    log_p += [log(p, 10)]
                    num_prime += 1
                elif leg == 0:
                    return p
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
            fudge = sum(log_p[i] for i, p in enumerate(prime) if p < min_prime)//4
            thresh -= fudge

            smooth, used_prime, partial = [], set(), {}
            num_smooth, num_used_prime, num_partial, num_poly, root_A = 0, 0, 0, 0, isqrt(root_2n // x_max)

            while num_smooth <= num_used_prime:
                # find an integer value A such that:
                # A is =~ sqrt(2*n) // x_max
                # A is a perfect square
                # sqrt(A) is prime, and n is a quadratic residue mod sqrt(A)
                while True:
                    root_A = nextprime(root_A)
                    leg = legendre(n, root_A)
                    if leg == 1:
                        break
                    elif leg == 0:
                        return root_A

                A = root_A**2

                # solve for an adequate B
                # B*B is a quadratic residue mod n, such that B*B-A*C = n
                # this is unsolvable if n is not a quadratic residue mod sqrt(A)
                b = mod_sqrt(n, root_A)
                B = (b + (n - b*b) * modinv(b + b, root_A)) % A
                C = (B*B - n) // A        # B*B-A*C = n <=> C = (B*B-n)//A

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
                    a, b, k = (mpz(((mod_root[i] - B) * inv_A) % p),
                               mpz(((p - mod_root[i] - B) * inv_A) % p),
                               0)
                    while k < x_max:
                        if k+a < x_max:
                            sums[k+a] += logp
                        if k+b < x_max:
                            sums[k+b] += logp
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
                        sieve_val = (A*x + 2*B)*x + C
                        if sieve_val < 0:
                            vec, sieve_val = {-1}, -sieve_val

                        for p in prime:
                            while sieve_val % p == 0:
                                if p in vec:
                                    sqr += [p]  # track perfect sqr facs to avoid sqrting something huge at the end
                                vec ^= {p}
                                sieve_val = mpz(sieve_val // p)
                        if sieve_val == 1:  # smooth
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
                        else:
                            partial[sieve_val] = (vec, (sqr, A*x+B, root_A))      # save partial for later pairing
                    i += 1

                num_smooth, num_used_prime = len(smooth), len(used_prime)

            used_prime = sorted(list(used_prime))

            # set up bit fields for gaussian elimination
            masks, mask, bitfields = [], 1, [0]*num_used_prime
            for vec, _ in smooth:
                masks += [mask]
                i = 0
                for p in used_prime:
                    if p in vec:
                        bitfields[i] |= mask
                    i += 1
                mask <<= 1

            # row echelon form
            offset = 0
            null_cols = []
            for col in xrange(num_smooth):
                pivot = bitfields[col-offset] & masks[col] == 0  # This occasionally throws IndexErrors.
                # TODO: figure out why it throws errors and fix it.
                for row in xrange(col+1-offset, num_used_prime):
                    if bitfields[row] & masks[col]:
                        if pivot:
                            bitfields[col-offset], bitfields[row], pivot = bitfields[row], bitfields[col-offset], False
                        else:
                            bitfields[row] ^= bitfields[col-offset]
                if pivot:
                    null_cols += [col]
                    offset += 1

            # reduced row echelon form
            for row in xrange(num_used_prime):
                mask = bitfields[row] & -bitfields[row]        # lowest set bit
                for up_row in xrange(row):
                    if bitfields[up_row] & mask:
                        bitfields[up_row] ^= bitfields[row]

            # check for non-trivial congruencies
            # TODO: if none exist, check combinations of null space columns...
            # if _still_ none exist, sieve more values
            for col in null_cols:
                all_vec, (lh, rh, rA) = smooth[col]
                lhs = lh   # sieved values (left hand side)
                rhs = [rh]  # sieved values - n (right hand side)
                rAs = [rA]  # root_As (cofactor of lhs)
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
                if 1 < factor < n:
                    return factor

        except IndexError:
            pass

        bound *= 1.2

# Note that the multiprocing incurs relatively significant overhead.
# Only call this if n is proving difficult to factor.
def multifactor(n, methods=(pollardRho_brent, pollard_pm1, williams_pp1,
                ecm, mpqs, fermat), verbose=False):
    def factory(method, n, output):
        output.put((method(n), str(method).split()[1]))
    factors = mpQueue()
    procs = [Process(target=factory, args=(m, n, factors)) for m in methods]
    for p in procs:
        p.start()
    (f, g) = factors.get()
    for p in procs:
        p.terminate()
    if verbose:
        names = {"pollardRho_brent": "prb",
                 "pollard_pm1": "p-1",
                 "williams_pp1": "p+1"}
        if g in names:
            name = names[g]
        else:
            name = g
        print("\033[1;31m" + name + "\033[;m", end=' ')
        stdout.flush()
    return f

'''
Obtains a complete factorization of n, yielding the prime factors as they are
obtained. If the user explicitly specifies a splitting method, use that method.
Otherwise,
1.  Pull out small factors with trial division.
2.  Do a few rounds of Pollard's Rho algorithm.
    TODO: a few rounds of ECM by itself?
    TODO: a certain amount of P-1?
3.  Launch multifactor on the remainder.  Multifactor has enough overhead that
    we want to be fairly sure that rho isn't likely to yield new factors soon. 
    The default value of rho_rounds=42000 seems good for that but is probably
    overkill.
'''

def primefac(n, trial_limit=1000, rho_rounds=42000, verbose=False,
             methods=(pollardRho_brent, pollard_pm1, williams_pp1, ecm, mpqs,
                      fermat)):
    if n < 2:
        return
    if isprime(n):
        yield n
        return

    factors, nroot = [], isqrt(n)
    # Note that we remove factors of 2 whether the user wants to or not.
    for p in primegen():
        if n % p == 0:
            while n % p == 0:
                yield p
                n //= p
            nroot = isqrt(n)
            if isprime(n):
                yield n
                return
        if p > nroot:
            if n != 1:
                yield n
            return
        if p >= trial_limit:
            break
    if isprime(n):
        yield n
        return

    if rho_rounds == "inf":
        factors = [n]
        while len(factors) != 0:
            n = min(factors)
            factors.remove(n)
            f = pollardRho_brent(n)
            if isprime(f):
                yield f
            else:
                factors.append(f)
            n //= f
            if isprime(n):
                yield n
            else:
                factors.append(n)
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
                while g == 1:
                    if rhocount >= rho_rounds:
                        raise Exception
                    rhocount += 1
                    x = (x**2 + c) % n
                    y = (y**2 + c) % n
                    y = (y**2 + c) % n
                    g = gcd(x-y, n)
            # We now have a nontrivial factor g of n.  If we took too long to get here, we're actually at the except statement.
            if isprime(g):
                yield g
            else:
                factors.append(g)
            n //= g
            if isprime(n):
                yield n
            else:
                factors.append(n)
        except Exception:
            difficult.append(n)  # Factoring n took too long.  We'll have multifactor chug on it.

    factors = difficult
    while len(factors) != 0:
        n = min(factors)
        factors.remove(n)
        f = multifactor(n, methods=methods, verbose=verbose)
        if isprime(f):
            yield f
        else:
            factors.append(f)
        n //= f
        if isprime(n):
            yield n
        else:
            factors.append(n)

def factorint(n, trial_limit=1000, rho_rounds=42000, methods=(pollardRho_brent, pollard_pm1, williams_pp1, ecm, mpqs, fermat)):
    out = {}
    for p in primefac(n, trial_limit=trial_limit, rho_rounds=rho_rounds, methods=methods):
        out[p] = out.get(p, 0) + 1
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
                             the multiple-polynomial quadratic sieve,
                             and fermat's factorization method.
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
"""  # TODO performance, credits

def rpn(instr):
    stack = []
    for token in instr.split():
        if set(token).issubset("1234567890"):
            stack.append(int(token))
        elif len(token) > 1 and token[0] == '-' and set(token[1:]).issubset("1234567890"):
            stack.append(int(token))
        elif token in ('+', '-', '*', '/', '%', '**', 'x', 'xx'):   # binary operators
            b = stack.pop()
            a = stack.pop()
            if token == '+':
                res = a + b
            elif token == '-':
                res = a - b
            elif token == '*':
                res = a * b
            elif token == 'x':
                res = a * b
            elif token == '/':
                res = a / b
            elif token == '%':
                res = a % b
            elif token == '**':
                res = a ** b
            elif token == 'xx':
                res = a ** b
            stack.append(res)
        elif token in ('!', '#', 'p!'):                             # unary operators
            a = stack.pop()
            if token == '!':
                res = listprod(xrange(1, a+1))
            elif token == '#':
                res = listprod(primes(a+1))
            elif token == 'p!':
                res = listprod(primes(a+1))
            stack.append(res)
        else:
            raise Exception("Failed to evaluate RPN expression: not sure what to do with '{t}'.".format(t=token))
    return [mpz(i) for i in stack]

def main(argv):
    if len(argv) == 1:
        sysexit(usage)
    rpx, tr, rr, veb, su = [], 1000, 42000, False, False
    ms = {"prb": pollardRho_brent,
          "p-1": pollard_pm1,
          "p+1": williams_pp1,
          "ecm": ecm,
          "mpqs": mpqs,
          "fermat": fermat}
    methods = (pollardRho_brent, pollard_pm1, williams_pp1, ecm, mpqs, fermat)
    try:
        for arg in argv[1:]:
            if arg in ("-v", "--verbose"):
                veb = True
            elif arg in ("-s", "--summary"):
                su = True
            elif arg in ("-vs", "-sv"):
                veb, su = True, True
            elif arg[:3] == "-t=":
                tr = "inf" if arg[3:] == "inf" else int(arg[3:])    # Maximum number for trial division
            elif arg[:3] == "-r=":
                rr = "inf" if arg[3:] == "inf" else int(arg[3:])    # Number of rho rounds before multifactor
            elif arg[:3] == "-m=":  # methods = tuple(ms[x] for x in arg[3:].split(',') if x in ms)
                methods = []
                for x in arg[3:].split(','):
                    if x in ms:
                        if x in ("p-1", "p+1", "mpqs") and ms[x] in methods:
                            continue
                        methods.append(ms[x])
            else:
                rpx.append(arg)
        nums = rpn(' '.join(rpx))
    except:
        sysexit("Error while parsing arguments")
    if su:
        print()
    for n in nums:
        print("%d: " % n, end='')
        f = {}
        for p in primefac(n, trial_limit=(n if tr == "inf" else tr),
                            rho_rounds=rr, verbose=veb, methods=methods):
            f[p] = f.get(p, 0) + 1
            print(p, end=' ')
            stdout.flush()
            assert isprime(p) and n % p == 0, (n, p)
        print()
        if su:
            print("Z%d  = " % len(str(n)), end='')
            outstr = ""
            for p in sorted(f):
                if f[p] == 1:
                    outstr += "P%d x " % len(str(p))
                else:
                    outstr += "P%d^%d x " % (len(str(p)), f[p])
            outstr = outstr[:-2] + " = "
            for p in sorted(f):
                if f[p] == 1:
                    outstr += " %d x" % p
                else:
                    outstr += " %d^%d x" % (p, f[p])
            print(outstr[:-2])
            print()

'''
main(['p', '-s',
'1489576198567193874913874619387459183543154617315437135656']) only test
'''

# TODO timeout?
if __name__ == "__main__":
    from sys import argv as arguments, stdout, exit as sysexit
    main(arguments)

'''
Fun examples:
primefac -v 1489576198567193874913874619387459183543154617315437135656
  On my system, the factor race is a bit unpredicatble on this number.
  prb, ecm, p-1, and mpqs all show up reasonably often.
primefac -v 12956921851257164598146425167654345673426523793463
  Z50 = P14 x P17 x P20 =
      24007127617807 x 28050585032291527 x 19240648901716863967.
  p-1 gets the P14 and p+1 gets the rest.
primefac -v 38 ! 1 +
    -->  Z45 = P23 x P23 = 14029308060317546154181 x 37280713718589679646221
  The MPQS (almost always) gets this one.
    Depending on the system running things,
        this can take from 10 seconds to 3 minutes.
'''
