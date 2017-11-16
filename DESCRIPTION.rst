For further updates, see the `primefac`__ package.
==================================================
__ https://pypi.python.org/pypi/primefac


primefac version 1.1
======================

This is a module and command-line utility for factoring integers.  As a module, we provide a primality test, several functions for extracting a non-trivial factor of an integer, and a generator that yields all of a number's prime factors (with multiplicity).  As a command-line utility, this project aims to replace GNU's ``factor`` command with a more versatile utility --- in particular, this utility can operate on arbitrarily large numbers, uses multiple cores in parallel, uses better algorithms, handles input in reverse Polish notation, and can be tweaked via command-line flags.  Specifically,

 - GNU's ``factor`` command won't factor anything greater than 2\ :sup:`127`\ -1.  primefac handles arbitrarily large integers.  If available, ``gmpy`` or ``gmpy2`` is imported, with the latter taking precedence over the former.
 - GNU's ``factor`` command uses Pollard's rho algorithm.  While this extracts small factors quickly, large factors take a while to find.  primefac uses, among other things, the elliptic curve method, which is far more efficient at extracting large factors.
 - GNU's ``factor`` command is a single-threaded application.  primefac uses by default five threads to take advantage of the multiple cores typically available on modern machines.  Each of these threads uses a different algorithm to factor the number:

   - One thread runs Brent's variation on Pollard's rho algorithm.  This is good for extracting smallish factors quickly.
   - One thread runs the two-stage version of Pollard's *p*\ -1 method.  This is good at finding factors *p* for which *p*\ -1 is a product of small primes.
   - One thread runs Williams' *p*\ +1 method.  This is good at finding factors *p* for which *p*\ +1 is a product of small primes.
   - One thread runs the elliptic curve method.  This is a bit slower than Pollard's rho algorithm when the factors extracted are small, but it has significantly better performance on difficult factors.
   - One thread runs the multiple-polynomial quadratic sieve.  This is the best algorithm for factoring "hard" numbers short of the horrifically complex general number field sieve.  However, it's (relatively speaking) more than a little slow when the numbers are small, and the time it takes depends only on the size of the number being factored rather than the size of the factors being extracted as with Pollard's rho algorithm and the elliptic curve method, so we use the preceding algorithms to handle those.

 - We also extend the utility by interpreting the command-line arguments as an expression in reverse Polish notation and factoring the numbers remaining on the evaluation stack when interpretation is complete.  For example, the command::

    python -m primefac 24 ! 1 - 38 ! 1 +

  will factor the numbers 24! - 1 = 620448401733239439359999 and 38! + 1 = 523022617466601111760007224100074291200000001.


Module Usage
============
The primary functions are ``isprime`` and ``primefac``, but we define a number of helper functions along the way.

.. code:: python

    gcd(a, b)

Computes the greatest common divisor of the integers ``a`` and ``b``.

.. code:: python

    isqrt(n)

Computes the greatest integer whose square does not exceed the non-negative integer ``n``.

.. code:: python

    introot(n, r=2)

For non-negative ``n``, returns the greatest integer less than or equal to the ``r``\ :sup:`th`\  root of ``n``.

For negative ``n``, returns the least integer greater than or equal to the ``r``\ :sup:`th`\  root of ``n``, or ``None`` if ``r`` is even.

.. code:: python

    primegen()

Non-terminating generator.  Yields the prime numbers.  It amounts to a recursive Sieve of Eratosthenes.  Memory usage is on the order of the square root of the most-recently-yielded prime.  See `this Programming Praxis post`__ for more about the algorithm.

__ http://programmingpraxis.com/2015/07/31/incremental-sieve-of-eratosthenes/

.. code:: python

    primes(n)

Returns a list of the primes strictly less than ``n``.

.. code:: python

    listprod(l)

Returns the product of the elements of ``l``, which can be any iterable (but should obviously terminate; e.g., ``listprod(primegen())`` would be a bad idea).

.. code:: python

    nextprime(n)

Determines, with some semblance of efficiency, the least prime number strictly greater than ``n``.

.. code:: python

    sprp(n, a, s=None, d=None)

Checks ``n`` for primality using the Strong Probable Primality Test to base ``a``.  If present, ``s`` and ``d`` should be the first and second items, respectively, of the tuple returned by the function ``pfactor(n)``.  We use this as a helper function for ``isprime``.

.. code:: python

    pfactor(n)

Helper function for ``sprp``.  Returns the tuple ``(x,y)`` where ``n - 1 == (2 ** x) * y`` and ``y`` is odd.  We have this bit separated out so that we don't waste time recomputing ``s`` and ``d`` for each base when we want to check ``n`` against multiple bases.

.. code:: python

    jacobi(a, p)

Computes the Jacobi symbol ``(a|p)``, where ``p`` is a positive odd number.  This is used in ``isprime``.

.. code:: python

    chain(n, u1, v1, u2, v2, d, q, m)

Helper function for ``isprime``.

.. code:: python

    isprime(n, tb=(3,5,7,11), eb=(2,), mrb=())

The main primality test.  It's an implementation of the BPSW test (Baillie-Pomerance-Selfridge-Wagstaff) with some prefiltes for speed and is deterministic for all numbers less than 2\ :sup:`64` --- in fact, while infinitely many false positives are conjectured to exist, no false positives are currently known.  The prefilters consist of trial division against 2 and the elements of the tuple ``tb``, checking whether ``n`` is square, and Euler's primality test to the bases in the tuple ``eb``.  If the number is less than 3825123056546413051, we use the Miller-Rabin test on a set of bases for which the test is known to be deterministic over this range.

.. code:: python

    ilog(x, b)

Returns the greatest integer ``l`` such that  ``b**l <= x``.

.. code:: python

    ispower(n)

Returns the largest integer that, when squared/cubed/etc, yields ``n``, or 0 if no such integer exists.  Note that the power to which this number is raised will be prime.

.. code:: python

    pollardRho_brent(n)

Brent's improvement on Pollard's rho algorithm.  Returns ``n`` if ``n`` is prime; otherwise, we keep chugging until we find a factor of ``n`` strictly between ``1`` and ``n``.

.. code:: python

    pollard_pm1(n, B1=100, B2=1000)

Pollard's *p*\ +1 algorithm, two-phase version.  Returns ``n`` if ``n`` is prime; otherwise, we keep chugging until we find a factor of ``n`` strictly between ``1`` and ``n``.

.. code:: python

    mlucas(v, a, n)

Helper function for ``williams_pp1``.  Multiplies along a Lucas sequence modulo ``n``.

.. code:: python

    williams_pp1(n)

Williams' *p*\ +1 algorithm.  Returns ``n`` if ``n`` is prime; otherwise, we keep chugging until we find a factor of ``n`` strictly between ``1`` and ``n``.

.. code:: python

    ecadd(p1, p2, p0, n)

Helper function for ``ecm``.  Adds two points ``p1`` and ``p2`` given point ``P0 = P1-P2`` modulo ``n``.

.. code:: python

    ecdub(p, A, n)

Helper function for ``ecm``.  Doubles point ``p`` on ``A`` modulo ``n``.

.. code:: python

    ecmul(m, p, A, n)

Helper function for ``ecm``.  Multiplies point ``p`` by ``m`` on curve ``A`` modulo ``n``.

.. code:: python

    ecm(n, B1=10, B2=20)

Factors ``n`` using the elliptic curve method, using Montgomery curves and an algorithm analogous to the two-phase variant of Pollard's *p*-1 method.  Returns ``n`` if ``n`` is prime; otherwise, we keep chugging until we find a factor of ``n`` strictly between ``1`` and ``n``.  For more details see `these`_ `two`_ Programming Praxis posts.

.. _these: http://programmingpraxis.com/2010/04/23/modern-elliptic-curve-factorization-part-1/
.. _two: http://programmingpraxis.com/2010/04/27/modern-elliptic-curve-factorization-part-2/

.. code:: python

    legendre(a,p), legendre1(a,p), legendre2(a,p)

Functions to comptue the Legendre symbol ``(a|p)``.  The return value isn't meaningful if ``p`` is composite.  We have three functions for this becaues of the details of the corresponding function in ``gmpy`` and how it's accessed.

.. code:: python

    mod_sqrt(n, p)

Computes a square root of ``n`` modulo the prime number ``p``.  The return value is not meaningful if ``n`` has no square root modulo ``p`` or if ``p`` is composite.

.. code:: python

    modinv(a, m)

Computes a multiplicative inverse of ``a`` modulo ``m``.  The return value is not meaningful if ``gcd(a,m) != 1``.

.. code:: python

    mpqs(n)

Factors ``n`` using the multiple polynomial quadratic sieve.  Returns ``n`` if ``n`` is prime; otherwise, we keep chugging until we find a factor of ``n`` strictly between ``1`` and ``n``.  This function was copied mostly verbatim from `this stackexchange post`__.

__ https://codegolf.stackexchange.com/questions/8629/9088#9088

.. code:: python

    multifactor(n, methods=(pollardRho_brent, pollard_pm1, williams_pp1, ecm, mpqs), verbose=False)

Runs several factoring algorithms on ``n`` simultaneously by loading them into their own threads via the ``multiprocessing`` module.  When one function returns, everything is killed off and that value gets returned.

.. code:: python

    primefac(n, trial_limit=1000, rho_rounds=42000, verbose=False,
             methods=(pollardRho_brent, pollard_pm1, williams_pp1, ecm, mpqs))

Generator.  Yields the prime factors of ``n``, with multiplicity.

.. code:: python

    factorint(n, trial_limit=1000, rho_rounds=42000,
              methods=(pollardRho_brent, pollard_pm1, williams_pp1, ecm, mpqs))

This collates ``primefac``'s output into a dict with the prime factors as the keys and their multiplicities as the data.  For example, ``factorint(5040)`` returns ``{2:4, 3:2, 5:1, 7:1}``.

.. code:: python

    rpn(instr)

Evaluates the string ``instr`` as an expression in reverse Polish notation.


Dependencies
------------

This package imports items from ``multiprocessing``, ``random``, ``itertools``, and ``math``.  These are all in the Python standard library.

We attempt to import items from ``gmpy2`` (or, failing that, ``gmpy``), but these packages are not necessary: the GMPY functions that would be imported are implemented natively if the import fails.


Command-Line Usage
==================

.. code:: sh

    python -m primefac [-vs] [-v|--verbose] [-s|--summary] [-t=NUM] [-r=NUM]
                    [-m=[prb][,p-1][,p+1][,ecm][,mpqs]] rpn

``rpn`` is an expression in revese Polish notation and is evaluated using integer arithmetic.  Each number that remains on the stack after evaluation is then factored.

``-t`` sets the trial division limit; the default value is 1000.  Use ``-t=inf`` to use trial division exclusively.

``-r`` sets the number of rounds of Pollard's rho algorithm to try before calling a factor "difficult".  The default value is 42,000.  Use ``-r=inf`` to use Pollard's rho exclusively once the trial division is completed.

If verbosity is invoked, we indicate in the output which algorithm produced which factors during the multifactor phase.

If the ``-s`` (or ``--summary``) flag is absent, then output is identical to the output of the GNU ``factor`` command, except possibly for the order of the factors and, if verbosity has been turned on, the annotations indicating which algorithm produced which factors.

If the ``-s`` (or ``--summary``) flag is present, then output is modified by adding a single newline between each item's output, before the first item, and after the last item.  Each item's output is also modified by printing a second line of data summarizing the results by describing the number of decimal digits in the input, the number of decimal digits in each prime factor, and the factors' multiplicities.  For example::

    >>> user@computer:~$ python -m primefac  -sv   24 ! 1 -   7 !
    >>> 
    >>> 620448401733239439359999: ecm 991459181683 625793187653
    >>> Z24  =  P12 x P12  =  625793187653 x 991459181683
    >>> 
    >>> 5040: 2 2 2 2 3 3 5 7
    >>> Z4  =  P1^4 x P1^2 x P1 x P1  =  2^4 x 3^2 x 5 x 7
    >>> 
    >>> user@computer:~$

Note that the primes in the summary lines are listed in strictly-increasing order, regardless of the order in which they were found.

The ``-v`` and ``-s`` flags may be combined into a single flag in either order --- i.e., into ``-vs`` or ``-sv``.

The `-m=` flag controls the functions used during the ``multifactor`` phase.  The options are ``prb``, ``p-1``, ``p+1``, ``ecm``, and ``mpqs``, representing Pollard's rho, Pollard's *p*\ -1, Williams' *p*\ +1, the elliptic curve method, and the multiple polynomial quadratic sieve, respectively.  The options must be separated by commas.  The options can be repeated: if ``prb`` is listed twice, for example, then ``multifactor`` will run two instances of ``pollardRho_brent`` simultaneously.  In the case of ``prb`` and ``ecm``, this decreases the expectation value of the time to find a factor, whereas the other three algorithms (*p*\ -1, *p*\ +1, and MPQS) have no randomized component so that running duplicate instances of these three algorithms confers no benefit.  We therefore ignore repeated listings of the latter three methods: for example, calling

.. code:: sh

    python -m primefac -m=prb,prb,ecm,ecm,ecm,mpqs,mpqs 38 ! 1 +

will run during the multifactor phase two instances of Pollard's rho, three instances of the elliptic curve method, and one instance of the MQPS.  Invoking more methods than you have cores available is unlikely to confer any benefit.


RPN
---

The acceptable binary operators are ``+``, ``-``, ``*``, ``/``, ``%``, and ``**``.  They all have the same meaning as they do in Python source code --- i.e., they are addition, subtraction, multiplication, integer division, remainder, and exponentiation, respectively.  The acceptable unary operators are ``!`` and ``#``.  They are the factorial and primorial, respectively.  To avoid triggering the shell's special characters, there are three aliases: ``x`` for ``*``, ``xx`` for ``**``, and ``p!`` for ``#``.  You may also enclose the RPN expression in quotes if this helps avoid interpretation problems with your shell.


What's New in v1.1
==================

Bugfixes:

 - In version 1.0.0, when neither ``gmpy`` nor ``gmpy2`` could be imported, ``legendre`` was not defined properly and errors were thrown.  This is fixed in version 1.1.

New features:

 - A new function ``factorint`` is added with the same argument structure as the ``primefac`` generator, minus the ``verbose`` option.  This collates ``primefac``'s output into a dict with the prime factors as the keys and their multiplicities as the data.  For example, ``factorint(5040)`` returns ``{2:4, 3:2, 5:1, 7:1}``.
