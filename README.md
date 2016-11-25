primefac fork
===============

This is fork of [primefac](https://pypi.python.org/pypi/primefac) Module.

* use a **fast** function, in `modinv`
* use [ecpy](https://github.com/elliptic-shiho/ecpy/) for faster ECM
* Implement fermat factorization

## License
This software released under the MIT License. 

---

Below is the original documentation for the primefac module. 

---

## primefac version 1.1

This is a module and command-line utility for factoring integers.  As a module, we provide a primality test, several functions for extracting a non-trivial factor of an integer, and a generator that yields all of a number's prime factors (with multiplicity).  As a command-line utility, this project aims to replace GNU's ``factor`` command with a more versatile utility --- in particular, this utility can operate on arbitrarily large numbers, uses multiple cores in parallel, uses better algorithms, handles input in reverse Polish notation, and can be tweaked via command-line flags.


## What's New in v1.1

Bugfixes:

 - In version 1.0.0, when neither ``gmpy`` nor ``gmpy2`` could be imported, ``legendre`` was not defined properly and errors were thrown.  This is fixed in version 1.1.

New features:

 - A new function ``factorint`` is added with the same argument structure as the ``primefac`` generator, minus the ``verbose`` option.  This collates ``primefac``'s output into a dict with the prime factors as the keys and their multiplicities as the data.  For example, ``factorint(5040)`` returns ``{2:4, 3:2, 5:1, 7:1}``.

---
