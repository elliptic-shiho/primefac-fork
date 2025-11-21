primefac fork 
===============

This is fork of [primefac](https://pypi.python.org/pypi/primefac) Module. 

**Caution: This software is forked from primefac v1.1 and has not kept up with the latest version of primefac.**

* use a **fast** function, in `modinv`
* Implement fermat factorization
* Use [factordb module](https://github.com/ryosan-470/factordb-pycli) for large number

## Requirements
* Python Python 3.x
  - To faster factorization, primefac-fork uses *gmpy2* so you need to install it (*gmpy* will install in below installation process). 
  - To factor with collective intelligence, primefac-fork uses *factordb-pycli* module. also this package will install in below. 
    + You can execute primefac-fork without factordb as follows: `python -m primefac -m=prb,p-1,p+1,ecm,mpqs,fermat YOUR_NUMBER_HERE`

## License
This software released under the MIT License. 

## Installation

To install directly from this repository:

``pip install git+https://github.com/elliptic-shiho/primefac-fork@master``

