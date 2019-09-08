#! /usr/bin/env python

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'DESCRIPTION.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='primefac',
    version='1.1.1-2',
    description='Module and command-line utility for factoring integers into '
    'primes.  Formerly called pyfac.',
    long_description=long_description,
    url='https://pypi.python.org/pypi/primefac',
    author='lucasbrown.cit',
    author_email='lucasbrown.cit@gmail.com',
    license='MIT',
    keywords='number numbers integer integers factoring factorization primes '
    'prime numbers math mathematics pollard\'s rho pollard\'s p-1 williams\' '
    'p+1 elliptic curve method ecm ecf multiple polynomial quadratic sieve '
    'mpqs',
    packages=find_packages(),
    py_modules=['primefac', '_primefac'],
    install_requires=['six', 'gmpy2', 'factordb-pycli'],

    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Software Development :: Libraries',
        'Topic :: Utilities',
    ],
)
