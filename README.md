# prime-sieve

Efficient implementations of the sieves of [Atkin](https://en.wikipedia.org/wiki/Sieve_of_Atkin) and [Eratosthenes](https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes) in C. The algorithm used is based on the size of the input since the sieve of Erastosthenes is faster for N < 10<sup>8</sup>. Also contains a relatively optimized implementation of a segmented version of the sieve of Eratosthenes.

# Usage
Usage is of the form

    ./prime_sieve [option] [args]

where <i>option</i> is either or 0 or 1. 
lists the primes under 100.

# Reference

A.O.L Atkin, D.J.Bernstein; [Prime Sieves using Binary Quadratic Forms](http://www.ams.org/journals/mcom/2004-73-246/S0025-5718-03-01501-1/S0025-5718-03-01501-1.pdf); *Mathematics of Computation*, 73-246: 1023-30
