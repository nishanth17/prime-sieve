# prime-sieve

Efficient implementations of the sieves of [Atkin](https://en.wikipedia.org/wiki/Sieve_of_Atkin) and [Eratosthenes](https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes) in C. The algorithm used is based on the size of the input since the sieve of Erastosthenes is faster for N < 10<sup>7</sup>.

# Usage

    ./prime_sieve 100
    
lists the primes under 100.

# Reference

A.O.L Atkin, D.J.Bernstein; [Prime Sieves using Binary Quadratic Forms](http://www.ams.org/journals/mcom/2004-73-246/S0025-5718-03-01501-1/S0025-5718-03-01501-1.pdf); *Mathematics of Computation*, 73-246: 1023-30
