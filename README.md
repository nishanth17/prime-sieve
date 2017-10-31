# prime-sieve

Efficient implementations of the sieves of [Atkin](https://en.wikipedia.org/wiki/Sieve_of_Atkin) and [Eratosthenes](https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes) in C. The algorithm used is based on the size of the input since the sieve of Erastosthenes is faster for N < 10<sup>8</sup>. Also contains a relatively optimized implementation of a segmented version of the sieve of Eratosthenes.

# Usage
This program requires AVX enabled and OpenMP. Usage is of the form

    ./prime_sieve [option] [args]

where `option` is either or `0` or `1`. The `0` option lists all the primes below a single argument using either the sieve of Eratosthenes or Atkin so

    ./prime_sieve 0 100

lists the primes under 100. The `1` option lists all the primes between 2 arguments (both inclusive) using the segmented sieve of Eratosthenes so

    ./prime_sieve 1 100 200

lists the primes between 100 and 200.


# Benchmarks
Benchmarks were computed using Macbook Pro with a core i7-7820HQ (2.9 GHz) processor and 16GB of RAM.


| Benchmarks (in secs)   | 10<sup>8</sup> | 10<sup>9</sup> | 3.3 x 10<sup>9</sup> | 10<sup>10</sup> | 2 x 10<sup>10</sup> | 4 x 10<sup>10</sup> |
|:--------------------   |---------------:|---------------:|---------------------:|----------------:|--------------------:|--------------------:|
| Atkin                  | 0.14           | 1.29           | 4.65                 | Overflow        | Overflow            | Overflow            |
| Segmented              | 0.18           | 1.77           | 6.09                 | 18.22           | 39.72               | 86.26               |
| Eratosthenes           | 0.32           | 4.66           | 15.76                | 47.89           | 100.03              | 228.13              |


# References
A.O.L Atkin, D.J.Bernstein; [Prime Sieves using Binary Quadratic Forms](http://www.ams.org/journals/mcom/2004-73-246/S0025-5718-03-01501-1/S0025-5718-03-01501-1.pdf); *Mathematics of Computation*, 73-246: 1023-30
