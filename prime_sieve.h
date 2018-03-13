#ifndef _primes_h
#define _primes_h

long num_primes_below(long n);

long sieve_of_eratosthenes(long n, long* primes);

long sieve_of_atkin(long n, long* primes);

long prime_sieve(long n, long* primes);

long segmented_sieve(long lo, long hi, long* primes);

#endif

