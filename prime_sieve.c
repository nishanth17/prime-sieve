#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <string.h> 
#include "prime_sieve.h"

/* COMPILER FLAGS: gcc -O3 -fopenmp -march=native -funroll-all-loops -Wa,-q prime_sieve.c -o prime_sieve
 
 NOTE: -> The -Wa,-q flag links gcc to the clang assembler since Apple's native llvm-gcc assembler doesn't
          support OpenMP instructions.
    -> The sieve of Atkin doesn't work for integers past MAX (~4 * 10^9 i.e. the limit of 32-bit integers).
 */
 
#define ERAT_MAX 79000000l
#define ATKIN_MAX 25000000000l

#define ERAT_SMALL_SEG_SIZE 65536l
#define ERAT_LARGE_SEG_SIZE 4194304l 
#define ATKIN_SEG_SIZE 1048576l

#define c2int(lo, k) (lo + (k << 1))
#define set_sieve(sieve, k) (sieve[k >> 3] |= 1 << (k & 7))
#define check_sieve(sieve, d) (((sieve[d >> 4] >> ((d >> 1) & 7)) & 1) == 0)

long DFG1[128][3] = { { 1, 0, 1 }, { 1, 0, 11 }, { 1, 0, 19 },
    { 1, 0, 29 }, { 1, 2, 15 }, { 1, 3, 5 }, { 1, 3, 25 }, { 1, 5, 9 },
    { 1, 5, 21 }, { 1, 7, 15 }, { 1, 8, 15 }, { 1, 10, 9 },
    { 1, 10, 21 }, { 1, 12, 5 }, { 1, 12, 25 }, { 1, 13, 15 },
    { 13, 1, 3 }, { 13, 1, 27 }, { 13, 4, 3 }, { 13, 4, 27 },
    { 13, 6, 7 }, { 13, 6, 13 }, { 13, 6, 17 }, { 13, 6, 23 },
    { 13, 9, 7 }, { 13, 9, 13 }, { 13, 9, 17 }, { 13, 9, 23 },
    { 13, 11, 3 }, { 13, 11, 27 }, { 13, 14, 3 }, { 13, 14, 27 },
    { 17, 2, 1 }, { 17, 2, 11 }, { 17, 2, 19 }, { 17, 2, 29 },
    { 17, 7, 1 }, { 17, 7, 11 }, { 17, 7, 19 }, { 17, 7, 29 },
    { 17, 8, 1 }, { 17, 8, 11 }, { 17, 8, 19 }, { 17, 8, 29 },
    { 17, 13, 1 }, { 17, 13, 11 }, { 17, 13, 19 }, { 17, 13, 29 },
    { 29, 1, 5 }, { 29, 1, 25 }, { 29, 4, 5 }, { 29, 4, 25 },
    { 29, 5, 7 }, { 29, 5, 13 }, { 29, 5, 17 }, { 29, 5, 23 },
    { 29, 10, 7 }, { 29, 10, 13 }, { 29, 10, 17 }, { 29, 10, 23 },
    { 29, 11, 5 }, { 29, 11, 25 }, { 29, 14, 5 }, { 29, 14, 25 },
    { 37, 2, 9 }, { 37, 2, 21 }, { 37, 3, 1 }, { 37, 3, 11 },
    { 37, 3, 19 }, { 37, 3, 29 }, { 37, 7, 9 }, { 37, 7, 21 },
    { 37, 8, 9 }, { 37, 8, 21 }, { 37, 12, 1 }, { 37, 12, 11 },
    { 37, 12, 19 }, { 37, 12, 29 }, { 37, 13, 9 }, { 37, 13, 21 },
    { 41, 2, 5 }, { 41, 2, 25 }, { 41, 5, 1 }, { 41, 5, 11 },
    { 41, 5, 19 }, { 41, 5, 29 }, { 41, 7, 5 }, { 41, 7, 25 },
    { 41, 8, 5 }, { 41, 8, 25 }, { 41, 10, 1 }, { 41, 10, 11 },
    { 41, 10, 19 }, { 41, 10, 29 }, { 41, 13, 5 }, { 41, 13, 25 },
    { 49, 0, 7 }, { 49, 0, 13 }, { 49, 0, 17 }, { 49, 0, 23 },
    { 49, 1, 15 }, { 49, 4, 15 }, { 49, 5, 3 }, { 49, 5, 27 },
    { 49, 6, 5 }, { 49, 6, 25 }, { 49, 9, 5 }, { 49, 9, 25 },
    { 49, 10, 3 }, { 49, 10, 27 }, { 49, 11, 15 }, { 49, 14, 15 },
    { 53, 1, 7 }, { 53, 1, 13 }, { 53, 1, 17 }, { 53, 1, 23 },
    { 53, 4, 7 }, { 53, 4, 13 }, { 53, 4, 17 }, { 53, 4, 23 },
    { 53, 11, 7 }, { 53, 11, 13 }, { 53, 11, 17 }, { 53, 11, 23 },
    { 53, 14, 7 }, { 53, 14, 13 }, { 53, 14, 17 }, { 53, 14, 23 } };

long DFG2[48][3] = { { 7, 1, 2 }, { 7, 1, 8 }, { 7, 1, 22 },
    { 7, 1, 28 }, { 7, 3, 10 }, { 7, 3, 20 }, { 7, 7, 10 },
    { 7, 7, 20 }, { 7, 9, 2 }, { 7, 9, 8 }, { 7, 9, 22 }, { 7, 9, 28 },
    { 19, 1, 4 }, { 19, 1, 14 }, { 19, 1, 16 }, { 19, 1, 26 },
    { 19, 5, 2 }, { 19, 5, 8 }, { 19, 5, 22 }, { 19, 5, 28 },
    { 19, 9, 4 }, { 19, 9, 14 }, { 19, 9, 16 }, { 19, 9, 26 },
    { 31, 3, 2 }, { 31, 3, 8 }, { 31, 3, 22 }, { 31, 3, 28 },
    { 31, 5, 4 }, { 31, 5, 14 }, { 31, 5, 16 }, { 31, 5, 26 },
    { 31, 7, 2 }, { 31, 7, 8 }, { 31, 7, 22 }, { 31, 7, 28 },
    { 43, 1, 10 }, { 43, 1, 20 }, { 43, 3, 4 }, { 43, 3, 14 },
    { 43, 3, 16 }, { 43, 3, 26 }, { 43, 7, 4 }, { 43, 7, 14 },
    { 43, 7, 16 }, { 43, 7, 26 }, { 43, 9, 10 }, { 43, 9, 20 } };

long DFG3[96][3] = { { 11, 0, 7 }, { 11, 0, 13 }, { 11, 0, 17 },
    { 11, 0, 23 }, { 11, 2, 1 }, { 11, 2, 11 }, { 11, 2, 19 },
    { 11, 2, 29 }, { 11, 3, 4 }, { 11, 3, 14 }, { 11, 3, 16 },
    { 11, 3, 26 }, { 11, 5, 2 }, { 11, 5, 8 }, { 11, 5, 22 },
    { 11, 5, 28 }, { 11, 7, 4 }, { 11, 7, 14 }, { 11, 7, 16 },
    { 11, 7, 26 }, { 11, 8, 1 }, { 11, 8, 11 }, { 11, 8, 19 },
    { 11, 8, 29 }, { 23, 1, 10 }, { 23, 1, 20 }, { 23, 2, 7 },
    { 23, 2, 13 }, { 23, 2, 17 }, { 23, 2, 23 }, { 23, 3, 2 },
    { 23, 3, 8 }, { 23, 3, 22 }, { 23, 3, 28 }, { 23, 4, 5 },
    { 23, 4, 25 }, { 23, 6, 5 }, { 23, 6, 25 }, { 23, 7, 2 },
    { 23, 7, 8 }, { 23, 7, 22 }, { 23, 7, 28 }, { 23, 8, 7 },
    { 23, 8, 13 }, { 23, 8, 17 }, { 23, 8, 23 }, { 23, 9, 10 },
    { 23, 9, 20 }, { 47, 1, 4 }, { 47, 1, 14 }, { 47, 1, 16 },
    { 47, 1, 26 }, { 47, 2, 5 }, { 47, 2, 25 }, { 47, 3, 10 },
    { 47, 3, 20 }, { 47, 4, 1 }, { 47, 4, 11 }, { 47, 4, 19 },
    { 47, 4, 29 }, { 47, 6, 1 }, { 47, 6, 11 }, { 47, 6, 19 },
    { 47, 6, 29 }, { 47, 7, 10 }, { 47, 7, 20 }, { 47, 8, 5 },
    { 47, 8, 25 }, { 47, 9, 4 }, { 47, 9, 14 }, { 47, 9, 16 },
    { 47, 9, 26 }, { 59, 0, 1 }, { 59, 0, 11 }, { 59, 0, 19 },
    { 59, 0, 29 }, { 59, 1, 2 }, { 59, 1, 8 }, { 59, 1, 22 },
    { 59, 1, 28 }, { 59, 4, 7 }, { 59, 4, 13 }, { 59, 4, 17 },
    { 59, 4, 23 }, { 59, 5, 4 }, { 59, 5, 14 }, { 59, 5, 16 },
    { 59, 5, 26 }, { 59, 6, 7 }, { 59, 6, 13 }, { 59, 6, 17 },
    { 59, 6, 23 }, { 59, 9, 2 }, { 59, 9, 8 }, { 59, 9, 22 },
    { 59, 9, 28 } };

long dAll[16] = { 1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59 };
long U60[17] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59 };

/* Returns a rough estimate of the number of primes below n. Based on the
 Taylor expansion of the first 2 terms of the logarithmic integral */
long num_primes_below(long n) {
    long u = n + 32;
    double lu = log(u);
    long num = (long) (u/lu + u/(lu*lu) * 1.5f);
    return num;
}

/* Returns a rough estimate of the number of primes below between lo and hi */
long num_primes_between(long lo, long hi) {
    return num_primes_below(hi) - num_primes_below(lo);
}

/* The extended GCD algorithm. Returns y such that ax + by = gcd(a, b) */
long extended_gcd(long a, long b) {
    long r = 0, s = 1;
    
    while (b != 0) {
        long c = a / b, d;
        d = a;
        a = b;
        b = d % b;
        d = r;
        r = s;
        s = d - c * s;
    }
    
    return r;
}

/* Binary search to find the index of the least element greater than 'x' in 'arr' */
long binary_search(long x, long* arr, long len) {
    if (x < arr[0]) {
        return 0;
    } else if (x > arr[len - 1]) {
        return len;
    }
    
    long l = 0, r = len - 1;
    while (l <= r) {
        long m = l + ((r - l) >> 1);
        if (arr[m] == x) {
            return m;
        } else if (arr[m] < x) {
            l = m + 1;
        } else {
            r = m - 1;
        }
    }
    
    return l;
}

/* Sieve of Eratosthenes w/ wheel factorization.
 * Note that this function returns the number of primes below 'n' and
 * populates the 'primes' array. */
long sieve_of_eratosthenes(long n, long* primes) {
    omp_set_num_threads(omp_get_num_procs());
    long nn = (((n >> 1) - 1) >> 5) + 1;
    long* sieve = (long*) calloc(nn, sizeof(long));
    long i, p, m, m2, ind, mm;
    
    #pragma omp parallel for
    for (i = 0; i < nn; i++)
        sieve[i] = 0xffffffff;
    
    sieve[nn - 1] = (1l << (((n >> 1) - 1) & 31)) - 1;
    long sq = ((long) sqrt(n) - 3) >> 1;
    for (p = 0; p <= sq; p++) {
        if ((sieve[p >> 5] & (1 << (p & 31))) != 0) {
            m = (p << 1) + 3; m2 = m << 1;
            
            #pragma omp parallel for private(ind)
            for (mm = m * m; mm <= n; mm += m2) {
                ind = (mm - 3) >> 1;
                sieve[ind >> 5] &= ~(1l << (ind & 31));
            }
        }
    }
    
    primes[0] = 2;
    long pos = 1;
    for (i = 0; i < nn; i++) {
        for (p = 0; p < 32; p++) {
            if ((sieve[i] & (1 << p)) != 0l)
                primes[pos++] = (((i << 5) | p) << 1) + 3;
        }
    }

    free(sieve);
    return pos;
}

/* Algorithm 3.1 : d < 60, f <= 15, g <= 30 such that d is congruent to 4f^2 + g^2(mod 60);
 * return all triples (x, y, k) with x > y > 0 and L <= k <= L + B such that
 * 4x^2 + y^2 = 60k + d */
void enum1(long f, long g, long d, long L, long B, long** segs) {
    long x = f, y0 = g;
    long temp = L + B;
    long k0 = ((f << 2) * f + g * g - d) / 60;
    
    while (k0 < temp) {
        k0 += (x << 1) + 15;
        x += 15;
    }
    
    while (1) {
        x -= 15;
        k0 -= (x << 1) + 15;
        if (x <= 0)
            break;
        
        while (k0 < L) {
            k0 += y0 + 15;
            y0 += 30;
        }
        
        long k = k0, y = y0;
        while (k < temp) {
            segs[d][(k - L) >> 5] ^= 1 << ((k - L) & 31);
            k += y + 15;
            y += 30;
        }
    }
}

/* Algorithm 3.2 : d < 60, f <= 10, g <= 30 such that d is congruent to 3f^2 + g^2(mod 60);
 * return all triples (x, y, k) with x > y > 0 and L <= k <= L + B such that
 * 3x^2 + y^2 = 60k + d */
void enum2(long f, long g, long d, long L, long B, long** segs) {
    long x = f, y0 = g;
    long temp = L + B;
    long k0 = (3 * f * f + g * g - d) / 60;
    
    while (k0 < temp) {
        k0 += x + 5;
        x += 10;
    }
    
    while (1) {
        x -= 10;
        k0 -= x + 5;
        if (x <= 0)
            break;
        
        while (k0 < L) {
            k0 += y0 + 15;
            y0 += 30;
        }
        
        long k = k0, y = y0;
        while (k < temp) {
            segs[d][(k - L) >> 5] ^= 1 << ((k - L) & 31);
            k += y + 15;
            y += 30;
        }
    }
}

/* Algorithm 3.3 : d < 60, f <= 10, g <= 30 such that d is congruent to 3f^2 - g^2(mod 60);
 * return all triples (x, y, k) with x > y > 0 and L <= k <= L + B such that
 * 3x^2 - y^2 = 60k + d */
void enum3(long f, long g, long d, long L, long B, long** segs) {
    long x = f, y0 = g;
    long temp = L + B;
    long k0 = (3 * f * f - g * g - d) / 60;
    
    while (1) {
        while (k0 >= temp) {
            if (x <= y0)
                return;
            k0 -= y0 + 15;
            y0 += 30;
        }
        
        long k = k0, y = y0;
        while (k >= L && y < x) {
            segs[d][(k - L) >> 5] ^= 1 << ((k - L) & 31);
            k -= y + 15;
            y += 30;
        }
        
        k0 += x + 5;
        x += 10;
    }
}

/* An implementation of a segmented version sieve of the Atkin in the paper referenced above.
 * This is far more efficient than Wikipedia's implementation of the same. Note that this
 function returns the number of primes below 'n' and populates the 'primes' array. */
long sieve_of_atkin(long n, long* primes) {
    omp_set_num_threads(omp_get_num_procs());
    long pos = 0, i, len;
    for (i = 0; i < 17; i++) {
        if (n < U60[i])
            return pos;
        primes[pos++] = U60[i];
    }
    
    long B = ATKIN_SEG_SIZE;
    long* base_primes = (long*) calloc(num_primes_below(sqrt(n)), sizeof(long));
    len = sieve_of_eratosthenes((long) sqrt(n), base_primes);
    long** segs = (long**) malloc(sizeof(long*) * 60);
    long p, p2, b, x, d, limit, j, lim1, base;

    // Precompute modular inverses
    long* px = (long*) malloc(len * sizeof(long));
    #pragma omp parallel for private(x, p2)
    for (i = 0; i < len; i++) {
        p2 = base_primes[i] * base_primes[i];
        x = -extended_gcd(p2, 60);
        if (x < 0)
            x += p2;
        px[i] = x;
    }
  
    long lim = n / 60l, L;
    long size = (B >> 5) + 1;
    for (i = 0; i < 16; i++)
        segs[dAll[i]] = (long*) calloc(size, sizeof(long));
    
    for (L = 1; L <= lim; L += B) {
        limit = 60 * (L + B);
        
        #pragma omp parallel for
        for (i = 0; i < 16; i++) {
            memset(segs[dAll[i]], 0, size * sizeof(long));
        }
        
        // Sieve off squarefree numbers
        #pragma omp parallel for
        for (i = 0; i < 128; i++)
            enum1(DFG1[i][1], DFG1[i][2], DFG1[i][0], L, B, segs);
        for (i = 0; i < 48; i++)
            enum2(DFG2[i][1], DFG2[i][2], DFG2[i][0], L, B, segs);
        for (i = 0; i < 96; i++)
            enum3(DFG3[i][1], DFG3[i][2], DFG3[i][0], L, B, segs);
        
        // Sieve off squareful numbers in parallel
        #pragma omp parallel for private(j, d, x, p, p2, b, i)
        for (i = 0; i < len; i++) {
            p = base_primes[i];
            p2 = p * p;
            if (p > 6 && p2 < limit) {
                b = px[i];
                for (j = 0; j < 16; j++) {
                    d = dAll[j];
                    x = (((b * L) % p2) * 60) + (b * d);
                    if (x > p2)
                        x %= p2;
                    for (x = x; x < B; x += p2) {
                        segs[d][x >> 5] &= ~(1 << (x & 31));
                    }
                }
            }
        }
        
        // Enumerate primes
        lim1 = (B >> 5) + 1;
        for (j = 0; j < lim1; j++) {
            for (x = 0; x < 32; x++) {
                base = 60 * (L + x + (j << 5));
                for (i = 0; i < 16; i++) {
                    d = dAll[i];
                    if (base + d > n)
                        goto exit;
                    if ((int)(segs[d][j] << 31 - x) < 0)
                        primes[pos++] = base + d;
                }
            }
        }
        
    exit:
        continue;
    }
    
    free(segs); free(base_primes); free(px);
    return pos;
}

/* A relativly simple segmented sieve of Eratosthenes to compute the primes between 'lo' and 'hi'
 * (both inclusive) using a wheel modulo 6. Returns the number of primes between them and
 * populates the 'primes' array. This can probably be optimized further with a wheel modulo
 * some primorial. */
long segmented_sieve(long lo, long hi, long* primes) {
    if (hi < lo) {
        return 0;
    }
    
    long max_prime = (long) sqrt(hi);
    long pos = 0, k = 0, lo_1 = 0, hi_1 = 0, i = 0, p = 0, n = 0, d = 0, end = 0;
    long* base_primes = (long*) calloc(num_primes_below(max_prime), sizeof(long));
    long num_base_primes = prime_sieve(max_prime, base_primes);
    
    // Include primes below sqrt(hi) if necessary
    if (lo < max_prime) {
        long bin_pos = binary_search(lo, base_primes, num_base_primes);
        for (k = bin_pos; k < num_base_primes; k++) {
            primes[pos++] = base_primes[k];
        }
        lo = max_prime;
    }
    
    // Compute segment and interval size
    long delta = (hi - lo > ERAT_LARGE_SEG_SIZE) ? ERAT_LARGE_SEG_SIZE : ERAT_SMALL_SEG_SIZE;
    long l = (delta >> 4) + 1, int_size = l << 3;
    byte* sieve = calloc(l, sizeof(byte));
    
    for (lo_1 = lo, hi_1 = lo + delta; lo_1 <= hi; lo_1 = hi_1 + 1, hi_1 += delta + 1) {
        // Re-zero sieve bytes when necessary
        if (lo_1 != lo) {
            memset(sieve, 0, l);
        }
        
        lo_1 = ((lo_1 & 1) == 0) ? lo_1 + 1 : lo_1;
        end = (hi_1 < hi) ? hi_1 : hi;
        // Mark and sieve primes
        for (i = 2; i < num_base_primes; i++) {
            p = base_primes[i];
            k = (p - (lo_1 % p)) % p;
            k = ((k & 1) == 1) ? (k + p) >> 1 : k >> 1;
            
            // Find closest multiple of p that = 1 (mod 6)
            switch ((c2int(lo_1, k) / p) % 3) {
                case 0:
                    k += p;
                case 2: {
                    if (k < int_size && c2int(lo_1, k) < end)
                        set_sieve(sieve, k);
                    k += p;
                    break;
                }
            }
            
            long twice_p = p << 1, step = twice_p + p;
            for (; k < int_size && c2int(lo_1, k) < end; k += step) {
                // Test 1 (mod 6)
                set_sieve(sieve, k);
                // Test (5 mod 6)
                if (k + twice_p < int_size && c2int(lo_1, k + twice_p) < end)
                    set_sieve(sieve, k + twice_p);
            }
        }
        
        // Find primes - start by finding closest number to lo_1 that = 1 (mod 6)
        long mod6 = lo_1 % 6;
        n = lo_1;
        if (mod6 <= 5 && mod6 >= 2) {
            d = 5 - mod6;
            if (lo_1 + d <= end && check_sieve(sieve, d)) {
                primes[pos++] = lo_1 + d;
            }
            n += d + 2;
        } else if (mod6 == 0) {
            n += 1;
        }
        
        // Then check everything that = 1 or 5 (mod 6)
        for (; n <= end; n += 6) {
            d = n - lo_1;
            if (check_sieve(sieve, d))
                primes[pos++] = n;
            if (n + 4 <= end && check_sieve(sieve, d+4))
                primes[pos++] = n + 4;
        }
    }
    
    free(base_primes); free(sieve); 
    return pos;
}

/* This function returns the number of primes below 'n' and populates the
 * 'primes' array. */
long prime_sieve(long n, long* primes) {
    if (n < ERAT_MAX)
        return sieve_of_eratosthenes(n, primes);
    else if (n <= ATKIN_MAX) {
        return sieve_of_atkin(n, primes);
    } else {
        return segmented_sieve(2, n, primes);
    }
}

int main(int argc, char** argv) {
    struct timeval start, end;
    int option = atoi(argv[1]);
    
    if (option == 0) {
        long N, i;
        long* primes;
        N = atol(argv[2]);
        primes = (long*) calloc(num_primes_below(N), sizeof(long));
    
        gettimeofday(&start, NULL);
        long t = sieve_of_atkin(N, primes);
        gettimeofday(&end, NULL);
    
        printf("\nPrimes below %lu:\n", N);
        for (i = 0; i < t; i++) {
           printf("%lu\n", primes[i]);
        }
        printf("\nNumber of primes below %lu: %lu\n", N, t);
    } else if (option == 1) {
        long lo, hi, i;
        long* primes;
        lo = atol(argv[2]);
        hi = atol(argv[3]);
        
        if (lo > hi) {
            printf("\nInvalid arguments\n");
            return 0;
        } else {
            long* primes = (long*) calloc(num_primes_between(lo, hi), sizeof(uint64_t));
            
            gettimeofday(&start, NULL);
            long t = segmented_sieve(lo, hi, primes);
            gettimeofday(&end, NULL);
            
            printf("\nPrimes primes between %lu and %lu:\n", lo, hi);
            for (i = 0; i < t; i++) {
               printf("%lu\n", primes[i]);
            }
            printf("Number of primes between %lu and %lu: %lu\n", lo, hi, t);
        }
    } else {
        printf("\nInvalid option\n");
        return 0;
    }
    
    double time = (double) ((end.tv_sec  - start.tv_sec) * 1000000u +
                            end.tv_usec - start.tv_usec) / 1.e6;
    printf("Time: %F seconds\n\n", time);
    return 0;
} 
