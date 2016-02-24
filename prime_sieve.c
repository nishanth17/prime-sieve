#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <immintrin.h>

/*  COMPILER FLAGS: gcc -O3 -fopenmp -msse4.2 -march=native -funroll-all-loops -Wa,-q prime_sieve.c -o prime_sieve
 
    NOTE: -> The -Wa,-q flag links gcc to the clang assembler since Apple's native llvm-gcc assembler doesn't
             support AVX instructions.
          -> The sieve of Atkin doesn't work for integers past MAX (~4 * 10^9 i.e. the limit of 32-bit integers)
 
    BENCHMARKS   |   10^8    |    10^9   |  3.35*10^9
    --------------------------------------------------
    Eratosthenes |   0.33s   |    3.9s   |    13.2s
    Atkin        |   0.14s   |    1.4s   |     4.9s
 
    Reference: "Prime Sieves using Binary Quadratic Forms", A.O.L Atkin and D.J.Bernstein, Mathematics of 
               Computation, 73-246: 1023-30, 
               http://www.ams.org/journals/mcom/2004-73-246/S0025-5718-03-01501-1/S0025-5718-03-01501-1.pdf
 */

#define MAX 3350000000l
#define THRESHOLD 79000000l

typedef unsigned long long uint64_t;

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
   Taylor expansion of the first 2 terms of the logarihtmic integral */
long num_primes_below(long n) {
    long u = n + 32;
    double lu = log(u);
    long num = (long) (u/lu + u/(lu*lu) * 1.5f);
    return num;
}

/* The extnded GCD algorithm. Returns a y such that ax + b = gcd(a, b) */
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

/* Sieve of Eratosthenes w/ wheel factorization.
 * Note that this function returns the number of primes below 'n' and
 * populates the 'primes' array. */
long sieve_of_eratosthenes(long n, long* primes) {
    omp_set_num_threads(omp_get_num_procs());
    long nn = (((n >> 1) - 1) >> 5) + 1;
    long* arr = (long*) calloc(nn, sizeof(long));
    long i, p, m, m2, ind, mm;

    #pragma omp parallel for
    for (i = 0; i < nn; i++)
        arr[i] = 0xffffffff;
    
    arr[nn - 1] = (1l << (((n >> 1) - 1) & 31)) - 1;
    long sq = ((long) sqrt(n) - 3) >> 1;
    
    for (p = 0; p <= sq; p++) {
        if ((arr[p >> 5] & (1 << (p & 31))) != 0) {
            m = (p << 1) + 3; m2 = m << 1;
            
            #pragma omp parallel for private(ind) 
            for (mm = m * m; mm <= n; mm += m2) {
                ind = (mm - 3) >> 1;
                arr[ind >> 5] &= ~(1l << (ind & 31));
            }
        }
    }
    
    primes[0] = 2;
    long pos = 1;
    for (i = 0; i < nn; i++) {
        for (p = 0; p < 32; p++) {
            if ((arr[i] & (1 << p)) != 0l)
                primes[pos++] = (((i << 5) | p) << 1) + 3;
        }
    }
    
    return pos;
}

/* Algorithm 4.1 : d < 60, f <= 15, g <= 30 such that d is congruent to 4f^2 + g^2(mod 60);
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

/* Algorithm 4.2 : d < 60, f <= 10, g <= 30 such that d is congruent to 3f^2 + g^2(mod 60);
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

/* Algorithm 4.3 : d < 60, f <= 10, g <= 30 such that d is congruent to 3f^2 - g^2(mod 60);
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

/* An implementation of the sieve of Atkin in the paper referenced above.  
 * This is far more efficient than Wikipedia's implementation of the same. 
 * Note that this function returns the number of primes below 'n' and
 * populates the 'primes' array. */
long sieve_of_atkin(long n, long* primes) {
    omp_set_num_threads(omp_get_num_procs());
    long pos = 0, i, len;
    for (i = 0; i < 17; i++) {
        if (n < U60[i])
            return pos;
        primes[pos++] = U60[i];
    }
    
    long B = 60 * ((long) sqrt(n));
    long* base_primes = (long*)calloc(num_primes_below(sqrt(n)), sizeof(long));
    len = sieve_of_eratosthenes((long) sqrt(n), base_primes);
    long** segs = (long**) malloc(sizeof(long*) * 60);
    
    long lim = n / 60l, L;
    long size = (B >> 5) + 1;
    long p, p2, b, x, d, limit, j, lim1, base, s2 = (size >> 2) << 2;;
    __m256i zero = _mm256_setzero_si256();
    
    for (i = 0; i < 16; i++)
        segs[dAll[i]] = (long*) calloc(size, sizeof(long));

    for (L = 1; L <= lim; L += B) {
        limit = 60 * (L + B);
        
        #pragma omp parallel for private(i, j)
        for (i = 0; i < 16; i++) {
            for (j = 0; j < s2; j += 4)
               _mm256_store_si256((__m256i*)&(segs[dAll[i]][j]), zero);
            for (j = s2; j < size; j++)
                segs[dAll[i]][j] = 0;
        }
        
        // Sieve off squarefree numbers
        #pragma omp parallel for private(i) 
        for (i = 0; i < 128; i++) 
            enum1(DFG1[i][1], DFG1[i][2], DFG1[i][0], L, B, segs);
        for (i = 0; i < 48; i++)
            enum2(DFG2[i][1], DFG2[i][2], DFG2[i][0], L, B, segs);
        for (i = 0; i < 96; i++)
            enum3(DFG3[i][1], DFG3[i][2], DFG3[i][0], L, B, segs);

        // Sieve off squareful numbers
        for (i = 0; i < len; i++) {
            p = base_primes[i];
            p2 = p * p;
            
            if (p2 > limit)
                break;
            if (p > 6) {
                b = -extended_gcd(p2, 60);
                if (b < 0)
                    b += p2;
                for (j = 0; j < 16; j++) {
                    d = dAll[j];
                    x = b * (60 * L + d);
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
    
    return pos;
}

/* Returns the number of primes below 'n' and populates the
 * 'primes' array. */
long prime_sieve(long n, long* primes) {
    if (n < THRESHOLD || n > MAX)
        return sieve_of_eratosthenes(n, primes);
    else
        return sieve_of_atkin(n, primes);
}


int main(int argc, char** argv) {
    struct timeval start, end;
    long N, i;
    long* primes;
    
    N = atoi(argv[1]);
    primes = (long*) calloc(num_primes_below(N), sizeof(long));
    
    gettimeofday(&start, NULL);
    long t = prime_sieve(N, primes);
    gettimeofday(&end, NULL);
    
    printf("\nPrimes below %lu:\n", N);
    for (i = 0; i < t; i++) {
        printf("%lu\n", primes[i]);
    }
    
    
    printf("\nNumber of primes below %lu: %lu\n", N, t);
    double time = (double) ((end.tv_sec  - start.tv_sec) * 1000000u +
                            end.tv_usec - start.tv_usec) / 1.e6;
    printf("Time: %F seconds\n\n", time);
    return 0;
} 