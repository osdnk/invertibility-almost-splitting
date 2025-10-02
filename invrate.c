// invrate.c
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/nmod_poly.h>

static uint64_t rdtsc_rng_state = 88172645463325252ull;

static inline uint64_t xorshift64(void){
    uint64_t x = rdtsc_rng_state;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    rdtsc_rng_state = x;
    return x;
}
static inline uint64_t rand_bounded(uint64_t bound_inclusive){
    // unbiased 0..bound_inclusive
    uint64_t m = bound_inclusive + 1;
    uint64_t lim = UINT64_MAX - (UINT64_MAX % m);
    uint64_t r;
    do { r = xorshift64(); } while (r > lim);
    return r % m;
}
static inline int64_t sample_in_range(int64_t B){
    // uniform in [-B, B]
    return (int64_t)rand_bounded((uint64_t)(2*B)) - B;
}

static int is_prime_ui(uint64_t q){
    fmpz_t n; fmpz_init(n); fmpz_set_ui(n, q);
    int r = fmpz_is_prime(n);   // correct FLINT function
    fmpz_clear(n);
    return r;
}

static int split_ok(uint64_t q, int f, int mode_afs){
    // f = 32. FS: q ≡ 1 (mod 2f). AFS: q ≡ 1 (mod f) but not mod 2f.
    uint64_t m = q - 1;
    if (mode_afs){
        if (m % (uint64_t)f != 0) return 0;
        if (m % (uint64_t)(2*f) == 0) return 0;
        return 1;
    }else{
        return (m % (uint64_t)(2*f) == 0);
    }
}

static uint64_t find_q_near_pow2(int k, int f, int mode_afs, uint64_t max_tries){
    if (k >= 63){ fprintf(stderr, "k too large for 64-bit.\n"); return 0; }
    uint64_t q = (uint64_t)1 << k;
    for (uint64_t t = 1; t <= max_tries; t++){
        q--;
        if (!is_prime_ui(q)) continue;
        if (split_ok(q, f, mode_afs)) return q;
    }
    return 0;
}

int main(int argc, char** argv){
    // Defaults
    const int deg = 16;       // phi(32)=32
    int f = deg * 2;               // Phi_32(x) = x^16 + 1
    int mode_afs = 1;         // 1 = AFS, 0 = FS
    int k = 15;               // search prime near 2^k
    uint64_t N = (uint64_t)pow(2, 30);     // trials per B
    int max_powB = 20;         // B in {2^0 .. 2^max_powB}
    uint64_t max_tries = 1000000;
    int use_inverse = 0;      // 0 = gcd check, 1 = inverse mod check

    // CLI: invrate [k] [N] [max_powB] [AFS=1|0] [seed]
    if (argc >= 2) k = atoi(argv[1]);
    if (argc >= 3) N = (uint64_t)pow(2, strtoull(argv[2], NULL, 10));
    if (argc >= 4) max_powB = atoi(argv[3]);
    if (argc >= 5) mode_afs = atoi(argv[4]);
    if (argc >= 6) rdtsc_rng_state = strtoull(argv[5], NULL, 10);

    uint64_t q = find_q_near_pow2(k, deg, mode_afs, max_tries);
    if (!q){
        fprintf(stderr, "No prime found near 2^%d under constraints.\n", k);
        return 1;
    }

    double expected = mode_afs
        ? 1.0 - pow(1.0 - 1.0 / ((double)q * (double)q), deg/2)
        : 1.0 - pow(1.0 - 1.0 / (double)q, deg);

    // Prepare phi(x) = x^16 + 1 mod q
    nmod_poly_t phi, c, g, inv;
    nmod_poly_init(phi, q);
    nmod_poly_init(c, q);
    nmod_poly_init(g, q);
    nmod_poly_init(inv, q);

    nmod_poly_zero(phi);
    nmod_poly_set_coeff_ui(phi, 0, 1);
    nmod_poly_set_coeff_ui(phi, deg, 1);

    // Pre-size c once to degree 15
    nmod_poly_fit_length(c, deg);
    _nmod_poly_set_length(c, deg);   // tell FLINT c has 16 coeff slots [0..15]

    printf("# k=%d, q=%" PRIu64 ", f=%d, split=%s, N=%" PRIu64 "\n",
           k, q, f, mode_afs ? "AFS" : "FS", N);
    printf("# Columns: B, invertible, N, noninvert_rate, expected_rate\n");

    for (int p = 0; p <= max_powB; p++){
        int64_t B = (int64_t)1 << p;
        uint64_t invertible = 0;

        for (uint64_t it = 0; it < N; it++){
            // Fill coefficients for c(x) with (a-b) ∈ [-2B, 2B], reduced mod q
            for (int i = 0; i < deg; i++){
                int64_t ai = sample_in_range(B);
                int64_t bi = sample_in_range(B);
                // printf(" // a[%d]=%" PRId64 ", b[%d]=%" PRId64 "\n", i, ai, i, bi);
                int64_t diff = ai - bi; // in [-2B..2B]
                uint64_t u;
                if (diff >= 0) u = (uint64_t)diff % q;
                else{
                    uint64_t t = (uint64_t)(-diff) % q;
                    u = t ? (q - t) : 0;
                }
                nmod_poly_set_coeff_ui(c, i, u);
            }
            // c length already set; no need to normalise each loop

            if (!use_inverse){
                nmod_poly_gcd(g, c, phi);
                if (nmod_poly_degree(g) == 0) invertible++;
            }else{
                if (nmod_poly_invmod(inv, c, phi)) invertible++;
            }
        }

        double rate = 1.0 - ((double)invertible / (double)N);
        printf("%" PRId64 ", %" PRIu64 ", %" PRIu64 ", %.10e, %.10e\n",
               B, invertible, N, rate, expected);
        fflush(stdout);
    }

    nmod_poly_clear(phi);
    nmod_poly_clear(c);
    nmod_poly_clear(g);
    nmod_poly_clear(inv);
    return 0;
}
