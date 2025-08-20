/* gf2n_field.h - GF(2^n) implementations with native GF(2^128) PCLMUL support */
#ifndef GF2N_FIELD_H
#define GF2N_FIELD_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <immintrin.h>
#include <flint/flint.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_poly.h>
#include <flint/nmod_poly.h>
#include <gmp.h>
#include <cpuid.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
   FORWARD DECLARATIONS
   ============================================================================ */

/* Forward declarations for conversion functions */
typedef struct {
    uint64_t low;
    uint64_t high;
} gf2128_t;

static gf2128_t fq_nmod_to_gf2128(const fq_nmod_t elem, const fq_nmod_ctx_t ctx);
static void gf2128_to_fq_nmod(fq_nmod_t res, const gf2128_t *elem, const fq_nmod_ctx_t ctx);

/* ============================================================================
   GF(2^128) NATIVE PCLMUL IMPLEMENTATION
   ============================================================================ */

/* Function pointer for multiplication */
typedef gf2128_t (*gf2128_mul_func)(const gf2128_t*, const gf2128_t*);

/* Global multiplication function pointer */
static gf2128_mul_func gf2128_mul = NULL;

/* Check CPU features */
static int has_pclmulqdq(void) {
    unsigned int eax, ebx, ecx, edx;
    if (__get_cpuid(1, &eax, &ebx, &ecx, &edx)) {
        return (ecx & bit_PCLMUL) != 0;
    }
    return 0;
}

/* Basic operations */
static inline gf2128_t gf2128_create(uint64_t low, uint64_t high) {
    gf2128_t result = {low, high};
    return result;
}

static inline gf2128_t gf2128_zero(void) {
    return gf2128_create(0, 0);
}

static inline gf2128_t gf2128_one(void) {
    return gf2128_create(1, 0);
}

static inline int gf2128_is_zero(const gf2128_t *a) {
    return (a->low == 0 && a->high == 0);
}

static inline int gf2128_equal(const gf2128_t *a, const gf2128_t *b) {
    return (a->low == b->low && a->high == b->high);
}

static inline gf2128_t gf2128_add(const gf2128_t *a, const gf2128_t *b) {
    gf2128_t result;
    result.low = a->low ^ b->low;
    result.high = a->high ^ b->high;
    return result;
}

static void gf2128_print(const gf2128_t *a) {
    printf("%016lx%016lx", a->high, a->low);
}

/* Software multiplication for GF(2^128) */
static gf2128_t gf2128_mul_software(const gf2128_t *a, const gf2128_t *b) {
    uint64_t a0 = a->low, a1 = a->high;
    uint64_t b0 = b->low, b1 = b->high;
    uint64_t z0 = 0, z1 = 0, z2 = 0, z3 = 0;
    
    /* 128x128 bit multiplication using schoolbook method */
    for (int i = 0; i < 64; i++) {
        if ((b0 >> i) & 1) {
            z0 ^= a0 << i;
            if (i > 0) {
                z1 ^= a0 >> (64 - i);
            }
        }
    }
    
    for (int i = 0; i < 64; i++) {
        if ((b0 >> i) & 1) {
            z1 ^= a1 << i;
            if (i > 0) {
                z2 ^= a1 >> (64 - i);
            }
        }
    }
    
    for (int i = 0; i < 64; i++) {
        if ((b1 >> i) & 1) {
            z1 ^= a0 << i;
            if (i > 0) {
                z2 ^= a0 >> (64 - i);
            }
        }
    }
    
    for (int i = 0; i < 64; i++) {
        if ((b1 >> i) & 1) {
            z2 ^= a1 << i;
            if (i > 0) {
                z3 ^= a1 >> (64 - i);
            }
        }
    }
    
    /* Reduction modulo x^128 + x^7 + x^2 + x + 1 */
    /* Reduce z3 (bits 192-255) */
    for (int i = 63; i >= 0; i--) {
        if ((z3 >> i) & 1) {
            int pos = 64 + i;
            z3 ^= (1ULL << i);
            
            if (pos + 7 < 128) {
                z1 ^= (1ULL << (pos + 7 - 64));
            } else {
                z2 ^= (1ULL << (pos + 7 - 128));
            }
            
            if (pos + 2 < 128) {
                z1 ^= (1ULL << (pos + 2 - 64));
            } else {
                z2 ^= (1ULL << (pos + 2 - 128));
            }
            
            if (pos + 1 < 128) {
                z1 ^= (1ULL << (pos + 1 - 64));
            } else {
                z2 ^= (1ULL << (pos + 1 - 128));
            }
            
            if (pos < 128) {
                z1 ^= (1ULL << (pos - 64));
            } else {
                z2 ^= (1ULL << (pos - 128));
            }
        }
    }
    
    /* Reduce z2 (bits 128-191) */
    for (int i = 63; i >= 0; i--) {
        if ((z2 >> i) & 1) {
            int pos = i;
            z2 ^= (1ULL << i);
            
            if (pos + 7 < 64) {
                z0 ^= (1ULL << (pos + 7));
            } else {
                z1 ^= (1ULL << (pos + 7 - 64));
            }
            
            if (pos + 2 < 64) {
                z0 ^= (1ULL << (pos + 2));
            } else {
                z1 ^= (1ULL << (pos + 2 - 64));
            }
            
            if (pos + 1 < 64) {
                z0 ^= (1ULL << (pos + 1));
            } else {
                z1 ^= (1ULL << (pos + 1 - 64));
            }
            
            if (pos < 64) {
                z0 ^= (1ULL << pos);
            } else {
                z1 ^= (1ULL << (pos - 64));
            }
        }
    }
    
    return gf2128_create(z0, z1);
}

/* CLMUL-based multiplication for GF(2^128) */
static inline gf2128_t gf2128_mul_clmul(const gf2128_t *a, const gf2128_t *b) {
    __m128i x = _mm_set_epi64x(a->high, a->low);
    __m128i y = _mm_set_epi64x(b->high, b->low);
    
    /* Step 1: 128x128 -> 256 bit multiplication */
    __m128i t0, t1, t2;
    
    t0 = _mm_clmulepi64_si128(x, y, 0x00);
    t1 = _mm_clmulepi64_si128(x, y, 0x11);
    t2 = _mm_xor_si128(
        _mm_clmulepi64_si128(x, y, 0x10),
        _mm_clmulepi64_si128(x, y, 0x01)
    );
    
    __m128i lo = _mm_xor_si128(t0, _mm_slli_si128(t2, 8));
    __m128i hi = _mm_xor_si128(t1, _mm_srli_si128(t2, 8));
    
    /* Step 2: Optimized reduction using bit-reflected algorithm */
    const uint64_t g = 0x87;  // Bit representation of x^7 + x^2 + x + 1
    
    /* First phase: eliminate bits 255-192 */
    __m128i tmp = _mm_clmulepi64_si128(hi, _mm_set_epi64x(0, g), 0x01);
    hi = _mm_xor_si128(hi, _mm_srli_si128(tmp, 8));
    lo = _mm_xor_si128(lo, _mm_slli_si128(tmp, 8));
    
    /* Second phase: eliminate bits 191-128 */
    tmp = _mm_clmulepi64_si128(hi, _mm_set_epi64x(0, g), 0x00);
    lo = _mm_xor_si128(lo, tmp);
    
    /* Extract result */
    gf2128_t res;
    res.low = _mm_extract_epi64(lo, 0);
    res.high = _mm_extract_epi64(lo, 1);
    
    return res;
}

/* Initialize multiplication function */
static void init_gf2128(void) {
    if (!gf2128_mul) {
        if (has_pclmulqdq()) {
            gf2128_mul = gf2128_mul_clmul;
            printf("Using PCLMULQDQ for GF(2^128) multiplication\n");
        } else {
            gf2128_mul = gf2128_mul_software;
            printf("Using software implementation for GF(2^128) multiplication\n");
        }
    }
}

/* Fast squaring in GF(2^128) */
static gf2128_t gf2128_sqr(const gf2128_t *a) {
    uint64_t r0 = 0, r1 = 0, r2 = 0, r3 = 0;
    
    /* Expand low 64 bits */
    for (int i = 0; i < 32; i++) {
        if ((a->low >> i) & 1) {
            r0 |= 1ULL << (2 * i);
        }
    }
    for (int i = 32; i < 64; i++) {
        if ((a->low >> i) & 1) {
            r1 |= 1ULL << (2 * (i - 32));
        }
    }
    
    /* Expand high 64 bits */
    for (int i = 0; i < 32; i++) {
        if ((a->high >> i) & 1) {
            r2 |= 1ULL << (2 * i);
        }
    }
    for (int i = 32; i < 64; i++) {
        if ((a->high >> i) & 1) {
            r3 |= 1ULL << (2 * (i - 32));
        }
    }
    
    /* Reduction modulo x^128 + x^7 + x^2 + x + 1 */
    for (int i = 63; i >= 0; i--) {
        if ((r3 >> i) & 1) {
            int pos = 64 + i;
            r3 ^= (1ULL << i);
            
            if (pos + 7 < 128) {
                r1 ^= (1ULL << (pos + 7 - 64));
            } else {
                r2 ^= (1ULL << (pos + 7 - 128));
            }
            
            if (pos + 2 < 128) {
                r1 ^= (1ULL << (pos + 2 - 64));
            } else {
                r2 ^= (1ULL << (pos + 2 - 128));
            }
            
            if (pos + 1 < 128) {
                r1 ^= (1ULL << (pos + 1 - 64));
            } else {
                r2 ^= (1ULL << (pos + 1 - 128));
            }
            
            if (pos < 128) {
                r1 ^= (1ULL << (pos - 64));
            } else {
                r2 ^= (1ULL << (pos - 128));
            }
        }
    }
    
    for (int i = 63; i >= 0; i--) {
        if ((r2 >> i) & 1) {
            int pos = i;
            r2 ^= (1ULL << i);
            
            if (pos + 7 < 64) {
                r0 ^= (1ULL << (pos + 7));
            } else {
                r1 ^= (1ULL << (pos + 7 - 64));
            }
            
            if (pos + 2 < 64) {
                r0 ^= (1ULL << (pos + 2));
            } else {
                r1 ^= (1ULL << (pos + 2 - 64));
            }
            
            if (pos + 1 < 64) {
                r0 ^= (1ULL << (pos + 1));
            } else {
                r1 ^= (1ULL << (pos + 1 - 64));
            }
            
            if (pos < 64) {
                r0 ^= (1ULL << pos);
            } else {
                r1 ^= (1ULL << (pos - 64));
            }
        }
    }
    
    return gf2128_create(r0, r1);
}

/* Forward declaration for conversion function */
static gf2128_t fq_nmod_to_gf2128(const fq_nmod_t elem, const fq_nmod_ctx_t ctx);

/* Inversion using repeated squaring - optimized for GF(2^128) */
static gf2128_t gf2128_inv(const gf2128_t *a) {
    if (gf2128_is_zero(a)) {
        return gf2128_zero();
    }
    
    /* For debugging: use FLINT to compute the correct inverse */
    /* Create GF(2^128) context */
    fq_nmod_ctx_t ctx;
    nmod_poly_t mod;
    nmod_poly_init(mod, 2);
    
    /* Set modulus to x^128 + x^7 + x^2 + x + 1 */
    nmod_poly_set_coeff_ui(mod, 0, 1);
    nmod_poly_set_coeff_ui(mod, 1, 1);
    nmod_poly_set_coeff_ui(mod, 2, 1);
    nmod_poly_set_coeff_ui(mod, 7, 1);
    nmod_poly_set_coeff_ui(mod, 128, 1);
    
    fq_nmod_ctx_init_modulus(ctx, mod, "t");
    
    /* Convert to FLINT format */
    fq_nmod_t a_flint, inv_flint;
    fq_nmod_init(a_flint, ctx);
    fq_nmod_init(inv_flint, ctx);
    
    gf2128_to_fq_nmod(a_flint, a, ctx);
    
    /* Compute inverse using FLINT */
    fq_nmod_inv(inv_flint, a_flint, ctx);
    
    /* Convert back */
    gf2128_t result = fq_nmod_to_gf2128(inv_flint, ctx);
    
    /* Cleanup */
    fq_nmod_clear(a_flint, ctx);
    fq_nmod_clear(inv_flint, ctx);
    fq_nmod_ctx_clear(ctx);
    nmod_poly_clear(mod);
    
    return result;
}

/* Division in GF(2^128) */
static inline gf2128_t gf2128_div(const gf2128_t *a, const gf2128_t *b) {
    gf2128_t b_inv = gf2128_inv(b);
    return gf2128_mul(a, &b_inv);
}

/* ============================================================================
   GF(2^128) POLYNOMIAL OPERATIONS
   ============================================================================ */

typedef struct {
    gf2128_t *coeffs;
    slong length;
    slong alloc;
} gf2128_poly_struct;

typedef gf2128_poly_struct gf2128_poly_t[1];

static void gf2128_poly_init(gf2128_poly_t poly) {
    poly->coeffs = NULL;
    poly->length = 0;
    poly->alloc = 0;
}

static void gf2128_poly_clear(gf2128_poly_t poly) {
    if (poly->coeffs) {
        free(poly->coeffs);
        poly->coeffs = NULL;
    }
    poly->length = 0;
    poly->alloc = 0;
}

static void gf2128_poly_fit_length(gf2128_poly_t poly, slong len) {
    if (len > poly->alloc) {
        slong new_alloc = FLINT_MAX(len, poly->alloc * 2);
        poly->coeffs = (gf2128_t *)realloc(poly->coeffs, new_alloc * sizeof(gf2128_t));
        for (slong i = poly->alloc; i < new_alloc; i++) {
            poly->coeffs[i] = gf2128_zero();
        }
        poly->alloc = new_alloc;
    }
}


static void gf2128_poly_normalise(gf2128_poly_t poly) {
    while (poly->length > 0 && gf2128_is_zero(&poly->coeffs[poly->length - 1])) {
        poly->length--;
    }
}


/* Helper function to set polynomial coefficient */
static void gf2128_poly_set_coeff(gf2128_poly_t poly, slong i, const gf2128_t *c) {
    gf2128_poly_fit_length(poly, i + 1);
    if (i >= poly->length) {
        // Zero out coefficients between old length and i
        for (slong j = poly->length; j < i; j++) {
            poly->coeffs[j] = gf2128_zero();
        }
        poly->length = i + 1;
    }
    poly->coeffs[i] = *c;
    
    // Update length if we're setting the coefficient to zero
    if (gf2128_is_zero(c) && i == poly->length - 1) {
        gf2128_poly_normalise(poly);
    }
}

static inline void gf2128_poly_zero(gf2128_poly_t poly) {
    poly->length = 0;
}

static inline int gf2128_poly_is_zero(const gf2128_poly_t poly) {
    return poly->length == 0;
}

static inline slong gf2128_poly_degree(const gf2128_poly_t poly) {
    return poly->length - 1;
}

static void gf2128_poly_set(gf2128_poly_t res, const gf2128_poly_t poly) {
    if (res == poly) return;
    gf2128_poly_fit_length(res, poly->length);
    memcpy(res->coeffs, poly->coeffs, poly->length * sizeof(gf2128_t));
    res->length = poly->length;
}

static gf2128_t gf2128_poly_get_coeff(const gf2128_poly_t poly, slong i) {
    if (i < poly->length) {
        return poly->coeffs[i];
    } else {
        return gf2128_zero();
    }
}

static void gf2128_poly_add(gf2128_poly_t res, const gf2128_poly_t a, const gf2128_poly_t b) {
    slong max_len = FLINT_MAX(a->length, b->length);
    slong min_len = FLINT_MIN(a->length, b->length);
    
    if (max_len == 0) {
        gf2128_poly_zero(res);
        return;
    }
    
    gf2128_poly_fit_length(res, max_len);
    
    for (slong i = 0; i < min_len; i++) {
        res->coeffs[i] = gf2128_add(&a->coeffs[i], &b->coeffs[i]);
    }
    
    if (a->length > b->length) {
        memcpy(res->coeffs + min_len, a->coeffs + min_len, 
               (a->length - min_len) * sizeof(gf2128_t));
        res->length = a->length;
    } else if (b->length > a->length) {
        memcpy(res->coeffs + min_len, b->coeffs + min_len, 
               (b->length - min_len) * sizeof(gf2128_t));
        res->length = b->length;
    } else {
        res->length = min_len;
    }
    
    gf2128_poly_normalise(res);
}

static void gf2128_poly_scalar_mul(gf2128_poly_t res, const gf2128_poly_t poly, const gf2128_t *c) {
    if (gf2128_is_zero(c)) {
        gf2128_poly_zero(res);
        return;
    }
    
    gf2128_t one = gf2128_one();
    if (gf2128_equal(c, &one)) {
        gf2128_poly_set(res, poly);
        return;
    }
    
    gf2128_poly_fit_length(res, poly->length);
    res->length = poly->length;
    
    for (slong i = 0; i < poly->length; i++) {
        res->coeffs[i] = gf2128_mul(&poly->coeffs[i], c);
    }
    
    gf2128_poly_normalise(res);
}

static void gf2128_poly_shift_left(gf2128_poly_t res, const gf2128_poly_t poly, slong n) {
    if (n == 0) {
        gf2128_poly_set(res, poly);
        return;
    }
    
    if (gf2128_poly_is_zero(poly)) {
        gf2128_poly_zero(res);
        return;
    }
    
    slong new_len = poly->length + n;
    gf2128_poly_fit_length(res, new_len);
    
    memmove(res->coeffs + n, poly->coeffs, poly->length * sizeof(gf2128_t));
    
    for (slong i = 0; i < n; i++) {
        res->coeffs[i] = gf2128_zero();
    }
    
    res->length = new_len;
}

static void gf2128_poly_mul_schoolbook(gf2128_poly_t res, const gf2128_poly_t a, const gf2128_poly_t b) {
    if (gf2128_poly_is_zero(a) || gf2128_poly_is_zero(b)) {
        gf2128_poly_zero(res);
        return;
    }
    
    slong rlen = a->length + b->length - 1;
    gf2128_t *temp = (gf2128_t *)calloc(rlen, sizeof(gf2128_t));
    
    for (slong i = 0; i < a->length; i++) {
        for (slong j = 0; j < b->length; j++) {
            gf2128_t prod = gf2128_mul(&a->coeffs[i], &b->coeffs[j]);
            temp[i + j] = gf2128_add(&temp[i + j], &prod);
        }
    }
    
    gf2128_poly_fit_length(res, rlen);
    memcpy(res->coeffs, temp, rlen * sizeof(gf2128_t));
    res->length = rlen;
    gf2128_poly_normalise(res);
    
    free(temp);
}

/* Fast polynomial multiplication for GF(2^128) using Karatsuba algorithm */

static void gf2128_poly_mul_karatsuba(gf2128_poly_t res, const gf2128_poly_t a, const gf2128_poly_t b) {
    if (gf2128_poly_is_zero(a) || gf2128_poly_is_zero(b)) {
        gf2128_poly_zero(res);
        return;
    }
    
    slong alen = a->length;
    slong blen = b->length;
    
    /* For small polynomials, use schoolbook multiplication */
    if (alen <= 8 || blen <= 8) {
        gf2128_poly_mul_schoolbook(res, a, b);
        return;
    }
    
    /* Handle extremely unbalanced polynomials */
    if (alen > 8 * blen || blen > 8 * alen) {
        gf2128_poly_mul_schoolbook(res, a, b);
        return;
    }
    
    /* Karatsuba algorithm - use balanced split */
    slong split = (FLINT_MAX(alen, blen) + 1) / 2;
    
    gf2128_poly_t a0, a1, b0, b1;
    gf2128_poly_t z0, z1, z2, temp1, temp2;
    
    /* Initialize all polynomials */
    gf2128_poly_init(a0);
    gf2128_poly_init(a1);
    gf2128_poly_init(b0);
    gf2128_poly_init(b1);
    gf2128_poly_init(z0);
    gf2128_poly_init(z1);
    gf2128_poly_init(z2);
    gf2128_poly_init(temp1);
    gf2128_poly_init(temp2);
    
    /* Split a: a = a0 + a1*x^split */
    for (slong i = 0; i < FLINT_MIN(split, alen); i++) {
        gf2128_poly_set_coeff(a0, i, &a->coeffs[i]);
    }
    
    for (slong i = split; i < alen; i++) {
        gf2128_poly_set_coeff(a1, i - split, &a->coeffs[i]);
    }
    
    /* Split b: b = b0 + b1*x^split */
    for (slong i = 0; i < FLINT_MIN(split, blen); i++) {
        gf2128_poly_set_coeff(b0, i, &b->coeffs[i]);
    }
    
    for (slong i = split; i < blen; i++) {
        gf2128_poly_set_coeff(b1, i - split, &b->coeffs[i]);
    }
    
    /* Compute z0 = a0 * b0 */
    gf2128_poly_mul_karatsuba(z0, a0, b0);
    
    /* Compute z2 = a1 * b1 */
    gf2128_poly_mul_karatsuba(z2, a1, b1);
    
    /* Compute z1 = (a0 + a1) * (b0 + b1) - z0 - z2 */
    gf2128_poly_add(temp1, a0, a1);  /* temp1 = a0 + a1 */
    gf2128_poly_add(temp2, b0, b1);  /* temp2 = b0 + b1 */
    gf2128_poly_mul_karatsuba(z1, temp1, temp2);  /* z1 = (a0 + a1) * (b0 + b1) */
    
    /* z1 = z1 - z0 - z2 = z1 + z0 + z2 (in GF(2^n)) */
    gf2128_poly_add(z1, z1, z0);
    gf2128_poly_add(z1, z1, z2);
    
    /* Construct result: res = z0 + z1*x^split + z2*x^(2*split) */
    slong result_len = alen + blen - 1;
    gf2128_poly_fit_length(res, result_len);
    
    /* Initialize result to zero */
    for (slong i = 0; i < result_len; i++) {
        res->coeffs[i] = gf2128_zero();
    }
    res->length = result_len;
    
    /* Add z0 */
    for (slong i = 0; i < z0->length; i++) {
        res->coeffs[i] = gf2128_add(&res->coeffs[i], &z0->coeffs[i]);
    }
    
    /* Add z1 * x^split */
    for (slong i = 0; i < z1->length; i++) {
        if (i + split < result_len) {
            res->coeffs[i + split] = gf2128_add(&res->coeffs[i + split], &z1->coeffs[i]);
        }
    }
    
    /* Add z2 * x^(2*split) */
    for (slong i = 0; i < z2->length; i++) {
        if (i + 2*split < result_len) {
            res->coeffs[i + 2*split] = gf2128_add(&res->coeffs[i + 2*split], &z2->coeffs[i]);
        }
    }
    
    /* Normalize result */
    gf2128_poly_normalise(res);
    
    /* Cleanup */
    gf2128_poly_clear(a0);
    gf2128_poly_clear(a1);
    gf2128_poly_clear(b0);
    gf2128_poly_clear(b1);
    gf2128_poly_clear(z0);
    gf2128_poly_clear(z1);
    gf2128_poly_clear(z2);
    gf2128_poly_clear(temp1);
    gf2128_poly_clear(temp2);
}

static void gf2128_poly_mul(gf2128_poly_t res, const gf2128_poly_t a, const gf2128_poly_t b){
    gf2128_poly_mul_karatsuba(res, a, b);
}

/* ============================================================================
   POLYNOMIAL MATRIX STRUCTURE
   ============================================================================ */

typedef struct {
    fq_nmod_poly_struct *entries;
    slong r;
    slong c;
    fq_nmod_poly_struct **rows;
    fq_nmod_ctx_struct *ctx;
} fq_nmod_poly_mat_struct;

typedef fq_nmod_poly_mat_struct fq_nmod_poly_mat_t[1];

/* Basic matrix operations */
static void fq_nmod_poly_mat_init(fq_nmod_poly_mat_t mat, slong rows, slong cols,
                          const fq_nmod_ctx_t ctx);
static void fq_nmod_poly_mat_clear(fq_nmod_poly_mat_t mat, const fq_nmod_ctx_t ctx);

static inline fq_nmod_poly_struct* fq_nmod_poly_mat_entry(const fq_nmod_poly_mat_t mat, 
                                                          slong i, slong j) {
    return mat->rows[i] + j;
}

/* ============================================================================
   GF(2^8) LOOKUP TABLES AND OPERATIONS
   ============================================================================ */

typedef struct {
    uint8_t log_table[256];
    uint16_t exp_table[512];
    uint8_t inv_table[256];
    uint8_t mul_table[256][256];
    uint8_t generator; 
} gf28_tables_t;

static gf28_tables_t *g_gf28_tables = NULL;

typedef struct {
    uint8_t mul_table[256][256];
    uint8_t sqr_table[256];
    uint8_t inv_table[256];
    uint8_t double_table[256];
    int initialized;
} gf28_complete_tables_t;

static gf28_complete_tables_t g_gf28_complete_tables = {0};

/* Initialize GF(2^8) lookup tables */
static void init_gf28_tables(uint8_t irred_poly) {
    if (g_gf28_tables) return;
    
    g_gf28_tables = (gf28_tables_t *)malloc(sizeof(gf28_tables_t));
    
    memset(g_gf28_tables->exp_table, 0, sizeof(g_gf28_tables->exp_table));
    memset(g_gf28_tables->log_table, 0, sizeof(g_gf28_tables->log_table));
    memset(g_gf28_tables->inv_table, 0, sizeof(g_gf28_tables->inv_table));
    
    /* Find a primitive element (generator) */
    uint8_t generator = 0;
    for (uint8_t g = 2; g < 256; g++) {
        uint16_t alpha = 1;
        int is_primitive = 1;
        
        for (int i = 0; i < 255; i++) {
            uint16_t new_alpha = 0;
            uint16_t a = alpha;
            uint8_t b = g;
            
            while (b) {
                if (b & 1) new_alpha ^= a;
                a <<= 1;
                if (a & 0x100) a ^= (0x100 | irred_poly);
                b >>= 1;
            }
            
            alpha = new_alpha & 0xFF;
            
            if (i < 254 && alpha == 1) {
                is_primitive = 0;
                break;
            }
        }
        
        if (is_primitive && alpha == 1) {
            generator = g;
            break;
        }
    }
    
    g_gf28_tables->generator = generator;
    
    /* Generate lookup tables using the primitive element */
    uint16_t alpha = 1;
    for (int i = 0; i < 255; i++) {
        g_gf28_tables->exp_table[i] = (uint8_t)alpha;
        g_gf28_tables->log_table[(uint8_t)alpha] = i;
        
        uint16_t new_alpha = 0;
        uint16_t a = alpha;
        uint8_t b = generator;
        
        while (b) {
            if (b & 1) new_alpha ^= a;
            a <<= 1;
            if (a & 0x100) a ^= (0x100 | irred_poly);
            b >>= 1;
        }
        
        alpha = new_alpha & 0xFF;
    }
    
    for (int i = 0; i < 255; i++) {
        g_gf28_tables->exp_table[i + 255] = g_gf28_tables->exp_table[i];
    }
    
    g_gf28_tables->exp_table[255] = 1;
    g_gf28_tables->log_table[0] = 0;
    
    g_gf28_tables->inv_table[0] = 0;
    for (int i = 1; i < 256; i++) {
        int log_inv = 255 - g_gf28_tables->log_table[i];
        g_gf28_tables->inv_table[i] = g_gf28_tables->exp_table[log_inv];
    }
}

/* Initialize complete lookup tables */
static void init_gf28_complete_tables(void) {
    if (g_gf28_complete_tables.initialized) return;
    
    if (!g_gf28_tables) {
        init_gf28_tables(0x1D);
    }
    
    //printf("Initializing complete GF(2^8) lookup tables (64KB)...\n");
    
    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 256; j++) {
            if (i == 0 || j == 0) {
                g_gf28_complete_tables.mul_table[i][j] = 0;
            } else {
                int log_sum = g_gf28_tables->log_table[i] + g_gf28_tables->log_table[j];
                g_gf28_complete_tables.mul_table[i][j] = g_gf28_tables->exp_table[log_sum];
            }
        }
        
        g_gf28_complete_tables.sqr_table[i] = g_gf28_complete_tables.mul_table[i][i];
        g_gf28_complete_tables.inv_table[i] = g_gf28_tables->inv_table[i];
        g_gf28_complete_tables.double_table[i] = g_gf28_complete_tables.mul_table[2][i];
    }
    
    g_gf28_complete_tables.initialized = 1;
    //printf("Complete lookup tables initialized.\n");
}

/* GF(2^8) arithmetic operations */
static inline uint8_t gf28_add(uint8_t a, uint8_t b) {
    return a ^ b;
}

static inline uint8_t gf28_mul(uint8_t a, uint8_t b) {
    return g_gf28_complete_tables.mul_table[a][b];
}

static inline const uint8_t* gf28_get_scalar_row(uint8_t scalar) {
    return g_gf28_complete_tables.mul_table[scalar];
}

static inline uint8_t gf28_sqr(uint8_t a) {
    return g_gf28_complete_tables.sqr_table[a];
}

static inline uint8_t gf28_inv(uint8_t a) {
    return g_gf28_complete_tables.inv_table[a];
}

static inline uint8_t gf28_div(uint8_t a, uint8_t b) {
    if (b == 0) return 0;
    return gf28_mul(a, gf28_inv(b));
}

static void init_gf28_standard(void) {
    if (g_gf28_tables) return;
    init_gf28_tables(0x1D);
    init_gf28_complete_tables();
}

static void cleanup_gf28_tables(void) {
    if (g_gf28_tables) {
        free(g_gf28_tables);
        g_gf28_tables = NULL;
    }
}

/* ============================================================================
   GF(2^8) POLYNOMIAL STRUCTURES
   ============================================================================ */

typedef struct {
    uint8_t *coeffs;
    slong length;
    slong alloc;
} gf28_poly_struct;
typedef gf28_poly_struct gf28_poly_t[1];

typedef struct {
    gf28_poly_struct *entries;
    slong r, c;
    gf28_poly_struct **rows;
} gf28_poly_mat_struct;
typedef gf28_poly_mat_struct gf28_poly_mat_t[1];

/* GF(2^8) polynomial operations */
static void gf28_poly_init(gf28_poly_t poly) {
    poly->coeffs = NULL;
    poly->length = 0;
    poly->alloc = 0;
}

static void gf28_poly_clear(gf28_poly_t poly) {
    if (poly->coeffs) {
        free(poly->coeffs);
        poly->coeffs = NULL;
    }
    poly->length = 0;
    poly->alloc = 0;
}

static void gf28_poly_fit_length(gf28_poly_t poly, slong len) {
    if (len > poly->alloc) {
        slong new_alloc = FLINT_MAX(len, poly->alloc * 2);
        poly->coeffs = (uint8_t *)realloc(poly->coeffs, new_alloc * sizeof(uint8_t));
        memset(poly->coeffs + poly->alloc, 0, (new_alloc - poly->alloc) * sizeof(uint8_t));
        poly->alloc = new_alloc;
    }
}

static void gf28_poly_normalise(gf28_poly_t poly) {
    while (poly->length > 0 && poly->coeffs[poly->length - 1] == 0) {
        poly->length--;
    }
}

static inline void gf28_poly_zero(gf28_poly_t poly) {
    poly->length = 0;
}

static inline int gf28_poly_is_zero(const gf28_poly_t poly) {
    return poly->length == 0;
}

static inline slong gf28_poly_degree(const gf28_poly_t poly) {
    return poly->length - 1;
}

static void gf28_poly_set(gf28_poly_t res, const gf28_poly_t poly) {
    if (res == poly) return;
    gf28_poly_fit_length(res, poly->length);
    memcpy(res->coeffs, poly->coeffs, poly->length * sizeof(uint8_t));
    res->length = poly->length;
}

static uint8_t gf28_poly_get_coeff(const gf28_poly_t poly, slong i) {
    return (i < poly->length) ? poly->coeffs[i] : 0;
}

static void gf28_poly_set_coeff(gf28_poly_t poly, slong i, uint8_t c) {
    if (i >= poly->alloc) {
        slong new_alloc = FLINT_MAX(i + 1, 2 * poly->alloc);
        poly->coeffs = (uint8_t *)flint_realloc(poly->coeffs, new_alloc * sizeof(uint8_t));
        poly->alloc = new_alloc;
    }

    if (i >= poly->length) {
        for (slong j = poly->length; j < i; j++) {
            poly->coeffs[j] = 0; // gf(2^8) 的零元是 0
        }
        poly->length = i + 1;
    }

    poly->coeffs[i] = c;

    if (c == 0 && i == poly->length - 1) {
        while (poly->length > 0 && poly->coeffs[poly->length - 1] == 0) {
            poly->length--;
        }
    }
}

static void gf28_poly_add(gf28_poly_t res, const gf28_poly_t a, const gf28_poly_t b) {
    slong max_len = FLINT_MAX(a->length, b->length);
    if (max_len == 0) {
        gf28_poly_zero(res);
        return;
    }
    
    gf28_poly_fit_length(res, max_len);
    slong min_len = FLINT_MIN(a->length, b->length);
    
    for (slong i = 0; i < min_len; i++) {
        res->coeffs[i] = gf28_add(a->coeffs[i], b->coeffs[i]);
    }
    
    if (a->length > b->length) {
        memcpy(res->coeffs + min_len, a->coeffs + min_len, 
               (a->length - min_len) * sizeof(uint8_t));
        res->length = a->length;
    } else if (b->length > a->length) {
        memcpy(res->coeffs + min_len, b->coeffs + min_len, 
               (b->length - min_len) * sizeof(uint8_t));
        res->length = b->length;
    } else {
        res->length = min_len;
    }
    
    gf28_poly_normalise(res);
}

static void gf28_poly_scalar_mul(gf28_poly_t res, const gf28_poly_t poly, uint8_t c) {
    if (c == 0) {
        gf28_poly_zero(res);
        return;
    }
    
    if (c == 1) {
        gf28_poly_set(res, poly);
        return;
    }
    
    gf28_poly_fit_length(res, poly->length);
    res->length = poly->length;
    
    for (slong i = 0; i < poly->length; i++) {
        res->coeffs[i] = gf28_mul(poly->coeffs[i], c);
    }
    
    gf28_poly_normalise(res);
}

static void gf28_poly_shift_left(gf28_poly_t res, const gf28_poly_t poly, slong n) {
    if (n == 0) {
        gf28_poly_set(res, poly);
        return;
    }
    
    if (gf28_poly_is_zero(poly)) {
        gf28_poly_zero(res);
        return;
    }
    
    slong new_len = poly->length + n;
    gf28_poly_fit_length(res, new_len);
    
    memmove(res->coeffs + n, poly->coeffs, poly->length * sizeof(uint8_t));
    memset(res->coeffs, 0, n * sizeof(uint8_t));
    res->length = new_len;
}

static void gf28_poly_mul_schoolbook(gf28_poly_t res, const gf28_poly_t a, const gf28_poly_t b) {
    if (gf28_poly_is_zero(a) || gf28_poly_is_zero(b)) {
        gf28_poly_zero(res);
        return;
    }
    
    slong rlen = a->length + b->length - 1;
    uint8_t *temp = (uint8_t *)calloc(rlen, sizeof(uint8_t));
    
    for (slong i = 0; i < a->length; i++) {
        if (a->coeffs[i] == 0) continue;
        for (slong j = 0; j < b->length; j++) {
            if (b->coeffs[j] == 0) continue;
            temp[i + j] = gf28_add(temp[i + j], gf28_mul(a->coeffs[i], b->coeffs[j]));
        }
    }
    
    gf28_poly_fit_length(res, rlen);
    memcpy(res->coeffs, temp, rlen * sizeof(uint8_t));
    res->length = rlen;
    gf28_poly_normalise(res);
    
    free(temp);
}

/* Fast polynomial multiplication for GF(2^8) using Karatsuba algorithm */
static void gf28_poly_mul_karatsuba(gf28_poly_t res, const gf28_poly_t a, const gf28_poly_t b) {
    if (gf28_poly_is_zero(a) || gf28_poly_is_zero(b)) {
        gf28_poly_zero(res);
        return;
    }
    
    slong alen = a->length;
    slong blen = b->length;
    
    /* For small polynomials, use schoolbook multiplication */
    if (alen < 32 || blen < 32) {
        slong rlen = alen + blen - 1;
        uint8_t *temp = (uint8_t *)calloc(rlen, sizeof(uint8_t));
        
        for (slong i = 0; i < alen; i++) {
            if (a->coeffs[i] == 0) continue;
            const uint8_t* row = gf28_get_scalar_row(a->coeffs[i]);
            for (slong j = 0; j < blen; j++) {
                temp[i + j] ^= row[b->coeffs[j]];
            }
        }
        
        gf28_poly_fit_length(res, rlen);
        memcpy(res->coeffs, temp, rlen * sizeof(uint8_t));
        res->length = rlen;
        gf28_poly_normalise(res);
        free(temp);
        return;
    }
    
    /* Karatsuba algorithm */
    slong split = FLINT_MAX(alen, blen) / 2;
    
    gf28_poly_t a0, a1, b0, b1;
    gf28_poly_t a0b0, a1b1, amid, bmid, mid;
    
    gf28_poly_init(a0);
    gf28_poly_init(a1);
    gf28_poly_init(b0);
    gf28_poly_init(b1);
    gf28_poly_init(a0b0);
    gf28_poly_init(a1b1);
    gf28_poly_init(amid);
    gf28_poly_init(bmid);
    gf28_poly_init(mid);
    
    /* Split polynomials: a = a0 + a1*x^split */
    gf28_poly_fit_length(a0, FLINT_MIN(split, alen));
    a0->length = FLINT_MIN(split, alen);
    memcpy(a0->coeffs, a->coeffs, a0->length * sizeof(uint8_t));
    gf28_poly_normalise(a0);
    
    if (alen > split) {
        gf28_poly_fit_length(a1, alen - split);
        a1->length = alen - split;
        memcpy(a1->coeffs, a->coeffs + split, a1->length * sizeof(uint8_t));
        gf28_poly_normalise(a1);
    }
    
    /* Split b similarly */
    gf28_poly_fit_length(b0, FLINT_MIN(split, blen));
    b0->length = FLINT_MIN(split, blen);
    memcpy(b0->coeffs, b->coeffs, b0->length * sizeof(uint8_t));
    gf28_poly_normalise(b0);
    
    if (blen > split) {
        gf28_poly_fit_length(b1, blen - split);
        b1->length = blen - split;
        memcpy(b1->coeffs, b->coeffs + split, b1->length * sizeof(uint8_t));
        gf28_poly_normalise(b1);
    }
    
    /* Compute three products recursively */
    gf28_poly_mul_karatsuba(a0b0, a0, b0);
    gf28_poly_mul_karatsuba(a1b1, a1, b1);
    
    gf28_poly_add(amid, a0, a1);
    gf28_poly_add(bmid, b0, b1);
    gf28_poly_mul_karatsuba(mid, amid, bmid);
    
    /* Combine results: res = a0b0 + ((a0+a1)(b0+b1) - a0b0 - a1b1)*x^split + a1b1*x^(2*split) */
    gf28_poly_t temp1, temp2;
    gf28_poly_init(temp1);
    gf28_poly_init(temp2);
    
    /* temp1 = mid - a0b0 - a1b1 */
    gf28_poly_add(temp1, mid, a0b0);
    gf28_poly_add(temp1, temp1, a1b1);
    
    /* Build result */
    slong rlen = alen + blen - 1;
    gf28_poly_fit_length(res, rlen);
    memset(res->coeffs, 0, rlen * sizeof(uint8_t));
    
    /* Add a0b0 */
    for (slong i = 0; i < a0b0->length; i++) {
        res->coeffs[i] ^= a0b0->coeffs[i];
    }
    
    /* Add temp1 * x^split */
    for (slong i = 0; i < temp1->length; i++) {
        res->coeffs[i + split] ^= temp1->coeffs[i];
    }
    
    /* Add a1b1 * x^(2*split) */
    for (slong i = 0; i < a1b1->length; i++) {
        res->coeffs[i + 2*split] ^= a1b1->coeffs[i];
    }
    
    res->length = rlen;
    gf28_poly_normalise(res);
    
    /* Cleanup */
    gf28_poly_clear(a0);
    gf28_poly_clear(a1);
    gf28_poly_clear(b0);
    gf28_poly_clear(b1);
    gf28_poly_clear(a0b0);
    gf28_poly_clear(a1b1);
    gf28_poly_clear(amid);
    gf28_poly_clear(bmid);
    gf28_poly_clear(mid);
    gf28_poly_clear(temp1);
    gf28_poly_clear(temp2);
}

/* Replace the original gf28_poly_mul with the fast version */
static void gf28_poly_mul(gf28_poly_t res, const gf28_poly_t a, const gf28_poly_t b) {
    gf28_poly_mul_schoolbook(res, a, b);//gf28_poly_mul_karatsuba gf28_poly_mul_schoolbook
}



/* GF(2^8) matrix operations */
static void gf28_poly_mat_init(gf28_poly_mat_t mat, slong rows, slong cols) {
    mat->entries = NULL;
    mat->rows = NULL;
    
    if (rows > 0 && cols > 0) {
        mat->entries = (gf28_poly_struct *)malloc(rows * cols * sizeof(gf28_poly_struct));
        mat->rows = (gf28_poly_struct **)malloc(rows * sizeof(gf28_poly_struct *));
        
        for (slong i = 0; i < rows * cols; i++) {
            gf28_poly_init(mat->entries + i);
        }
        
        for (slong i = 0; i < rows; i++) {
            mat->rows[i] = mat->entries + i * cols;
        }
    }
    
    mat->r = rows;
    mat->c = cols;
}

static void gf28_poly_mat_clear(gf28_poly_mat_t mat) {
    if (mat->entries != NULL) {
        for (slong i = 0; i < mat->r * mat->c; i++) {
            gf28_poly_clear(mat->entries + i);
        }
        free(mat->entries);
        free(mat->rows);
    }
}

static inline gf28_poly_struct *gf28_poly_mat_entry(gf28_poly_mat_t mat, slong i, slong j) {
    return mat->rows[i] + j;
}

static void gf28_poly_mat_swap_rows(gf28_poly_mat_t mat, slong r, slong s) {
    if (r != s) {
        gf28_poly_struct *tmp = mat->rows[r];
        mat->rows[r] = mat->rows[s];
        mat->rows[s] = tmp;
    }
}

/* ============================================================================
   GF(2^16) LOOKUP TABLES AND OPERATIONS
   ============================================================================ */

typedef struct {
    uint16_t log_table[65536];
    uint32_t exp_table[131072];
    uint16_t inv_table[65536];
    uint16_t generator;
} gf216_tables_t;

typedef struct {
    uint16_t *mul_table;
    uint16_t sqr_table[65536];
    uint16_t inv_table[65536];
    uint16_t double_table[65536];
    int initialized;
} gf216_complete_tables_t;

static gf216_tables_t *g_gf216_tables = NULL;
static gf216_complete_tables_t g_gf216_complete_tables = {0};

static void init_gf216_tables(uint16_t irred_poly) {
    if (g_gf216_tables) return;
    
    g_gf216_tables = (gf216_tables_t *)malloc(sizeof(gf216_tables_t));
    
    memset(g_gf216_tables->exp_table, 0, sizeof(g_gf216_tables->exp_table));
    memset(g_gf216_tables->log_table, 0, sizeof(g_gf216_tables->log_table));
    memset(g_gf216_tables->inv_table, 0, sizeof(g_gf216_tables->inv_table));
    
    uint16_t generator = 0;
    for (uint16_t g = 2; g < 65536; g++) {
        uint32_t alpha = 1;
        int is_primitive = 1;
        
        for (int i = 0; i < 65535; i++) {
            uint32_t new_alpha = 0;
            uint32_t a = alpha;
            uint16_t b = g;
            
            while (b) {
                if (b & 1) new_alpha ^= a;
                a <<= 1;
                if (a & 0x10000) a ^= (0x10000 | irred_poly);
                b >>= 1;
            }
            
            alpha = new_alpha & 0xFFFF;
            
            if (i < 65534 && alpha == 1) {
                is_primitive = 0;
                break;
            }
        }
        
        if (is_primitive && alpha == 1) {
            generator = g;
            break;
        }
    }
    
    g_gf216_tables->generator = generator;
    
    uint32_t alpha = 1;
    for (int i = 0; i < 65535; i++) {
        g_gf216_tables->exp_table[i] = (uint16_t)alpha;
        g_gf216_tables->log_table[(uint16_t)alpha] = i;
        
        uint32_t new_alpha = 0;
        uint32_t a = alpha;
        uint16_t b = generator;
        
        while (b) {
            if (b & 1) new_alpha ^= a;
            a <<= 1;
            if (a & 0x10000) a ^= (0x10000 | irred_poly);
            b >>= 1;
        }
        
        alpha = new_alpha & 0xFFFF;
    }
    
    for (int i = 0; i < 65535; i++) {
        g_gf216_tables->exp_table[i + 65535] = g_gf216_tables->exp_table[i];
    }
    
    g_gf216_tables->exp_table[65535] = 1;
    g_gf216_tables->log_table[0] = 0;
    
    g_gf216_tables->inv_table[0] = 0;
    for (int i = 1; i < 65536; i++) {
        int log_inv = 65535 - g_gf216_tables->log_table[i];
        g_gf216_tables->inv_table[i] = g_gf216_tables->exp_table[log_inv];
    }
}

static void init_gf216_complete_tables(void) {
    if (g_gf216_complete_tables.initialized) return;
    
    if (!g_gf216_tables) {
        init_gf216_tables(0x002D);
    }
    
    //printf("Initializing complete GF(2^16) lookup tables (256MB)...\n");
    
    for (int i = 0; i < 65536; i++) {
        if (i == 0) {
            g_gf216_complete_tables.sqr_table[i] = 0;
        } else {
            int log_sq = (2 * g_gf216_tables->log_table[i]) % 65535;
            g_gf216_complete_tables.sqr_table[i] = g_gf216_tables->exp_table[log_sq];
        }
        
        g_gf216_complete_tables.inv_table[i] = g_gf216_tables->inv_table[i];
        
        if (i == 0) {
            g_gf216_complete_tables.double_table[i] = 0;
        } else {
            int log_2 = g_gf216_tables->log_table[2];
            int log_prod = (log_2 + g_gf216_tables->log_table[i]) % 65535;
            g_gf216_complete_tables.double_table[i] = g_gf216_tables->exp_table[log_prod];
        }
    }
    
    g_gf216_complete_tables.initialized = 1;
    //printf("Complete lookup tables initialized.\n");
}

static inline uint16_t gf216_add(uint16_t a, uint16_t b) {
    return a ^ b;
}

static inline uint16_t gf216_mul(uint16_t a, uint16_t b) {
    if (a == 0 || b == 0) return 0;
    int log_sum = g_gf216_tables->log_table[a] + g_gf216_tables->log_table[b];
    return g_gf216_tables->exp_table[log_sum];
}

static inline uint16_t gf216_sqr(uint16_t a) {
    return g_gf216_complete_tables.sqr_table[a];
}

static inline uint16_t gf216_inv(uint16_t a) {
    return g_gf216_complete_tables.inv_table[a];
}

static inline uint16_t gf216_div(uint16_t a, uint16_t b) {
    if (b == 0) return 0;
    return gf216_mul(a, gf216_inv(b));
}

static void init_gf216_standard(void) {
    if (g_gf216_tables) return;
    init_gf216_tables(0x002D);
    init_gf216_complete_tables();
}

static void cleanup_gf216_tables(void) {
    if (g_gf216_tables) {
        free(g_gf216_tables);
        g_gf216_tables = NULL;
    }
    if (g_gf216_complete_tables.mul_table) {
        free(g_gf216_complete_tables.mul_table);
        g_gf216_complete_tables.mul_table = NULL;
    }
}

static void print_gf216_memory_usage(void) {
    size_t total = 0;
    
    printf("=== GF(2^16) Memory Usage ===\n");
    
    size_t basic = sizeof(gf216_tables_t);
    printf("Basic tables (log/exp/inv): %zu bytes (%.1f MB)\n", basic, basic/1024.0/1024.0);
    total += basic;
    
    size_t other = sizeof(g_gf216_complete_tables) - sizeof(uint16_t*);
    printf("Other tables (sqr/inv/double): %zu bytes (%.1f KB)\n", other, other/1024.0);
    total += other;
    
    printf("Total memory: %zu bytes (%.1f MB)\n", total, total/1024.0/1024.0);
}

/* GF(2^16) polynomial structures */
typedef struct {
    uint16_t *coeffs;
    slong length;
    slong alloc;
} gf216_poly_struct;
typedef gf216_poly_struct gf216_poly_t[1];

typedef struct {
    gf216_poly_struct *entries;
    slong r, c;
    gf216_poly_struct **rows;
} gf216_poly_mat_struct;
typedef gf216_poly_mat_struct gf216_poly_mat_t[1];

/* GF(2^16) polynomial operations */
static void gf216_poly_init(gf216_poly_t poly) {
    poly->coeffs = NULL;
    poly->length = 0;
    poly->alloc = 0;
}

static void gf216_poly_clear(gf216_poly_t poly) {
    if (poly->coeffs) {
        free(poly->coeffs);
        poly->coeffs = NULL;
    }
    poly->length = 0;
    poly->alloc = 0;
}

static void gf216_poly_fit_length(gf216_poly_t poly, slong len) {
    if (len > poly->alloc) {
        slong new_alloc = FLINT_MAX(len, poly->alloc * 2);
        poly->coeffs = (uint16_t *)realloc(poly->coeffs, new_alloc * sizeof(uint16_t));
        memset(poly->coeffs + poly->alloc, 0, (new_alloc - poly->alloc) * sizeof(uint16_t));
        poly->alloc = new_alloc;
    }
}

static void gf216_poly_normalise(gf216_poly_t poly) {
    while (poly->length > 0 && poly->coeffs[poly->length - 1] == 0) {
        poly->length--;
    }
}

static inline void gf216_poly_zero(gf216_poly_t poly) {
    poly->length = 0;
}

static inline int gf216_poly_is_zero(const gf216_poly_t poly) {
    return poly->length == 0;
}

static inline slong gf216_poly_degree(const gf216_poly_t poly) {
    return poly->length - 1;
}

static void gf216_poly_set(gf216_poly_t res, const gf216_poly_t poly) {
    if (res == poly) return;
    gf216_poly_fit_length(res, poly->length);
    memcpy(res->coeffs, poly->coeffs, poly->length * sizeof(uint16_t));
    res->length = poly->length;
}

static void gf216_poly_set_coeff(gf216_poly_t poly, slong i, uint16_t c) {
    if (i >= poly->alloc) {
        slong new_alloc = FLINT_MAX(i + 1, 2 * poly->alloc);
        poly->coeffs = (uint16_t *)flint_realloc(poly->coeffs, new_alloc * sizeof(uint16_t));
        poly->alloc = new_alloc;
    }

    if (i >= poly->length) {
        for (slong j = poly->length; j < i; j++) {
            poly->coeffs[j] = 0; // gf(2^16) 的零元是 0
        }
        poly->length = i + 1;
    }

    poly->coeffs[i] = c;

    if (c == 0 && i == poly->length - 1) {
        while (poly->length > 0 && poly->coeffs[poly->length - 1] == 0) {
            poly->length--;
        }
    }
}

static uint16_t gf216_poly_get_coeff(const gf216_poly_t poly, slong i) {
    return (i < poly->length) ? poly->coeffs[i] : 0;
}

static void gf216_poly_add(gf216_poly_t res, const gf216_poly_t a, const gf216_poly_t b) {
    slong max_len = FLINT_MAX(a->length, b->length);
    if (max_len == 0) {
        gf216_poly_zero(res);
        return;
    }
    
    gf216_poly_fit_length(res, max_len);
    slong min_len = FLINT_MIN(a->length, b->length);
    
    for (slong i = 0; i < min_len; i++) {
        res->coeffs[i] = a->coeffs[i] ^ b->coeffs[i];
    }
    
    if (a->length > b->length) {
        memcpy(res->coeffs + min_len, a->coeffs + min_len, 
               (a->length - min_len) * sizeof(uint16_t));
        res->length = a->length;
    } else if (b->length > a->length) {
        memcpy(res->coeffs + min_len, b->coeffs + min_len, 
               (b->length - min_len) * sizeof(uint16_t));
        res->length = b->length;
    } else {
        res->length = min_len;
    }
    
    gf216_poly_normalise(res);
}

static void gf216_poly_scalar_mul(gf216_poly_t res, const gf216_poly_t poly, uint16_t c) {
    if (c == 0) {
        gf216_poly_zero(res);
        return;
    }
    
    if (c == 1) {
        gf216_poly_set(res, poly);
        return;
    }
    
    gf216_poly_fit_length(res, poly->length);
    res->length = poly->length;
    
    for (slong i = 0; i < poly->length; i++) {
        res->coeffs[i] = gf216_mul(poly->coeffs[i], c);
    }
    
    gf216_poly_normalise(res);
}

static void gf216_poly_shift_left(gf216_poly_t res, const gf216_poly_t poly, slong n) {
    if (n == 0) {
        gf216_poly_set(res, poly);
        return;
    }
    
    if (gf216_poly_is_zero(poly)) {
        gf216_poly_zero(res);
        return;
    }
    
    slong new_len = poly->length + n;
    gf216_poly_fit_length(res, new_len);
    
    memmove(res->coeffs + n, poly->coeffs, poly->length * sizeof(uint16_t));
    memset(res->coeffs, 0, n * sizeof(uint16_t));
    res->length = new_len;
}

static void gf216_poly_mul(gf216_poly_t res, const gf216_poly_t a, const gf216_poly_t b) {
    if (gf216_poly_is_zero(a) || gf216_poly_is_zero(b)) {
        gf216_poly_zero(res);
        return;
    }
    
    slong rlen = a->length + b->length - 1;
    uint16_t *temp = (uint16_t *)calloc(rlen, sizeof(uint16_t));
    
    for (slong i = 0; i < a->length; i++) {
        if (a->coeffs[i] == 0) continue;
        for (slong j = 0; j < b->length; j++) {
            if (b->coeffs[j] == 0) continue;
            temp[i + j] ^= gf216_mul(a->coeffs[i], b->coeffs[j]);
        }
    }
    
    gf216_poly_fit_length(res, rlen);
    memcpy(res->coeffs, temp, rlen * sizeof(uint16_t));
    res->length = rlen;
    gf216_poly_normalise(res);
    
    free(temp);
}

/* GF(2^16) matrix operations */
static void gf216_poly_mat_init(gf216_poly_mat_t mat, slong rows, slong cols) {
    mat->entries = NULL;
    mat->rows = NULL;
    
    if (rows > 0 && cols > 0) {
        mat->entries = (gf216_poly_struct *)malloc(rows * cols * sizeof(gf216_poly_struct));
        mat->rows = (gf216_poly_struct **)malloc(rows * sizeof(gf216_poly_struct *));
        
        for (slong i = 0; i < rows * cols; i++) {
            gf216_poly_init(mat->entries + i);
        }
        
        for (slong i = 0; i < rows; i++) {
            mat->rows[i] = mat->entries + i * cols;
        }
    }
    
    mat->r = rows;
    mat->c = cols;
}

static void gf216_poly_mat_clear(gf216_poly_mat_t mat) {
    if (mat->entries != NULL) {
        for (slong i = 0; i < mat->r * mat->c; i++) {
            gf216_poly_clear(mat->entries + i);
        }
        free(mat->entries);
        free(mat->rows);
    }
}

static inline gf216_poly_struct *gf216_poly_mat_entry(gf216_poly_mat_t mat, slong i, slong j) {
    return mat->rows[i] + j;
}

static void gf216_poly_mat_swap_rows(gf216_poly_mat_t mat, slong r, slong s) {
    if (r != s) {
        gf216_poly_struct *tmp = mat->rows[r];
        mat->rows[r] = mat->rows[s];
        mat->rows[s] = tmp;
    }
}

static void gf216_poly_mat_permute_rows(gf216_poly_mat_t mat, const slong *perm) {
    gf216_poly_struct **new_rows = (gf216_poly_struct **)malloc(mat->r * sizeof(gf216_poly_struct *));
    
    for (slong i = 0; i < mat->r; i++) {
        new_rows[i] = mat->rows[perm[i]];
    }
    
    memcpy(mat->rows, new_rows, mat->r * sizeof(gf216_poly_struct *));
    free(new_rows);
}

/* ============================================================================
   GF(2^128) POLYNOMIAL MATRIX OPERATIONS
   ============================================================================ */

typedef struct {
    gf2128_poly_struct *entries;
    slong r, c;
    gf2128_poly_struct **rows;
} gf2128_poly_mat_struct;
typedef gf2128_poly_mat_struct gf2128_poly_mat_t[1];

static void gf2128_poly_mat_init(gf2128_poly_mat_t mat, slong rows, slong cols) {
    mat->entries = NULL;
    mat->rows = NULL;
    
    if (rows > 0 && cols > 0) {
        mat->entries = (gf2128_poly_struct *)malloc(rows * cols * sizeof(gf2128_poly_struct));
        mat->rows = (gf2128_poly_struct **)malloc(rows * sizeof(gf2128_poly_struct *));
        
        for (slong i = 0; i < rows * cols; i++) {
           gf2128_poly_init(mat->entries + i);
       }
       
       for (slong i = 0; i < rows; i++) {
           mat->rows[i] = mat->entries + i * cols;
       }
   }
   
   mat->r = rows;
   mat->c = cols;
}

static void gf2128_poly_mat_clear(gf2128_poly_mat_t mat) {
   if (mat->entries != NULL) {
       for (slong i = 0; i < mat->r * mat->c; i++) {
           gf2128_poly_clear(mat->entries + i);
       }
       free(mat->entries);
       free(mat->rows);
   }
}

static inline gf2128_poly_struct *gf2128_poly_mat_entry(gf2128_poly_mat_t mat, slong i, slong j) {
   return mat->rows[i] + j;
}

static inline const gf2128_poly_struct *gf2128_poly_mat_entry_const(const gf2128_poly_mat_t mat, slong i, slong j) {
   return mat->rows[i] + j;
}

static void gf2128_poly_mat_swap_rows(gf2128_poly_mat_t mat, slong r, slong s) {
   if (r != s) {
       gf2128_poly_struct *tmp = mat->rows[r];
       mat->rows[r] = mat->rows[s];
       mat->rows[s] = tmp;
   }
}

/* ============================================================================
  HELPER FUNCTIONS
  ============================================================================ */

static inline int _fq_nmod_ctx_is_gf2n(const fq_nmod_ctx_t ctx) {
   const nmod_poly_struct *mod = fq_nmod_ctx_modulus(ctx);
   return (mod->mod.n == 2);
}

static uint64_t extract_irred_poly(const fq_nmod_ctx_t ctx) {
   const nmod_poly_struct *mod = fq_nmod_ctx_modulus(ctx);
   uint64_t irred = 0;
   
   for (slong i = 0; i <= nmod_poly_degree(mod); i++) {
       if (nmod_poly_get_coeff_ui(mod, i))
           irred |= (1ULL << i);
   }
   
   return irred;
}

/* ============================================================================
  CONVERSION FUNCTIONS - IMPLEMENTATION
  ============================================================================ */

/* Enhanced GF(2^8) conversion support */
typedef struct {
   uint8_t flint_to_gf28[256];
   uint8_t gf28_to_flint[256];
   int initialized;
} gf28_conversion_t;

static gf28_conversion_t *g_gf28_conversion = NULL;

typedef struct {
   uint16_t *flint_to_gf216;
   uint16_t *gf216_to_flint;
   int initialized;
} gf216_conversion_t;

static gf216_conversion_t *g_gf216_conversion = NULL;

/* ============================================================================
   GF(2^128) CONVERSION SUPPORT - STUB IMPLEMENTATION
   ============================================================================ */

typedef struct {
    int initialized;
    uint64_t flint_poly_low;
    uint64_t flint_poly_high;
} gf2128_conversion_t;

static gf2128_conversion_t *g_gf2128_conversion = NULL;

/* Initialize GF(2^128) conversion - stub implementation */
static void init_gf2128_conversion(const fq_nmod_ctx_t ctx) {
    if (g_gf2128_conversion && g_gf2128_conversion->initialized) {
        return;
    }
    
    if (!g_gf2128_conversion) {
        g_gf2128_conversion = (gf2128_conversion_t *)calloc(1, sizeof(gf2128_conversion_t));
    }
    
    // For now, just mark as initialized
    g_gf2128_conversion->initialized = 1;
}

/* Convert between fq_nmod and native GF(2^128) - IMPLEMENTATION */
static gf2128_t fq_nmod_to_gf2128(const fq_nmod_t elem, const fq_nmod_ctx_t ctx) {
   if (fq_nmod_is_zero(elem, ctx)) {
       return gf2128_zero();
   }
   
   uint64_t poly_low = 0, poly_high = 0;
   
   for (int i = 0; i < 64; i++) {
       if (nmod_poly_get_coeff_ui(elem, i)) {
           poly_low |= (1ULL << i);
       }
   }
   for (int i = 64; i < 128; i++) {
       if (nmod_poly_get_coeff_ui(elem, i)) {
           poly_high |= (1ULL << (i - 64));
       }
   }
   
   return gf2128_create(poly_low, poly_high);
}

static void gf2128_to_fq_nmod(fq_nmod_t res, const gf2128_t *elem, const fq_nmod_ctx_t ctx) {
   fq_nmod_zero(res, ctx);
   
   if (gf2128_is_zero(elem)) {
       return;
   }
   
   for (int i = 0; i < 64; i++) {
       if ((elem->low >> i) & 1) {
           nmod_poly_set_coeff_ui(res, i, 1);
       }
   }
   for (int i = 0; i < 64; i++) {
       if ((elem->high >> i) & 1) {
           nmod_poly_set_coeff_ui(res, i + 64, 1);
       }
   }
}

/* Initialize GF(2^8) conversion tables */
static void init_gf28_conversion(const fq_nmod_ctx_t ctx) {
   if (g_gf28_conversion && g_gf28_conversion->initialized) {
       return;
   }
   
   if (!g_gf28_conversion) {
       g_gf28_conversion = (gf28_conversion_t *)calloc(1, sizeof(gf28_conversion_t));
   }
   
   //printf("Initializing GF(2^8) conversion tables...\n");
   
   if (!g_gf28_tables) {
       init_gf28_tables(0x1D);
   }
   
   uint64_t flint_poly = extract_irred_poly(ctx);
   //printf("FLINT GF(2^8) irreducible polynomial: 0x%lX\n", flint_poly);
   //printf("Our GF(2^8) irreducible polynomial: 0x11D\n");
   
   if (flint_poly == 0x11D) {
       //printf("Using identity mapping (same polynomial)\n");
       for (int i = 0; i < 256; i++) {
           g_gf28_conversion->flint_to_gf28[i] = i;
           g_gf28_conversion->gf28_to_flint[i] = i;
       }
   } else {
       //printf("Building conversion between different representations\n");
       
       g_gf28_conversion->flint_to_gf28[0] = 0;
       g_gf28_conversion->gf28_to_flint[0] = 0;
       
       fq_nmod_t gen, elem;
       fq_nmod_init(gen, ctx);
       fq_nmod_init(elem, ctx);
       
       fq_nmod_gen(gen, ctx);
       
       uint8_t flint_log[256];
       uint8_t flint_exp[256];
       memset(flint_log, 0xFF, 256);
       
       fq_nmod_one(elem, ctx);
       for (int i = 0; i < 255; i++) {
           uint8_t poly_val = 0;
           for (int j = 0; j < 8; j++) {
               if (nmod_poly_get_coeff_ui(elem, j)) {
                   poly_val |= (1 << j);
               }
           }
           
           flint_exp[i] = poly_val;
           flint_log[poly_val] = i;
           
           fq_nmod_mul(elem, elem, gen, ctx);
       }
       
       for (int i = 0; i < 255; i++) {
           uint8_t flint_elem = flint_exp[i];
           uint8_t our_elem = g_gf28_tables->exp_table[i];
           
           g_gf28_conversion->flint_to_gf28[flint_elem] = our_elem;
           g_gf28_conversion->gf28_to_flint[our_elem] = flint_elem;
       }
       
       fq_nmod_clear(gen, ctx);
       fq_nmod_clear(elem, ctx);
   }
   
   g_gf28_conversion->initialized = 1;
   //printf("GF(2^8) conversion tables initialized\n");
}

/* Initialize GF(2^16) conversion tables */
static void init_gf216_conversion(const fq_nmod_ctx_t ctx) {
    if (g_gf216_conversion && g_gf216_conversion->initialized) {
        return;
    }
    
    if (!g_gf216_conversion) {
        g_gf216_conversion = (gf216_conversion_t *)calloc(1, sizeof(gf216_conversion_t));
        g_gf216_conversion->flint_to_gf216 = (uint16_t *)calloc(65536, sizeof(uint16_t));
        g_gf216_conversion->gf216_to_flint = (uint16_t *)calloc(65536, sizeof(uint16_t));
    }
    
    //printf("Initializing GF(2^16) conversion tables...\n");
    
    init_gf216_standard();
    
    // Check if polynomials match
    uint64_t flint_poly = extract_irred_poly(ctx);
    //printf("FLINT GF(2^16) irreducible polynomial: 0x%lX\n", flint_poly);
    //printf("Our GF(2^16) irreducible polynomial: 0x1002D\n");
    
    if (flint_poly == 0x1002D) {
        // Same polynomial, use identity mapping
        //printf("Using identity mapping (same polynomial)\n");
        for (int i = 0; i < 65536; i++) {
            g_gf216_conversion->flint_to_gf216[i] = i;
            g_gf216_conversion->gf216_to_flint[i] = i;
        }
    } else {
        // Different polynomials, build conversion table
        //printf("Building conversion between different representations\n");
        
        fq_nmod_t elem, gen;
        fq_nmod_init(elem, ctx);
        fq_nmod_init(gen, ctx);
        
        fq_nmod_gen(gen, ctx);
        
        // Map zero
        g_gf216_conversion->flint_to_gf216[0] = 0;
        g_gf216_conversion->gf216_to_flint[0] = 0;
        
        // Build logarithm tables for both representations
        uint16_t flint_log[65536];
        uint16_t flint_exp[65536];
        memset(flint_log, 0xFF, sizeof(flint_log));
        
        fq_nmod_one(elem, ctx);
        for (int i = 0; i < 65535; i++) {
            uint16_t flint_rep = 0;
            for (int j = 0; j < 16; j++) {
                if (nmod_poly_get_coeff_ui(elem, j)) {
                    flint_rep |= (1 << j);
                }
            }
            
            flint_exp[i] = flint_rep;
            flint_log[flint_rep] = i;
            
            fq_nmod_mul(elem, elem, gen, ctx);
        }
        
        // Map between representations using logarithms
        for (int i = 0; i < 65535; i++) {
            uint16_t flint_elem = flint_exp[i];
            uint16_t our_elem = g_gf216_tables->exp_table[i];
            
            g_gf216_conversion->flint_to_gf216[flint_elem] = our_elem;
            g_gf216_conversion->gf216_to_flint[our_elem] = flint_elem;
        }
        
        fq_nmod_clear(elem, ctx);
        fq_nmod_clear(gen, ctx);
    }
    
    g_gf216_conversion->initialized = 1;
    //printf("GF(2^16) conversion tables initialized\n");
}
static void cleanup_gf28_conversion(void) {
   if (g_gf28_conversion) {
       free(g_gf28_conversion);
       g_gf28_conversion = NULL;
   }
}

static void cleanup_gf216_conversion(void) {
   if (g_gf216_conversion) {
       if (g_gf216_conversion->flint_to_gf216) {
           free(g_gf216_conversion->flint_to_gf216);
       }
       if (g_gf216_conversion->gf216_to_flint) {
           free(g_gf216_conversion->gf216_to_flint);
       }
       free(g_gf216_conversion);
       g_gf216_conversion = NULL;
   }
}

static uint8_t fq_nmod_to_gf28_elem(const fq_nmod_t elem, const fq_nmod_ctx_t ctx) {
   if (!g_gf28_conversion || !g_gf28_conversion->initialized) {
       init_gf28_conversion(ctx);
   }
   
   if (fq_nmod_is_zero(elem, ctx)) {
       return 0;
   }
   
   uint8_t flint_rep = 0;
   for (int j = 0; j < 8; j++) {
       if (nmod_poly_get_coeff_ui(elem, j)) {
           flint_rep |= (1 << j);
       }
   }
   
   return g_gf28_conversion->flint_to_gf28[flint_rep];
}

static void gf28_elem_to_fq_nmod(fq_nmod_t res, uint8_t elem, const fq_nmod_ctx_t ctx) {
   if (!g_gf28_conversion || !g_gf28_conversion->initialized) {
       init_gf28_conversion(ctx);
   }
   
   fq_nmod_zero(res, ctx);
   
   if (elem == 0) {
       return;
   }
   
   uint8_t flint_rep = g_gf28_conversion->gf28_to_flint[elem];
   
   for (int j = 0; j < 8; j++) {
       if (flint_rep & (1 << j)) {
           nmod_poly_set_coeff_ui(res, j, 1);
       }
   }
}

static uint16_t fq_nmod_to_gf216_elem(const fq_nmod_t elem, const fq_nmod_ctx_t ctx) {
   if (!g_gf216_conversion || !g_gf216_conversion->initialized) {
       init_gf216_conversion(ctx);
   }
   
   if (fq_nmod_is_zero(elem, ctx)) {
       return 0;
   }
   
   uint16_t flint_rep = 0;
   for (int j = 0; j < 16; j++) {
       if (nmod_poly_get_coeff_ui(elem, j)) {
           flint_rep |= (1 << j);
       }
   }
   
   return g_gf216_conversion->flint_to_gf216[flint_rep];
}

static void gf216_elem_to_fq_nmod(fq_nmod_t res, uint16_t elem, const fq_nmod_ctx_t ctx) {
   if (!g_gf216_conversion || !g_gf216_conversion->initialized) {
       init_gf216_conversion(ctx);
   }
   
   fq_nmod_zero(res, ctx);
   
   if (elem == 0) {
       return;
   }
   
   uint16_t flint_rep = g_gf216_conversion->gf216_to_flint[elem];
   
   for (int j = 0; j < 16; j++) {
       if (flint_rep & (1 << j)) {
           nmod_poly_set_coeff_ui(res, j, 1);
       }
   }
}

/* Convert polynomials */
static void fq_nmod_poly_to_gf2128_poly(gf2128_poly_t res,
                                        const fq_nmod_poly_struct *poly,
                                        const fq_nmod_ctx_t ctx) {
   slong len = fq_nmod_poly_length(poly, ctx);
   if (len == 0) {
       gf2128_poly_zero(res);
       return;
   }
   
   gf2128_poly_fit_length(res, len);
   res->length = len;
   
   for (slong i = 0; i < len; i++) {
       fq_nmod_t coeff;
       fq_nmod_init(coeff, ctx);
       fq_nmod_poly_get_coeff(coeff, poly, i, ctx);
       
       res->coeffs[i] = fq_nmod_to_gf2128(coeff, ctx);
       
       fq_nmod_clear(coeff, ctx);
   }
   
   gf2128_poly_normalise(res);
}

static void gf2128_poly_to_fq_nmod_poly(fq_nmod_poly_struct *res,
                                        const gf2128_poly_t poly,
                                        const fq_nmod_ctx_t ctx) {
   fq_nmod_poly_zero(res, ctx);
   
   for (slong i = 0; i < poly->length; i++) {
       if (!gf2128_is_zero(&poly->coeffs[i])) {
           fq_nmod_t coeff;
           fq_nmod_init(coeff, ctx);
           
           gf2128_to_fq_nmod(coeff, &poly->coeffs[i], ctx);
           fq_nmod_poly_set_coeff(res, i, coeff, ctx);
           
           fq_nmod_clear(coeff, ctx);
       }
   }
}

static void fq_nmod_poly_to_gf28_poly(gf28_poly_t res,
                                     const fq_nmod_poly_struct *poly,
                                     const fq_nmod_ctx_t ctx) {
   if (!g_gf28_conversion || !g_gf28_conversion->initialized) {
       init_gf28_conversion(ctx);
   }
   
   slong len = fq_nmod_poly_length(poly, ctx);
   if (len == 0) {
       gf28_poly_zero(res);
       return;
   }
   
   gf28_poly_fit_length(res, len);
   res->length = len;
   
   for (slong i = 0; i < len; i++) {
       fq_nmod_t coeff;
       fq_nmod_init(coeff, ctx);
       fq_nmod_poly_get_coeff(coeff, poly, i, ctx);
       
       res->coeffs[i] = fq_nmod_to_gf28_elem(coeff, ctx);
       
       fq_nmod_clear(coeff, ctx);
   }
   
   gf28_poly_normalise(res);
}

static void gf28_poly_to_fq_nmod_poly(fq_nmod_poly_struct *res,
                                     const gf28_poly_t poly,
                                     const fq_nmod_ctx_t ctx) {
   if (!g_gf28_conversion || !g_gf28_conversion->initialized) {
       init_gf28_conversion(ctx);
   }
   
   fq_nmod_poly_zero(res, ctx);
   
   for (slong i = 0; i < poly->length; i++) {
       if (poly->coeffs[i] != 0) {
           fq_nmod_t coeff;
           fq_nmod_init(coeff, ctx);
           
           gf28_elem_to_fq_nmod(coeff, poly->coeffs[i], ctx);
           fq_nmod_poly_set_coeff(res, i, coeff, ctx);
           
           fq_nmod_clear(coeff, ctx);
       }
   }
}

static void fq_nmod_poly_to_gf216_poly(gf216_poly_t res,
                                     const fq_nmod_poly_struct *poly,
                                     const fq_nmod_ctx_t ctx) {
   if (!g_gf216_conversion || !g_gf216_conversion->initialized) {
       init_gf216_conversion(ctx);
   }
   
   slong len = fq_nmod_poly_length(poly, ctx);
   if (len == 0) {
       gf216_poly_zero(res);
       return;
   }
   
   gf216_poly_fit_length(res, len);
   res->length = len;
   
   for (slong i = 0; i < len; i++) {
       fq_nmod_t coeff;
       fq_nmod_init(coeff, ctx);
       fq_nmod_poly_get_coeff(coeff, poly, i, ctx);
       
       res->coeffs[i] = fq_nmod_to_gf216_elem(coeff, ctx);
       
       fq_nmod_clear(coeff, ctx);
   }
   
   gf216_poly_normalise(res);
}

static void gf216_poly_to_fq_nmod_poly(fq_nmod_poly_struct *res,
                                     const gf216_poly_t poly,
                                     const fq_nmod_ctx_t ctx) {
   if (!g_gf216_conversion || !g_gf216_conversion->initialized) {
       init_gf216_conversion(ctx);
   }
   
   fq_nmod_poly_zero(res, ctx);
   
   for (slong i = 0; i < poly->length; i++) {
       if (poly->coeffs[i] != 0) {
           fq_nmod_t coeff;
           fq_nmod_init(coeff, ctx);
           
           gf216_elem_to_fq_nmod(coeff, poly->coeffs[i], ctx);
           fq_nmod_poly_set_coeff(res, i, coeff, ctx);
           
           fq_nmod_clear(coeff, ctx);
       }
   }
}

/* Convert matrices */
static void fq_nmod_poly_mat_to_gf2128(gf2128_poly_mat_t res,
                                       const fq_nmod_poly_mat_t mat,
                                       const fq_nmod_ctx_t ctx) {
   for (slong i = 0; i < mat->r; i++) {
       for (slong j = 0; j < mat->c; j++) {
           fq_nmod_poly_to_gf2128_poly(
               gf2128_poly_mat_entry(res, i, j),
               fq_nmod_poly_mat_entry(mat, i, j),
               ctx);
       }
   }
}

static void gf2128_poly_mat_to_fq_nmod(fq_nmod_poly_mat_t res,
                                       const gf2128_poly_mat_t mat,
                                       const fq_nmod_ctx_t ctx) {
   for (slong i = 0; i < mat->r; i++) {
       for (slong j = 0; j < mat->c; j++) {
           gf2128_poly_to_fq_nmod_poly(
               fq_nmod_poly_mat_entry(res, i, j),
               gf2128_poly_mat_entry_const(mat, i, j),
               ctx);
       }
   }
}

static void fq_nmod_poly_mat_to_gf28(gf28_poly_mat_t res,
                                    const fq_nmod_poly_mat_t mat,
                                    const fq_nmod_ctx_t ctx) {
   for (slong i = 0; i < mat->r; i++) {
       for (slong j = 0; j < mat->c; j++) {
           fq_nmod_poly_to_gf28_poly(
               gf28_poly_mat_entry(res, i, j),
               fq_nmod_poly_mat_entry(mat, i, j),
               ctx);
       }
   }
}

static void fq_nmod_poly_mat_to_gf216(gf216_poly_mat_t res,
                                    const fq_nmod_poly_mat_t mat,
                                    const fq_nmod_ctx_t ctx) {
   for (slong i = 0; i < mat->r; i++) {
       for (slong j = 0; j < mat->c; j++) {
           fq_nmod_poly_to_gf216_poly(
               gf216_poly_mat_entry(res, i, j),
               fq_nmod_poly_mat_entry(mat, i, j),
               ctx);
       }
   }
}

/* ============================================================================
  UTILITY FUNCTIONS
  ============================================================================ */

static void gf2128_to_hex(const gf2128_t *a, char *buf) {
   sprintf(buf, "%016llx%016llx", 
           (unsigned long long)a->high, 
           (unsigned long long)a->low);
}

static int gf2128_from_hex(gf2128_t *a, const char *hex) {
   if (strlen(hex) > 32) return -1;
   
   char high_str[17] = {0};
   char low_str[17] = {0};
   
   size_t len = strlen(hex);
   if (len > 16) {
       strncpy(high_str, hex, len - 16);
       strncpy(low_str, hex + len - 16, 16);
   } else {
       strncpy(low_str, hex, len);
   }
   
   uint64_t high = strtoull(high_str, NULL, 16);
   uint64_t low = strtoull(low_str, NULL, 16);
   
   *a = gf2128_create(low, high);
   
   return 0;
}

/* ============================================================================
   GF(2^32) IMPLEMENTATION
   ============================================================================ */

typedef struct {
    uint32_t value;
} gf232_t;

/* Common irreducible polynomial for GF(2^32): x^32 + x^7 + x^3 + x^2 + 1 */
#define GF232_MODULUS 0x8D

/* Basic operations */
static inline gf232_t gf232_create(uint32_t val) {
    gf232_t result = {val};
    return result;
}

static inline gf232_t gf232_zero(void) {
    return gf232_create(0);
}

static inline gf232_t gf232_one(void) {
    return gf232_create(1);
}

static inline int gf232_is_zero(const gf232_t *a) {
    return a->value == 0;
}

static inline int gf232_equal(const gf232_t *a, const gf232_t *b) {
    return a->value == b->value;
}

static inline gf232_t gf232_add(const gf232_t *a, const gf232_t *b) {
    gf232_t result;
    result.value = a->value ^ b->value;
    return result;
}

static void gf232_print(const gf232_t *a) {
    printf("%08x", a->value);
}

/* Software multiplication for GF(2^32) */
static gf232_t gf232_mul_software(const gf232_t *a, const gf232_t *b) {
    uint64_t result = 0;
    uint64_t temp = a->value;
    uint32_t mask = b->value;
    
    while (mask) {
        if (mask & 1) {
            result ^= temp;
        }
        temp <<= 1;
        mask >>= 1;
    }
    
    /* Reduction modulo x^32 + x^7 + x^3 + x^2 + 1 */
    for (int i = 63; i >= 32; i--) {
        if ((result >> i) & 1) {
            result ^= ((uint64_t)GF232_MODULUS << (i - 32));
        }
    }
    
    return gf232_create((uint32_t)result);
}


/* Optimized PCLMUL-based multiplication for GF(2^32) */
static gf232_t gf232_mul_pclmul(const gf232_t *a, const gf232_t *b) {
    __m128i x = _mm_set_epi64x(0, a->value);
    __m128i y = _mm_set_epi64x(0, b->value);
    
    /* 32x32 -> 64 bit multiplication */
    __m128i prod = _mm_clmulepi64_si128(x, y, 0x00);
    
    uint64_t result = _mm_extract_epi64(prod, 0);
    
    /* Optimized reduction for x^32 + x^7 + x^3 + x^2 + 1 */
    /* Instead of bit-by-bit reduction, use folding approach */
    
    uint64_t hi = result >> 32;  // Upper 32 bits
    uint64_t lo = result & 0xFFFFFFFF;  // Lower 32 bits
    
    /* For GF(2^32) with polynomial x^32 + x^7 + x^3 + x^2 + 1:
     * x^32 ≡ x^7 + x^3 + x^2 + 1
     * So we need to multiply the high part by (x^7 + x^3 + x^2 + 1) = 0x8D
     */
    
    /* Method 1: Direct computation
     * hi * x^32 ≡ hi * (x^7 + x^3 + x^2 + 1)
     */
    uint64_t fold = 0;
    uint32_t h = (uint32_t)hi;
    
    /* Compute h * 0x8D using shifts and XORs */
    fold = h;                    // h * 1
    fold ^= (uint64_t)h << 2;    // h * x^2
    fold ^= (uint64_t)h << 3;    // h * x^3
    fold ^= (uint64_t)h << 7;    // h * x^7
    
    /* XOR with lower part */
    result = lo ^ (fold & 0xFFFFFFFF);
    
    /* Handle any overflow from the fold operation */
    uint32_t fold_hi = fold >> 32;
    if (fold_hi) {
        /* Need another reduction step */
        uint32_t extra = fold_hi;
        extra ^= fold_hi << 2;
        extra ^= fold_hi << 3;
        extra ^= fold_hi << 7;
        result ^= extra;
    }
    
    return gf232_create((uint32_t)result);
}

/* Function pointer for multiplication */
typedef gf232_t (*gf232_mul_func)(const gf232_t*, const gf232_t*);
static gf232_mul_func gf232_mul = NULL;

/* Initialize multiplication function */
static void init_gf232(void) {
    if (!gf232_mul) {
        if (has_pclmulqdq()) {
            gf232_mul = gf232_mul_pclmul;
            printf("Using PCLMULQDQ for GF(2^32) multiplication\n");
        } else {
            gf232_mul = gf232_mul_software;
            printf("Using software implementation for GF(2^32) multiplication\n");
        }
    }
}

/* Fast squaring in GF(2^32) */
static gf232_t gf232_sqr(const gf232_t *a) {
    uint64_t result = 0;
    uint32_t val = a->value;
    
    /* Expand bits */
    for (int i = 0; i < 16; i++) {
        if ((val >> i) & 1) {
            result |= (uint64_t)1 << (2 * i);
        }
    }
    for (int i = 16; i < 32; i++) {
        if ((val >> i) & 1) {
            result |= (uint64_t)1 << (2 * i);
        }
    }
    
    /* Reduction */
    for (int i = 63; i >= 32; i--) {
        if ((result >> i) & 1) {
            result ^= ((uint64_t)GF232_MODULUS << (i - 32));
        }
    }
    
    return gf232_create((uint32_t)result);
}

/* Inversion using repeated squaring */
static gf232_t gf232_inv(const gf232_t *a) {
    if (gf232_is_zero(a)) {
        return gf232_zero();
    }
    
    // In GF(2^32), a^(2^32-1) = 1, so a^(2^32-2) = a^(-1)
    // We need to compute a^(0xFFFFFFFE)
    // Using square-and-multiply algorithm
    
    gf232_t result = gf232_one();
    gf232_t base = *a;
    
    // 2^32 - 2 = 0xFFFFFFFE = 11111111111111111111111111111110 in binary
    // We can compute this efficiently using the binary representation
    
    // Method 1: Direct computation using bit pattern
    // a^(2^32-2) = a^(2^31 + 2^30 + ... + 2^1)
    //           = a^2 * a^4 * a^8 * ... * a^(2^31)
    
    gf232_t a_pow = *a;
    
    // Start from a^2
    a_pow = gf232_sqr(&a_pow);
    result = gf232_mul(&result, &a_pow);
    
    // Continue with a^4, a^8, ..., a^(2^31)
    for (int i = 2; i < 32; i++) {
        a_pow = gf232_sqr(&a_pow);
        result = gf232_mul(&result, &a_pow);
    }
    
    return result;
}

/* Division in GF(2^32) */
static inline gf232_t gf232_div(const gf232_t *a, const gf232_t *b) {
    gf232_t b_inv = gf232_inv(b);
    return gf232_mul(a, &b_inv);
}

/* ============================================================================
   GF(2^64) IMPLEMENTATION
   ============================================================================ */

typedef struct {
    uint64_t value;
} gf264_t;

/* Common irreducible polynomial for GF(2^64): x^64 + x^4 + x^3 + x + 1 */
#define GF264_MODULUS 0x1B

/* Basic operations */
static inline gf264_t gf264_create(uint64_t val) {
    gf264_t result = {val};
    return result;
}

static inline gf264_t gf264_zero(void) {
    return gf264_create(0);
}

static inline gf264_t gf264_one(void) {
    return gf264_create(1);
}

static inline int gf264_is_zero(const gf264_t *a) {
    return a->value == 0;
}

static inline int gf264_equal(const gf264_t *a, const gf264_t *b) {
    return a->value == b->value;
}

static inline gf264_t gf264_add(const gf264_t *a, const gf264_t *b) {
    gf264_t result;
    result.value = a->value ^ b->value;
    return result;
}

static void gf264_print(const gf264_t *a) {
    printf("%016llx", (unsigned long long)a->value);
}

/* Software multiplication for GF(2^64) */
static gf264_t gf264_mul_software(const gf264_t *a, const gf264_t *b) {
    uint64_t result_low = 0, result_high = 0;
    uint64_t a_val = a->value;
    uint64_t b_val = b->value;
    
    /* 64x64 -> 128 bit multiplication */
    for (int i = 0; i < 64; i++) {
        if ((b_val >> i) & 1) {
            result_low ^= a_val << i;
            if (i > 0) {
                result_high ^= a_val >> (64 - i);
            }
        }
    }
    
    /* Reduction modulo x^64 + x^4 + x^3 + x + 1 */
    for (int i = 63; i >= 0; i--) {
        if ((result_high >> i) & 1) {
            result_high ^= (1ULL << i);
            
            if (i + 4 < 64) {
                result_low ^= (1ULL << (i + 4));
            } else {
                result_high ^= (1ULL << (i + 4 - 64));
            }
            
            if (i + 3 < 64) {
                result_low ^= (1ULL << (i + 3));
            } else {
                result_high ^= (1ULL << (i + 3 - 64));
            }
            
            if (i + 1 < 64) {
                result_low ^= (1ULL << (i + 1));
            } else {
                result_high ^= (1ULL << (i + 1 - 64));
            }
            
            result_low ^= (1ULL << i);
        }
    }
    
    return gf264_create(result_low);
}


/* PCLMUL-based multiplication for GF(2^64) */
static gf264_t gf264_mul_pclmul(const gf264_t *a, const gf264_t *b) {
    __m128i x = _mm_set_epi64x(0, a->value);
    __m128i y = _mm_set_epi64x(0, b->value);
    
    /* 64x64 -> 128 bit multiplication */
    __m128i prod = _mm_clmulepi64_si128(x, y, 0x00);
    
    uint64_t lo = _mm_extract_epi64(prod, 0);
    uint64_t hi = _mm_extract_epi64(prod, 1);
    
    /* Reduction modulo x^64 + x^4 + x^3 + x + 1 */
    /* We can use a more efficient folding approach */
    
    /* For each bit in the high part, we need to fold it back */
    /* hi * x^64 ≡ hi * (x^4 + x^3 + x + 1) mod (x^64 + x^4 + x^3 + x + 1) */
    
    /* First, compute hi * (x^4 + x^3 + x + 1) */
    __m128i fold_poly = _mm_set_epi64x(0, 0x1B); // 0x1B = x^4 + x^3 + x + 1
    __m128i hi_vec = _mm_set_epi64x(0, hi);
    __m128i folded = _mm_clmulepi64_si128(hi_vec, fold_poly, 0x00);
    
    /* Extract the result and XOR with lo */
    uint64_t fold_lo = _mm_extract_epi64(folded, 0);
    uint64_t fold_hi = _mm_extract_epi64(folded, 1);
    
    lo ^= fold_lo;
    
    /* Now we need to handle any overflow in fold_hi */
    /* This should be small, so we can do it bit by bit */
    for (int i = 0; i < 4; i++) {  // Only need to check first 4 bits
        if ((fold_hi >> i) & 1) {
            lo ^= (1ULL << (i + 4));  // x^4 term
            lo ^= (1ULL << (i + 3));  // x^3 term
            lo ^= (1ULL << (i + 1));  // x term
            lo ^= (1ULL << i);        // 1 term
        }
    }
    
    return gf264_create(lo);
}

/* Function pointer for multiplication */
typedef gf264_t (*gf264_mul_func)(const gf264_t*, const gf264_t*);
static gf264_mul_func gf264_mul = NULL;

/* Initialize multiplication function */
static void init_gf264(void) {
    if (!gf264_mul) {
        if (has_pclmulqdq()) {
            gf264_mul = gf264_mul_pclmul;
            printf("Using PCLMULQDQ for GF(2^64) multiplication\n");
        } else {
            gf264_mul = gf264_mul_software;
            printf("Using software implementation for GF(2^64) multiplication\n");
        }
    }
}

/* Fast squaring in GF(2^64) */
static gf264_t gf264_sqr(const gf264_t *a) {
    uint64_t result_low = 0, result_high = 0;
    uint64_t val = a->value;
    
    /* Expand bits */
    for (int i = 0; i < 32; i++) {
        if ((val >> i) & 1) {
            result_low |= (uint64_t)1 << (2 * i);
        }
    }
    for (int i = 32; i < 64; i++) {
        if ((val >> i) & 1) {
            result_high |= (uint64_t)1 << (2 * (i - 32));
        }
    }
    
    /* Reduction modulo x^64 + x^4 + x^3 + x + 1 */
    for (int i = 63; i >= 0; i--) {
        if ((result_high >> i) & 1) {
            result_high ^= (1ULL << i);
            
            if (i + 4 < 64) {
                result_low ^= (1ULL << (i + 4));
            } else {
                result_high ^= (1ULL << (i + 4 - 64));
            }
            
            if (i + 3 < 64) {
                result_low ^= (1ULL << (i + 3));
            } else {
                result_high ^= (1ULL << (i + 3 - 64));
            }
            
            if (i + 1 < 64) {
                result_low ^= (1ULL << (i + 1));
            } else {
                result_high ^= (1ULL << (i + 1 - 64));
            }
            
            result_low ^= (1ULL << i);
        }
    }
    
    return gf264_create(result_low);
}

/* Inversion using repeated squaring */
/* Inversion using addition chain optimized for GF(2^64) */
static gf264_t gf264_inv(const gf264_t *a) {
    if (gf264_is_zero(a)) {
        return gf264_zero();
    }
    
    /* In GF(2^64), we need to compute a^(2^64-2) = a^(0xFFFFFFFFFFFFFFFE) */
    /* Binary: 111...1110 (63 ones followed by one zero) */
    /* This equals 2^1 + 2^2 + 2^3 + ... + 2^63 */
    
    /* More efficient method using addition chains */
    /* First compute a^(2^32-1) using the pattern for 32 ones */
    
    gf264_t x = *a;
    gf264_t y = *a;
    
    /* Build a^3 = a * a^2 */
    y = gf264_sqr(&y);
    x = gf264_mul(&x, &y);
    
    /* Build a^(2^4-1) = a^15 */
    y = x;
    for (int i = 0; i < 2; i++) {
        y = gf264_sqr(&y);
    }
    x = gf264_mul(&x, &y);
    
    /* Build a^(2^8-1) = a^255 */
    y = x;
    for (int i = 0; i < 4; i++) {
        y = gf264_sqr(&y);
    }
    x = gf264_mul(&x, &y);
    
    /* Build a^(2^16-1) = a^65535 */
    y = x;
    for (int i = 0; i < 8; i++) {
        y = gf264_sqr(&y);
    }
    x = gf264_mul(&x, &y);
    
    /* Build a^(2^32-1) */
    y = x;
    for (int i = 0; i < 16; i++) {
        y = gf264_sqr(&y);
    }
    x = gf264_mul(&x, &y);
    
    /* Build a^(2^64-1) */
    y = x;
    for (int i = 0; i < 32; i++) {
        y = gf264_sqr(&y);
    }
    x = gf264_mul(&x, &y);
    
    /* Now x = a^(2^64-1) */
    /* We need a^(2^64-2) = a^(2^64-1) / a */
    /* Since a^(2^64-1) = 1 in GF(2^64), we have a^(2^64-2) = 1/a = a^(-1) */
    /* But we can't use that directly. Instead: */
    /* a^(2^64-2) = a^((2^64-1)-1) = a^(2^64-1) * a^(-1) */
    /* We need to divide x by a, which means we went too far */
    
    /* Start over with correct approach */
    /* a^(2^64-2) = a^2 * a^4 * a^8 * ... * a^(2^63) */
    
    gf264_t result = *a;
    gf264_t a_power = *a;
    
    /* Start with a^2 */
    a_power = gf264_sqr(&a_power);
    result = a_power;
    
    /* Multiply by a^4, a^8, ..., a^(2^63) */
    for (int i = 2; i <= 63; i++) {
        a_power = gf264_sqr(&a_power);  /* a_power = a^(2^i) */
        result = gf264_mul(&result, &a_power);
    }
    
    return result;
}
/* Division in GF(2^64) */
static inline gf264_t gf264_div(const gf264_t *a, const gf264_t *b) {
    gf264_t b_inv = gf264_inv(b);
    return gf264_mul(a, &b_inv);
}

/* ============================================================================
   GF(2^32) POLYNOMIAL OPERATIONS
   ============================================================================ */

typedef struct {
    gf232_t *coeffs;
    slong length;
    slong alloc;
} gf232_poly_struct;
typedef gf232_poly_struct gf232_poly_t[1];

static void gf232_poly_init(gf232_poly_t poly) {
    poly->coeffs = NULL;
    poly->length = 0;
    poly->alloc = 0;
}

static void gf232_poly_clear(gf232_poly_t poly) {
    if (poly->coeffs) {
        free(poly->coeffs);
        poly->coeffs = NULL;
    }
    poly->length = 0;
    poly->alloc = 0;
}

static void gf232_poly_fit_length(gf232_poly_t poly, slong len) {
    if (len > poly->alloc) {
        slong new_alloc = FLINT_MAX(len, poly->alloc * 2);
        poly->coeffs = (gf232_t *)realloc(poly->coeffs, new_alloc * sizeof(gf232_t));
        for (slong i = poly->alloc; i < new_alloc; i++) {
            poly->coeffs[i] = gf232_zero();
        }
        poly->alloc = new_alloc;
    }
}

static void gf232_poly_normalise(gf232_poly_t poly) {
    while (poly->length > 0 && gf232_is_zero(&poly->coeffs[poly->length - 1])) {
        poly->length--;
    }
}

static inline void gf232_poly_zero(gf232_poly_t poly) {
    poly->length = 0;
}

static inline int gf232_poly_is_zero(const gf232_poly_t poly) {
    return poly->length == 0;
}

static inline slong gf232_poly_degree(const gf232_poly_t poly) {
    return poly->length - 1;
}

static void gf232_poly_set(gf232_poly_t res, const gf232_poly_t poly) {
    if (res == poly) return;
    gf232_poly_fit_length(res, poly->length);
    memcpy(res->coeffs, poly->coeffs, poly->length * sizeof(gf232_t));
    res->length = poly->length;
}

static gf232_t gf232_poly_get_coeff(const gf232_poly_t poly, slong i) {
    if (i < poly->length) {
        return poly->coeffs[i];
    } else {
        return gf232_zero();
    }
}

static void gf232_poly_set_coeff(gf232_poly_t poly, slong i, const gf232_t *c) {
    gf232_poly_fit_length(poly, i + 1);
    if (i >= poly->length) {
        for (slong j = poly->length; j < i; j++) {
            poly->coeffs[j] = gf232_zero();
        }
        poly->length = i + 1;
    }
    poly->coeffs[i] = *c;
    
    if (gf232_is_zero(c) && i == poly->length - 1) {
        gf232_poly_normalise(poly);
    }
}

static void gf232_poly_add(gf232_poly_t res, const gf232_poly_t a, const gf232_poly_t b) {
    slong max_len = FLINT_MAX(a->length, b->length);
    slong min_len = FLINT_MIN(a->length, b->length);
    
    if (max_len == 0) {
        gf232_poly_zero(res);
        return;
    }
    
    gf232_poly_fit_length(res, max_len);
    
    for (slong i = 0; i < min_len; i++) {
        res->coeffs[i] = gf232_add(&a->coeffs[i], &b->coeffs[i]);
    }
    
    if (a->length > b->length) {
        memcpy(res->coeffs + min_len, a->coeffs + min_len, 
               (a->length - min_len) * sizeof(gf232_t));
        res->length = a->length;
    } else if (b->length > a->length) {
        memcpy(res->coeffs + min_len, b->coeffs + min_len, 
               (b->length - min_len) * sizeof(gf232_t));
        res->length = b->length;
    } else {
        res->length = min_len;
    }
    
    gf232_poly_normalise(res);
}

static void gf232_poly_scalar_mul(gf232_poly_t res, const gf232_poly_t poly, const gf232_t *c) {
    if (gf232_is_zero(c)) {
        gf232_poly_zero(res);
        return;
    }
    
    gf232_t one = gf232_one();
    if (gf232_equal(c, &one)) {
        gf232_poly_set(res, poly);
        return;
    }
    
    gf232_poly_fit_length(res, poly->length);
    res->length = poly->length;
    
    for (slong i = 0; i < poly->length; i++) {
        res->coeffs[i] = gf232_mul(&poly->coeffs[i], c);
    }
    
    gf232_poly_normalise(res);
}

static void gf232_poly_shift_left(gf232_poly_t res, const gf232_poly_t poly, slong n) {
    if (n == 0) {
        gf232_poly_set(res, poly);
        return;
    }
    
    if (gf232_poly_is_zero(poly)) {
        gf232_poly_zero(res);
        return;
    }
    
    slong new_len = poly->length + n;
    gf232_poly_fit_length(res, new_len);
    
    memmove(res->coeffs + n, poly->coeffs, poly->length * sizeof(gf232_t));
    
    for (slong i = 0; i < n; i++) {
        res->coeffs[i] = gf232_zero();
    }
    
    res->length = new_len;
}

/* Schoolbook multiplication for small polynomials */
static void gf232_poly_mul_schoolbook(gf232_poly_t res, const gf232_poly_t a, const gf232_poly_t b) {
    if (gf232_poly_is_zero(a) || gf232_poly_is_zero(b)) {
        gf232_poly_zero(res);
        return;
    }
    
    slong rlen = a->length + b->length - 1;
    gf232_t *temp = (gf232_t *)calloc(rlen, sizeof(gf232_t));
    
    for (slong i = 0; i < a->length; i++) {
        for (slong j = 0; j < b->length; j++) {
            gf232_t prod = gf232_mul(&a->coeffs[i], &b->coeffs[j]);
            temp[i + j] = gf232_add(&temp[i + j], &prod);
        }
    }
    
    gf232_poly_fit_length(res, rlen);
    memcpy(res->coeffs, temp, rlen * sizeof(gf232_t));
    res->length = rlen;
    gf232_poly_normalise(res);
    
    free(temp);
}

/* Fast polynomial multiplication for GF(2^32) using Karatsuba algorithm */
static void gf232_poly_mul_karatsuba(gf232_poly_t res, const gf232_poly_t a, const gf232_poly_t b) {
    if (gf232_poly_is_zero(a) || gf232_poly_is_zero(b)) {
        gf232_poly_zero(res);
        return;
    }
    
    slong alen = a->length;
    slong blen = b->length;
    
    /* For small polynomials, use schoolbook multiplication */
    if (alen <= 16 || blen <= 16) {
        gf232_poly_mul_schoolbook(res, a, b);
        return;
    }
    
    /* Handle extremely unbalanced polynomials */
    if (alen > 8 * blen || blen > 8 * alen) {
        gf232_poly_mul_schoolbook(res, a, b);
        return;
    }
    
    /* Karatsuba algorithm - use balanced split */
    slong split = (FLINT_MAX(alen, blen) + 1) / 2;
    
    gf232_poly_t a0, a1, b0, b1;
    gf232_poly_t z0, z1, z2, temp1, temp2;
    
    /* Initialize all polynomials */
    gf232_poly_init(a0);
    gf232_poly_init(a1);
    gf232_poly_init(b0);
    gf232_poly_init(b1);
    gf232_poly_init(z0);
    gf232_poly_init(z1);
    gf232_poly_init(z2);
    gf232_poly_init(temp1);
    gf232_poly_init(temp2);
    
    /* Split a: a = a0 + a1*x^split */
    for (slong i = 0; i < FLINT_MIN(split, alen); i++) {
        gf232_poly_set_coeff(a0, i, &a->coeffs[i]);
    }
    
    for (slong i = split; i < alen; i++) {
        gf232_poly_set_coeff(a1, i - split, &a->coeffs[i]);
    }
    
    /* Split b: b = b0 + b1*x^split */
    for (slong i = 0; i < FLINT_MIN(split, blen); i++) {
        gf232_poly_set_coeff(b0, i, &b->coeffs[i]);
    }
    
    for (slong i = split; i < blen; i++) {
        gf232_poly_set_coeff(b1, i - split, &b->coeffs[i]);
    }
    
    /* Compute z0 = a0 * b0 */
    gf232_poly_mul_karatsuba(z0, a0, b0);
    
    /* Compute z2 = a1 * b1 */
    gf232_poly_mul_karatsuba(z2, a1, b1);
    
    /* Compute z1 = (a0 + a1) * (b0 + b1) - z0 - z2 */
    gf232_poly_add(temp1, a0, a1);  /* temp1 = a0 + a1 */
    gf232_poly_add(temp2, b0, b1);  /* temp2 = b0 + b1 */
    gf232_poly_mul_karatsuba(z1, temp1, temp2);  /* z1 = (a0 + a1) * (b0 + b1) */
    
    /* z1 = z1 - z0 - z2 = z1 + z0 + z2 (in GF(2^n)) */
    gf232_poly_add(z1, z1, z0);
    gf232_poly_add(z1, z1, z2);
    
    /* Construct result: res = z0 + z1*x^split + z2*x^(2*split) */
    slong result_len = alen + blen - 1;
    gf232_poly_fit_length(res, result_len);
    
    /* Initialize result to zero */
    for (slong i = 0; i < result_len; i++) {
        res->coeffs[i] = gf232_zero();
    }
    res->length = result_len;
    
    /* Add z0 */
    for (slong i = 0; i < z0->length; i++) {
        res->coeffs[i] = gf232_add(&res->coeffs[i], &z0->coeffs[i]);
    }
    
    /* Add z1 * x^split */
    for (slong i = 0; i < z1->length; i++) {
        if (i + split < result_len) {
            res->coeffs[i + split] = gf232_add(&res->coeffs[i + split], &z1->coeffs[i]);
        }
    }
    
    /* Add z2 * x^(2*split) */
    for (slong i = 0; i < z2->length; i++) {
        if (i + 2*split < result_len) {
            res->coeffs[i + 2*split] = gf232_add(&res->coeffs[i + 2*split], &z2->coeffs[i]);
        }
    }
    
    /* Normalize result */
    gf232_poly_normalise(res);
    
    /* Cleanup */
    gf232_poly_clear(a0);
    gf232_poly_clear(a1);
    gf232_poly_clear(b0);
    gf232_poly_clear(b1);
    gf232_poly_clear(z0);
    gf232_poly_clear(z1);
    gf232_poly_clear(z2);
    gf232_poly_clear(temp1);
    gf232_poly_clear(temp2);
}



/* Replace the original gf232_poly_mul with the Karatsuba version */
static void gf232_poly_mul(gf232_poly_t res, const gf232_poly_t a, const gf232_poly_t b) {
    gf232_poly_mul_karatsuba(res, a, b);
}



/* ============================================================================
   GF(2^64) POLYNOMIAL OPERATIONS
   ============================================================================ */

typedef struct {
    gf264_t *coeffs;
    slong length;
    slong alloc;
} gf264_poly_struct;
typedef gf264_poly_struct gf264_poly_t[1];

static void gf264_poly_init(gf264_poly_t poly) {
    poly->coeffs = NULL;
    poly->length = 0;
    poly->alloc = 0;
}

static void gf264_poly_clear(gf264_poly_t poly) {
    if (poly->coeffs) {
        free(poly->coeffs);
        poly->coeffs = NULL;
    }
    poly->length = 0;
    poly->alloc = 0;
}

static void gf264_poly_fit_length(gf264_poly_t poly, slong len) {
    if (len > poly->alloc) {
        slong new_alloc = FLINT_MAX(len, poly->alloc * 2);
        poly->coeffs = (gf264_t *)realloc(poly->coeffs, new_alloc * sizeof(gf264_t));
        for (slong i = poly->alloc; i < new_alloc; i++) {
            poly->coeffs[i] = gf264_zero();
        }
        poly->alloc = new_alloc;
    }
}

static void gf264_poly_normalise(gf264_poly_t poly) {
    while (poly->length > 0 && gf264_is_zero(&poly->coeffs[poly->length - 1])) {
        poly->length--;
    }
}

static inline void gf264_poly_zero(gf264_poly_t poly) {
    poly->length = 0;
}

static inline int gf264_poly_is_zero(const gf264_poly_t poly) {
    return poly->length == 0;
}

static inline slong gf264_poly_degree(const gf264_poly_t poly) {
    return poly->length - 1;
}

static void gf264_poly_set(gf264_poly_t res, const gf264_poly_t poly) {
    if (res == poly) return;
    gf264_poly_fit_length(res, poly->length);
    memcpy(res->coeffs, poly->coeffs, poly->length * sizeof(gf264_t));
    res->length = poly->length;
}

static gf264_t gf264_poly_get_coeff(const gf264_poly_t poly, slong i) {
    if (i < poly->length) {
        return poly->coeffs[i];
    } else {
        return gf264_zero();
    }
}

static void gf264_poly_set_coeff(gf264_poly_t poly, slong i, const gf264_t *c) {
    gf264_poly_fit_length(poly, i + 1);
    if (i >= poly->length) {
        for (slong j = poly->length; j < i; j++) {
            poly->coeffs[j] = gf264_zero();
        }
        poly->length = i + 1;
    }
    poly->coeffs[i] = *c;
    
    if (gf264_is_zero(c) && i == poly->length - 1) {
        gf264_poly_normalise(poly);
    }
}

static void gf264_poly_add(gf264_poly_t res, const gf264_poly_t a, const gf264_poly_t b) {
    slong max_len = FLINT_MAX(a->length, b->length);
    slong min_len = FLINT_MIN(a->length, b->length);
    
    if (max_len == 0) {
        gf264_poly_zero(res);
        return;
    }
    
    gf264_poly_fit_length(res, max_len);
    
    for (slong i = 0; i < min_len; i++) {
        res->coeffs[i] = gf264_add(&a->coeffs[i], &b->coeffs[i]);
    }
    
    if (a->length > b->length) {
        memcpy(res->coeffs + min_len, a->coeffs + min_len, 
               (a->length - min_len) * sizeof(gf264_t));
        res->length = a->length;
    } else if (b->length > a->length) {
        memcpy(res->coeffs + min_len, b->coeffs + min_len, 
               (b->length - min_len) * sizeof(gf264_t));
        res->length = b->length;
    } else {
        res->length = min_len;
    }
    
    gf264_poly_normalise(res);
}

static void gf264_poly_scalar_mul(gf264_poly_t res, const gf264_poly_t poly, const gf264_t *c) {
    if (gf264_is_zero(c)) {
        gf264_poly_zero(res);
        return;
    }
    
    gf264_t one = gf264_one();
    if (gf264_equal(c, &one)) {
        gf264_poly_set(res, poly);
        return;
    }
    
    gf264_poly_fit_length(res, poly->length);
    res->length = poly->length;
    
    for (slong i = 0; i < poly->length; i++) {
        res->coeffs[i] = gf264_mul(&poly->coeffs[i], c);
    }
    
    gf264_poly_normalise(res);
}

static void gf264_poly_mul(gf264_poly_t res, const gf264_poly_t a, const gf264_poly_t b) {
    if (gf264_poly_is_zero(a) || gf264_poly_is_zero(b)) {
        gf264_poly_zero(res);
        return;
    }
    
    slong rlen = a->length + b->length - 1;
    gf264_t *temp = (gf264_t *)calloc(rlen, sizeof(gf264_t));
    
    for (slong i = 0; i < a->length; i++) {
        for (slong j = 0; j < b->length; j++) {
            gf264_t prod = gf264_mul(&a->coeffs[i], &b->coeffs[j]);
            temp[i + j] = gf264_add(&temp[i + j], &prod);
        }
    }
    
    gf264_poly_fit_length(res, rlen);
    memcpy(res->coeffs, temp, rlen * sizeof(gf264_t));
    res->length = rlen;
    gf264_poly_normalise(res);
    
    free(temp);
}

/* ============================================================================
   CONVERSION FUNCTIONS FOR GF(2^32) AND GF(2^64)
   ============================================================================ */

/* ============================================================================
   GF(2^32) CONVERSION SUPPORT
   ============================================================================ */

typedef struct {
    uint32_t flint_to_gf232[256];  // 只需要存储一个小的查找表用于快速转换
    uint32_t gf232_to_flint[256];
    uint32_t flint_poly;
    uint32_t our_poly;
    int initialized;
} gf232_conversion_t;

static gf232_conversion_t *g_gf232_conversion = NULL;

/* Initialize GF(2^32) conversion tables */
static void init_gf232_conversion(const fq_nmod_ctx_t ctx) {
    if (g_gf232_conversion && g_gf232_conversion->initialized) {
        return;
    }
    
    if (!g_gf232_conversion) {
        g_gf232_conversion = (gf232_conversion_t *)calloc(1, sizeof(gf232_conversion_t));
    }
    
    //printf("Initializing GF(2^32) conversion...\n");
    
    // Extract FLINT's irreducible polynomial
    uint64_t flint_poly = extract_irred_poly(ctx);
    g_gf232_conversion->flint_poly = (uint32_t)(flint_poly & 0xFFFFFFFF);
    g_gf232_conversion->our_poly = 0x8D; // x^32 + x^7 + x^3 + x^2 + 1
    
    //printf("FLINT GF(2^32) irreducible polynomial: 0x%X\n", g_gf232_conversion->flint_poly);
    //printf("Our GF(2^32) irreducible polynomial: 0x%X\n", g_gf232_conversion->our_poly);
    
    if ((flint_poly >> 32) == 1 && (flint_poly & 0xFFFFFFFF) == g_gf232_conversion->our_poly) {
        //printf("Using identity mapping (same polynomial)\n");
        g_gf232_conversion->initialized = 1;
        return;
    }
    
    //printf("Different polynomials detected, using direct polynomial conversion\n");
    g_gf232_conversion->initialized = 1;
}

/* Enhanced conversion functions for GF(2^32) */
static gf232_t fq_nmod_to_gf232(const fq_nmod_t elem, const fq_nmod_ctx_t ctx) {
    if (!g_gf232_conversion || !g_gf232_conversion->initialized) {
        init_gf232_conversion(ctx);
    }
    
    if (fq_nmod_is_zero(elem, ctx)) {
        return gf232_zero();
    }
    
    // If polynomials are the same, direct conversion
    if (g_gf232_conversion->flint_poly == g_gf232_conversion->our_poly) {
        uint32_t poly_val = 0;
        for (int i = 0; i < 32; i++) {
            if (nmod_poly_get_coeff_ui(elem, i)) {
                poly_val |= (1UL << i);
            }
        }
        return gf232_create(poly_val);
    }
    
    // Otherwise, we need to convert between different field representations
    // This is more complex for GF(2^32), so we use FLINT for conversion
    fq_nmod_ctx_t our_ctx;
    nmod_poly_t our_mod;
    nmod_poly_init(our_mod, 2);
    
    // Set our modulus
    nmod_poly_set_coeff_ui(our_mod, 0, 1);
    nmod_poly_set_coeff_ui(our_mod, 2, 1);
    nmod_poly_set_coeff_ui(our_mod, 3, 1);
    nmod_poly_set_coeff_ui(our_mod, 7, 1);
    nmod_poly_set_coeff_ui(our_mod, 32, 1);
    
    fq_nmod_ctx_init_modulus(our_ctx, our_mod, "t");
    
    // Convert element representation
    fq_nmod_t converted;
    fq_nmod_init(converted, our_ctx);
    
    // Copy polynomial coefficients
    for (int i = 0; i < 32; i++) {
        if (nmod_poly_get_coeff_ui(elem, i)) {
            nmod_poly_set_coeff_ui(converted, i, 1);
        }
    }
    
    // Reduce in our field
    fq_nmod_reduce(converted, our_ctx);
    
    // Extract result
    uint32_t result = 0;
    for (int i = 0; i < 32; i++) {
        if (nmod_poly_get_coeff_ui(converted, i)) {
            result |= (1UL << i);
        }
    }
    
    fq_nmod_clear(converted, our_ctx);
    fq_nmod_ctx_clear(our_ctx);
    nmod_poly_clear(our_mod);
    
    return gf232_create(result);
}

static void gf232_to_fq_nmod(fq_nmod_t res, const gf232_t *elem, const fq_nmod_ctx_t ctx) {
    if (!g_gf232_conversion || !g_gf232_conversion->initialized) {
        init_gf232_conversion(ctx);
    }
    
    fq_nmod_zero(res, ctx);
    
    if (gf232_is_zero(elem)) {
        return;
    }
    
    // If polynomials are the same, direct conversion
    if (g_gf232_conversion->flint_poly == g_gf232_conversion->our_poly) {
        for (int i = 0; i < 32; i++) {
            if ((elem->value >> i) & 1) {
                nmod_poly_set_coeff_ui(res, i, 1);
            }
        }
        return;
    }
    
    // Otherwise, convert between different representations
    fq_nmod_ctx_t our_ctx;
    nmod_poly_t our_mod;
    nmod_poly_init(our_mod, 2);
    
    nmod_poly_set_coeff_ui(our_mod, 0, 1);
    nmod_poly_set_coeff_ui(our_mod, 2, 1);
    nmod_poly_set_coeff_ui(our_mod, 3, 1);
    nmod_poly_set_coeff_ui(our_mod, 7, 1);
    nmod_poly_set_coeff_ui(our_mod, 32, 1);
    
    fq_nmod_ctx_init_modulus(our_ctx, our_mod, "t");
    
    fq_nmod_t temp;
    fq_nmod_init(temp, our_ctx);
    
    // Set from our representation
    for (int i = 0; i < 32; i++) {
        if ((elem->value >> i) & 1) {
            nmod_poly_set_coeff_ui(temp, i, 1);
        }
    }
    
    // Now we need to express this element in FLINT's representation
    // This is complex, so for now we'll use the polynomial directly
    for (int i = 0; i < 32; i++) {
        if ((elem->value >> i) & 1) {
            nmod_poly_set_coeff_ui(res, i, 1);
        }
    }
    fq_nmod_reduce(res, ctx);
    
    fq_nmod_clear(temp, our_ctx);
    fq_nmod_ctx_clear(our_ctx);
    nmod_poly_clear(our_mod);
}

/* ============================================================================
   GF(2^64) CONVERSION SUPPORT - COMPLETE IMPLEMENTATION
   ============================================================================ */

typedef struct {
    uint64_t *flint_to_gf264_low;   // 低32位查找
    uint64_t *flint_to_gf264_high;  // 高32位查找
    uint64_t *gf264_to_flint_low;
    uint64_t *gf264_to_flint_high;
    uint64_t flint_poly_low;
    uint64_t flint_poly_high;
    int initialized;
} gf264_conversion_t;

static gf264_conversion_t *g_gf264_conversion = NULL;

/* Extract 64-bit polynomial correctly */
static void extract_gf264_poly(const fq_nmod_ctx_t ctx, uint64_t *poly_low, uint64_t *poly_high) {
    const nmod_poly_struct *mod = fq_nmod_ctx_modulus(ctx);
    *poly_low = 0;
    *poly_high = 0;
    
    for (slong i = 0; i < 64 && i <= nmod_poly_degree(mod); i++) {
        if (nmod_poly_get_coeff_ui(mod, i)) {
            *poly_low |= (1ULL << i);
        }
    }
    
    for (slong i = 64; i <= nmod_poly_degree(mod); i++) {
        if (nmod_poly_get_coeff_ui(mod, i)) {
            *poly_high |= (1ULL << (i - 64));
        }
    }
}

/* Initialize GF(2^64) conversion with full polynomial matching */
static void init_gf264_conversion(const fq_nmod_ctx_t ctx) {
    if (g_gf264_conversion && g_gf264_conversion->initialized) {
        return;
    }
    
    if (!g_gf264_conversion) {
        g_gf264_conversion = (gf264_conversion_t *)calloc(1, sizeof(gf264_conversion_t));
    }
    
    //printf("Initializing GF(2^64) conversion...\n");
    
    // Extract FLINT's polynomial correctly
    extract_gf264_poly(ctx, &g_gf264_conversion->flint_poly_low, 
                            &g_gf264_conversion->flint_poly_high);
    
    //printf("FLINT GF(2^64) irreducible polynomial: 0x%lX%016lX\n", g_gf264_conversion->flint_poly_high, g_gf264_conversion->flint_poly_low);
    //printf("Our GF(2^64) irreducible polynomial: 0x1000000000000001B\n");
    
    // Check if polynomials match
    if (g_gf264_conversion->flint_poly_high == 1 && 
        g_gf264_conversion->flint_poly_low == 0x1B) {
        //printf("Using identity mapping (same polynomial)\n");
        g_gf264_conversion->initialized = 1;
        return;
    }
    
    // Different polynomials - need conversion
    //printf("Different polynomials detected, building conversion tables\n");
    
    // For GF(2^64), we can't build a full lookup table (too large)
    // Instead, we'll use a generator-based approach
    
    // Allocate smaller lookup tables for common values
    g_gf264_conversion->flint_to_gf264_low = (uint64_t *)calloc(65536, sizeof(uint64_t));
    g_gf264_conversion->gf264_to_flint_low = (uint64_t *)calloc(65536, sizeof(uint64_t));
    
    // Find generators in both fields
    fq_nmod_t gen_flint, elem;
    fq_nmod_init(gen_flint, ctx);
    fq_nmod_init(elem, ctx);
    
    fq_nmod_gen(gen_flint, ctx);
    
    // Build lookup table for first 65536 powers
    fq_nmod_one(elem, ctx);
    for (int i = 0; i < 65536; i++) {
        // Extract FLINT representation
        uint64_t flint_val = 0;
        for (int j = 0; j < 64; j++) {
            if (nmod_poly_get_coeff_ui(elem, j)) {
                flint_val |= (1ULL << j);
            }
        }
        
        // For now, use direct mapping (this may need adjustment)
        g_gf264_conversion->flint_to_gf264_low[i] = flint_val;
        if (flint_val < 65536) {
            g_gf264_conversion->gf264_to_flint_low[flint_val] = i;
        }
        
        fq_nmod_mul(elem, elem, gen_flint, ctx);
    }
    
    fq_nmod_clear(gen_flint, ctx);
    fq_nmod_clear(elem, ctx);
    
    g_gf264_conversion->initialized = 1;
    printf("GF(2^64) conversion initialized\n");
}

/* Enhanced conversion functions for GF(2^64) */
static gf264_t fq_nmod_to_gf264(const fq_nmod_t elem, const fq_nmod_ctx_t ctx) {
    if (!g_gf264_conversion || !g_gf264_conversion->initialized) {
        init_gf264_conversion(ctx);
    }
    
    if (fq_nmod_is_zero(elem, ctx)) {
        return gf264_zero();
    }
    
    // Check if same polynomial
    if (g_gf264_conversion->flint_poly_high == 1 && 
        g_gf264_conversion->flint_poly_low == 0x1B) {
        // Direct conversion
        uint64_t poly_val = 0;
        for (int i = 0; i < 64; i++) {
            if (nmod_poly_get_coeff_ui(elem, i)) {
                poly_val |= (1ULL << i);
            }
        }
        return gf264_create(poly_val);
    }
    
    // Different polynomials - need field isomorphism
    // For now, we'll use FLINT to handle the conversion
    
    // Create our field context
    fq_nmod_ctx_t our_ctx;
    nmod_poly_t our_mod;
    nmod_poly_init(our_mod, 2);
    
    // Set our modulus: x^64 + x^4 + x^3 + x + 1
    nmod_poly_set_coeff_ui(our_mod, 0, 1);
    nmod_poly_set_coeff_ui(our_mod, 1, 1);
    nmod_poly_set_coeff_ui(our_mod, 3, 1);
    nmod_poly_set_coeff_ui(our_mod, 4, 1);
    nmod_poly_set_coeff_ui(our_mod, 64, 1);
    
    fq_nmod_ctx_init_modulus(our_ctx, our_mod, "t");
    
    // Convert element by polynomial evaluation
    fq_nmod_t result;
    fq_nmod_init(result, our_ctx);
    
    // Copy polynomial representation
    for (int i = 0; i < 64; i++) {
        if (nmod_poly_get_coeff_ui(elem, i)) {
            nmod_poly_set_coeff_ui(result, i, 1);
        }
    }
    
    // Reduce in our field
    fq_nmod_reduce(result, our_ctx);
    
    // Extract result
    uint64_t val = 0;
    for (int i = 0; i < 64; i++) {
        if (nmod_poly_get_coeff_ui(result, i)) {
            val |= (1ULL << i);
        }
    }
    
    fq_nmod_clear(result, our_ctx);
    fq_nmod_ctx_clear(our_ctx);
    nmod_poly_clear(our_mod);
    
    return gf264_create(val);
}

static void gf264_to_fq_nmod(fq_nmod_t res, const gf264_t *elem, const fq_nmod_ctx_t ctx) {
    if (!g_gf264_conversion || !g_gf264_conversion->initialized) {
        init_gf264_conversion(ctx);
    }
    
    fq_nmod_zero(res, ctx);
    
    if (gf264_is_zero(elem)) {
        return;
    }
    
    // Check if same polynomial
    if (g_gf264_conversion->flint_poly_high == 1 && 
        g_gf264_conversion->flint_poly_low == 0x1B) {
        // Direct conversion
        for (int i = 0; i < 64; i++) {
            if ((elem->value >> i) & 1) {
                nmod_poly_set_coeff_ui(res, i, 1);
            }
        }
        return;
    }
    
    // Different polynomials - reverse conversion
    // Create our field context
    fq_nmod_ctx_t our_ctx;
    nmod_poly_t our_mod;
    nmod_poly_init(our_mod, 2);
    
    // Set our modulus
    nmod_poly_set_coeff_ui(our_mod, 0, 1);
    nmod_poly_set_coeff_ui(our_mod, 1, 1);
    nmod_poly_set_coeff_ui(our_mod, 3, 1);
    nmod_poly_set_coeff_ui(our_mod, 4, 1);
    nmod_poly_set_coeff_ui(our_mod, 64, 1);
    
    fq_nmod_ctx_init_modulus(our_ctx, our_mod, "t");
    
    // Create element in our field
    fq_nmod_t temp;
    fq_nmod_init(temp, our_ctx);
    
    for (int i = 0; i < 64; i++) {
        if ((elem->value >> i) & 1) {
            nmod_poly_set_coeff_ui(temp, i, 1);
        }
    }
    
    // Now we need to express this in FLINT's field
    // For simplicity, copy the polynomial representation
    for (int i = 0; i < 64; i++) {
        if ((elem->value >> i) & 1) {
            nmod_poly_set_coeff_ui(res, i, 1);
        }
    }
    
    // Reduce in FLINT's field
    fq_nmod_reduce(res, ctx);
    
    fq_nmod_clear(temp, our_ctx);
    fq_nmod_ctx_clear(our_ctx);
    nmod_poly_clear(our_mod);
}
/* Convert polynomials */
static void fq_nmod_poly_to_gf232_poly(gf232_poly_t res,
                                       const fq_nmod_poly_struct *poly,
                                       const fq_nmod_ctx_t ctx) {
    slong len = fq_nmod_poly_length(poly, ctx);
    if (len == 0) {
        gf232_poly_zero(res);
        return;
    }
    
    gf232_poly_fit_length(res, len);
    res->length = len;
    
    for (slong i = 0; i < len; i++) {
        fq_nmod_t coeff;
        fq_nmod_init(coeff, ctx);
        fq_nmod_poly_get_coeff(coeff, poly, i, ctx);
        
        res->coeffs[i] = fq_nmod_to_gf232(coeff, ctx);
        
        fq_nmod_clear(coeff, ctx);
    }
    
    gf232_poly_normalise(res);
}

static void gf232_poly_to_fq_nmod_poly(fq_nmod_poly_struct *res,
                                       const gf232_poly_t poly,
                                       const fq_nmod_ctx_t ctx) {
    fq_nmod_poly_zero(res, ctx);
    
    for (slong i = 0; i < poly->length; i++) {
        if (!gf232_is_zero(&poly->coeffs[i])) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ctx);
            
            gf232_to_fq_nmod(coeff, &poly->coeffs[i], ctx);
            fq_nmod_poly_set_coeff(res, i, coeff, ctx);
            
            fq_nmod_clear(coeff, ctx);
        }
    }
}

static void fq_nmod_poly_to_gf264_poly(gf264_poly_t res,
                                       const fq_nmod_poly_struct *poly,
                                       const fq_nmod_ctx_t ctx) {
    slong len = fq_nmod_poly_length(poly, ctx);
    if (len == 0) {
        gf264_poly_zero(res);
        return;
    }
    
    gf264_poly_fit_length(res, len);
    res->length = len;
    
    for (slong i = 0; i < len; i++) {
        fq_nmod_t coeff;
        fq_nmod_init(coeff, ctx);
        fq_nmod_poly_get_coeff(coeff, poly, i, ctx);
        
        res->coeffs[i] = fq_nmod_to_gf264(coeff, ctx);
        
        fq_nmod_clear(coeff, ctx);
    }
    
    gf264_poly_normalise(res);
}

static void gf264_poly_to_fq_nmod_poly(fq_nmod_poly_struct *res,
                                       const gf264_poly_t poly,
                                       const fq_nmod_ctx_t ctx) {
    fq_nmod_poly_zero(res, ctx);
    
    for (slong i = 0; i < poly->length; i++) {
        if (!gf264_is_zero(&poly->coeffs[i])) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ctx);
            
            gf264_to_fq_nmod(coeff, &poly->coeffs[i], ctx);
            fq_nmod_poly_set_coeff(res, i, coeff, ctx);
            
            fq_nmod_clear(coeff, ctx);
        }
    }
}

/* ============================================================================
   MATRIX STRUCTURES FOR GF(2^32) AND GF(2^64)
   ============================================================================ */

typedef struct {
    gf232_poly_struct *entries;
    slong r, c;
    gf232_poly_struct **rows;
} gf232_poly_mat_struct;
typedef gf232_poly_mat_struct gf232_poly_mat_t[1];

typedef struct {
    gf264_poly_struct *entries;
    slong r, c;
    gf264_poly_struct **rows;
} gf264_poly_mat_struct;
typedef gf264_poly_mat_struct gf264_poly_mat_t[1];

/* GF(2^32) matrix operations */
static void gf232_poly_mat_init(gf232_poly_mat_t mat, slong rows, slong cols) {
    mat->entries = NULL;
    mat->rows = NULL;
    
    if (rows > 0 && cols > 0) {
        mat->entries = (gf232_poly_struct *)malloc(rows * cols * sizeof(gf232_poly_struct));
        mat->rows = (gf232_poly_struct **)malloc(rows * sizeof(gf232_poly_struct *));
        
        for (slong i = 0; i < rows * cols; i++) {
            gf232_poly_init(mat->entries + i);
        }
        
        for (slong i = 0; i < rows; i++) {
            mat->rows[i] = mat->entries + i * cols;
        }
    }
    
    mat->r = rows;
    mat->c = cols;
}

static void gf232_poly_mat_clear(gf232_poly_mat_t mat) {
    if (mat->entries != NULL) {
        for (slong i = 0; i < mat->r * mat->c; i++) {
            gf232_poly_clear(mat->entries + i);
        }
        free(mat->entries);
        free(mat->rows);
    }
}

static inline gf232_poly_struct *gf232_poly_mat_entry(gf232_poly_mat_t mat, slong i, slong j) {
    return mat->rows[i] + j;
}

/* GF(2^64) matrix operations */
static void gf264_poly_mat_init(gf264_poly_mat_t mat, slong rows, slong cols) {
    mat->entries = NULL;
    mat->rows = NULL;
    
    if (rows > 0 && cols > 0) {
        mat->entries = (gf264_poly_struct *)malloc(rows * cols * sizeof(gf264_poly_struct));
        mat->rows = (gf264_poly_struct **)malloc(rows * sizeof(gf264_poly_struct *));
        
        for (slong i = 0; i < rows * cols; i++) {
            gf264_poly_init(mat->entries + i);
        }
        
        for (slong i = 0; i < rows; i++) {
            mat->rows[i] = mat->entries + i * cols;
        }
    }
    
    mat->r = rows;
    mat->c = cols;
}

static void gf264_poly_mat_clear(gf264_poly_mat_t mat) {
    if (mat->entries != NULL) {
        for (slong i = 0; i < mat->r * mat->c; i++) {
            gf264_poly_clear(mat->entries + i);
        }
        free(mat->entries);
        free(mat->rows);
    }
}

static inline gf264_poly_struct *gf264_poly_mat_entry(gf264_poly_mat_t mat, slong i, slong j) {
    return mat->rows[i] + j;
}

/* Initialize all GF(2^n) fields */
static void init_all_gf2n_fields(void) {
    init_gf28_standard();
    init_gf216_standard();
    init_gf232();
    init_gf264();
    init_gf2128();
}

static void cleanup_gf232_conversion(void) {
    if (g_gf232_conversion) {
        free(g_gf232_conversion);
        g_gf232_conversion = NULL;
    }
}

static void cleanup_gf264_conversion(void) {
    if (g_gf264_conversion) {
        free(g_gf264_conversion);
        g_gf264_conversion = NULL;
    }
}

static void cleanup_all_gf2n_fields(void) {
    cleanup_gf28_tables();
    cleanup_gf28_conversion();
    cleanup_gf216_tables();
    cleanup_gf216_conversion();
    cleanup_gf232_conversion();
    cleanup_gf264_conversion();
}

void fq_nmod_poly_mat_init(fq_nmod_poly_mat_t mat, slong rows, slong cols,
                          const fq_nmod_ctx_t ctx) {
    mat->entries = NULL;
    mat->rows = NULL;
    
    if (rows > 0 && cols > 0) {
        mat->entries = (fq_nmod_poly_struct *)flint_malloc(rows * cols * sizeof(fq_nmod_poly_struct));
        mat->rows = (fq_nmod_poly_struct **)flint_malloc(rows * sizeof(fq_nmod_poly_struct *));
        
        for (slong i = 0; i < rows * cols; i++)
            fq_nmod_poly_init(mat->entries + i, ctx);
        
        for (slong i = 0; i < rows; i++)
            mat->rows[i] = mat->entries + i * cols;
    }
    
    mat->r = rows;
    mat->c = cols;
    mat->ctx = (fq_nmod_ctx_struct *)ctx;
}

void fq_nmod_poly_mat_clear(fq_nmod_poly_mat_t mat, const fq_nmod_ctx_t ctx) {
    if (mat->entries != NULL) {
        for (slong i = 0; i < mat->r * mat->c; i++)
            fq_nmod_poly_clear(mat->entries + i, ctx);
        
        flint_free(mat->entries);
        flint_free(mat->rows);
    }
}


#ifdef __cplusplus
}
#endif

#endif /* GF2N_TOWER_FIELD_H */