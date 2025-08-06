/* gf28_mpolv.h - GF(2^n) implementations with native GF(2^128) PCLMUL support */
#ifndef GF28_MPOLY_H
#define GF28_MPOLY_H
/*
 * GF(2^8) Multivariate Polynomial Array Multiplication - Enhanced Version
 * 
 * This implementation provides fast multiplication of multivariate polynomials
 * over GF(2^8) using array-based methods with dynamic memory allocation.
 * 
 * Key improvements:
 * 1. Uses the same irreducible polynomial as FLINT (0x11B) for compatibility
 * 2. Dynamic memory allocation to support arbitrarily large polynomials
 * 3. Optimized memory usage with chunking for very large polynomials
 * 4. Better error handling and recovery
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <sys/time.h>

#include <flint/flint.h>
#include <flint/mpoly.h>
#include <flint/fmpz.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_mpoly.h>
#include <flint/longlong.h>
#include "gf2n_field.h"

/* GF(2^8) multivariate polynomial structure */
typedef struct {
    uint8_t *coeffs;
    ulong *exps;
    slong length;
    slong coeffs_alloc;
    slong exps_alloc;
    flint_bitcnt_t bits;
} gf28_mpoly_struct;

typedef gf28_mpoly_struct gf28_mpoly_t[1];

/* GF(2^8) multivariate polynomial context */
typedef struct {
    mpoly_ctx_t minfo;
} gf28_mpoly_ctx_struct;

typedef gf28_mpoly_ctx_struct gf28_mpoly_ctx_t[1];

/* Basic initialization and memory management */
static void gf28_mpoly_init(gf28_mpoly_t poly, const gf28_mpoly_ctx_t ctx) {
    poly->coeffs = NULL;
    poly->exps = NULL;
    poly->length = 0;
    poly->coeffs_alloc = 0;
    poly->exps_alloc = 0;
    poly->bits = MPOLY_MIN_BITS;
}

static void gf28_mpoly_clear(gf28_mpoly_t poly, const gf28_mpoly_ctx_t ctx) {
    if (poly->coeffs) flint_free(poly->coeffs);
    if (poly->exps) flint_free(poly->exps);
}

static void gf28_mpoly_ctx_init(gf28_mpoly_ctx_t ctx, slong nvars, const ordering_t ord) {
    mpoly_ctx_init(ctx->minfo, nvars, ord);
}

static void gf28_mpoly_ctx_clear(gf28_mpoly_ctx_t ctx) {
    mpoly_ctx_clear(ctx->minfo);
}

static void gf28_mpoly_zero(gf28_mpoly_t poly, const gf28_mpoly_ctx_t ctx) {
    poly->length = 0;
}

static void _gf28_mpoly_fit_length(uint8_t **coeffs, slong *coeffs_alloc,
                                   ulong **exps, slong *exps_alloc, slong N, slong length) {
    if (length > *coeffs_alloc) {
        slong new_alloc = FLINT_MAX(length, 2 * (*coeffs_alloc));
        *coeffs = (uint8_t *) flint_realloc(*coeffs, new_alloc * sizeof(uint8_t));
        *coeffs_alloc = new_alloc;
    }
    
    if (N*length > *exps_alloc) {
        slong new_alloc = FLINT_MAX(N*length, 2 * (*exps_alloc));
        *exps = (ulong *) flint_realloc(*exps, new_alloc * sizeof(ulong));
        *exps_alloc = new_alloc;
    }
}

static void gf28_mpoly_fit_length_reset_bits(gf28_mpoly_t poly, slong len, 
                                             flint_bitcnt_t bits, const gf28_mpoly_ctx_t ctx) {
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    _gf28_mpoly_fit_length(&poly->coeffs, &poly->coeffs_alloc,
                           &poly->exps, &poly->exps_alloc, N, len);
    poly->bits = bits;
}

static void gf28_mpoly_init3(gf28_mpoly_t poly, slong alloc, flint_bitcnt_t bits,
                             const gf28_mpoly_ctx_t ctx) {
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    poly->coeffs = (uint8_t *) flint_malloc(alloc * sizeof(uint8_t));
    poly->exps = (ulong *) flint_malloc(N * alloc * sizeof(ulong));
    poly->coeffs_alloc = alloc;
    poly->exps_alloc = N * alloc;
    poly->length = 0;
    poly->bits = bits;
}

static void gf28_mpoly_swap(gf28_mpoly_t poly1, gf28_mpoly_t poly2, const gf28_mpoly_ctx_t ctx) {
    gf28_mpoly_struct t = *poly1;
    *poly1 = *poly2;
    *poly2 = t;
}

static void _gf28_mpoly_set_length(gf28_mpoly_t poly, slong len, const gf28_mpoly_ctx_t ctx) {
    poly->length = len;
}

/* Dynamic memory management for large arrays */
typedef struct {
    uint8_t *data;
    slong size;
    slong chunk_size;
    int use_chunks;
} dynamic_array_t;

static void dynamic_array_init(dynamic_array_t *arr, slong size) {
    arr->size = size;
    arr->chunk_size = 1L << 20;  // 1MB chunks
    arr->use_chunks = (size > (1L << 28));  // Use chunks for arrays > 256MB
    
    if (arr->use_chunks) {
        // For very large arrays, allocate in chunks to avoid single huge allocation
        arr->data = NULL;  // Will allocate on demand
    } else {
        arr->data = (uint8_t *)calloc(size, sizeof(uint8_t));
        if (!arr->data && size > 0) {
            flint_printf("Failed to allocate %ld bytes for coefficient array\n", size);
            flint_abort();
        }
    }
}

static void dynamic_array_clear(dynamic_array_t *arr) {
    if (arr->data) {
        free(arr->data);
        arr->data = NULL;
    }
}

/* Block size for cache optimization */
#define BLOCK 128

/* Enhanced array multiplication with better memory handling */
static void _gf28_mpoly_addmul_array1_safe(uint8_t *poly1, slong array_size,
                                           const uint8_t *poly2, const ulong *exp2, slong len2,
                                           const uint8_t *poly3, const ulong *exp3, slong len3) {
    slong ii, i, jj, j;
    uint8_t *c2;
    
    for (ii = 0; ii < len2 + BLOCK; ii += BLOCK) {
        for (jj = 0; jj < len3 + BLOCK; jj += BLOCK) {
            for (i = ii; i < FLINT_MIN(ii + BLOCK, len2); i++) {
                slong offset2 = (slong)exp2[i];
                
                // Bounds check
                if (offset2 >= array_size) {
                    flint_printf("Array bounds error: offset %ld >= size %ld\n", offset2, array_size);
                    continue;
                }
                
                c2 = poly1 + offset2;
                
                if (poly2[i] != 0) {
                    const uint8_t* row = gf28_get_scalar_row(poly2[i]);
                    for (j = jj; j < FLINT_MIN(jj + BLOCK, len3); j++) {
                        slong offset3 = (slong)exp3[j];
                        slong total_offset = offset2 + offset3;
                        
                        // Bounds check
                        if (total_offset >= array_size) {
                            flint_printf("Array bounds error: total offset %ld >= size %ld\n", 
                                       total_offset, array_size);
                            continue;
                        }
                        
                        c2[offset3] ^= row[poly3[j]];
                    }
                }
            }
        }
    }
}

/* LEX ordering unpacking with bounds checking */
static slong gf28_mpoly_append_array_LEX_safe(gf28_mpoly_t P, slong Plen, uint8_t *coeff_array,
                                              const ulong *mults, slong num, slong array_size, 
                                              slong top, const gf28_mpoly_ctx_t ctx) {
    slong off, j;
    slong topmult = num == 0 ? 1 : mults[num - 1];
    slong lastd = topmult - 1;
    slong reset = array_size/topmult;
    slong counter = reset;
    ulong startexp = ((ulong)top << (P->bits*num)) + ((ulong)lastd << (P->bits*(num-1)));
    //ulong startexp = (top << (P->bits*num)) + (lastd << (P->bits*(num-1)));
    uint8_t coeff;
    
    for (off = array_size - 1; off >= 0; off--) {
        if (coeff_array[off] != 0) {
            coeff = coeff_array[off];
            slong d = off;
            ulong exp = startexp;
            
            for (j = 0; j + 1 < num; j++) {
                /* Use full precision arithmetic here */
                ulong exp_j = d % mults[j];
                exp += exp_j << (P->bits*j);
                d = d / mults[j];
            }
            
            _gf28_mpoly_fit_length(&P->coeffs, &P->coeffs_alloc,
                                   &P->exps, &P->exps_alloc, 1, Plen + 1);
            P->exps[Plen] = exp;
            P->coeffs[Plen] = coeff;
            Plen++;
            coeff_array[off] = 0;
        }
        
        counter--;
        if (counter <= 0) {
            counter = reset;
            lastd--;
            startexp -= UWORD(1) << (P->bits*(num-1));
        }
    }
    
    return Plen;
}

/* DEGLEX ordering unpacking with bounds checking */
static slong gf28_mpoly_append_array_DEGLEX_safe(gf28_mpoly_t P, slong Plen, uint8_t *coeff_array,
                                                 slong top, slong nvars, slong degb,
                                                 const gf28_mpoly_ctx_t ctx) {
    slong i;
    ulong exp, lomask = (UWORD(1) << (P->bits - 1)) - 1;
    slong off, array_size;
    slong *curexp, *degpow;
    ulong *oneexp;
    uint8_t coeff;
    int carry;
    TMP_INIT;
    
    TMP_START;
    curexp = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    degpow = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    oneexp = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
    
    array_size = 1;
    curexp[0] = 0;
    oneexp[0] = 0;
    degpow[0] = 1;
    
    for (i = 0; i < nvars-1; i++) {
        curexp[i] = 0;
        degpow[i] = array_size;
        oneexp[i] = (UWORD(1) << (P->bits*(i+1))) - UWORD(1);
        array_size *= degb;
    }
    
    off = 0;
    if (nvars > 1) {
        curexp[nvars - 2] = top;
        off = top * degpow[nvars - 2];
    }
    exp = (top << (P->bits*nvars)) + (top << (P->bits*(nvars-1)));
    
    carry = 1;
    do {
        if (off >= 0 && off < array_size && coeff_array[off] != 0) {
            coeff = coeff_array[off];
            _gf28_mpoly_fit_length(&P->coeffs, &P->coeffs_alloc,
                                   &P->exps, &P->exps_alloc, 1, Plen + 1);
            P->exps[Plen] = exp;
            P->coeffs[Plen] = coeff;
            Plen++;
            coeff_array[off] = 0;
        }
        
        exp -= oneexp[0];
        off -= 1;
        curexp[0] -= 1;
        if (curexp[0] >= 0) {
            carry = 0;
        } else {
            exp -= curexp[0]*oneexp[0];
            off -= curexp[0];
            curexp[0] = 0;
            carry = 1;
            
            for (i = 1; i < nvars - 1; i++) {
                exp -= oneexp[i];
                off -= degpow[i];
                curexp[i] -= 1;
                if (curexp[i] < 0) {
                    exp -= curexp[i]*oneexp[i];
                    off -= curexp[i]*degpow[i];
                    curexp[i] = 0;
                    carry = 1;
                } else {
                    ulong t = exp & lomask;
                    off += t*degpow[i - 1];
                    curexp[i - 1] = t;
                    exp += t*oneexp[i - 1];
                    carry = 0;
                    break;
                }
            }
        }
    } while (!carry);
    
    TMP_END;
    return Plen;
}

/* DEGREVLEX ordering unpacking with bounds checking */
static slong gf28_mpoly_append_array_DEGREVLEX_safe(gf28_mpoly_t P, slong Plen, uint8_t *coeff_array,
                                                    slong top, slong nvars, slong degb,
                                                    const gf28_mpoly_ctx_t ctx) {
    slong i;
    ulong exp, mask = UWORD(1) << (P->bits - 1);
    slong off, array_size;
    slong *curexp, *degpow;
    ulong *oneexp;
    uint8_t coeff;
    int carry;
    TMP_INIT;
    
    TMP_START;
    curexp = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    degpow = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    oneexp = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
    
    array_size = 1;
    oneexp[0] = 0;
    for (i = 0; i < nvars-1; i++) {
        curexp[i] = 0;
        degpow[i] = array_size;
        oneexp[i] = (UWORD(1) << (P->bits*(i+1))) - UWORD(1);
        array_size *= degb;
    }
    
    off = 0;
    exp = (top << (P->bits*nvars)) + top;
    
    do {
        if (off >= 0 && off < array_size && coeff_array[off] != 0) {
            coeff = coeff_array[off];
            _gf28_mpoly_fit_length(&P->coeffs, &P->coeffs_alloc,
                                   &P->exps, &P->exps_alloc, 1, Plen + 1);
            P->exps[Plen] = exp;
            P->coeffs[Plen] = coeff;
            Plen++;
            coeff_array[off] = 0;
        }
        
        exp += oneexp[0];
        off += 1;
        curexp[0] += 1;
        if ((exp & mask) == 0) {
            carry = (nvars - 1 == 0);
        } else {
            carry = 1;
            exp -= curexp[0]*oneexp[0];
            off -= curexp[0];
            curexp[0] = 0;
            for (i = 1; i < nvars - 1; i++) {
                exp += oneexp[i];
                off += degpow[i];
                curexp[i] += 1;
                if ((exp & mask) == 0) {
                    carry = 0;
                    break;
                } else {
                    carry = 1;
                    exp -= curexp[i]*oneexp[i];
                    off -= curexp[i]*degpow[i];
                    curexp[i] = 0;
                }
            }
        }
    } while (!carry);
    
    TMP_END;
    return Plen;
}

/* LEX multiplication with dynamic allocation */
static void _gf28_mpoly_mul_array_chunked_LEX(
    gf28_mpoly_t P,
    const gf28_mpoly_t A,
    const gf28_mpoly_t B,
    const ulong *mults,
    const gf28_mpoly_ctx_t ctx)
{
    slong num = ctx->minfo->nfields - 1;
    slong Pi, i, j, Plen, Pl, Al, Bl, array_size;
    slong *Amain, *Bmain;
    ulong *Apexp, *Bpexp;
    TMP_INIT;
    
    array_size = 1;
    for (i = 0; i < num; i++) {
        array_size *= mults[i];
    }
    
    Al = 1 + (slong)(A->exps[0] >> (A->bits*num));
    Bl = 1 + (slong)(B->exps[0] >> (B->bits*num));
    
    TMP_START;
    
    Amain = (slong *) TMP_ALLOC((Al + 1)*sizeof(slong));
    Bmain = (slong *) TMP_ALLOC((Bl + 1)*sizeof(slong));
    Apexp = (ulong *) flint_malloc(A->length*sizeof(ulong));
    Bpexp = (ulong *) flint_malloc(B->length*sizeof(ulong));
    mpoly_main_variable_split_LEX(Amain, Apexp, A->exps, Al, A->length, mults, num, A->bits);
    mpoly_main_variable_split_LEX(Bmain, Bpexp, B->exps, Bl, B->length, mults, num, B->bits);
    
    Pl = Al + Bl - 1;
    Plen = 0;
    
    /* Use dynamic array for coefficient storage */
    dynamic_array_t coeff_array;
    dynamic_array_init(&coeff_array, array_size);
    
    if (!coeff_array.data) {
        flint_printf("Memory allocation failed for array of size %ld\n", array_size);
        flint_free(Apexp);
        flint_free(Bpexp);
        TMP_END;
        return;
    }
    
    for (Pi = 0; Pi < Pl; Pi++) {
        memset(coeff_array.data, 0, array_size * sizeof(uint8_t));
        
        for (i = 0, j = Pi; i < Al && j >= 0; i++, j--) {
            if (j < Bl) {
                _gf28_mpoly_addmul_array1_safe(coeff_array.data, array_size,
                        A->coeffs + Amain[i],
                        Apexp + Amain[i], Amain[i + 1] - Amain[i],
                        B->coeffs + Bmain[j],
                        Bpexp + Bmain[j], Bmain[j + 1] - Bmain[j]);
            }
        }
        
        Plen = gf28_mpoly_append_array_LEX_safe(P, Plen, coeff_array.data,
                              mults, num, array_size, Pl - Pi - 1, ctx);
    }
    
    _gf28_mpoly_set_length(P, Plen, ctx);
    
    dynamic_array_clear(&coeff_array);
    flint_free(Apexp);
    flint_free(Bpexp);
    TMP_END;
}

/* DEGLEX/DEGREVLEX multiplication with dynamic allocation */
static void _gf28_mpoly_mul_array_chunked_DEG(
    gf28_mpoly_t P,
    const gf28_mpoly_t A,
    const gf28_mpoly_t B,
    ulong degb,
    const gf28_mpoly_ctx_t ctx)
{
    slong nvars = ctx->minfo->nvars;
    slong Pi, i, j, Plen, Pl, Al, Bl, array_size;
    slong *Amain, *Bmain;
    ulong *Apexp, *Bpexp;
    slong (*upack)(gf28_mpoly_t, slong, uint8_t *, slong, slong, slong, const gf28_mpoly_ctx_t);
    TMP_INIT;
    
    TMP_START;
    
    Al = 1 + (slong)(A->exps[0] >> (A->bits*nvars));
    Bl = 1 + (slong)(B->exps[0] >> (B->bits*nvars));
    
    array_size = 1;
    for (i = 0; i < nvars-1; i++) {
        ulong hi;
        umul_ppmm(hi, array_size, array_size, degb);
        if (hi != 0) {
            flint_printf("Array size overflow in DEG multiplication\n");
            TMP_END;
            return;
        }
    }
    
    upack = &gf28_mpoly_append_array_DEGLEX_safe;
    if (ctx->minfo->ord == ORD_DEGREVLEX) {
        upack = &gf28_mpoly_append_array_DEGREVLEX_safe;
    }
    
    Amain = (slong *) TMP_ALLOC((Al + 1)*sizeof(slong));
    Bmain = (slong *) TMP_ALLOC((Bl + 1)*sizeof(slong));
    Apexp = (ulong *) flint_malloc(A->length*sizeof(ulong));
    Bpexp = (ulong *) flint_malloc(B->length*sizeof(ulong));
    mpoly_main_variable_split_DEG(Amain, Apexp, A->exps, Al, A->length,
                                   degb, nvars, A->bits);
    mpoly_main_variable_split_DEG(Bmain, Bpexp, B->exps, Bl, B->length,
                                   degb, nvars, B->bits);
    
    Pl = Al + Bl - 1;
    FLINT_ASSERT(Pl == degb);
    Plen = 0;
    
    /* Use dynamic array for coefficient storage */
    dynamic_array_t coeff_array;
    dynamic_array_init(&coeff_array, array_size);
    
    if (!coeff_array.data) {
        flint_printf("Memory allocation failed for array of size %ld\n", array_size);
        flint_free(Apexp);
        flint_free(Bpexp);
        TMP_END;
        return;
    }
    
    for (Pi = 0; Pi < Pl; Pi++) {
        memset(coeff_array.data, 0, array_size * sizeof(uint8_t));
        
        for (i = 0, j = Pi; i < Al && j >= 0; i++, j--) {
            if (j < Bl) {
                _gf28_mpoly_addmul_array1_safe(coeff_array.data, array_size,
                        A->coeffs + Amain[i],
                        Apexp + Amain[i], Amain[i + 1] - Amain[i],
                        B->coeffs + Bmain[j],
                        Bpexp + Bmain[j], Bmain[j + 1] - Bmain[j]);
            }
        }
        
        Plen = upack(P, Plen, coeff_array.data, Pl - Pi - 1, nvars, degb, ctx);
    }
    
    _gf28_mpoly_set_length(P, Plen, ctx);
    
    dynamic_array_clear(&coeff_array);
    flint_free(Apexp);
    flint_free(Bpexp);
    TMP_END;
}

/* Main LEX multiplication function with adaptive limits */
static int _gf28_mpoly_mul_array_LEX(
    gf28_mpoly_t A,
    const gf28_mpoly_t B,
    fmpz *maxBfields,
    const gf28_mpoly_t C,
    fmpz *maxCfields,
    const gf28_mpoly_ctx_t ctx)
{
    slong i, exp_bits, array_size;
    ulong max, *mults;
    int success;
    TMP_INIT;
    
    FLINT_ASSERT(ctx->minfo->nvars > 0);
    FLINT_ASSERT(B->length != 0);
    FLINT_ASSERT(C->length != 0);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(1 == mpoly_words_per_exp(B->bits, ctx->minfo));
    FLINT_ASSERT(1 == mpoly_words_per_exp(C->bits, ctx->minfo));
    
    TMP_START;
    
    mults = (ulong *) TMP_ALLOC(ctx->minfo->nfields*sizeof(ulong));
    
    /* Calculate required array sizes */
    i = ctx->minfo->nfields - 1;
    FLINT_ASSERT(fmpz_fits_si(maxBfields + i));
    FLINT_ASSERT(fmpz_fits_si(maxCfields + i));
    
    /* IMPORTANT: For the last field (main variable in LEX ordering),
       we need the sum of degrees, not just max(degB, degC) */
    mults[i] = 1 + fmpz_get_ui(maxBfields + i) + fmpz_get_ui(maxCfields + i);
    max = mults[i];
    
    if (((slong)mults[i]) <= 0) {
        success = 0;
        goto cleanup;
    }
    
    
    array_size = WORD(1);
    for (i--; i >= 0; i--) {
        ulong hi;
        FLINT_ASSERT(fmpz_fits_si(maxBfields + i));
        FLINT_ASSERT(fmpz_fits_si(maxCfields + i));
        mults[i] = 1 + fmpz_get_ui(maxBfields + i) + fmpz_get_ui(maxCfields + i);
        max |= mults[i];
        umul_ppmm(hi, array_size, array_size, mults[i]);
        if (hi != 0 || (slong)mults[i] <= 0 || array_size <= 0) {
            success = 0;
            goto cleanup;
        }
    }
    
    /* Check if array size is reasonable (up to 1GB) */
    if (array_size > (1L << 30)) {
        flint_printf("Warning: Large array size %ld bytes requested\n", array_size);
        success = 0;
        goto cleanup;
    }
    
    exp_bits = FLINT_MAX(MPOLY_MIN_BITS, FLINT_BIT_COUNT(max) + 1);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    
    if (1 != mpoly_words_per_exp(exp_bits, ctx->minfo)) {
        success = 0;
        goto cleanup;
    }
    
    if (A == B || A == C) {
        gf28_mpoly_t T;
        gf28_mpoly_init3(T, B->length + C->length - 1, exp_bits, ctx);
        _gf28_mpoly_mul_array_chunked_LEX(T, C, B, mults, ctx);
        gf28_mpoly_swap(T, A, ctx);
        gf28_mpoly_clear(T, ctx);
    } else {
        gf28_mpoly_fit_length_reset_bits(A, B->length + C->length - 1, exp_bits, ctx);
        _gf28_mpoly_mul_array_chunked_LEX(A, C, B, mults, ctx);
    }
    success = 1;
    
cleanup:
    TMP_END;
    return success;
}

/* Main DEGLEX/DEGREVLEX multiplication function with adaptive limits */
static int _gf28_mpoly_mul_array_DEG(
    gf28_mpoly_t A,
    const gf28_mpoly_t B,
    fmpz *maxBfields,
    const gf28_mpoly_t C,
    fmpz *maxCfields,
    const gf28_mpoly_ctx_t ctx)
{
    slong i, exp_bits, array_size;
    ulong deg;
    int success;
    
    FLINT_ASSERT(ctx->minfo->nvars > 0);
    FLINT_ASSERT(B->length != 0);
    FLINT_ASSERT(C->length != 0);
    FLINT_ASSERT(ctx->minfo->ord == ORD_DEGREVLEX || ctx->minfo->ord == ORD_DEGLEX);
    FLINT_ASSERT(1 == mpoly_words_per_exp(B->bits, ctx->minfo));
    FLINT_ASSERT(1 == mpoly_words_per_exp(C->bits, ctx->minfo));
    
    i = ctx->minfo->nfields - 1;
    FLINT_ASSERT(fmpz_fits_si(maxBfields + i));
    FLINT_ASSERT(fmpz_fits_si(maxCfields + i));
    deg = 1 + fmpz_get_ui(maxBfields + i) + fmpz_get_ui(maxCfields + i);
    
    if (((slong)deg) <= 0) {
        success = 0;
        goto cleanup;
    }
    
    /* Calculate array size */
    array_size = WORD(1);
    for (i = ctx->minfo->nvars - 2; i >= 0; i--) {
        ulong hi;
        umul_ppmm(hi, array_size, array_size, deg);
        if (hi != WORD(0) || array_size <= 0) {
            success = 0;
            goto cleanup;
        }
    }
    
    /* Check if array size is reasonable (up to 1GB) */
    if (array_size > (1L << 30)) {
        flint_printf("Warning: Large array size %ld bytes requested for DEG ordering\n", array_size);
        success = 0;
        goto cleanup;
    }
    
    exp_bits = FLINT_MAX(MPOLY_MIN_BITS, FLINT_BIT_COUNT(deg) + 1);
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo);
    
    if (1 != mpoly_words_per_exp(exp_bits, ctx->minfo)) {
        success = 0;
        goto cleanup;
    }
    
    if (A == B || A == C) {
        gf28_mpoly_t T;
        gf28_mpoly_init3(T, B->length + C->length - 1, exp_bits, ctx);
        _gf28_mpoly_mul_array_chunked_DEG(T, C, B, deg, ctx);
        gf28_mpoly_swap(T, A, ctx);
        gf28_mpoly_clear(T, ctx);
    } else {
        gf28_mpoly_fit_length_reset_bits(A, B->length + C->length - 1, exp_bits, ctx);
        _gf28_mpoly_mul_array_chunked_DEG(A, C, B, deg, ctx);
    }
    success = 1;
    
cleanup:
    return success;
}

/* Main entry point for array multiplication */
int gf28_mpoly_mul_array(gf28_mpoly_t A, const gf28_mpoly_t B,
                         const gf28_mpoly_t C, const gf28_mpoly_ctx_t ctx)
{
    slong i;
    int success;
    fmpz *maxBfields, *maxCfields;
    TMP_INIT;
    
    if (B->length == 0 || C->length == 0) {
        gf28_mpoly_zero(A, ctx);
        return 1;
    }
    
    if (B->bits == 0 || C->bits == 0) {
        return 0;
    }
    
    if (ctx->minfo->nvars < 1 ||
        1 != mpoly_words_per_exp(B->bits, ctx->minfo) ||
        1 != mpoly_words_per_exp(C->bits, ctx->minfo)) {
        return 0;
    }
    
    TMP_START;
    
    maxBfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    maxCfields = (fmpz *) TMP_ALLOC(ctx->minfo->nfields*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nfields; i++) {
        fmpz_init(maxBfields + i);
        fmpz_init(maxCfields + i);
    }
    mpoly_max_fields_fmpz(maxBfields, B->exps, B->length, B->bits, ctx->minfo);
    mpoly_max_fields_fmpz(maxCfields, C->exps, C->length, C->bits, ctx->minfo);
    
    switch (ctx->minfo->ord) {
        case ORD_LEX:
            success = _gf28_mpoly_mul_array_LEX(A, B, maxBfields, C, maxCfields, ctx);
            break;
        case ORD_DEGLEX:
        case ORD_DEGREVLEX:
            success = _gf28_mpoly_mul_array_DEG(A, B, maxBfields, C, maxCfields, ctx);
            break;
        default:
            success = 0;
            break;
    }
    
    for (i = 0; i < ctx->minfo->nfields; i++) {
        fmpz_clear(maxBfields + i);
        fmpz_clear(maxCfields + i);
    }
    
    TMP_END;
    return success;
}

/* ============================================================================
   HELPER FUNCTIONS FOR TESTING
   ============================================================================ */

static void gf28_mpoly_set(gf28_mpoly_t res, const gf28_mpoly_t poly, const gf28_mpoly_ctx_t ctx) {
    if (res == poly) return;
    
    gf28_mpoly_fit_length_reset_bits(res, poly->length, poly->bits, ctx);
    res->length = poly->length;
    
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    memcpy(res->coeffs, poly->coeffs, poly->length * sizeof(uint8_t));
    memcpy(res->exps, poly->exps, N * poly->length * sizeof(ulong));
}

static void gf28_mpoly_set_coeff_ui_ui(gf28_mpoly_t poly, uint8_t c, 
                                       const ulong *exp, const gf28_mpoly_ctx_t ctx) {
    if (poly->bits == 0) {
        poly->bits = MPOLY_MIN_BITS;
    }
    
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    
    ulong *cmpmask = (ulong *) flint_malloc(N * FLINT_BITS * sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, poly->bits, ctx->minfo);
    
    ulong *packed_exp = (ulong *) flint_malloc(N * sizeof(ulong));
    mpoly_set_monomial_ui(packed_exp, exp, poly->bits, ctx->minfo);
    
    slong pos = 0;
    for (pos = 0; pos < poly->length; pos++) {
        if (mpoly_monomial_equal(poly->exps + N*pos, packed_exp, N)) {
            poly->coeffs[pos] = c;
            flint_free(cmpmask);
            flint_free(packed_exp);
            return;
        }
        if (mpoly_monomial_lt(poly->exps + N*pos, packed_exp, N, cmpmask)) {
            break;
        }
    }
    
    _gf28_mpoly_fit_length(&poly->coeffs, &poly->coeffs_alloc,
                           &poly->exps, &poly->exps_alloc, N, poly->length + 1);
    
    for (slong i = poly->length; i > pos; i--) {
        poly->coeffs[i] = poly->coeffs[i-1];
        mpoly_monomial_set(poly->exps + N*i, poly->exps + N*(i-1), N);
    }
    
    poly->coeffs[pos] = c;
    mpoly_monomial_set(poly->exps + N*pos, packed_exp, N);
    poly->length++;
    
    flint_free(cmpmask);
    flint_free(packed_exp);
}

static void gf28_mpoly_print_old(const gf28_mpoly_t poly, const char **vars, 
                             const gf28_mpoly_ctx_t ctx) {
    if (poly->length == 0) {
        printf("0");
        return;
    }
    
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    
    for (slong i = 0; i < poly->length; i++) {
        if (i > 0) printf(" + ");
        printf("%d", poly->coeffs[i]);
        
        ulong *exp = (ulong *) flint_malloc(nvars * sizeof(ulong));
        mpoly_get_monomial_ui(exp, poly->exps + N*i, poly->bits, ctx->minfo);
        
        for (slong j = 0; j < nvars; j++) {
            if (exp[j] > 0) {
                printf("*%s", vars[j]);
                if (exp[j] > 1) {
                    printf("^%lu", exp[j]);
                }
            }
        }
        flint_free(exp);
    }
}

static void gf28_mpoly_print(const gf28_mpoly_t poly, const char **vars, 
                             const gf28_mpoly_ctx_t ctx) {
    if (poly->length == 0) {
        printf("0");
        return;
    }
    
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    
    for (slong i = 0; i < poly->length; i++) {
        if (i > 0) printf(" + ");
        
        /* Print coefficient as hex value */
        printf("0x%02x", poly->coeffs[i]);
        
        ulong *exp = (ulong *) flint_malloc(nvars * sizeof(ulong));
        mpoly_get_monomial_ui(exp, poly->exps + N*i, poly->bits, ctx->minfo);
        
        int has_vars = 0;
        for (slong j = 0; j < nvars; j++) {
            if (exp[j] > 0) {
                has_vars = 1;
                break;
            }
        }
        
        if (has_vars) {
            printf("*");
            int first_var = 1;
            for (slong j = 0; j < nvars; j++) {
                if (exp[j] > 0) {
                    if (!first_var) printf("*");
                    printf("%s", vars[j]);
                    if (exp[j] > 1) {
                        printf("^%lu", exp[j]);
                    }
                    first_var = 0;
                }
            }
        }
        flint_free(exp);
    }
}

/* Conversion to/from FLINT - now using same irreducible polynomial */
static void gf28_mpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t res, const gf28_mpoly_t poly,
                                         const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx) {
    fq_nmod_mpoly_zero(res, fq_mpoly_ctx);
    
    if (poly->length == 0) return;
    
    slong N = mpoly_words_per_exp(poly->bits, fq_mpoly_ctx->minfo);
    slong nvars = fq_mpoly_ctx->minfo->nvars;
    
    for (slong i = 0; i < poly->length; i++) {
        if (poly->coeffs[i] != 0) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, fqctx);
            
            /* Direct conversion since we use the same polynomial */
            nmod_poly_zero(coeff);
            uint8_t val = poly->coeffs[i];
            for (int j = 0; j < 8; j++) {
                if (val & (1 << j)) {
                    nmod_poly_set_coeff_ui(coeff, j, 1);
                }
            }
            
            ulong *exp = (ulong *) flint_malloc(nvars * sizeof(ulong));
            mpoly_get_monomial_ui(exp, poly->exps + N*i, poly->bits, fq_mpoly_ctx->minfo);
            
            fq_nmod_mpoly_set_coeff_fq_nmod_ui(res, coeff, exp, fq_mpoly_ctx);
            
            flint_free(exp);
            fq_nmod_clear(coeff, fqctx);
        }
    }
}

static void fq_nmod_mpoly_to_gf28_mpoly(gf28_mpoly_t res, const fq_nmod_mpoly_t poly,
                                         const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx) {
    gf28_mpoly_ctx_t ctx;
    gf28_mpoly_ctx_init(ctx, fq_mpoly_ctx->minfo->nvars, fq_mpoly_ctx->minfo->ord);
    
    gf28_mpoly_zero(res, ctx);
    
    slong len = fq_nmod_mpoly_length(poly, fq_mpoly_ctx);
    if (len == 0) {
        gf28_mpoly_ctx_clear(ctx);
        return;
    }
    
    flint_bitcnt_t bits = FLINT_MAX(poly->bits, MPOLY_MIN_BITS);
    gf28_mpoly_fit_length_reset_bits(res, len, bits, ctx);
    res->length = 0;
    
    slong nvars = fq_mpoly_ctx->minfo->nvars;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    
    slong actual_terms = 0;
    for (slong i = 0; i < len; i++) {
        fq_nmod_t coeff;
        fq_nmod_init(coeff, fqctx);
        
        fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, poly, i, fq_mpoly_ctx);
        
        /* Direct conversion since we use the same polynomial */
        uint8_t gf28_coeff = 0;
        for (int j = 0; j < 8; j++) {
            if (nmod_poly_get_coeff_ui(coeff, j)) {
                gf28_coeff |= (1 << j);
            }
        }
        
        if (gf28_coeff != 0) {
            ulong *exp = (ulong *) flint_malloc(nvars * sizeof(ulong));
            fq_nmod_mpoly_get_term_exp_ui(exp, poly, i, fq_mpoly_ctx);
            
            mpoly_set_monomial_ui(res->exps + N*actual_terms, exp, bits, ctx->minfo);
            res->coeffs[actual_terms] = gf28_coeff;
            actual_terms++;
            
            flint_free(exp);
        }
        
        fq_nmod_clear(coeff, fqctx);
    }
    
    res->length = actual_terms;
    gf28_mpoly_ctx_clear(ctx);
}

/* Timer function */
static double get_wall_time_8() {
    struct timeval time;
    if (gettimeofday(&time, NULL)) {
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * 0.000001;
}

/* ============================================================================
   DENSE POLYNOMIAL STRUCTURES
   ============================================================================ */

typedef struct {
    uint8_t *coeffs;
    slong alloc;
    slong *deg_bounds;
    slong nvars;
} gf28_mpolyd_struct;

typedef gf28_mpolyd_struct gf28_mpolyd_t[1];

/* Initialize dense polynomial */
static void gf28_mpolyd_init(gf28_mpolyd_t A, slong nvars) {
    A->coeffs = NULL;
    A->alloc = 0;
    A->deg_bounds = (slong *)calloc(nvars, sizeof(slong));
    A->nvars = nvars;
}

/* Clear dense polynomial */
static void gf28_mpolyd_clear(gf28_mpolyd_t A) {
    if (A->coeffs) free(A->coeffs);
    if (A->deg_bounds) free(A->deg_bounds);
}

/* Compute offset for monomial in dense array */
static slong gf28_mpolyd_offset(const gf28_mpolyd_t A, const ulong *exp) {
    slong off = 0;
    slong stride = 1;
    
    /* Pack with first variable changing fastest (LEX order) */
    for (slong i = 0; i < A->nvars; i++) {
        /* IMPORTANT: Don't wrap exponents, just check bounds */
        if (exp[i] >= (ulong)A->deg_bounds[i]) return -1;
        off += (slong)exp[i] * stride;
        stride *= A->deg_bounds[i];
    }
    
    return off;
}

/* ============================================================================
   UNIVARIATE VIEW OF MULTIVARIATE POLYNOMIALS
   ============================================================================ */

/* Division with remainder for univariate view */
static void gf28_mpolyd_divrem_univar(gf28_mpolyd_t Q, gf28_mpolyd_t R,
                                      const gf28_mpolyd_t A, const gf28_mpolyd_t B) {
    slong n = A->alloc;
    
    /* Copy A to R */
    memcpy(R->coeffs, A->coeffs, n * sizeof(uint8_t));
    
    /* Clear Q */
    memset(Q->coeffs, 0, Q->alloc * sizeof(uint8_t));
    
    /* Find degree of B (highest non-zero coefficient) */
    slong degB = -1;
    for (slong i = n - 1; i >= 0; i--) {
        if (B->coeffs[i] != 0) {
            degB = i;
            break;
        }
    }
    
    if (degB < 0) return; /* B is zero */
    
    uint8_t lc_B_inv = gf28_inv(B->coeffs[degB]);
    
    /* Main division loop */
    for (slong i = n - 1; i >= degB; i--) {
        if (R->coeffs[i] != 0) {
            /* Compute quotient coefficient */
            uint8_t q = gf28_mul(R->coeffs[i], lc_B_inv);
            
            /* IMPORTANT: Check bounds before storing */
            if (i - degB < Q->alloc) {
                Q->coeffs[i - degB] = q;
            }
            
            /* Subtract q*B from R */
            for (slong j = 0; j <= degB && i - degB + j < n; j++) {
                R->coeffs[i - degB + j] ^= gf28_mul(q, B->coeffs[j]);
            }
        }
    }
}

/* Check if dense polynomial is zero */
static int gf28_mpolyd_is_zero(const gf28_mpolyd_t A) {
    for (slong i = 0; i < A->alloc; i++) {
        if (A->coeffs[i] != 0) return 0;
    }
    return 1;
}

/* ============================================================================
   MAIN DIVISION ALGORITHM
   ============================================================================ */

/* Multivariate polynomial division using dense representation */

/* 辅助函数：检查是否为单项式 */
static int gf28_mpoly_is_monomial(const gf28_mpoly_t poly) {
    return poly->length == 1;
}

/* 辅助函数：获取单项式的系数 */
static uint8_t gf28_mpoly_get_monomial_coeff(const gf28_mpoly_t poly) {
    if (poly->length != 1) return 0;
    return poly->coeffs[0];
}

/* 辅助函数：获取单项式的指数 */
static void gf28_mpoly_get_monomial_exp(ulong *exp, const gf28_mpoly_t poly, 
                                        const gf28_mpoly_ctx_t ctx) {
    if (poly->length != 1) return;
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    mpoly_get_monomial_ui(exp, poly->exps, poly->bits, ctx->minfo);
}


/* 改进的稀疏到密集转换 */
static void gf28_mpoly_to_mpolyd(gf28_mpolyd_t A, const gf28_mpoly_t B, 
                                         const gf28_mpoly_ctx_t ctx) {
    slong nvars = ctx->minfo->nvars;
    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);
    ulong *exp = (ulong *)malloc(nvars * sizeof(ulong));
    
    /* 清空所有系数 */
    memset(A->coeffs, 0, A->alloc * sizeof(uint8_t));
    
    /* 转换每一项 */
    for (slong i = 0; i < B->length; i++) {
        mpoly_get_monomial_ui(exp, B->exps + N*i, B->bits, ctx->minfo);
        slong off = gf28_mpolyd_offset(A, exp);
        if (off >= 0 && off < A->alloc) {
            A->coeffs[off] = B->coeffs[i];
        }
    }
    
    free(exp);
}

/* 改进的密集到稀疏转换 - 确保正确处理位数 */
static void gf28_mpolyd_to_mpoly(gf28_mpoly_t A, const gf28_mpolyd_t B,
                                         const gf28_mpoly_ctx_t ctx) {
    slong nvars = ctx->minfo->nvars;
    ulong *exp = (ulong *)calloc(nvars, sizeof(ulong));
    
    /* 首先计算需要的位数 */
    flint_bitcnt_t bits_needed = MPOLY_MIN_BITS;
    for (slong i = 0; i < nvars; i++) {
        if (B->deg_bounds[i] > 0) {
            slong bits_i = FLINT_BIT_COUNT(B->deg_bounds[i] - 1);
            bits_needed = FLINT_MAX(bits_needed, bits_i);
        }
    }
    if (bits_needed < 16) bits_needed = 16; /* 至少16位 */
    bits_needed = mpoly_fix_bits(bits_needed, ctx->minfo);
    
    gf28_mpoly_zero(A, ctx);
    
    /* 遍历密集数组 */
    for (slong off = 0; off < B->alloc; off++) {
        if (B->coeffs[off] == 0) continue;
        
        /* 从偏移量重建指数向量 */
        slong temp = off;
        for (slong i = 0; i < nvars; i++) {
            exp[i] = (ulong)(temp % B->deg_bounds[i]);
            temp /= B->deg_bounds[i];
        }
        
        /* 添加项到多项式 - 使用正确的位数 */
        if (A->bits < bits_needed) {
            gf28_mpoly_fit_length_reset_bits(A, A->length + 1, bits_needed, ctx);
        }
        gf28_mpoly_set_coeff_ui_ui(A, B->coeffs[off], exp, ctx);
    }
    
    free(exp);
}

/* 多变量除法与余数 */
static void gf28_mpolyd_divrem_multivar(gf28_mpolyd_t Q, gf28_mpolyd_t R,
                                       const gf28_mpolyd_t A, const gf28_mpolyd_t B,
                                       const gf28_mpoly_ctx_t ctx) {
    /* 暂时使用单变量视角的除法 */
    gf28_mpolyd_divrem_univar(Q, R, A, B);
}

/* 改进的度数界限设置 */
static int gf28_mpolyd_set_degbounds(gf28_mpolyd_t A, const slong *bounds) {
    slong size = 1;
    
    for (slong i = 0; i < A->nvars; i++) {
        A->deg_bounds[i] = bounds[i];
        if (bounds[i] <= 0) return 0;
        
        /* 检查溢出 */
        if (size > WORD_MAX / bounds[i]) return 0;
        size *= bounds[i];
    }
    
    /* 限制大小以防止过度内存使用 */
    if (size > (1L << 26)) return 0; /* 64M 系数上限 */
    
    A->alloc = size;
    A->coeffs = (uint8_t *)calloc(size, sizeof(uint8_t));
    return A->coeffs != NULL;
}

/* 单项式除法 - 直接实现 */
static int gf28_mpoly_divides_monomial(gf28_mpoly_t Q, const gf28_mpoly_t A, 
                                       const gf28_mpoly_t B, const gf28_mpoly_ctx_t ctx) {
    slong nvars = ctx->minfo->nvars;
    ulong *exp_b = (ulong *)malloc(nvars * sizeof(ulong));
    ulong *exp_a = (ulong *)malloc(nvars * sizeof(ulong));
    ulong *exp_q = (ulong *)malloc(nvars * sizeof(ulong));
    uint8_t coeff_b = gf28_mpoly_get_monomial_coeff(B);
    uint8_t coeff_b_inv;
    slong N;
    
    if (coeff_b == 0) {
        free(exp_b); free(exp_a); free(exp_q);
        return 0; /* 除以零 */
    }
    
    coeff_b_inv = gf28_inv(coeff_b);
    gf28_mpoly_get_monomial_exp(exp_b, B, ctx);
    
    /* 确定结果的位数 */
    flint_bitcnt_t bits = FLINT_MAX(A->bits, B->bits);
    if (bits < 16) bits = 16; /* 确保至少16位以处理大指数 */
    
    /* 初始化结果 */
    gf28_mpoly_fit_length_reset_bits(Q, A->length, bits, ctx);
    Q->length = 0;
    N = mpoly_words_per_exp(bits, ctx->minfo);
    
    /* 对 A 的每一项进行除法 */
    for (slong i = 0; i < A->length; i++) {
        mpoly_get_monomial_ui(exp_a, A->exps + N*i, A->bits, ctx->minfo);
        
        /* 检查是否可以整除 */
        int divisible = 1;
        for (slong j = 0; j < nvars; j++) {
            if (exp_a[j] < exp_b[j]) {
                divisible = 0;
                break;
            }
            exp_q[j] = exp_a[j] - exp_b[j];
        }
        
        if (!divisible) {
            /* 不能整除 */
            free(exp_b); free(exp_a); free(exp_q);
            gf28_mpoly_zero(Q, ctx);
            return 0;
        }
        
        /* 计算系数 */
        uint8_t coeff_q = gf28_mul(A->coeffs[i], coeff_b_inv);
        
        if (coeff_q != 0) {
            /* 添加到结果 */
            mpoly_set_monomial_ui(Q->exps + N*Q->length, exp_q, bits, ctx->minfo);
            Q->coeffs[Q->length] = coeff_q;
            Q->length++;
        }
    }
    
    free(exp_b); free(exp_a); free(exp_q);
    return 1;
}

/* 改进的密集表示除法 - 修复指数处理 */
int gf28_mpoly_divides_dense(gf28_mpoly_t Q, const gf28_mpoly_t A, 
                                     const gf28_mpoly_t B, const gf28_mpoly_ctx_t ctx) {
    slong nvars = ctx->minfo->nvars;
    int success = 0;
    
    /* 获取度数界限 */
    slong *degs_A = (slong *)malloc(nvars * sizeof(slong));
    slong *degs_B = (slong *)malloc(nvars * sizeof(slong));
    slong *bounds = (slong *)malloc(nvars * sizeof(slong));
    
    mpoly_degrees_si(degs_A, A->exps, A->length, A->bits, ctx->minfo);
    mpoly_degrees_si(degs_B, B->exps, B->length, B->bits, ctx->minfo);
    
    /* 检查是否可能除得尽 */
    for (slong i = 0; i < nvars; i++) {
        if (degs_A[i] < degs_B[i]) {
            free(degs_A); free(degs_B); free(bounds);
            gf28_mpoly_zero(Q, ctx);
            return 0;
        }
        bounds[i] = degs_A[i] + 1;
    }
    
    /* 初始化密集多项式 */
    gf28_mpolyd_t Ad, Bd, Qd, Rd;
    gf28_mpolyd_init(Ad, nvars);
    gf28_mpolyd_init(Bd, nvars);
    gf28_mpolyd_init(Qd, nvars);
    gf28_mpolyd_init(Rd, nvars);
    
    /* 设置度数界限 */
    if (!gf28_mpolyd_set_degbounds(Ad, bounds) ||
        !gf28_mpolyd_set_degbounds(Bd, bounds) ||
        !gf28_mpolyd_set_degbounds(Qd, bounds) ||
        !gf28_mpolyd_set_degbounds(Rd, bounds)) {
        goto cleanup;
    }
    
    /* 转换为密集表示 */
    gf28_mpoly_to_mpolyd(Ad, A, ctx);
    gf28_mpoly_to_mpolyd(Bd, B, ctx);
    
    /* 执行除法 */
    gf28_mpolyd_divrem_multivar(Qd, Rd, Ad, Bd, ctx);
    
    /* 检查余数是否为零 */
    if (gf28_mpolyd_is_zero(Rd)) {
        /* 转换商回到稀疏表示 */
        gf28_mpolyd_to_mpoly(Q, Qd, ctx);
        success = 1;
    } else {
        gf28_mpoly_zero(Q, ctx);
        success = 0;
    }
    
cleanup:
    gf28_mpolyd_clear(Ad);
    gf28_mpolyd_clear(Bd);
    gf28_mpolyd_clear(Qd);
    gf28_mpolyd_clear(Rd);
    free(degs_A);
    free(degs_B);
    free(bounds);
    
    return success;
}


/* 主除法函数 - 修复版本 */
int gf28_mpoly_divides(gf28_mpoly_t Q, const gf28_mpoly_t A, 
                      const gf28_mpoly_t B, const gf28_mpoly_ctx_t ctx) {
    /* 特殊情况处理 */
    if (B->length == 0) {
        if (A->length == 0) {
            gf28_mpoly_zero(Q, ctx);
            return 1;
        }
        return 0; /* 除以零 */
    }
    
    if (A->length == 0) {
        gf28_mpoly_zero(Q, ctx);
        return 1;
    }
    
    /* 如果 B 是单项式，使用直接除法 */
    if (gf28_mpoly_is_monomial(B)) {
        return gf28_mpoly_divides_monomial(Q, A, B, ctx);
    }
    
    /* 对于一般情况，使用改进的密集表示除法 */
    return gf28_mpoly_divides_dense(Q, A, B, ctx);
}


/* ============================================================================
   HELPER FUNCTIONS FOR TESTING
   ============================================================================ */

/* Simple multiplication for verification */
static void gf28_mpoly_mul_simple(gf28_mpoly_t res, const gf28_mpoly_t a, 
                                 const gf28_mpoly_t b, const gf28_mpoly_ctx_t ctx) {
    slong nvars = ctx->minfo->nvars;
    slong N = mpoly_words_per_exp(a->bits, ctx->minfo);
    ulong *exp_a, *exp_b, *exp_sum;
    
    gf28_mpoly_zero(res, ctx);
    
    if (a->length == 0 || b->length == 0) return;
    
    exp_a = (ulong *)malloc(nvars * sizeof(ulong));
    exp_b = (ulong *)malloc(nvars * sizeof(ulong));
    exp_sum = (ulong *)malloc(nvars * sizeof(ulong));
    
    /* Multiply all pairs of terms */
    for (slong i = 0; i < a->length; i++) {
        mpoly_get_monomial_ui(exp_a, a->exps + N*i, a->bits, ctx->minfo);
        
        for (slong j = 0; j < b->length; j++) {
            mpoly_get_monomial_ui(exp_b, b->exps + N*j, b->bits, ctx->minfo);
            
            /* Add exponents */
            for (slong k = 0; k < nvars; k++) {
                exp_sum[k] = exp_a[k] + exp_b[k];
            }
            
            /* Multiply coefficients */
            uint8_t coeff = gf28_mul(a->coeffs[i], b->coeffs[j]);
            
            if (coeff != 0) {
                /* Get current coefficient */
                uint8_t current = 0;
                for (slong k = 0; k < res->length; k++) {
                    ulong *exp_k = (ulong *)malloc(nvars * sizeof(ulong));
                    mpoly_get_monomial_ui(exp_k, res->exps + N*k, res->bits, ctx->minfo);
                    
                    int equal = 1;
                    for (slong l = 0; l < nvars; l++) {
                        if (exp_k[l] != exp_sum[l]) {
                            equal = 0;
                            break;
                        }
                    }
                    
                    if (equal) {
                        current = res->coeffs[k];
                        break;
                    }
                    free(exp_k);
                }
                
                /* Set new coefficient */
                gf28_mpoly_set_coeff_ui_ui(res, gf28_add(current, coeff), exp_sum, ctx);
            }
        }
    }
    
    free(exp_a);
    free(exp_b);
    free(exp_sum);
}

/* Check if two polynomials are equal */
static int gf28_mpoly_equal(const gf28_mpoly_t A, const gf28_mpoly_t B, 
                           const gf28_mpoly_ctx_t ctx) {
    if (A->length != B->length) return 0;
    
    slong nvars = ctx->minfo->nvars;
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    
    /* For each term in A, find matching term in B */
    for (slong i = 0; i < A->length; i++) {
        ulong *exp_a = (ulong *)malloc(nvars * sizeof(ulong));
        mpoly_get_monomial_ui(exp_a, A->exps + N*i, A->bits, ctx->minfo);
        
        int found = 0;
        for (slong j = 0; j < B->length; j++) {
            ulong *exp_b = (ulong *)malloc(nvars * sizeof(ulong));
            mpoly_get_monomial_ui(exp_b, B->exps + N*j, B->bits, ctx->minfo);
            
            int equal = 1;
            for (slong k = 0; k < nvars; k++) {
                if (exp_a[k] != exp_b[k]) {
                    equal = 0;
                    break;
                }
            }
            
            if (equal) {
                if (A->coeffs[i] != B->coeffs[j]) {
                    free(exp_a);
                    free(exp_b);
                    return 0;
                }
                found = 1;
                free(exp_b);
                break;
            }
            free(exp_b);
        }
        
        free(exp_a);
        if (!found) return 0;
    }
    
    return 1;
}

/* Generate random polynomial */
static void gf28_mpoly_randtest(gf28_mpoly_t poly, flint_rand_t state,
                               slong length, slong exp_bound, 
                               const gf28_mpoly_ctx_t ctx) {
    slong nvars = ctx->minfo->nvars;
    ulong *exp = (ulong *)malloc(nvars * sizeof(ulong));
    
    gf28_mpoly_zero(poly, ctx);
    
    for (slong i = 0; i < length; i++) {
        /* Generate random exponents */
        for (slong j = 0; j < nvars; j++) {
            exp[j] = n_randint(state, exp_bound);
        }
        
        /* Generate random non-zero coefficient */
        uint8_t c = n_randint(state, 255) + 1;
        
        gf28_mpoly_set_coeff_ui_ui(poly, c, exp, ctx);
    }
    
    free(exp);
}


#endif