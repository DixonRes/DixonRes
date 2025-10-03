/* gf2128_mpoly.c - GF(2^128) Multivariate Polynomial Implementation */

#include "gf2128_mpoly.h"

/* ============================================================================
   BASIC INITIALIZATION AND MEMORY MANAGEMENT
   ============================================================================ */

void gf2128_mpoly_init(gf2128_mpoly_t poly, const gf2128_mpoly_ctx_t ctx) {
    poly->coeffs = NULL;
    poly->exps = NULL;
    poly->length = 0;
    poly->coeffs_alloc = 0;
    poly->exps_alloc = 0;
    poly->bits = MPOLY_MIN_BITS;
}

void gf2128_mpoly_clear(gf2128_mpoly_t poly, const gf2128_mpoly_ctx_t ctx) {
    if (poly->coeffs) flint_free(poly->coeffs);
    if (poly->exps) flint_free(poly->exps);
}

void gf2128_mpoly_ctx_init(gf2128_mpoly_ctx_t ctx, slong nvars, const ordering_t ord) {
    mpoly_ctx_init(ctx->minfo, nvars, ord);
}

void gf2128_mpoly_ctx_clear(gf2128_mpoly_ctx_t ctx) {
    mpoly_ctx_clear(ctx->minfo);
}

void gf2128_mpoly_zero(gf2128_mpoly_t poly, const gf2128_mpoly_ctx_t ctx) {
    poly->length = 0;
}

void _gf2128_mpoly_fit_length(gf2128_t **coeffs, slong *coeffs_alloc,
                             ulong **exps, slong *exps_alloc, slong N, slong length) {
    if (length > *coeffs_alloc) {
        slong new_alloc = FLINT_MAX(length, 2 * (*coeffs_alloc));
        *coeffs = (gf2128_t *) flint_realloc(*coeffs, new_alloc * sizeof(gf2128_t));
        *coeffs_alloc = new_alloc;
    }
    
    if (N*length > *exps_alloc) {
        slong new_alloc = FLINT_MAX(N*length, 2 * (*exps_alloc));
        *exps = (ulong *) flint_realloc(*exps, new_alloc * sizeof(ulong));
        *exps_alloc = new_alloc;
    }
}

void gf2128_mpoly_fit_length_reset_bits(gf2128_mpoly_t poly, slong len, 
                                       flint_bitcnt_t bits, const gf2128_mpoly_ctx_t ctx) {
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    _gf2128_mpoly_fit_length(&poly->coeffs, &poly->coeffs_alloc,
                             &poly->exps, &poly->exps_alloc, N, len);
    poly->bits = bits;
}

void gf2128_mpoly_init3(gf2128_mpoly_t poly, slong alloc, flint_bitcnt_t bits,
                       const gf2128_mpoly_ctx_t ctx) {
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    poly->coeffs = (gf2128_t *) flint_malloc(alloc * sizeof(gf2128_t));
    poly->exps = (ulong *) flint_malloc(N * alloc * sizeof(ulong));
    poly->coeffs_alloc = alloc;
    poly->exps_alloc = N * alloc;
    poly->length = 0;
    poly->bits = bits;
}

void gf2128_mpoly_swap(gf2128_mpoly_t poly1, gf2128_mpoly_t poly2, const gf2128_mpoly_ctx_t ctx) {
    gf2128_mpoly_struct t = *poly1;
    *poly1 = *poly2;
    *poly2 = t;
}

void _gf2128_mpoly_set_length(gf2128_mpoly_t poly, slong len, const gf2128_mpoly_ctx_t ctx) {
    poly->length = len;
}

/* ============================================================================
   DYNAMIC MEMORY MANAGEMENT FOR LARGE ARRAYS
   ============================================================================ */

void gf2128_dynamic_array_init(gf2128_dynamic_array_t *arr, slong size) {
    arr->size = size;
    arr->chunk_size = 1L << 16;  // 64K elements chunks (1MB per chunk)
    arr->use_chunks = (size > (1L << 22));  // Use chunks for arrays > 64M elements
    
    if (arr->use_chunks) {
        arr->data = NULL;  // Will allocate on demand
    } else {
        arr->data = (gf2128_t *)calloc(size, sizeof(gf2128_t));
        if (!arr->data && size > 0) {
            flint_printf("Failed to allocate %ld GF(2^128) elements\n", size);
            flint_abort();
        }
    }
}

void gf2128_dynamic_array_clear(gf2128_dynamic_array_t *arr) {
    if (arr->data) {
        free(arr->data);
        arr->data = NULL;
    }
}

/* ============================================================================
   ARRAY MULTIPLICATION OPERATIONS
   ============================================================================ */

/* Array multiplication with bounds checking */
void _gf2128_mpoly_addmul_array1_safe(gf2128_t *poly1, slong array_size,
                                      const gf2128_t *poly2, const ulong *exp2, slong len2,
                                      const gf2128_t *poly3, const ulong *exp3, slong len3) {
    slong ii, i, jj, j;
    gf2128_t *c2;
    
    for (ii = 0; ii < len2 + BLOCK; ii += BLOCK) {
        for (jj = 0; jj < len3 + BLOCK; jj += BLOCK) {
            for (i = ii; i < FLINT_MIN(ii + BLOCK, len2); i++) {
                slong offset2 = (slong)exp2[i];
                
                if (offset2 >= array_size) {
                    flint_printf("Array bounds error: offset %ld >= size %ld\n", offset2, array_size);
                    continue;
                }
                
                c2 = poly1 + offset2;
                
                if (!gf2128_is_zero(&poly2[i])) {
                    for (j = jj; j < FLINT_MIN(jj + BLOCK, len3); j++) {
                        slong offset3 = (slong)exp3[j];
                        slong total_offset = offset2 + offset3;
                        
                        if (total_offset >= array_size) {
                            flint_printf("Array bounds error: total offset %ld >= size %ld\n", 
                                       total_offset, array_size);
                            continue;
                        }
                        
                        gf2128_t prod = gf2128_mul(&poly2[i], &poly3[j]);
                        c2[offset3] = gf2128_add(&c2[offset3], &prod);
                    }
                }
            }
        }
    }
}

/* LEX ordering unpacking */
slong gf2128_mpoly_append_array_LEX_safe(gf2128_mpoly_t P, slong Plen, gf2128_t *coeff_array,
                                         const ulong *mults, slong num, slong array_size, 
                                         slong top, const gf2128_mpoly_ctx_t ctx) {
    slong off, j;
    slong topmult = num == 0 ? 1 : mults[num - 1];
    slong lastd = topmult - 1;
    slong reset = array_size/topmult;
    slong counter = reset;
    ulong startexp = ((ulong)top << (P->bits*num)) + ((ulong)lastd << (P->bits*(num-1)));
    gf2128_t coeff;
    
    for (off = array_size - 1; off >= 0; off--) {
        if (!gf2128_is_zero(&coeff_array[off])) {
            coeff = coeff_array[off];
            slong d = off;
            ulong exp = startexp;
            
            for (j = 0; j + 1 < num; j++) {
                ulong exp_j = d % mults[j];
                exp += exp_j << (P->bits*j);
                d = d / mults[j];
            }
            
            _gf2128_mpoly_fit_length(&P->coeffs, &P->coeffs_alloc,
                                     &P->exps, &P->exps_alloc, 1, Plen + 1);
            P->exps[Plen] = exp;
            P->coeffs[Plen] = coeff;
            Plen++;
            coeff_array[off] = gf2128_zero();
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

/* DEGLEX ordering unpacking */
slong gf2128_mpoly_append_array_DEGLEX_safe(gf2128_mpoly_t P, slong Plen, gf2128_t *coeff_array,
                                            slong top, slong nvars, slong degb,
                                            const gf2128_mpoly_ctx_t ctx) {
    slong i;
    ulong exp, lomask = (UWORD(1) << (P->bits - 1)) - 1;
    slong off, array_size;
    slong *curexp, *degpow;
    ulong *oneexp;
    gf2128_t coeff;
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
        if (off >= 0 && off < array_size && !gf2128_is_zero(&coeff_array[off])) {
            coeff = coeff_array[off];
            _gf2128_mpoly_fit_length(&P->coeffs, &P->coeffs_alloc,
                                     &P->exps, &P->exps_alloc, 1, Plen + 1);
            P->exps[Plen] = exp;
            P->coeffs[Plen] = coeff;
            Plen++;
            coeff_array[off] = gf2128_zero();
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

/* DEGREVLEX ordering unpacking */
slong gf2128_mpoly_append_array_DEGREVLEX_safe(gf2128_mpoly_t P, slong Plen, gf2128_t *coeff_array,
                                               slong top, slong nvars, slong degb,
                                               const gf2128_mpoly_ctx_t ctx) {
    slong i;
    ulong exp, mask = UWORD(1) << (P->bits - 1);
    slong off, array_size;
    slong *curexp, *degpow;
    ulong *oneexp;
    gf2128_t coeff;
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
        if (off >= 0 && off < array_size && !gf2128_is_zero(&coeff_array[off])) {
            coeff = coeff_array[off];
            _gf2128_mpoly_fit_length(&P->coeffs, &P->coeffs_alloc,
                                     &P->exps, &P->exps_alloc, 1, Plen + 1);
            P->exps[Plen] = exp;
            P->coeffs[Plen] = coeff;
            Plen++;
            coeff_array[off] = gf2128_zero();
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
void _gf2128_mpoly_mul_array_chunked_LEX(gf2128_mpoly_t P,
                                         const gf2128_mpoly_t A,
                                         const gf2128_mpoly_t B,
                                         const ulong *mults,
                                         const gf2128_mpoly_ctx_t ctx)
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
    
    gf2128_dynamic_array_t coeff_array;
    gf2128_dynamic_array_init(&coeff_array, array_size);
    
    if (!coeff_array.data) {
        flint_printf("Memory allocation failed for array of size %ld\n", array_size);
        flint_free(Apexp);
        flint_free(Bpexp);
        TMP_END;
        return;
    }
    
    for (Pi = 0; Pi < Pl; Pi++) {
        // Clear array
        for (slong k = 0; k < array_size; k++) {
            coeff_array.data[k] = gf2128_zero();
        }
        
        for (i = 0, j = Pi; i < Al && j >= 0; i++, j--) {
            if (j < Bl) {
                _gf2128_mpoly_addmul_array1_safe(coeff_array.data, array_size,
                        A->coeffs + Amain[i],
                        Apexp + Amain[i], Amain[i + 1] - Amain[i],
                        B->coeffs + Bmain[j],
                        Bpexp + Bmain[j], Bmain[j + 1] - Bmain[j]);
            }
        }
        
        Plen = gf2128_mpoly_append_array_LEX_safe(P, Plen, coeff_array.data,
                              mults, num, array_size, Pl - Pi - 1, ctx);
    }
    
    _gf2128_mpoly_set_length(P, Plen, ctx);
    
    gf2128_dynamic_array_clear(&coeff_array);
    flint_free(Apexp);
    flint_free(Bpexp);
    TMP_END;
}

/* DEGLEX/DEGREVLEX multiplication with dynamic allocation */
void _gf2128_mpoly_mul_array_chunked_DEG(gf2128_mpoly_t P,
                                         const gf2128_mpoly_t A,
                                         const gf2128_mpoly_t B,
                                         ulong degb,
                                         const gf2128_mpoly_ctx_t ctx)
{
    slong nvars = ctx->minfo->nvars;
    slong Pi, i, j, Plen, Pl, Al, Bl, array_size;
    slong *Amain, *Bmain;
    ulong *Apexp, *Bpexp;
    slong (*upack)(gf2128_mpoly_t, slong, gf2128_t *, slong, slong, slong, const gf2128_mpoly_ctx_t);
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
    
    upack = &gf2128_mpoly_append_array_DEGLEX_safe;
    if (ctx->minfo->ord == ORD_DEGREVLEX) {
        upack = &gf2128_mpoly_append_array_DEGREVLEX_safe;
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
    
    gf2128_dynamic_array_t coeff_array;
    gf2128_dynamic_array_init(&coeff_array, array_size);
    
    if (!coeff_array.data) {
        flint_printf("Memory allocation failed for array of size %ld\n", array_size);
        flint_free(Apexp);
        flint_free(Bpexp);
        TMP_END;
        return;
    }
    
    for (Pi = 0; Pi < Pl; Pi++) {
        // Clear array
        for (slong k = 0; k < array_size; k++) {
            coeff_array.data[k] = gf2128_zero();
        }
        
        for (i = 0, j = Pi; i < Al && j >= 0; i++, j--) {
            if (j < Bl) {
                _gf2128_mpoly_addmul_array1_safe(coeff_array.data, array_size,
                        A->coeffs + Amain[i],
                        Apexp + Amain[i], Amain[i + 1] - Amain[i],
                        B->coeffs + Bmain[j],
                        Bpexp + Bmain[j], Bmain[j + 1] - Bmain[j]);
            }
        }
        
        Plen = upack(P, Plen, coeff_array.data, Pl - Pi - 1, nvars, degb, ctx);
    }
    
    _gf2128_mpoly_set_length(P, Plen, ctx);
    
    gf2128_dynamic_array_clear(&coeff_array);
    flint_free(Apexp);
    flint_free(Bpexp);
    TMP_END;
}

/* Main LEX multiplication function */
int _gf2128_mpoly_mul_array_LEX(gf2128_mpoly_t A,
                               const gf2128_mpoly_t B,
                               fmpz *maxBfields,
                               const gf2128_mpoly_t C,
                               fmpz *maxCfields,
                               const gf2128_mpoly_ctx_t ctx)
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
    
    i = ctx->minfo->nfields - 1;
    FLINT_ASSERT(fmpz_fits_si(maxBfields + i));
    FLINT_ASSERT(fmpz_fits_si(maxCfields + i));
    
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
    
    /* Check if array size is reasonable (limit for GF(2^128) is smaller) */
    if (array_size > (1L << 26)) {  // 64M elements
        flint_printf("Warning: Large array size %ld elements requested\n", array_size);
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
        gf2128_mpoly_t T;
        gf2128_mpoly_init3(T, B->length + C->length - 1, exp_bits, ctx);
        _gf2128_mpoly_mul_array_chunked_LEX(T, C, B, mults, ctx);
        gf2128_mpoly_swap(T, A, ctx);
        gf2128_mpoly_clear(T, ctx);
    } else {
        gf2128_mpoly_fit_length_reset_bits(A, B->length + C->length - 1, exp_bits, ctx);
        _gf2128_mpoly_mul_array_chunked_LEX(A, C, B, mults, ctx);
    }
    success = 1;
    
cleanup:
    TMP_END;
    return success;
}

/* Main DEGLEX/DEGREVLEX multiplication function */
int _gf2128_mpoly_mul_array_DEG(gf2128_mpoly_t A,
                               const gf2128_mpoly_t B,
                               fmpz *maxBfields,
                               const gf2128_mpoly_t C,
                               fmpz *maxCfields,
                               const gf2128_mpoly_ctx_t ctx)
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
    
    array_size = WORD(1);
    for (i = ctx->minfo->nvars - 2; i >= 0; i--) {
        ulong hi;
        umul_ppmm(hi, array_size, array_size, deg);
        if (hi != WORD(0) || array_size <= 0) {
            success = 0;
            goto cleanup;
        }
    }
    
    /* Check if array size is reasonable */
    if (array_size > (1L << 26)) {  // 64M elements
        flint_printf("Warning: Large array size %ld elements requested for DEG ordering\n", array_size);
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
        gf2128_mpoly_t T;
        gf2128_mpoly_init3(T, B->length + C->length - 1, exp_bits, ctx);
        _gf2128_mpoly_mul_array_chunked_DEG(T, C, B, deg, ctx);
        gf2128_mpoly_swap(T, A, ctx);
        gf2128_mpoly_clear(T, ctx);
    } else {
        gf2128_mpoly_fit_length_reset_bits(A, B->length + C->length - 1, exp_bits, ctx);
        _gf2128_mpoly_mul_array_chunked_DEG(A, C, B, deg, ctx);
    }
    success = 1;
    
cleanup:
    return success;
}

/* Main entry point for array multiplication */
int gf2128_mpoly_mul_array(gf2128_mpoly_t A, const gf2128_mpoly_t B,
                          const gf2128_mpoly_t C, const gf2128_mpoly_ctx_t ctx)
{
    slong i;
    int success;
    fmpz *maxBfields, *maxCfields;
    TMP_INIT;
    
    if (B->length == 0 || C->length == 0) {
        gf2128_mpoly_zero(A, ctx);
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
            success = _gf2128_mpoly_mul_array_LEX(A, B, maxBfields, C, maxCfields, ctx);
            break;
        case ORD_DEGLEX:
        case ORD_DEGREVLEX:
            success = _gf2128_mpoly_mul_array_DEG(A, B, maxBfields, C, maxCfields, ctx);
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
   HELPER FUNCTIONS
   ============================================================================ */

void gf2128_mpoly_set(gf2128_mpoly_t res, const gf2128_mpoly_t poly, const gf2128_mpoly_ctx_t ctx) {
    if (res == poly) return;
    
    gf2128_mpoly_fit_length_reset_bits(res, poly->length, poly->bits, ctx);
    res->length = poly->length;
    
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    memcpy(res->coeffs, poly->coeffs, poly->length * sizeof(gf2128_t));
    memcpy(res->exps, poly->exps, N * poly->length * sizeof(ulong));
}

void gf2128_mpoly_set_coeff_ui_ui(gf2128_mpoly_t poly, const gf2128_t *c, 
                                  const ulong *exp, const gf2128_mpoly_ctx_t ctx) {
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
            poly->coeffs[pos] = *c;
            flint_free(cmpmask);
            flint_free(packed_exp);
            return;
        }
        if (mpoly_monomial_lt(poly->exps + N*pos, packed_exp, N, cmpmask)) {
            break;
        }
    }
    
    _gf2128_mpoly_fit_length(&poly->coeffs, &poly->coeffs_alloc,
                             &poly->exps, &poly->exps_alloc, N, poly->length + 1);
    
    for (slong i = poly->length; i > pos; i--) {
        poly->coeffs[i] = poly->coeffs[i-1];
        mpoly_monomial_set(poly->exps + N*i, poly->exps + N*(i-1), N);
    }
    
    poly->coeffs[pos] = *c;
    mpoly_monomial_set(poly->exps + N*pos, packed_exp, N);
    poly->length++;
    
    flint_free(cmpmask);
    flint_free(packed_exp);
}

void gf2128_mpoly_print(const gf2128_mpoly_t poly, const char **vars, 
                       const gf2128_mpoly_ctx_t ctx) {
    if (poly->length == 0) {
        printf("0");
        return;
    }
    
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    
    for (slong i = 0; i < poly->length; i++) {
        if (i > 0) printf(" + ");
        
        /* Print coefficient as hex value */
        printf("0x");
        gf2128_print(&poly->coeffs[i]);
        
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

/* ============================================================================
   DIVISION SUPPORT STRUCTURES AND FUNCTIONS
   ============================================================================ */

/* Initialize dense polynomial */
void gf2128_mpolyd_init(gf2128_mpolyd_t A, slong nvars) {
    A->coeffs = NULL;
    A->alloc = 0;
    A->deg_bounds = (slong *)calloc(nvars, sizeof(slong));
    A->nvars = nvars;
}

/* Clear dense polynomial */
void gf2128_mpolyd_clear(gf2128_mpolyd_t A) {
    if (A->coeffs) free(A->coeffs);
    if (A->deg_bounds) free(A->deg_bounds);
}

/* Compute offset for monomial in dense array */
slong gf2128_mpolyd_offset(const gf2128_mpolyd_t A, const ulong *exp) {
    slong off = 0;
    slong stride = 1;
    
    for (slong i = 0; i < A->nvars; i++) {
        if (exp[i] >= (ulong)A->deg_bounds[i]) return -1;
        off += (slong)exp[i] * stride;
        stride *= A->deg_bounds[i];
    }
    
    return off;
}

/* Division with remainder for univariate view */
void gf2128_mpolyd_divrem_univar(gf2128_mpolyd_t Q, gf2128_mpolyd_t R,
                                 const gf2128_mpolyd_t A, const gf2128_mpolyd_t B) {
    slong n = A->alloc;
    
    /* Copy A to R */
    memcpy(R->coeffs, A->coeffs, n * sizeof(gf2128_t));
    
    /* Clear Q */
    for (slong i = 0; i < Q->alloc; i++) {
        Q->coeffs[i] = gf2128_zero();
    }
    
    /* Find degree of B */
    slong degB = -1;
    for (slong i = n - 1; i >= 0; i--) {
        if (!gf2128_is_zero(&B->coeffs[i])) {
            degB = i;
            break;
        }
    }
    
    if (degB < 0) return;
    
    gf2128_t lc_B_inv = gf2128_inv(&B->coeffs[degB]);
    
    /* Main division loop */
    for (slong i = n - 1; i >= degB; i--) {
        if (!gf2128_is_zero(&R->coeffs[i])) {
            gf2128_t q = gf2128_mul(&R->coeffs[i], &lc_B_inv);
            
            if (i - degB < Q->alloc) {
                Q->coeffs[i - degB] = q;
            }
            
            for (slong j = 0; j <= degB && i - degB + j < n; j++) {
                gf2128_t prod = gf2128_mul(&q, &B->coeffs[j]);
                R->coeffs[i - degB + j] = gf2128_add(&R->coeffs[i - degB + j], &prod);
            }
        }
    }
}

/* Check if dense polynomial is zero */
int gf2128_mpolyd_is_zero(const gf2128_mpolyd_t A) {
    for (slong i = 0; i < A->alloc; i++) {
        if (!gf2128_is_zero(&A->coeffs[i])) return 0;
    }
    return 1;
}

/* Check if polynomial is monomial */
int gf2128_mpoly_is_monomial(const gf2128_mpoly_t poly) {
    return poly->length == 1;
}

/* Get monomial coefficient */
gf2128_t gf2128_mpoly_get_monomial_coeff(const gf2128_mpoly_t poly) {
    if (poly->length != 1) return gf2128_zero();
    return poly->coeffs[0];
}

/* Get monomial exponent */
void gf2128_mpoly_get_monomial_exp(ulong *exp, const gf2128_mpoly_t poly, 
                                   const gf2128_mpoly_ctx_t ctx) {
    if (poly->length != 1) return;
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    mpoly_get_monomial_ui(exp, poly->exps, poly->bits, ctx->minfo);
}

/* Sparse to dense conversion */
void gf2128_mpoly_to_mpolyd(gf2128_mpolyd_t A, const gf2128_mpoly_t B, 
                            const gf2128_mpoly_ctx_t ctx) {
    slong nvars = ctx->minfo->nvars;
    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);
    ulong *exp = (ulong *)malloc(nvars * sizeof(ulong));
    
    /* Clear all coefficients */
    for (slong i = 0; i < A->alloc; i++) {
        A->coeffs[i] = gf2128_zero();
    }
    
    /* Convert each term */
    for (slong i = 0; i < B->length; i++) {
        mpoly_get_monomial_ui(exp, B->exps + N*i, B->bits, ctx->minfo);
        slong off = gf2128_mpolyd_offset(A, exp);
        if (off >= 0 && off < A->alloc) {
            A->coeffs[off] = B->coeffs[i];
        }
    }
    
    free(exp);
}

/* Dense to sparse conversion */
void gf2128_mpolyd_to_mpoly(gf2128_mpoly_t A, const gf2128_mpolyd_t B,
                            const gf2128_mpoly_ctx_t ctx) {
    slong nvars = ctx->minfo->nvars;
    ulong *exp = (ulong *)calloc(nvars, sizeof(ulong));
    
    /* Calculate needed bits */
    flint_bitcnt_t bits_needed = MPOLY_MIN_BITS;
    for (slong i = 0; i < nvars; i++) {
        if (B->deg_bounds[i] > 0) {
            slong bits_i = FLINT_BIT_COUNT(B->deg_bounds[i] - 1);
            bits_needed = FLINT_MAX(bits_needed, bits_i);
        }
    }
    if (bits_needed < 16) bits_needed = 16;
    bits_needed = mpoly_fix_bits(bits_needed, ctx->minfo);
    
    gf2128_mpoly_zero(A, ctx);
    
    /* Traverse dense array */
    for (slong off = 0; off < B->alloc; off++) {
        if (gf2128_is_zero(&B->coeffs[off])) continue;
        
        /* Reconstruct exponent vector from offset */
        slong temp = off;
        for (slong i = 0; i < nvars; i++) {
            exp[i] = (ulong)(temp % B->deg_bounds[i]);
            temp /= B->deg_bounds[i];
        }
        
        /* Add term to polynomial */
        if (A->bits < bits_needed) {
            gf2128_mpoly_fit_length_reset_bits(A, A->length + 1, bits_needed, ctx);
        }
        gf2128_mpoly_set_coeff_ui_ui(A, &B->coeffs[off], exp, ctx);
    }
    
    free(exp);
}

/* Multivariate division with remainder */
void gf2128_mpolyd_divrem_multivar(gf2128_mpolyd_t Q, gf2128_mpolyd_t R,
                                   const gf2128_mpolyd_t A, const gf2128_mpolyd_t B,
                                   const gf2128_mpoly_ctx_t ctx) {
    /* For now, use univariate view division */
    gf2128_mpolyd_divrem_univar(Q, R, A, B);
}

/* Set degree bounds */
int gf2128_mpolyd_set_degbounds(gf2128_mpolyd_t A, const slong *bounds) {
    slong size = 1;
    
    for (slong i = 0; i < A->nvars; i++) {
        A->deg_bounds[i] = bounds[i];
        if (bounds[i] <= 0) return 0;
        
        /* Check overflow */
        if (size > WORD_MAX / bounds[i]) return 0;
        size *= bounds[i];
    }
    
    /* Limit size for GF(2^128) - smaller than GF(2^8) due to larger element size */
    if (size > (1L << 22)) return 0; /* 4M coefficients max */
    
    A->alloc = size;
    A->coeffs = (gf2128_t *)calloc(size, sizeof(gf2128_t));
    return A->coeffs != NULL;
}

/* Monomial division */
int gf2128_mpoly_divides_monomial(gf2128_mpoly_t Q, const gf2128_mpoly_t A, 
                                  const gf2128_mpoly_t B, const gf2128_mpoly_ctx_t ctx) {
    slong nvars = ctx->minfo->nvars;
    ulong *exp_b = (ulong *)malloc(nvars * sizeof(ulong));
    ulong *exp_a = (ulong *)malloc(nvars * sizeof(ulong));
    ulong *exp_q = (ulong *)malloc(nvars * sizeof(ulong));
    gf2128_t coeff_b = gf2128_mpoly_get_monomial_coeff(B);
    gf2128_t coeff_b_inv;
    slong N;
    
    if (gf2128_is_zero(&coeff_b)) {
        free(exp_b); free(exp_a); free(exp_q);
        return 0;
    }
    
    coeff_b_inv = gf2128_inv(&coeff_b);
    gf2128_mpoly_get_monomial_exp(exp_b, B, ctx);
    
    /* Determine result bits */
    flint_bitcnt_t bits = FLINT_MAX(A->bits, B->bits);
    if (bits < 16) bits = 16;
    
    /* Initialize result */
    gf2128_mpoly_fit_length_reset_bits(Q, A->length, bits, ctx);
    Q->length = 0;
    N = mpoly_words_per_exp(bits, ctx->minfo);
    
    /* Divide each term of A */
    for (slong i = 0; i < A->length; i++) {
        mpoly_get_monomial_ui(exp_a, A->exps + N*i, A->bits, ctx->minfo);
        
        /* Check divisibility */
        int divisible = 1;
        for (slong j = 0; j < nvars; j++) {
            if (exp_a[j] < exp_b[j]) {
                divisible = 0;
                break;
            }
            exp_q[j] = exp_a[j] - exp_b[j];
        }
        
        if (!divisible) {
            free(exp_b); free(exp_a); free(exp_q);
            gf2128_mpoly_zero(Q, ctx);
            return 0;
        }
        
        /* Compute coefficient */
        gf2128_t coeff_q = gf2128_mul(&A->coeffs[i], &coeff_b_inv);
        
        if (!gf2128_is_zero(&coeff_q)) {
            /* Add to result */
            mpoly_set_monomial_ui(Q->exps + N*Q->length, exp_q, bits, ctx->minfo);
            Q->coeffs[Q->length] = coeff_q;
            Q->length++;
        }
    }
    
    free(exp_b); free(exp_a); free(exp_q);
    return 1;
}

/* Dense representation division */
int gf2128_mpoly_divides_dense(gf2128_mpoly_t Q, const gf2128_mpoly_t A, 
                               const gf2128_mpoly_t B, const gf2128_mpoly_ctx_t ctx) {
    slong nvars = ctx->minfo->nvars;
    int success = 0;
    
    /* Get degree bounds */
    slong *degs_A = (slong *)malloc(nvars * sizeof(slong));
    slong *degs_B = (slong *)malloc(nvars * sizeof(slong));
    slong *bounds = (slong *)malloc(nvars * sizeof(slong));
    
    mpoly_degrees_si(degs_A, A->exps, A->length, A->bits, ctx->minfo);
    mpoly_degrees_si(degs_B, B->exps, B->length, B->bits, ctx->minfo);
    
    /* Check if division is possible */
    for (slong i = 0; i < nvars; i++) {
        if (degs_A[i] < degs_B[i]) {
            free(degs_A); free(degs_B); free(bounds);
            gf2128_mpoly_zero(Q, ctx);
            return 0;
        }
        bounds[i] = degs_A[i] + 1;
    }
    
    /* Initialize dense polynomials */
    gf2128_mpolyd_t Ad, Bd, Qd, Rd;
    gf2128_mpolyd_init(Ad, nvars);
    gf2128_mpolyd_init(Bd, nvars);
    gf2128_mpolyd_init(Qd, nvars);
    gf2128_mpolyd_init(Rd, nvars);
    
    /* Set degree bounds */
    if (!gf2128_mpolyd_set_degbounds(Ad, bounds) ||
        !gf2128_mpolyd_set_degbounds(Bd, bounds) ||
        !gf2128_mpolyd_set_degbounds(Qd, bounds) ||
        !gf2128_mpolyd_set_degbounds(Rd, bounds)) {
        goto cleanup;
    }
    
    /* Convert to dense representation */
    gf2128_mpoly_to_mpolyd(Ad, A, ctx);
    gf2128_mpoly_to_mpolyd(Bd, B, ctx);
    
    /* Perform division */
    gf2128_mpolyd_divrem_multivar(Qd, Rd, Ad, Bd, ctx);
    
    /* Check if remainder is zero */
    if (gf2128_mpolyd_is_zero(Rd)) {
        /* Convert quotient back to sparse representation */
        gf2128_mpolyd_to_mpoly(Q, Qd, ctx);
        success = 1;
    } else {
        gf2128_mpoly_zero(Q, ctx);
        success = 0;
    }
    
cleanup:
    gf2128_mpolyd_clear(Ad);
    gf2128_mpolyd_clear(Bd);
    gf2128_mpolyd_clear(Qd);
    gf2128_mpolyd_clear(Rd);
    free(degs_A);
    free(degs_B);
    free(bounds);
    
    return success;
}

/* Main division function */
int gf2128_mpoly_divides(gf2128_mpoly_t Q, const gf2128_mpoly_t A, 
                         const gf2128_mpoly_t B, const gf2128_mpoly_ctx_t ctx) {
    /* Special cases */
    if (B->length == 0) {
        if (A->length == 0) {
            gf2128_mpoly_zero(Q, ctx);
            return 1;
        }
        return 0; /* Division by zero */
    }
    
    if (A->length == 0) {
        gf2128_mpoly_zero(Q, ctx);
        return 1;
    }
    
    /* If B is monomial, use direct division */
    if (gf2128_mpoly_is_monomial(B)) {
        return gf2128_mpoly_divides_monomial(Q, A, B, ctx);
    }
    
    /* For general case, use dense representation division */
    return gf2128_mpoly_divides_dense(Q, A, B, ctx);
}

/* ============================================================================
   ADDITIONAL HELPER FUNCTIONS
   ============================================================================ */

/* Simple multiplication for verification */
void gf2128_mpoly_mul_simple(gf2128_mpoly_t res, const gf2128_mpoly_t a, 
                             const gf2128_mpoly_t b, const gf2128_mpoly_ctx_t ctx) {
    slong nvars = ctx->minfo->nvars;
    slong N = mpoly_words_per_exp(a->bits, ctx->minfo);
    ulong *exp_a, *exp_b, *exp_sum;
    
    gf2128_mpoly_zero(res, ctx);
    
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
            gf2128_t coeff = gf2128_mul(&a->coeffs[i], &b->coeffs[j]);
            
            if (!gf2128_is_zero(&coeff)) {
                /* Get current coefficient */
                gf2128_t current = gf2128_zero();
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
                gf2128_t new_coeff = gf2128_add(&current, &coeff);
                gf2128_mpoly_set_coeff_ui_ui(res, &new_coeff, exp_sum, ctx);
            }
        }
    }
    
    free(exp_a);
    free(exp_b);
    free(exp_sum);
}

/* Check if two polynomials are equal */
int gf2128_mpoly_equal(const gf2128_mpoly_t A, const gf2128_mpoly_t B, 
                       const gf2128_mpoly_ctx_t ctx) {
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
                if (!gf2128_equal(&A->coeffs[i], &B->coeffs[j])) {
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
void gf2128_mpoly_randtest(gf2128_mpoly_t poly, flint_rand_t state,
                           slong length, slong exp_bound, 
                           const gf2128_mpoly_ctx_t ctx) {
    slong nvars = ctx->minfo->nvars;
    ulong *exp = (ulong *)malloc(nvars * sizeof(ulong));
    
    gf2128_mpoly_zero(poly, ctx);
    
    for (slong i = 0; i < length; i++) {
        /* Generate random exponents */
        for (slong j = 0; j < nvars; j++) {
            exp[j] = n_randint(state, exp_bound);
        }
        
        /* Generate random non-zero coefficient */
        gf2128_t c;
        do {
            c.low = n_randtest(state);
            c.high = n_randtest(state);
        } while (gf2128_is_zero(&c));
        
        gf2128_mpoly_set_coeff_ui_ui(poly, &c, exp, ctx);
    }
    
    free(exp);
}

/* Conversion functions for FLINT compatibility */
void gf2128_mpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t res, const gf2128_mpoly_t poly,
                                   const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx) {
    fq_nmod_mpoly_zero(res, fq_mpoly_ctx);
    
    if (poly->length == 0) return;
    
    slong N = mpoly_words_per_exp(poly->bits, fq_mpoly_ctx->minfo);
    slong nvars = fq_mpoly_ctx->minfo->nvars;
    
    for (slong i = 0; i < poly->length; i++) {
        if (!gf2128_is_zero(&poly->coeffs[i])) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, fqctx);
            
            /* Convert GF(2^128) element to FLINT */
            gf2128_to_fq_nmod(coeff, &poly->coeffs[i], fqctx);
            
            ulong *exp = (ulong *) flint_malloc(nvars * sizeof(ulong));
            mpoly_get_monomial_ui(exp, poly->exps + N*i, poly->bits, fq_mpoly_ctx->minfo);
            
            fq_nmod_mpoly_set_coeff_fq_nmod_ui(res, coeff, exp, fq_mpoly_ctx);
            
            flint_free(exp);
            fq_nmod_clear(coeff, fqctx);
        }
    }
}

void fq_nmod_mpoly_to_gf2128_mpoly(gf2128_mpoly_t res, const fq_nmod_mpoly_t poly,
                                   const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx) {
    gf2128_mpoly_ctx_t ctx;
    gf2128_mpoly_ctx_init(ctx, fq_mpoly_ctx->minfo->nvars, fq_mpoly_ctx->minfo->ord);
    
    gf2128_mpoly_zero(res, ctx);
    
    slong len = fq_nmod_mpoly_length(poly, fq_mpoly_ctx);
    if (len == 0) {
        gf2128_mpoly_ctx_clear(ctx);
        return;
    }
    
    flint_bitcnt_t bits = FLINT_MAX(poly->bits, MPOLY_MIN_BITS);
    gf2128_mpoly_fit_length_reset_bits(res, len, bits, ctx);
    res->length = 0;
    
    slong nvars = fq_mpoly_ctx->minfo->nvars;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    
    slong actual_terms = 0;
    for (slong i = 0; i < len; i++) {
        fq_nmod_t coeff;
        fq_nmod_init(coeff, fqctx);
        
        fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, poly, i, fq_mpoly_ctx);
        
        /* Convert FLINT element to GF(2^128) */
        gf2128_t gf2128_coeff = fq_nmod_to_gf2128(coeff, fqctx);
        
        if (!gf2128_is_zero(&gf2128_coeff)) {
            ulong *exp = (ulong *) flint_malloc(nvars * sizeof(ulong));
            fq_nmod_mpoly_get_term_exp_ui(exp, poly, i, fq_mpoly_ctx);
            
            mpoly_set_monomial_ui(res->exps + N*actual_terms, exp, bits, ctx->minfo);
            res->coeffs[actual_terms] = gf2128_coeff;
            actual_terms++;
            
            flint_free(exp);
        }
        
        fq_nmod_clear(coeff, fqctx);
    }
    
    res->length = actual_terms;
    gf2128_mpoly_ctx_clear(ctx);
}

/* Timer function */
double get_wall_time_128(void) {
    struct timeval time;
    if (gettimeofday(&time, NULL)) {
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * 0.000001;
}