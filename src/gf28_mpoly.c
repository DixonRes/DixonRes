/* gf28_mpoly.c - GF(2^8) Multivariate Polynomial Implementation */

#include "gf28_mpoly.h"
#include <emmintrin.h>  // For SSE2 instructions

/* ============================================================================
   BASIC INITIALIZATION AND MEMORY MANAGEMENT
   ============================================================================ */

void gf28_mpoly_init(gf28_mpoly_t poly, const gf28_mpoly_ctx_t ctx) {
    poly->coeffs = NULL;
    poly->exps = NULL;
    poly->length = 0;
    poly->coeffs_alloc = 0;
    poly->exps_alloc = 0;
    poly->bits = MPOLY_MIN_BITS;
}

void gf28_mpoly_clear(gf28_mpoly_t poly, const gf28_mpoly_ctx_t ctx) {
    if (poly->coeffs) flint_free(poly->coeffs);
    if (poly->exps) flint_free(poly->exps);
}

void gf28_mpoly_ctx_init(gf28_mpoly_ctx_t ctx, slong nvars, const ordering_t ord) {
    mpoly_ctx_init(ctx->minfo, nvars, ord);
}

void gf28_mpoly_ctx_clear(gf28_mpoly_ctx_t ctx) {
    mpoly_ctx_clear(ctx->minfo);
}

void gf28_mpoly_zero(gf28_mpoly_t poly, const gf28_mpoly_ctx_t ctx) {
    poly->length = 0;
}

void _gf28_mpoly_fit_length(uint8_t **coeffs, slong *coeffs_alloc,
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

void gf28_mpoly_fit_length_reset_bits(gf28_mpoly_t poly, slong len, 
                                     flint_bitcnt_t bits, const gf28_mpoly_ctx_t ctx) {
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    _gf28_mpoly_fit_length(&poly->coeffs, &poly->coeffs_alloc,
                           &poly->exps, &poly->exps_alloc, N, len);
    poly->bits = bits;
}

void gf28_mpoly_init3(gf28_mpoly_t poly, slong alloc, flint_bitcnt_t bits,
                     const gf28_mpoly_ctx_t ctx) {
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    poly->coeffs = (uint8_t *) flint_malloc(alloc * sizeof(uint8_t));
    poly->exps = (ulong *) flint_malloc(N * alloc * sizeof(ulong));
    poly->coeffs_alloc = alloc;
    poly->exps_alloc = N * alloc;
    poly->length = 0;
    poly->bits = bits;
}

void gf28_mpoly_swap(gf28_mpoly_t poly1, gf28_mpoly_t poly2, const gf28_mpoly_ctx_t ctx) {
    gf28_mpoly_struct t = *poly1;
    *poly1 = *poly2;
    *poly2 = t;
}

void _gf28_mpoly_set_length(gf28_mpoly_t poly, slong len, const gf28_mpoly_ctx_t ctx) {
    poly->length = len;
}

/* ============================================================================
   HELPER FUNCTIONS
   ============================================================================ */

void gf28_mpoly_set(gf28_mpoly_t res, const gf28_mpoly_t poly, const gf28_mpoly_ctx_t ctx) {
    if (res == poly) return;
    
    gf28_mpoly_fit_length_reset_bits(res, poly->length, poly->bits, ctx);
    res->length = poly->length;
    
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    memcpy(res->coeffs, poly->coeffs, poly->length * sizeof(uint8_t));
    memcpy(res->exps, poly->exps, N * poly->length * sizeof(ulong));
}

void gf28_mpoly_set_coeff_ui_ui(gf28_mpoly_t poly, uint8_t c, 
                                const ulong *exp, const gf28_mpoly_ctx_t ctx) {
    if (poly->bits == 0) {
        poly->bits = MPOLY_MIN_BITS;
    }
    
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    
    ulong *cmpmask = (ulong *) flint_malloc(N * FLINT_BITS * sizeof(ulong));
    mpoly_get_cmpmask(cmpmask, N, poly->bits, ctx->minfo);
    
    ulong *packed_exp = (ulong *) flint_malloc(N * sizeof(ulong));
    mpoly_set_monomial_ui(packed_exp, exp, poly->bits, ctx->minfo);
    
    /* Find position */
    slong pos = 0;
    for (pos = 0; pos < poly->length; pos++) {
        if (mpoly_monomial_equal(poly->exps + N*pos, packed_exp, N)) {
            /* Update existing coefficient */
            if (c == 0) {
                /* Remove this term */
                for (slong i = pos; i < poly->length - 1; i++) {
                    poly->coeffs[i] = poly->coeffs[i + 1];
                    mpoly_monomial_set(poly->exps + N*i, poly->exps + N*(i + 1), N);
                }
                poly->length--;
            } else {
                poly->coeffs[pos] = c;
            }
            flint_free(cmpmask);
            flint_free(packed_exp);
            return;
        }
        if (mpoly_monomial_lt(poly->exps + N*pos, packed_exp, N, cmpmask)) {
            break;
        }
    }
    
    /* Insert new term if c != 0 */
    if (c != 0) {
        _gf28_mpoly_fit_length(&poly->coeffs, &poly->coeffs_alloc,
                               &poly->exps, &poly->exps_alloc, N, poly->length + 1);
        
        for (slong i = poly->length; i > pos; i--) {
            poly->coeffs[i] = poly->coeffs[i-1];
            mpoly_monomial_set(poly->exps + N*i, poly->exps + N*(i-1), N);
        }
        
        poly->coeffs[pos] = c;
        mpoly_monomial_set(poly->exps + N*pos, packed_exp, N);
        poly->length++;
    }
    
    flint_free(cmpmask);
    flint_free(packed_exp);
}

void gf28_mpoly_print(const gf28_mpoly_t poly, const char **vars, 
                     const gf28_mpoly_ctx_t ctx) {
    if (poly->length == 0) {
        printf("0");
        return;
    }
    
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    
    for (slong i = 0; i < poly->length; i++) {
        if (i > 0) printf(" + ");
        
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

/* ============================================================================
   ARRAY MULTIPLICATION
   ============================================================================ */

void _gf28_mpoly_addmul_array1_safe(uint8_t *poly1, slong array_size,
                                   const uint8_t *poly2, const ulong *exp2, slong len2,
                                   const uint8_t *poly3, const ulong *exp3, slong len3) {
    slong ii, i, jj, j;
    uint8_t *c2;
    
    for (ii = 0; ii < len2 + BLOCK; ii += BLOCK) {
        for (jj = 0; jj < len3 + BLOCK; jj += BLOCK) {
            for (i = ii; i < FLINT_MIN(ii + BLOCK, len2); i++) {
                slong offset2 = (slong)exp2[i];
                
                if (offset2 >= array_size) {
                    continue;
                }
                
                c2 = poly1 + offset2;
                
                if (poly2[i] != 0) {
                    const uint8_t* mul_row = gf28_get_scalar_row(poly2[i]);
                    
                    for (j = jj; j < FLINT_MIN(jj + BLOCK, len3); j++) {
                        slong offset3 = (slong)exp3[j];
                        slong total_offset = offset2 + offset3;
                        
                        if (total_offset >= array_size) {
                            continue;
                        }
                        
                        c2[offset3] ^= mul_row[poly3[j]];
                    }
                }
            }
        }
    }
}

slong gf28_mpoly_append_array_LEX_safe(gf28_mpoly_t P, slong Plen, uint8_t *coeff_array,
                                      const ulong *mults, slong num, slong array_size, 
                                      slong top, const gf28_mpoly_ctx_t ctx) {
    slong off, j;
    slong topmult = num == 0 ? 1 : mults[num - 1];
    slong lastd = topmult - 1;
    slong reset = array_size/topmult;
    slong counter = reset;
    ulong startexp = ((ulong)top << (P->bits*num)) + ((ulong)lastd << (P->bits*(num-1)));
    uint8_t coeff;
    
    for (off = array_size - 1; off >= 0; off--) {
        if (coeff_array[off] != 0) {
            coeff = coeff_array[off];
            slong d = off;
            ulong exp = startexp;
            
            for (j = 0; j + 1 < num; j++) {
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

void gf28_dynamic_array_init(gf28_dynamic_array_t *arr, slong size) {
    arr->size = size;
    arr->data = (uint8_t *)calloc(size, sizeof(uint8_t));
    if (!arr->data && size > 0) {
        flint_printf("Failed to allocate %ld GF(2^8) elements\n", size);
        flint_abort();
    }
}

void gf28_dynamic_array_clear(gf28_dynamic_array_t *arr) {
    if (arr->data) {
        free(arr->data);
        arr->data = NULL;
    }
}

void _gf28_mpoly_mul_array_chunked_LEX(gf28_mpoly_t P,
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
    
    gf28_dynamic_array_t coeff_array;
    gf28_dynamic_array_init(&coeff_array, array_size);
    
    if (!coeff_array.data) {
        flint_printf("Memory allocation failed for array of size %ld\n", array_size);
        flint_free(Apexp);
        flint_free(Bpexp);
        TMP_END;
        return;
    }
    
    for (Pi = 0; Pi < Pl; Pi++) {
        memset(coeff_array.data, 0, array_size);
        
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
    
    gf28_dynamic_array_clear(&coeff_array);
    flint_free(Apexp);
    flint_free(Bpexp);
    TMP_END;
}

int _gf28_mpoly_mul_array_LEX(gf28_mpoly_t A,
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
    
    if (array_size > (1L << 28)) {
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

double get_wall_time_8(void) {
    struct timeval time;
    if (gettimeofday(&time, NULL)) {
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * 0.000001;
}

void fq_nmod_mpoly_to_gf28_mpoly(gf28_mpoly_t res, const fq_nmod_mpoly_t poly,
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
        
        uint8_t gf28_coeff = fq_nmod_to_gf28_elem(coeff, fqctx);
        
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

/* ============================================================================
   INTELLIGENT MULTIPLICATION DISPATCH
   ============================================================================ */

/* Check if multivariate polynomial is effectively univariate */
int gf28_mpoly_is_univariate(const gf28_mpoly_t poly, const gf28_mpoly_ctx_t ctx, 
                            slong *main_var) {
    if (poly->length == 0) return 1;
    
    slong nvars = ctx->minfo->nvars;
    if (nvars == 1) {
        *main_var = 0;
        return 1;
    }
    
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    ulong *exp = (ulong *)malloc(nvars * sizeof(ulong));
    
    int found_var = -1;
    
    for (slong i = 0; i < poly->length; i++) {
        mpoly_get_monomial_ui(exp, poly->exps + N*i, poly->bits, ctx->minfo);
        
        for (slong j = 0; j < nvars; j++) {
            if (exp[j] > 0) {
                if (found_var == -1) {
                    found_var = j;  // First variable found
                } else if (found_var != j) {
                    // Found a different variable, not univariate
                    free(exp);
                    return 0;
                }
            }
        }
    }
    
    free(exp);
    *main_var = found_var >= 0 ? found_var : 0;
    return 1;
}

/* Get degree of univariate polynomial */
slong gf28_mpoly_univariate_degree(const gf28_mpoly_t poly, const gf28_mpoly_ctx_t ctx,
                                  slong main_var) {
    if (poly->length == 0) return -1;
    
    slong nvars = ctx->minfo->nvars;
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    ulong *exp = (ulong *)malloc(nvars * sizeof(ulong));
    
    slong max_degree = 0;
    
    for (slong i = 0; i < poly->length; i++) {
        mpoly_get_monomial_ui(exp, poly->exps + N*i, poly->bits, ctx->minfo);
        if ((slong)exp[main_var] > max_degree) {
            max_degree = (slong)exp[main_var];
        }
    }
    
    free(exp);
    return max_degree;
}

/* Convert gf28_mpoly to gf28_poly (for univariate case) */
void gf28_mpoly_to_gf28_poly_univar(gf28_poly_t res, const gf28_mpoly_t poly,
                                   const gf28_mpoly_ctx_t ctx, slong main_var) {
    gf28_poly_zero(res);
    
    if (poly->length == 0) return;
    
    slong nvars = ctx->minfo->nvars;
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    ulong *exp = (ulong *)malloc(nvars * sizeof(ulong));
    
    for (slong i = 0; i < poly->length; i++) {
        mpoly_get_monomial_ui(exp, poly->exps + N*i, poly->bits, ctx->minfo);
        gf28_poly_set_coeff(res, exp[main_var], poly->coeffs[i]);
    }
    
    free(exp);
}

/* Convert gf28_poly to gf28_mpoly (for univariate case) */
void gf28_poly_to_gf28_mpoly_univar(gf28_mpoly_t res, const gf28_poly_t poly,
                                   const gf28_mpoly_ctx_t ctx, slong main_var) {
    gf28_mpoly_zero(res, ctx);
    
    if (poly->length == 0) return;
    
    /* Count non-zero terms first */
    slong nonzero_count = 0;
    for (slong i = 0; i < poly->length; i++) {
        if (poly->coeffs[i] != 0) {
            nonzero_count++;
        }
    }
    
    if (nonzero_count == 0) return;
    
    /* Calculate required bits for the maximum degree */
    slong max_degree = poly->length - 1;
    flint_bitcnt_t required_bits = FLINT_BIT_COUNT(max_degree) + 1;
    required_bits = FLINT_MAX(required_bits, MPOLY_MIN_BITS);
    required_bits = mpoly_fix_bits(required_bits, ctx->minfo);
    
    /* Pre-allocate exactly what we need */
    gf28_mpoly_fit_length_reset_bits(res, nonzero_count, required_bits, ctx);
    
    slong nvars = ctx->minfo->nvars;
    slong N = mpoly_words_per_exp(required_bits, ctx->minfo);
    
    /* Build polynomial directly in reverse order (highest degree first for LEX ordering) */
    slong term_idx = 0;
    for (slong i = poly->length - 1; i >= 0; i--) {
        if (poly->coeffs[i] != 0) {
            /* Set coefficient */
            res->coeffs[term_idx] = poly->coeffs[i];
            
            /* Set exponent directly - much faster than using mpoly_set_monomial_ui */
            memset(res->exps + N * term_idx, 0, N * sizeof(ulong));
            
            /* For LEX ordering with multiple variables, we need to pack the exponent correctly */
            if (required_bits <= FLINT_BITS) {
                /* Single word case - pack degree into the appropriate field */
                ulong exp_word = 0;
                
                /* Calculate bit position for main_var in LEX ordering */
                slong bits_per_var = required_bits;
                slong bit_pos = (nvars - 1 - main_var) * bits_per_var;
                
                if (bit_pos + bits_per_var <= FLINT_BITS) {
                    exp_word = ((ulong)i) << bit_pos;
                    res->exps[N * term_idx] = exp_word;
                } else {
                    /* Fall back to safe method for complex packing */
                    ulong *temp_exp = (ulong *)calloc(nvars, sizeof(ulong));
                    temp_exp[main_var] = i;
                    mpoly_set_monomial_ui(res->exps + N * term_idx, temp_exp, required_bits, ctx->minfo);
                    free(temp_exp);
                }
            } else {
                /* Multi-word case - use safe method */
                ulong *temp_exp = (ulong *)calloc(nvars, sizeof(ulong));
                temp_exp[main_var] = i;
                mpoly_set_monomial_ui(res->exps + N * term_idx, temp_exp, required_bits, ctx->minfo);
                free(temp_exp);
            }
            
            term_idx++;
        }
    }
    
    res->length = nonzero_count;
}

/* Multiply using FLINT fq_nmod_poly_mul with conversions */
void gf28_mpoly_mul_flint_univar(gf28_mpoly_t res, const gf28_mpoly_t a, 
                                const gf28_mpoly_t b, const gf28_mpoly_ctx_t ctx,
                                slong main_var) {
    /* Initialize FLINT context for GF(2^8) */
    nmod_poly_t modulus;
    fq_nmod_ctx_t fq_ctx;
    
    nmod_poly_init(modulus, 2);
    nmod_poly_set_coeff_ui(modulus, 0, 1);  // x^8 + x^4 + x^3 + x + 1
    nmod_poly_set_coeff_ui(modulus, 2, 1);
    nmod_poly_set_coeff_ui(modulus, 3, 1);
    nmod_poly_set_coeff_ui(modulus, 4, 1);
    nmod_poly_set_coeff_ui(modulus, 8, 1);
    
    fq_nmod_ctx_init_modulus(fq_ctx, modulus, "x");
    
    /* Convert to GF(2^8) polynomials */
    gf28_poly_t poly_a, poly_b, poly_res;
    gf28_poly_init(poly_a);
    gf28_poly_init(poly_b);
    gf28_poly_init(poly_res);
    
    gf28_mpoly_to_gf28_poly_univar(poly_a, a, ctx, main_var);
    gf28_mpoly_to_gf28_poly_univar(poly_b, b, ctx, main_var);
    
    /* Convert to FLINT format */
    fq_nmod_poly_t fa, fb, fresult;
    fq_nmod_poly_init(fa, fq_ctx);
    fq_nmod_poly_init(fb, fq_ctx);
    fq_nmod_poly_init(fresult, fq_ctx);
    
    gf28_poly_to_fq_nmod_poly(fa, poly_a, fq_ctx);
    gf28_poly_to_fq_nmod_poly(fb, poly_b, fq_ctx);
    
    /* Multiply using FLINT */
    fq_nmod_poly_mul(fresult, fa, fb, fq_ctx);
    
    /* Convert back */
    fq_nmod_poly_to_gf28_poly(poly_res, fresult, fq_ctx);
    gf28_poly_to_gf28_mpoly_univar(res, poly_res, ctx, main_var);
    
    /* Cleanup */
    gf28_poly_clear(poly_a);
    gf28_poly_clear(poly_b);
    gf28_poly_clear(poly_res);
    fq_nmod_poly_clear(fa, fq_ctx);
    fq_nmod_poly_clear(fb, fq_ctx);
    fq_nmod_poly_clear(fresult, fq_ctx);
    fq_nmod_ctx_clear(fq_ctx);
    nmod_poly_clear(modulus);
}

/* ============================================================================
   DIVISION IMPLEMENTATION WITH DEBUG
   ============================================================================ */

void gf28_mpolyd_init(gf28_mpolyd_t A, slong nvars) {
    A->coeffs = NULL;
    A->alloc = 0;
    A->deg_bounds = (slong *)calloc(nvars, sizeof(slong));
    A->nvars = nvars;
}

void gf28_mpolyd_clear(gf28_mpolyd_t A) {
    if (A->coeffs) free(A->coeffs);
    if (A->deg_bounds) free(A->deg_bounds);
}

slong gf28_mpolyd_offset(const gf28_mpolyd_t A, const ulong *exp) {
    slong off = 0;
    slong stride = 1;
    
    for (slong i = 0; i < A->nvars; i++) {
        if (exp[i] >= (ulong)A->deg_bounds[i]) return -1;
        off += (slong)exp[i] * stride;
        stride *= A->deg_bounds[i];
    }
    
    return off;
}

void batch_xor_sse2(uint8_t* dst, const uint8_t* src, slong len) {
    slong i = 0;
    
    // Process 16-byte aligned parts
    for (; i + 16 <= len; i += 16) {
        __m128i a = _mm_loadu_si128((__m128i*)(dst + i));
        __m128i b = _mm_loadu_si128((__m128i*)(src + i));
        __m128i result = _mm_xor_si128(a, b);
        _mm_storeu_si128((__m128i*)(dst + i), result);
    }
    
    // Process remaining bytes
    for (; i < len; i++) {
        dst[i] ^= src[i];
    }
}

void gf28_mpolyd_divrem_univar(gf28_mpolyd_t Q, gf28_mpolyd_t R,
                              const gf28_mpolyd_t A, const gf28_mpolyd_t B) {
    slong n = A->alloc;
    
    #if DEBUG_DIVISION
    printf("\n[DEBUG] gf28_mpolyd_divrem_univar: n=%ld\n", n);
    #endif
    
    memcpy(R->coeffs, A->coeffs, n * sizeof(uint8_t));
    memset(Q->coeffs, 0, Q->alloc * sizeof(uint8_t));
    
    slong degB = -1;
    for (slong i = n - 1; i >= 0; i--) {
        if (B->coeffs[i] != 0) {
            degB = i;
            break;
        }
    }
    
    if (degB < 0) {
        #if DEBUG_DIVISION
        printf("[DEBUG] B is zero polynomial\n");
        #endif
        return;
    }
    
    #if DEBUG_DIVISION
    printf("[DEBUG] degB = %ld, lc(B) = 0x%02x\n", degB, B->coeffs[degB]);
    #endif
    
    uint8_t lc_B_inv = gf28_inv(B->coeffs[degB]);
    
    #if DEBUG_DIVISION
    printf("[DEBUG] lc(B)^{-1} = 0x%02x\n", lc_B_inv);
    printf("[DEBUG] Check: lc(B) * lc(B)^{-1} = 0x%02x\n", gf28_mul(B->coeffs[degB], lc_B_inv));
    #endif
    
    // Dynamic allocation to avoid stack overflow
    uint8_t (*B_multiples)[degB + 1] = malloc(256 * (degB + 1) * sizeof(uint8_t));
    if (!B_multiples) {
        printf("[ERROR] Failed to allocate memory for B_multiples\n");
        return;
    }
    
    // Precomputation
    for (int k = 0; k < 256; k++) {
        const uint8_t* row = gf28_get_scalar_row(k);
        for (slong j = 0; j <= degB; j++) {
            B_multiples[k][j] = row[B->coeffs[j]];
        }
    }
    
    for (slong i = n - 1; i >= degB; i--) {
        if (R->coeffs[i] != 0) {
            uint8_t q = gf28_mul(R->coeffs[i], lc_B_inv);
            
            if (i - degB < Q->alloc) {  // Add boundary check
                Q->coeffs[i - degB] = q;
            }
            
            // Use precomputed multiples directly
            const uint8_t* B_times_q = B_multiples[q];
            
            // Ensure no boundary violations
            slong len = FLINT_MIN(degB + 1, n - (i - degB));
            
            if (0 && len >= 16) {
                batch_xor_sse2(R->coeffs + i - degB, B_times_q, len);
            } else {
                // Otherwise use standard loop
                for (slong j = 0; j < len; j++) {
                    R->coeffs[i - degB + j] ^= B_times_q[j];
                }
            }
        }
    }
    
    // Free memory
    free(B_multiples);
    
    #if DEBUG_DIVISION
    printf("[DEBUG] Division complete. Checking remainder...\n");
    slong degR = -1;
    for (slong i = n - 1; i >= 0; i--) {
        if (R->coeffs[i] != 0) {
            degR = i;
            break;
        }
    }
    printf("[DEBUG] degR = %ld (should be < %ld)\n", degR, degB);
    #endif
}

int gf28_mpolyd_is_zero(const gf28_mpolyd_t A) {
    for (slong i = 0; i < A->alloc; i++) {
        if (A->coeffs[i] != 0) return 0;
    }
    return 1;
}

int gf28_mpoly_is_monomial(const gf28_mpoly_t poly) {
    return poly->length == 1;
}

uint8_t gf28_mpoly_get_monomial_coeff(const gf28_mpoly_t poly) {
    if (poly->length != 1) return 0;
    return poly->coeffs[0];
}

void gf28_mpoly_get_monomial_exp(ulong *exp, const gf28_mpoly_t poly, 
                                const gf28_mpoly_ctx_t ctx) {
    if (poly->length != 1) return;
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    mpoly_get_monomial_ui(exp, poly->exps, poly->bits, ctx->minfo);
}

void gf28_mpoly_to_mpolyd(gf28_mpolyd_t A, const gf28_mpoly_t B, 
                         const gf28_mpoly_ctx_t ctx) {
    slong nvars = ctx->minfo->nvars;
    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);
    ulong *exp = (ulong *)malloc(nvars * sizeof(ulong));
    
    memset(A->coeffs, 0, A->alloc * sizeof(uint8_t));
    
    #if DEBUG_DIVISION
    printf("\n[DEBUG] gf28_mpoly_to_mpolyd: B has %ld terms\n", B->length);
    #endif
    
    for (slong i = 0; i < B->length; i++) {
        mpoly_get_monomial_ui(exp, B->exps + N*i, B->bits, ctx->minfo);
        slong off = gf28_mpolyd_offset(A, exp);
        if (off >= 0 && off < A->alloc) {
            A->coeffs[off] = B->coeffs[i];
            #if DEBUG_DIVISION
            if (i < 5 || i >= B->length - 5) {
                printf("[DEBUG] Term %ld: coeff=0x%02x, offset=%ld\n", i, B->coeffs[i], off);
            }
            #endif
        }
    }
    
    free(exp);
}

void gf28_mpolyd_to_mpoly_fast(gf28_mpoly_t A, const gf28_mpolyd_t B,
                              const gf28_mpoly_ctx_t ctx) {
    slong nvars = ctx->minfo->nvars;
    
    /* Pre-calculate required bits */
    flint_bitcnt_t bits_needed = MPOLY_MIN_BITS;
    for (slong i = 0; i < nvars; i++) {
        if (B->deg_bounds[i] > 0) {
            slong bits_i = FLINT_BIT_COUNT(B->deg_bounds[i] - 1);
            bits_needed = FLINT_MAX(bits_needed, bits_i);
        }
    }
    if (bits_needed < 16) bits_needed = 16;
    bits_needed = mpoly_fix_bits(bits_needed, ctx->minfo);
    
    slong N = mpoly_words_per_exp(bits_needed, ctx->minfo);
    
    /* First pass: count nonzero terms */
    slong nonzero_count = 0;
    for (slong off = 0; off < B->alloc; off++) {
        if (B->coeffs[off] != 0) {
            nonzero_count++;
        }
    }
    
    if (nonzero_count == 0) {
        gf28_mpoly_zero(A, ctx);
        return;
    }
    
    /* Pre-allocate exact space needed */
    gf28_mpoly_fit_length_reset_bits(A, nonzero_count, bits_needed, ctx);
    A->length = 0;
    
    /* Temporary arrays for batch processing */
    ulong *exp_buffer = (ulong *)malloc(nvars * sizeof(ulong));
    ulong *packed_exp_buffer = (ulong *)malloc(N * sizeof(ulong));
    
    /* Build terms directly in reverse offset order for LEX ordering */
    for (slong off = B->alloc - 1; off >= 0; off--) {
        if (B->coeffs[off] != 0) {
            /* Convert offset to exponent vector */
            slong temp = off;
            for (slong i = 0; i < nvars; i++) {
                exp_buffer[i] = (ulong)(temp % B->deg_bounds[i]);
                temp /= B->deg_bounds[i];
            }
            
            /* Pack exponent directly */
            mpoly_set_monomial_ui(packed_exp_buffer, exp_buffer, bits_needed, ctx->minfo);
            
            /* Store directly without searching/inserting */
            A->coeffs[A->length] = B->coeffs[off];
            mpoly_monomial_set(A->exps + N * A->length, packed_exp_buffer, N);
            A->length++;
        }
    }
    
    free(exp_buffer);
    free(packed_exp_buffer);
    
    #if DEBUG_DIVISION
    printf("[DEBUG] gf28_mpolyd_to_mpoly_fast: built %ld terms\n", A->length);
    #endif
}

void gf28_mpolyd_to_mpoly_univariate(gf28_mpoly_t A, const gf28_mpolyd_t B,
                                    const gf28_mpoly_ctx_t ctx) {
    if (ctx->minfo->nvars != 1) {
        gf28_mpolyd_to_mpoly_fast(A, B, ctx);
        return;
    }
    
    /* For univariate, we can be extremely fast */
    slong degree_bound = B->deg_bounds[0];
    slong nonzero_count = 0;
    
    /* Count and find highest degree */
    slong highest_degree = -1;
    for (slong i = 0; i < degree_bound; i++) {
        if (B->coeffs[i] != 0) {
            nonzero_count++;
            highest_degree = i;
        }
    }
    
    if (nonzero_count == 0) {
        gf28_mpoly_zero(A, ctx);
        return;
    }
    
    /* Calculate required bits */
    flint_bitcnt_t bits_needed = FLINT_BIT_COUNT(highest_degree) + 1;
    bits_needed = FLINT_MAX(bits_needed, MPOLY_MIN_BITS);
    bits_needed = mpoly_fix_bits(bits_needed, ctx->minfo);
    
    gf28_mpoly_fit_length_reset_bits(A, nonzero_count, bits_needed, ctx);
    
    A->length = 0;
    
    /* Process from highest to lowest degree for LEX order */
    for (slong deg = highest_degree; deg >= 0; deg--) {
        if (B->coeffs[deg] != 0) {
            A->coeffs[A->length] = B->coeffs[deg];
            
            /* For univariate, exponent packing is trivial */
            if (mpoly_words_per_exp(bits_needed, ctx->minfo) == 1) {
                A->exps[A->length] = (ulong)deg;
            } else {
                ulong exp[1] = {(ulong)deg};
                mpoly_set_monomial_ui(A->exps + A->length * mpoly_words_per_exp(bits_needed, ctx->minfo), 
                                     exp, bits_needed, ctx->minfo);
            }
            
            A->length++;
        }
    }
    
    printf("[UNIVARIATE] Ultra-fast conversion completed: %ld terms\n", A->length);
}

/* Main optimized function with automatic dispatch */
void gf28_mpolyd_to_mpoly(gf28_mpoly_t A, const gf28_mpolyd_t B,
                         const gf28_mpoly_ctx_t ctx) {
    if (ctx->minfo->nvars == 1) {
        gf28_mpolyd_to_mpoly_univariate(A, B, ctx);
    } else {
        gf28_mpolyd_to_mpoly_fast(A, B, ctx);
    }
}

void gf28_mpolyd_divrem_multivar(gf28_mpolyd_t Q, gf28_mpolyd_t R,
                                const gf28_mpolyd_t A, const gf28_mpolyd_t B,
                                const gf28_mpoly_ctx_t ctx) {
    gf28_mpolyd_divrem_univar(Q, R, A, B);
}

int gf28_mpolyd_set_degbounds(gf28_mpolyd_t A, const slong *bounds) {
    slong size = 1;
    
    for (slong i = 0; i < A->nvars; i++) {
        A->deg_bounds[i] = bounds[i];
        if (bounds[i] <= 0) return 0;
        
        if (size > WORD_MAX / bounds[i]) return 0;
        size *= bounds[i];
    }
    
    if (size > (1L << 28)) return 0;
    
    A->alloc = size;
    A->coeffs = (uint8_t *)calloc(size, sizeof(uint8_t));
    
    #if DEBUG_DIVISION
    printf("[DEBUG] gf28_mpolyd_set_degbounds: allocated %ld coefficients\n", size);
    #endif
    
    return A->coeffs != NULL;
}

int gf28_mpoly_divides_monomial(gf28_mpoly_t Q, const gf28_mpoly_t A, 
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
        return 0;
    }
    
    coeff_b_inv = gf28_inv(coeff_b);
    gf28_mpoly_get_monomial_exp(exp_b, B, ctx);
    
    flint_bitcnt_t bits = FLINT_MAX(A->bits, B->bits);
    if (bits < 16) bits = 16;
    
    gf28_mpoly_fit_length_reset_bits(Q, A->length, bits, ctx);
    Q->length = 0;
    N = mpoly_words_per_exp(bits, ctx->minfo);
    
    for (slong i = 0; i < A->length; i++) {
        mpoly_get_monomial_ui(exp_a, A->exps + N*i, A->bits, ctx->minfo);
        
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
            gf28_mpoly_zero(Q, ctx);
            return 0;
        }
        
        uint8_t coeff_q = gf28_mul(A->coeffs[i], coeff_b_inv);
        
        if (coeff_q != 0) {
            mpoly_set_monomial_ui(Q->exps + N*Q->length, exp_q, bits, ctx->minfo);
            Q->coeffs[Q->length] = coeff_q;
            Q->length++;
        }
    }
    
    free(exp_b); free(exp_a); free(exp_q);
    return 1;
}

int gf28_mpoly_divides_dense(gf28_mpoly_t Q, const gf28_mpoly_t A, 
                            const gf28_mpoly_t B, const gf28_mpoly_ctx_t ctx) {
    slong nvars = ctx->minfo->nvars;
    int success = 0;
    
    #if DEBUG_DIVISION
    printf("\n[DEBUG] gf28_mpoly_divides_dense: A has %ld terms, B has %ld terms\n", 
           A->length, B->length);
    #endif
    
    slong *degs_A = (slong *)malloc(nvars * sizeof(slong));
    slong *degs_B = (slong *)malloc(nvars * sizeof(slong));
    slong *bounds = (slong *)malloc(nvars * sizeof(slong));
    
    mpoly_degrees_si(degs_A, A->exps, A->length, A->bits, ctx->minfo);
    mpoly_degrees_si(degs_B, B->exps, B->length, B->bits, ctx->minfo);
    
    #if DEBUG_DIVISION
    printf("[DEBUG] Degrees: A=(");
    for (slong i = 0; i < nvars; i++) printf("%ld%s", degs_A[i], i < nvars-1 ? "," : "");
    printf("), B=(");
    for (slong i = 0; i < nvars; i++) printf("%ld%s", degs_B[i], i < nvars-1 ? "," : "");
    printf(")\n");
    #endif
    
    slong main_var_a, main_var_b;
    int is_univar_a = gf28_mpoly_is_univariate(A, ctx, &main_var_a);
    int is_univar_b = gf28_mpoly_is_univariate(B, ctx, &main_var_b);
    
    for (slong i = 0; i < nvars; i++) {
        if (degs_A[i] < degs_B[i]) {
            #if DEBUG_DIVISION
            printf("[DEBUG] Division impossible: deg_A[%ld] < deg_B[%ld]\n", i, i);
            #endif
            free(degs_A); free(degs_B); free(bounds);
            gf28_mpoly_zero(Q, ctx);
            return 0;
        }
        bounds[i] = degs_A[i] + 1;
    }
    
    gf28_mpolyd_t Ad, Bd, Qd, Rd;
    gf28_mpolyd_init(Ad, nvars);
    gf28_mpolyd_init(Bd, nvars);
    gf28_mpolyd_init(Qd, nvars);
    gf28_mpolyd_init(Rd, nvars);
    
    if (!gf28_mpolyd_set_degbounds(Ad, bounds) ||
        !gf28_mpolyd_set_degbounds(Bd, bounds) ||
        !gf28_mpolyd_set_degbounds(Qd, bounds) ||
        !gf28_mpolyd_set_degbounds(Rd, bounds)) {
        #if DEBUG_DIVISION
        printf("[DEBUG] Failed to allocate dense arrays\n");
        #endif
        goto cleanup;
    }
    
    gf28_mpoly_to_mpolyd(Ad, A, ctx);
    gf28_mpoly_to_mpolyd(Bd, B, ctx);
    
    /* TIMING START */
    double start_time = get_wall_time_8();
    printf("begin univar divrem\n");
    
    gf28_mpolyd_divrem_univar(Qd, Rd, Ad, Bd);
    
    printf("end univar divrem\n");
    double end_time = get_wall_time_8();
    printf("Univariate division time: %.6f seconds\n", end_time - start_time);
    /* TIMING END */
    
    if (gf28_mpolyd_is_zero(Rd)) {
        #if DEBUG_DIVISION
        printf("[DEBUG] Remainder is zero, division successful\n");
        #endif            
        if (is_univar_a && is_univar_b && main_var_a == main_var_b) {
            printf("gf28_mpolyd_to_mpoly_univariate\n");
            gf28_mpolyd_to_mpoly_univariate(Q, Qd, ctx);
        } else {
            gf28_mpolyd_to_mpoly_fast(Q, Qd, ctx);
        }
        success = 1;
    } else {
        #if DEBUG_DIVISION
        printf("[DEBUG] Remainder is nonzero, division failed\n");
        #endif
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

int gf28_mpoly_divides_flint_univar(gf28_mpoly_t Q, const gf28_mpoly_t A, 
                                   const gf28_mpoly_t B, const gf28_mpoly_ctx_t ctx,
                                   slong main_var) {
    #if DEBUG_DIVISION
    printf("[DEBUG] Using FLINT univariate division for variable %ld\n", main_var);
    printf("[DEBUG] A degree: %ld, B degree: %ld\n", 
           gf28_mpoly_univariate_degree(A, ctx, main_var),
           gf28_mpoly_univariate_degree(B, ctx, main_var));
    #endif
    
    /* Initialize FLINT context for GF(2^8) */
    nmod_poly_t modulus;
    fq_nmod_ctx_t fq_ctx;
    
    nmod_poly_init(modulus, 2);
    nmod_poly_set_coeff_ui(modulus, 0, 1);  // x^8 + x^4 + x^3 + x + 1
    nmod_poly_set_coeff_ui(modulus, 2, 1);
    nmod_poly_set_coeff_ui(modulus, 3, 1);
    nmod_poly_set_coeff_ui(modulus, 4, 1);
    nmod_poly_set_coeff_ui(modulus, 8, 1);
    
    fq_nmod_ctx_init_modulus(fq_ctx, modulus, "x");
    
    /* Initialize conversion if needed */
    init_gf28_conversion(fq_ctx);
    
    /* Convert to GF(2^8) polynomials */
    gf28_poly_t poly_a, poly_b, poly_q;
    gf28_poly_init(poly_a);
    gf28_poly_init(poly_b);
    gf28_poly_init(poly_q);
    
    gf28_mpoly_to_gf28_poly_univar(poly_a, A, ctx, main_var);
    gf28_mpoly_to_gf28_poly_univar(poly_b, B, ctx, main_var);
    
    #if DEBUG_DIVISION
    printf("[DEBUG] Converted to univariate: poly_a length=%ld, poly_b length=%ld\n",
           poly_a->length, poly_b->length);
    #endif
    
    /* Convert to FLINT format */
    fq_nmod_poly_t fa, fb, fq;
    fq_nmod_poly_init(fa, fq_ctx);
    fq_nmod_poly_init(fb, fq_ctx);
    fq_nmod_poly_init(fq, fq_ctx);
    
    gf28_poly_to_fq_nmod_poly(fa, poly_a, fq_ctx);
    gf28_poly_to_fq_nmod_poly(fb, poly_b, fq_ctx);
    
    #if DEBUG_DIVISION
    printf("[DEBUG] FLINT polynomials: fa degree=%ld, fb degree=%ld\n",
           fq_nmod_poly_degree(fa, fq_ctx), fq_nmod_poly_degree(fb, fq_ctx));
    #endif
    
    /* Use FLINT's optimized division */
    int divides = fq_nmod_poly_divides(fq, fa, fb, fq_ctx);
    
    #if DEBUG_DIVISION
    printf("[DEBUG] FLINT division result: %s\n", divides ? "exact" : "not exact");
    if (divides) {
        printf("[DEBUG] Quotient degree: %ld\n", fq_nmod_poly_degree(fq, fq_ctx));
    }
    #endif
    
    if (divides) {
        /* Convert result back */
        fq_nmod_poly_to_gf28_poly(poly_q, fq, fq_ctx);
        gf28_poly_to_gf28_mpoly_univar(Q, poly_q, ctx, main_var);
        
        #if DEBUG_DIVISION
        printf("[DEBUG] Final quotient has %ld terms\n", Q->length);
        #endif
    } else {
        gf28_mpoly_zero(Q, ctx);
    }
    
    /* Cleanup */
    gf28_poly_clear(poly_a);
    gf28_poly_clear(poly_b);
    gf28_poly_clear(poly_q);
    fq_nmod_poly_clear(fa, fq_ctx);
    fq_nmod_poly_clear(fb, fq_ctx);
    fq_nmod_poly_clear(fq, fq_ctx);
    fq_nmod_ctx_clear(fq_ctx);
    nmod_poly_clear(modulus);
    
    return divides;
}

/* Intelligent multiplication function */
int gf28_mpoly_mul(gf28_mpoly_t res, const gf28_mpoly_t a, const gf28_mpoly_t b, 
                   const gf28_mpoly_ctx_t ctx) {
    /* Handle zero polynomials */
    if (a->length == 0 || b->length == 0) {
        gf28_mpoly_zero(res, ctx);
        return 1;
    }
    
    /* Check if both polynomials are effectively univariate */
    slong main_var_a, main_var_b;
    int is_univar_a = gf28_mpoly_is_univariate(a, ctx, &main_var_a);
    int is_univar_b = gf28_mpoly_is_univariate(b, ctx, &main_var_b);

    if (is_univar_a && is_univar_b && main_var_a == main_var_b) {
        /* Both are univariate in the same variable */
        slong deg_a = gf28_mpoly_univariate_degree(a, ctx, main_var_a);
        slong deg_b = gf28_mpoly_univariate_degree(b, ctx, main_var_b);
        slong max_deg = FLINT_MAX(deg_a, deg_b);
        
        /* Convert to univariate polynomials */
        gf28_poly_t poly_a, poly_b, poly_res;
        gf28_poly_init(poly_a);
        gf28_poly_init(poly_b);
        gf28_poly_init(poly_res);
        
        gf28_mpoly_to_gf28_poly_univar(poly_a, a, ctx, main_var_a);
        gf28_mpoly_to_gf28_poly_univar(poly_b, b, ctx, main_var_b);
        
        if (max_deg < 10) {
            gf28_poly_mul_schoolbook(poly_res, poly_a, poly_b);
        } else if (max_deg <= 10000) {
            gf28_poly_mul_karatsuba(poly_res, poly_a, poly_b);
        } else {
            gf28_poly_clear(poly_a);
            gf28_poly_clear(poly_b);
            gf28_poly_clear(poly_res);
            
            gf28_mpoly_mul_flint_univar(res, a, b, ctx, main_var_a);
            return 1;
        }
        
        /* Convert result back to multivariate */
        gf28_poly_to_gf28_mpoly_univar(res, poly_res, ctx, main_var_a);
        gf28_poly_clear(poly_a);
        gf28_poly_clear(poly_b);
        gf28_poly_clear(poly_res);
        return 1;
    } else {
        /* Multivariate case: use array multiplication */
        return gf28_mpoly_mul_array(res, a, b, ctx);
    }
}
    
int gf28_mpoly_divides(gf28_mpoly_t Q, const gf28_mpoly_t A, 
                      const gf28_mpoly_t B, const gf28_mpoly_ctx_t ctx) {
    #if DEBUG_DIVISION
    printf("\n[DEBUG] === Starting gf28_mpoly_divides ===\n");
    printf("[DEBUG] A: %ld terms, B: %ld terms\n", A->length, B->length);
    #endif
    
    /* Handle zero polynomials */
    if (B->length == 0) {
        if (A->length == 0) {
            gf28_mpoly_zero(Q, ctx);
            #if DEBUG_DIVISION
            printf("[DEBUG] Both polynomials are zero: exact division\n");
            #endif
            return 1;
        }
        #if DEBUG_DIVISION
        printf("[DEBUG] Division by zero: not exact\n");
        #endif
        return 0;
    }
    
    if (A->length == 0) {
        gf28_mpoly_zero(Q, ctx);
        #if DEBUG_DIVISION
        printf("[DEBUG] Zero dividend: exact division (quotient = 0)\n");
        #endif
        return 1;
    }
    
    /* Check if both polynomials are effectively univariate */
    slong main_var_a, main_var_b;
    int is_univar_a = gf28_mpoly_is_univariate(A, ctx, &main_var_a);
    int is_univar_b = gf28_mpoly_is_univariate(B, ctx, &main_var_b);
    
    #if DEBUG_DIVISION
    printf("[DEBUG] Univariate check: A=%s (var %ld), B=%s (var %ld)\n",
           is_univar_a ? "yes" : "no", main_var_a,
           is_univar_b ? "yes" : "no", main_var_b);
    #endif
    
    if (0 && is_univar_a && is_univar_b && main_var_a == main_var_b) {
        /* Both are univariate in the same variable - use FLINT optimization */
        slong deg_a = gf28_mpoly_univariate_degree(A, ctx, main_var_a);
        slong deg_b = gf28_mpoly_univariate_degree(B, ctx, main_var_b);
        
        #if DEBUG_DIVISION
        printf("[DEBUG] Using FLINT path: degrees A=%ld, B=%ld\n", deg_a, deg_b);
        #endif
        
        if (deg_a < deg_b) {
            #if DEBUG_DIVISION
            printf("[DEBUG] Degree of A < degree of B: division not exact\n");
            #endif
            gf28_mpoly_zero(Q, ctx);
            return 0;
        }
        
        return gf28_mpoly_divides_flint_univar(Q, A, B, ctx, main_var_a);
    }
    
    /* Handle monomial divisors */
    if (gf28_mpoly_is_monomial(B)) {
        #if DEBUG_DIVISION
        printf("[DEBUG] Using monomial division\n");
        #endif
        return gf28_mpoly_divides_monomial(Q, A, B, ctx);
    }
    
    /* Fallback to dense multivariate division */
    #if DEBUG_DIVISION
    printf("[DEBUG] Falling back to dense multivariate division\n");
    #endif
    printf("gf28_mpoly_divides_dense\n");
    return gf28_mpoly_divides_dense(Q, A, B, ctx);
}

/* ============================================================================
   ADDITIONAL FUNCTIONS
   ============================================================================ */

void gf28_mpoly_mul_simple(gf28_mpoly_t res, const gf28_mpoly_t a, 
                          const gf28_mpoly_t b, const gf28_mpoly_ctx_t ctx) {
    slong nvars = ctx->minfo->nvars;
    slong N = mpoly_words_per_exp(a->bits, ctx->minfo);
    ulong *exp_a, *exp_b, *exp_sum;
    
    gf28_mpoly_zero(res, ctx);
    
    if (a->length == 0 || b->length == 0) return;
    
    exp_a = (ulong *)malloc(nvars * sizeof(ulong));
    exp_b = (ulong *)malloc(nvars * sizeof(ulong));
    exp_sum = (ulong *)malloc(nvars * sizeof(ulong));
    
    for (slong i = 0; i < a->length; i++) {
        mpoly_get_monomial_ui(exp_a, a->exps + N*i, a->bits, ctx->minfo);
        
        for (slong j = 0; j < b->length; j++) {
            mpoly_get_monomial_ui(exp_b, b->exps + N*j, b->bits, ctx->minfo);
            
            for (slong k = 0; k < nvars; k++) {
                exp_sum[k] = exp_a[k] + exp_b[k];
            }
            
            uint8_t coeff = gf28_mul(a->coeffs[i], b->coeffs[j]);
            
            if (coeff != 0) {
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
                        free(exp_k);
                        break;
                    }
                    free(exp_k);
                }
                
                uint8_t new_coeff = gf28_add(current, coeff);
                gf28_mpoly_set_coeff_ui_ui(res, new_coeff, exp_sum, ctx);
            }
        }
    }
    
    free(exp_a);
    free(exp_b);
    free(exp_sum);
}

int gf28_mpoly_equal(const gf28_mpoly_t A, const gf28_mpoly_t B, 
                    const gf28_mpoly_ctx_t ctx) {
    if (A->length != B->length) return 0;
    
    slong nvars = ctx->minfo->nvars;
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    
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

void gf28_mpoly_randtest(gf28_mpoly_t poly, flint_rand_t state,
                        slong length, slong exp_bound, 
                        const gf28_mpoly_ctx_t ctx) {
    slong nvars = ctx->minfo->nvars;
    ulong *exp = (ulong *)malloc(nvars * sizeof(ulong));
    
    gf28_mpoly_zero(poly, ctx);
    
    for (slong i = 0; i < length; i++) {
        for (slong j = 0; j < nvars; j++) {
            exp[j] = n_randint(state, exp_bound);
        }
        
        uint8_t c = n_randint(state, 255) + 1;
        
        gf28_mpoly_set_coeff_ui_ui(poly, c, exp, ctx);
    }
    
    free(exp);
}

void gf28_mpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t res, const gf28_mpoly_t poly,
                                const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx) {
    fq_nmod_mpoly_zero(res, fq_mpoly_ctx);
    
    if (poly->length == 0) return;
    
    slong N = mpoly_words_per_exp(poly->bits, fq_mpoly_ctx->minfo);
    slong nvars = fq_mpoly_ctx->minfo->nvars;
    
    for (slong i = 0; i < poly->length; i++) {
        if (poly->coeffs[i] != 0) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, fqctx);
            
            gf28_elem_to_fq_nmod(coeff, poly->coeffs[i], fqctx);
            
            ulong *exp = (ulong *) flint_malloc(nvars * sizeof(ulong));
            mpoly_get_monomial_ui(exp, poly->exps + N*i, poly->bits, fq_mpoly_ctx->minfo);
            
            fq_nmod_mpoly_set_coeff_fq_nmod_ui(res, coeff, exp, fq_mpoly_ctx);
            
            flint_free(exp);
            fq_nmod_clear(coeff, fqctx);
        }
    }
}