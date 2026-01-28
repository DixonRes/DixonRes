/* gf2n_mpoly.c - Unified GF(2^n) Multivariate Polynomial Implementation */

#include "gf2n_mpoly.h"
#include <emmintrin.h>  // For SSE2 instructions

/* ============================================================================
   MACRO TEMPLATE SYSTEM FOR CORE ALGORITHMS
   ============================================================================ */

/*
 * IMPORTANT: Dense structures (gfXX_mpolyd_t) are ALREADY DEFINED in gf2n_mpoly.h
 * This file only provides the IMPLEMENTATIONS, not the type definitions.
 */

#define DEFINE_MPOLY_ARRAY_MUL(FIELD, COEFF_T, IS_ZERO, ADD_OP, MUL_OP, ARRAY_LIMIT) \
\
/* Array multiplication helper */ \
static void _##FIELD##_mpoly_addmul_array1_safe( \
    COEFF_T *poly1, slong array_size, \
    const COEFF_T *poly2, const ulong *exp2, slong len2, \
    const COEFF_T *poly3, const ulong *exp3, slong len3) \
{ \
    for (slong ii = 0; ii < len2 + BLOCK; ii += BLOCK) { \
        for (slong jj = 0; jj < len3 + BLOCK; jj += BLOCK) { \
            for (slong i = ii; i < FLINT_MIN(ii + BLOCK, len2); i++) { \
                slong offset2 = (slong)exp2[i]; \
                if (offset2 >= array_size) continue; \
                \
                COEFF_T *c2 = poly1 + offset2; \
                \
                if (!(IS_ZERO(poly2[i]))) { \
                    for (slong j = jj; j < FLINT_MIN(jj + BLOCK, len3); j++) { \
                        slong offset3 = (slong)exp3[j]; \
                        slong total_offset = offset2 + offset3; \
                        if (total_offset >= array_size) continue; \
                        \
                        COEFF_T prod = MUL_OP(poly2[i], poly3[j]); \
                        c2[offset3] = ADD_OP(c2[offset3], prod); \
                    } \
                } \
            } \
        } \
    } \
} \
\
/* LEX ordering unpacking */ \
static slong FIELD##_mpoly_append_array_LEX_safe( \
    FIELD##_mpoly_t P, slong Plen, COEFF_T *coeff_array, \
    const ulong *mults, slong num, slong array_size, \
    slong top, const FIELD##_mpoly_ctx_t ctx) \
{ \
    slong topmult = num == 0 ? 1 : mults[num - 1]; \
    slong lastd = topmult - 1; \
    slong reset = array_size / topmult; \
    slong counter = reset; \
    ulong startexp = ((ulong)top << (P->bits * num)) + \
                     ((ulong)lastd << (P->bits * (num - 1))); \
    \
    for (slong off = array_size - 1; off >= 0; off--) { \
        if (!(IS_ZERO(coeff_array[off]))) { \
            COEFF_T coeff = coeff_array[off]; \
            slong d = off; \
            ulong exp = startexp; \
            \
            for (slong j = 0; j + 1 < num; j++) { \
                ulong exp_j = d % mults[j]; \
                exp += exp_j << (P->bits * j); \
                d = d / mults[j]; \
            } \
            \
            if (Plen >= P->coeffs_alloc) { \
                slong new_alloc = FLINT_MAX(Plen + 1, 2 * P->coeffs_alloc); \
                P->coeffs = (COEFF_T *)flint_realloc(P->coeffs, \
                    new_alloc * sizeof(COEFF_T)); \
                P->coeffs_alloc = new_alloc; \
            } \
            if (Plen >= P->exps_alloc) { \
                slong new_alloc = FLINT_MAX(Plen + 1, 2 * P->exps_alloc); \
                P->exps = (ulong *)flint_realloc(P->exps, \
                    new_alloc * sizeof(ulong)); \
                P->exps_alloc = new_alloc; \
            } \
            \
            P->exps[Plen] = exp; \
            P->coeffs[Plen] = coeff; \
            Plen++; \
        } \
        \
        counter--; \
        if (counter <= 0) { \
            counter = reset; \
            lastd--; \
            startexp -= UWORD(1) << (P->bits * (num - 1)); \
        } \
    } \
    \
    return Plen; \
} \
\
/* Chunked LEX multiplication */ \
static void _##FIELD##_mpoly_mul_array_chunked_LEX( \
    FIELD##_mpoly_t P, \
    const FIELD##_mpoly_t A, \
    const FIELD##_mpoly_t B, \
    const ulong *mults, \
    const FIELD##_mpoly_ctx_t ctx) \
{ \
    slong num = ctx->minfo->nfields - 1; \
    slong array_size = 1; \
    for (slong i = 0; i < num; i++) { \
        array_size *= mults[i]; \
    } \
    \
    slong Al = 1 + (slong)(A->exps[0] >> (A->bits * num)); \
    slong Bl = 1 + (slong)(B->exps[0] >> (B->bits * num)); \
    \
    slong *Amain = (slong *)flint_malloc((Al + 1) * sizeof(slong)); \
    slong *Bmain = (slong *)flint_malloc((Bl + 1) * sizeof(slong)); \
    ulong *Apexp = (ulong *)flint_malloc(A->length * sizeof(ulong)); \
    ulong *Bpexp = (ulong *)flint_malloc(B->length * sizeof(ulong)); \
    \
    mpoly_main_variable_split_LEX(Amain, Apexp, A->exps, Al, \
        A->length, mults, num, A->bits); \
    mpoly_main_variable_split_LEX(Bmain, Bpexp, B->exps, Bl, \
        B->length, mults, num, B->bits); \
    \
    slong Pl = Al + Bl - 1; \
    slong Plen = 0; \
    \
    COEFF_T *coeff_array = (COEFF_T *)calloc(array_size, sizeof(COEFF_T)); \
    if (!coeff_array) { \
        flint_free(Amain); flint_free(Bmain); \
        flint_free(Apexp); flint_free(Bpexp); \
        return; \
    } \
    \
    for (slong Pi = 0; Pi < Pl; Pi++) { \
        memset(coeff_array, 0, array_size * sizeof(COEFF_T)); \
        \
        for (slong i = 0, j = Pi; i < Al && j >= 0; i++, j--) { \
            if (j < Bl) { \
                _##FIELD##_mpoly_addmul_array1_safe(coeff_array, array_size, \
                    A->coeffs + Amain[i], Apexp + Amain[i], \
                    Amain[i + 1] - Amain[i], \
                    B->coeffs + Bmain[j], Bpexp + Bmain[j], \
                    Bmain[j + 1] - Bmain[j]); \
            } \
        } \
        \
        Plen = FIELD##_mpoly_append_array_LEX_safe(P, Plen, coeff_array, \
            mults, num, array_size, Pl - Pi - 1, ctx); \
    } \
    \
    P->length = Plen; \
    \
    free(coeff_array); \
    flint_free(Amain); flint_free(Bmain); \
    flint_free(Apexp); flint_free(Bpexp); \
} \
\
/* Main LEX multiplication */ \
static int _##FIELD##_mpoly_mul_array_LEX( \
    FIELD##_mpoly_t A, \
    const FIELD##_mpoly_t B, fmpz *maxBfields, \
    const FIELD##_mpoly_t C, fmpz *maxCfields, \
    const FIELD##_mpoly_ctx_t ctx) \
{ \
    FLINT_ASSERT(ctx->minfo->nvars > 0); \
    FLINT_ASSERT(B->length != 0 && C->length != 0); \
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX); \
    \
    slong i = ctx->minfo->nfields - 1; \
    FLINT_ASSERT(fmpz_fits_si(maxBfields + i)); \
    FLINT_ASSERT(fmpz_fits_si(maxCfields + i)); \
    \
    ulong *mults = (ulong *)flint_malloc(ctx->minfo->nfields * sizeof(ulong)); \
    mults[i] = 1 + fmpz_get_ui(maxBfields + i) + fmpz_get_ui(maxCfields + i); \
    ulong max = mults[i]; \
    \
    if (((slong)mults[i]) <= 0) { \
        flint_free(mults); \
        return 0; \
    } \
    \
    slong array_size = WORD(1); \
    for (i--; i >= 0; i--) { \
        mults[i] = 1 + fmpz_get_ui(maxBfields + i) + \
                   fmpz_get_ui(maxCfields + i); \
        max |= mults[i]; \
        ulong hi; \
        umul_ppmm(hi, array_size, array_size, mults[i]); \
        if (hi != 0 || (slong)mults[i] <= 0 || array_size <= 0) { \
            flint_free(mults); \
            return 0; \
        } \
    } \
    \
    if (array_size > ARRAY_LIMIT) { \
        flint_free(mults); \
        return 0; \
    } \
    \
    flint_bitcnt_t exp_bits = FLINT_MAX(MPOLY_MIN_BITS, \
        FLINT_BIT_COUNT(max) + 1); \
    exp_bits = mpoly_fix_bits(exp_bits, ctx->minfo); \
    \
    if (1 != mpoly_words_per_exp(exp_bits, ctx->minfo)) { \
        flint_free(mults); \
        return 0; \
    } \
    \
    if (A == B || A == C) { \
        FIELD##_mpoly_t T; \
        FIELD##_mpoly_init(T, ctx); \
        T->bits = exp_bits; \
        _##FIELD##_mpoly_mul_array_chunked_LEX(T, C, B, mults, ctx); \
        /* Swap using memcpy to handle array types */ \
        FIELD##_mpoly_struct temp_swap = *T; \
        *T = *A; \
        *A = temp_swap; \
        FIELD##_mpoly_clear(T, ctx); \
    } else { \
        A->bits = exp_bits; \
        _##FIELD##_mpoly_mul_array_chunked_LEX(A, C, B, mults, ctx); \
    } \
    \
    flint_free(mults); \
    return 1; \
}

/* ============================================================================
   MACRO INVOCATIONS FOR EACH FIELD
   ============================================================================ */

/* GF(2^8): Simple XOR operations */
#define GF28_IS_ZERO(x) ((x) == 0)
#define GF28_ADD(a, b) ((a) ^ (b))
#define GF28_MUL(a, b) (gf28_mul((a), (b)))

DEFINE_MPOLY_ARRAY_MUL(gf28, uint8_t, GF28_IS_ZERO, GF28_ADD, GF28_MUL, (1L << 28))

/* GF(2^16): Simple XOR operations */
#define GF216_IS_ZERO(x) ((x) == 0)
#define GF216_ADD(a, b) ((a) ^ (b))
#define GF216_MUL(a, b) (gf216_mul((a), (b)))

DEFINE_MPOLY_ARRAY_MUL(gf216, uint16_t, GF216_IS_ZERO, GF216_ADD, GF216_MUL, (1L << 27))

/* GF(2^32): Struct operations */
#define GF232_IS_ZERO(x) (gf232_is_zero(&(x)))
#define GF232_ADD(a, b) (gf232_add(&(a), &(b)))
#define GF232_MUL(a, b) (gf232_mul(&(a), &(b)))

DEFINE_MPOLY_ARRAY_MUL(gf232, gf232_t, GF232_IS_ZERO, GF232_ADD, GF232_MUL, (1L << 27))

/* GF(2^64): Struct operations */
#define GF264_IS_ZERO(x) (gf264_is_zero(&(x)))
#define GF264_ADD(a, b) (gf264_add(&(a), &(b)))
#define GF264_MUL(a, b) (gf264_mul(&(a), &(b)))

DEFINE_MPOLY_ARRAY_MUL(gf264, gf264_t, GF264_IS_ZERO, GF264_ADD, GF264_MUL, (1L << 26))

/* GF(2^128): Struct operations with smaller limit */
#define GF2128_IS_ZERO(x) (gf2128_is_zero(&(x)))
#define GF2128_ADD(a, b) (gf2128_add(&(a), &(b)))
#define GF2128_MUL(a, b) (gf2128_mul(&(a), &(b)))

DEFINE_MPOLY_ARRAY_MUL(gf2128, gf2128_t, GF2128_IS_ZERO, GF2128_ADD, GF2128_MUL, (1L << 26))

/* ============================================================================
   DIVISION IMPLEMENTATION (Uses existing structures from header)
   ============================================================================ */

#define DEFINE_MPOLY_DIVISION_IMPL(FIELD, COEFF_T, IS_ZERO, ADD_OP, INV_OP, MUL_OP, SET_COEFF_CALL) \
\
static void FIELD##_mpolyd_init(FIELD##_mpolyd_t A, slong nvars) { \
    A->coeffs = NULL; \
    A->alloc = 0; \
    A->deg_bounds = (slong *)calloc(nvars, sizeof(slong)); \
    A->nvars = nvars; \
} \
\
static void FIELD##_mpolyd_clear(FIELD##_mpolyd_t A) { \
    if (A->coeffs) free(A->coeffs); \
    if (A->deg_bounds) free(A->deg_bounds); \
} \
\
static slong FIELD##_mpolyd_offset(const FIELD##_mpolyd_t A, const ulong *exp) { \
    slong off = 0, stride = 1; \
    for (slong i = 0; i < A->nvars; i++) { \
        if (exp[i] >= (ulong)A->deg_bounds[i]) return -1; \
        off += (slong)exp[i] * stride; \
        stride *= A->deg_bounds[i]; \
    } \
    return off; \
} \
\
static int FIELD##_mpolyd_set_degbounds(FIELD##_mpolyd_t A, const slong *bounds) { \
    slong size = 1; \
    for (slong i = 0; i < A->nvars; i++) { \
        A->deg_bounds[i] = bounds[i]; \
        if (bounds[i] <= 0) return 0; \
        if (size > WORD_MAX / bounds[i]) return 0; \
        size *= bounds[i]; \
    } \
    if (size > (1L << 26)) return 0; \
    A->alloc = size; \
    A->coeffs = (COEFF_T *)calloc(size, sizeof(COEFF_T)); \
    return A->coeffs != NULL; \
} \
\
static void FIELD##_mpoly_to_mpolyd(FIELD##_mpolyd_t A, \
    const FIELD##_mpoly_t B, const FIELD##_mpoly_ctx_t ctx) \
{ \
    slong nvars = ctx->minfo->nvars; \
    slong N = mpoly_words_per_exp(B->bits, ctx->minfo); \
    ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong)); \
    \
    memset(A->coeffs, 0, A->alloc * sizeof(COEFF_T)); \
    \
    for (slong i = 0; i < B->length; i++) { \
        mpoly_get_monomial_ui(exp, B->exps + N*i, B->bits, ctx->minfo); \
        slong off = FIELD##_mpolyd_offset(A, exp); \
        if (off >= 0 && off < A->alloc) { \
            A->coeffs[off] = B->coeffs[i]; \
        } \
    } \
    \
    flint_free(exp); \
} \
\
static void FIELD##_mpolyd_to_mpoly(FIELD##_mpoly_t A, \
    const FIELD##_mpolyd_t B, const FIELD##_mpoly_ctx_t ctx) \
{ \
    slong nvars = ctx->minfo->nvars; \
    ulong *exp = (ulong *)calloc(nvars, sizeof(ulong)); \
    \
    flint_bitcnt_t bits_needed = MPOLY_MIN_BITS; \
    for (slong i = 0; i < nvars; i++) { \
        if (B->deg_bounds[i] > 0) { \
            slong bits_i = FLINT_BIT_COUNT(B->deg_bounds[i] - 1); \
            bits_needed = FLINT_MAX(bits_needed, bits_i); \
        } \
    } \
    bits_needed = mpoly_fix_bits(bits_needed, ctx->minfo); \
    \
    FIELD##_mpoly_zero(A, ctx); \
    A->bits = bits_needed; \
    \
    for (slong off = B->alloc - 1; off >= 0; off--) { \
        if (IS_ZERO(B->coeffs[off])) continue; \
        \
        slong temp = off; \
        for (slong i = 0; i < nvars; i++) { \
            exp[i] = (ulong)(temp % B->deg_bounds[i]); \
            temp /= B->deg_bounds[i]; \
        } \
        \
        SET_COEFF_CALL(A, B->coeffs[off], exp, ctx); \
    } \
    \
    free(exp); \
} \
\
static void FIELD##_mpolyd_divrem_univar(FIELD##_mpolyd_t Q, \
    FIELD##_mpolyd_t R, const FIELD##_mpolyd_t A, \
    const FIELD##_mpolyd_t B) \
{ \
    slong n = A->alloc; \
    memcpy(R->coeffs, A->coeffs, n * sizeof(COEFF_T)); \
    memset(Q->coeffs, 0, Q->alloc * sizeof(COEFF_T)); \
    \
    slong degB = -1; \
    for (slong i = n - 1; i >= 0; i--) { \
        if (!(IS_ZERO(B->coeffs[i]))) { \
            degB = i; \
            break; \
        } \
    } \
    if (degB < 0) return; \
    \
    COEFF_T lc_B_inv = INV_OP(B->coeffs[degB]); \
    \
    for (slong i = n - 1; i >= degB; i--) { \
        if (!(IS_ZERO(R->coeffs[i]))) { \
            COEFF_T q = MUL_OP(R->coeffs[i], lc_B_inv); \
            \
            if (i - degB < Q->alloc) { \
                Q->coeffs[i - degB] = q; \
            } \
            \
            for (slong j = 0; j <= degB && i - degB + j < n; j++) { \
                COEFF_T prod = MUL_OP(q, B->coeffs[j]); \
                R->coeffs[i - degB + j] = ADD_OP(R->coeffs[i - degB + j], prod); \
            } \
        } \
    } \
} \
\
static int FIELD##_mpolyd_is_zero(const FIELD##_mpolyd_t A) { \
    for (slong i = 0; i < A->alloc; i++) { \
        if (!(IS_ZERO(A->coeffs[i]))) return 0; \
    } \
    return 1; \
} \
\
static int FIELD##_mpoly_divides_dense(FIELD##_mpoly_t Q, \
    const FIELD##_mpoly_t A, const FIELD##_mpoly_t B, \
    const FIELD##_mpoly_ctx_t ctx) \
{ \
    slong nvars = ctx->minfo->nvars; \
    slong *degs_A = (slong *)flint_malloc(nvars * sizeof(slong)); \
    slong *degs_B = (slong *)flint_malloc(nvars * sizeof(slong)); \
    slong *bounds = (slong *)flint_malloc(nvars * sizeof(slong)); \
    \
    mpoly_degrees_si(degs_A, A->exps, A->length, A->bits, ctx->minfo); \
    mpoly_degrees_si(degs_B, B->exps, B->length, B->bits, ctx->minfo); \
    \
    for (slong i = 0; i < nvars; i++) { \
        if (degs_A[i] < degs_B[i]) { \
            flint_free(degs_A); flint_free(degs_B); flint_free(bounds); \
            FIELD##_mpoly_zero(Q, ctx); \
            return 0; \
        } \
        bounds[i] = degs_A[i] + 1; \
    } \
    \
    FIELD##_mpolyd_t Ad, Bd, Qd, Rd; \
    FIELD##_mpolyd_init(Ad, nvars); \
    FIELD##_mpolyd_init(Bd, nvars); \
    FIELD##_mpolyd_init(Qd, nvars); \
    FIELD##_mpolyd_init(Rd, nvars); \
    \
    int success = 0; \
    if (!FIELD##_mpolyd_set_degbounds(Ad, bounds) || \
        !FIELD##_mpolyd_set_degbounds(Bd, bounds) || \
        !FIELD##_mpolyd_set_degbounds(Qd, bounds) || \
        !FIELD##_mpolyd_set_degbounds(Rd, bounds)) { \
        goto cleanup; \
    } \
    \
    FIELD##_mpoly_to_mpolyd(Ad, A, ctx); \
    FIELD##_mpoly_to_mpolyd(Bd, B, ctx); \
    FIELD##_mpolyd_divrem_univar(Qd, Rd, Ad, Bd); \
    \
    if (FIELD##_mpolyd_is_zero(Rd)) { \
        FIELD##_mpolyd_to_mpoly(Q, Qd, ctx); \
        success = 1; \
    } else { \
        FIELD##_mpoly_zero(Q, ctx); \
        success = 0; \
    } \
    \
cleanup: \
    FIELD##_mpolyd_clear(Ad); \
    FIELD##_mpolyd_clear(Bd); \
    FIELD##_mpolyd_clear(Qd); \
    FIELD##_mpolyd_clear(Rd); \
    flint_free(degs_A); flint_free(degs_B); flint_free(bounds); \
    return success; \
}

/* Instantiate division for each field with correct set_coeff call */
DEFINE_MPOLY_DIVISION_IMPL(gf28, uint8_t, GF28_IS_ZERO, GF28_ADD, gf28_inv, GF28_MUL, gf28_mpoly_set_coeff_ui_ui)
DEFINE_MPOLY_DIVISION_IMPL(gf216, uint16_t, GF216_IS_ZERO, GF216_ADD, gf216_inv, GF216_MUL, gf216_mpoly_set_coeff_ui_ui)

/* For struct types, need wrapper macros */
#define GF232_INV(x) (gf232_inv(&(x)))
#define GF232_SET_COEFF(A, coeff, exp, ctx) gf232_mpoly_set_coeff_ui_ui(A, &(coeff), exp, ctx)
DEFINE_MPOLY_DIVISION_IMPL(gf232, gf232_t, GF232_IS_ZERO, GF232_ADD, GF232_INV, GF232_MUL, GF232_SET_COEFF)

#define GF264_INV(x) (gf264_inv(&(x)))
#define GF264_SET_COEFF(A, coeff, exp, ctx) gf264_mpoly_set_coeff_ui_ui(A, &(coeff), exp, ctx)
DEFINE_MPOLY_DIVISION_IMPL(gf264, gf264_t, GF264_IS_ZERO, GF264_ADD, GF264_INV, GF264_MUL, GF264_SET_COEFF)

#define GF2128_INV(x) (gf2128_inv(&(x)))
#define GF2128_SET_COEFF(A, coeff, exp, ctx) gf2128_mpoly_set_coeff_ui_ui(A, &(coeff), exp, ctx)
DEFINE_MPOLY_DIVISION_IMPL(gf2128, gf2128_t, GF2128_IS_ZERO, GF2128_ADD, GF2128_INV, GF2128_MUL, GF2128_SET_COEFF)

/* ============================================================================
   PUBLIC API IMPLEMENTATIONS - BASIC OPERATIONS
   ============================================================================ */

double get_wall_time(void) {
    struct timeval time;
    if (gettimeofday(&time, NULL)) return 0;
    return (double)time.tv_sec + (double)time.tv_usec * 0.000001;
}

/* Macro for basic operations - identical for all fields */
#define DEFINE_BASIC_OPS(FIELD, COEFF_T, IS_ZERO_CHECK) \
\
void FIELD##_mpoly_init(FIELD##_mpoly_t poly, const FIELD##_mpoly_ctx_t ctx) { \
    poly->coeffs = NULL; \
    poly->exps = NULL; \
    poly->length = 0; \
    poly->coeffs_alloc = 0; \
    poly->exps_alloc = 0; \
    poly->bits = MPOLY_MIN_BITS; \
} \
\
void FIELD##_mpoly_clear(FIELD##_mpoly_t poly, const FIELD##_mpoly_ctx_t ctx) { \
    if (poly->coeffs) flint_free(poly->coeffs); \
    if (poly->exps) flint_free(poly->exps); \
} \
\
void FIELD##_mpoly_ctx_init(FIELD##_mpoly_ctx_t ctx, slong nvars, const ordering_t ord) { \
    mpoly_ctx_init(ctx->minfo, nvars, ord); \
} \
\
void FIELD##_mpoly_ctx_clear(FIELD##_mpoly_ctx_t ctx) { \
    mpoly_ctx_clear(ctx->minfo); \
} \
\
void FIELD##_mpoly_zero(FIELD##_mpoly_t poly, const FIELD##_mpoly_ctx_t ctx) { \
    poly->length = 0; \
} \
\
void FIELD##_mpoly_set(FIELD##_mpoly_t res, const FIELD##_mpoly_t poly, \
    const FIELD##_mpoly_ctx_t ctx) \
{ \
    if (res == poly) return; \
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo); \
    if (res->coeffs_alloc < poly->length) { \
        res->coeffs = (COEFF_T *)flint_realloc(res->coeffs, \
            poly->length * sizeof(COEFF_T)); \
        res->coeffs_alloc = poly->length; \
    } \
    if (res->exps_alloc < N * poly->length) { \
        res->exps = (ulong *)flint_realloc(res->exps, \
            N * poly->length * sizeof(ulong)); \
        res->exps_alloc = N * poly->length; \
    } \
    memcpy(res->coeffs, poly->coeffs, poly->length * sizeof(COEFF_T)); \
    memcpy(res->exps, poly->exps, N * poly->length * sizeof(ulong)); \
    res->length = poly->length; \
    res->bits = poly->bits; \
} \
\
int FIELD##_mpoly_equal(const FIELD##_mpoly_t A, const FIELD##_mpoly_t B, \
    const FIELD##_mpoly_ctx_t ctx) \
{ \
    if (A->length != B->length) return 0; \
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo); \
    for (slong i = 0; i < A->length; i++) { \
        if (IS_ZERO_CHECK(A->coeffs[i], B->coeffs[i])) return 0; \
        if (!mpoly_monomial_equal(A->exps + N*i, B->exps + N*i, N)) return 0; \
    } \
    return 1; \
} \
\
int FIELD##_mpoly_mul_array(FIELD##_mpoly_t A, const FIELD##_mpoly_t B, \
    const FIELD##_mpoly_t C, const FIELD##_mpoly_ctx_t ctx) \
{ \
    if (B->length == 0 || C->length == 0) { \
        FIELD##_mpoly_zero(A, ctx); \
        return 1; \
    } \
    \
    fmpz *maxBfields = (fmpz *)flint_malloc(ctx->minfo->nfields * sizeof(fmpz)); \
    fmpz *maxCfields = (fmpz *)flint_malloc(ctx->minfo->nfields * sizeof(fmpz)); \
    for (slong i = 0; i < ctx->minfo->nfields; i++) { \
        fmpz_init(maxBfields + i); \
        fmpz_init(maxCfields + i); \
    } \
    \
    mpoly_max_fields_fmpz(maxBfields, B->exps, B->length, B->bits, ctx->minfo); \
    mpoly_max_fields_fmpz(maxCfields, C->exps, C->length, C->bits, ctx->minfo); \
    \
    int success = 0; \
    switch (ctx->minfo->ord) { \
        case ORD_LEX: \
            success = _##FIELD##_mpoly_mul_array_LEX(A, B, maxBfields, \
                C, maxCfields, ctx); \
            break; \
        default: \
            success = 0; \
    } \
    \
    for (slong i = 0; i < ctx->minfo->nfields; i++) { \
        fmpz_clear(maxBfields + i); \
        fmpz_clear(maxCfields + i); \
    } \
    flint_free(maxBfields); \
    flint_free(maxCfields); \
    \
    return success; \
} \
\
int FIELD##_mpoly_mul(FIELD##_mpoly_t res, const FIELD##_mpoly_t a, \
    const FIELD##_mpoly_t b, const FIELD##_mpoly_ctx_t ctx) \
{ \
    return FIELD##_mpoly_mul_array(res, a, b, ctx); \
} \
\
int FIELD##_mpoly_divides(FIELD##_mpoly_t Q, const FIELD##_mpoly_t A, \
    const FIELD##_mpoly_t B, const FIELD##_mpoly_ctx_t ctx) \
{ \
    if (B->length == 0) return A->length == 0; \
    if (A->length == 0) { \
        FIELD##_mpoly_zero(Q, ctx); \
        return 1; \
    } \
    return FIELD##_mpoly_divides_dense(Q, A, B, ctx); \
}

/* Instantiate basic operations with appropriate equality checks */
#define GF28_EQ_CHECK(a, b) ((a) != (b))
#define GF216_EQ_CHECK(a, b) ((a) != (b))
#define GF232_EQ_CHECK(a, b) (!gf232_equal(&(a), &(b)))
#define GF264_EQ_CHECK(a, b) (!gf264_equal(&(a), &(b)))
#define GF2128_EQ_CHECK(a, b) (!gf2128_equal(&(a), &(b)))

DEFINE_BASIC_OPS(gf28, uint8_t, GF28_EQ_CHECK)
DEFINE_BASIC_OPS(gf216, uint16_t, GF216_EQ_CHECK)
DEFINE_BASIC_OPS(gf232, gf232_t, GF232_EQ_CHECK)
DEFINE_BASIC_OPS(gf264, gf264_t, GF264_EQ_CHECK)
DEFINE_BASIC_OPS(gf2128, gf2128_t, GF2128_EQ_CHECK)

/* ============================================================================
   COEFFICIENT SETTING - Field-specific implementations
   ============================================================================ */

void gf28_mpoly_set_coeff_ui_ui(gf28_mpoly_t poly, uint8_t c,
    const ulong *exp, const gf28_mpoly_ctx_t ctx)
{
    if (poly->bits == 0) poly->bits = MPOLY_MIN_BITS;
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    ulong *packed_exp = (ulong *)flint_malloc(N * sizeof(ulong));
    mpoly_set_monomial_ui(packed_exp, exp, poly->bits, ctx->minfo);
    
    for (slong pos = 0; pos < poly->length; pos++) {
        if (mpoly_monomial_equal(poly->exps + N*pos, packed_exp, N)) {
            if (c == 0) {
                for (slong i = pos; i < poly->length - 1; i++) {
                    poly->coeffs[i] = poly->coeffs[i + 1];
                    mpoly_monomial_set(poly->exps + N*i, poly->exps + N*(i+1), N);
                }
                poly->length--;
            } else {
                poly->coeffs[pos] = c;
            }
            flint_free(packed_exp);
            return;
        }
    }
    
    if (c != 0) {
        if (poly->length >= poly->coeffs_alloc) {
            slong new_alloc = FLINT_MAX(poly->length + 1, 2 * poly->coeffs_alloc);
            poly->coeffs = (uint8_t *)flint_realloc(poly->coeffs,
                new_alloc * sizeof(uint8_t));
            poly->coeffs_alloc = new_alloc;
        }
        if (N * poly->length >= poly->exps_alloc) {
            slong new_alloc = FLINT_MAX(N * (poly->length + 1), 2 * poly->exps_alloc);
            poly->exps = (ulong *)flint_realloc(poly->exps,
                new_alloc * sizeof(ulong));
            poly->exps_alloc = new_alloc;
        }
        poly->coeffs[poly->length] = c;
        mpoly_monomial_set(poly->exps + N*poly->length, packed_exp, N);
        poly->length++;
    }
    
    flint_free(packed_exp);
}

void gf216_mpoly_set_coeff_ui_ui(gf216_mpoly_t poly, uint16_t c,
    const ulong *exp, const gf216_mpoly_ctx_t ctx)
{
    if (poly->bits == 0) poly->bits = MPOLY_MIN_BITS;
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    ulong *packed_exp = (ulong *)flint_malloc(N * sizeof(ulong));
    mpoly_set_monomial_ui(packed_exp, exp, poly->bits, ctx->minfo);
    
    for (slong pos = 0; pos < poly->length; pos++) {
        if (mpoly_monomial_equal(poly->exps + N*pos, packed_exp, N)) {
            if (c == 0) {
                for (slong i = pos; i < poly->length - 1; i++) {
                    poly->coeffs[i] = poly->coeffs[i + 1];
                    mpoly_monomial_set(poly->exps + N*i, poly->exps + N*(i+1), N);
                }
                poly->length--;
            } else {
                poly->coeffs[pos] = c;
            }
            flint_free(packed_exp);
            return;
        }
    }
    
    if (c != 0) {
        if (poly->length >= poly->coeffs_alloc) {
            slong new_alloc = FLINT_MAX(poly->length + 1, 2 * poly->coeffs_alloc);
            poly->coeffs = (uint16_t *)flint_realloc(poly->coeffs,
                new_alloc * sizeof(uint16_t));
            poly->coeffs_alloc = new_alloc;
        }
        if (N * poly->length >= poly->exps_alloc) {
            slong new_alloc = FLINT_MAX(N * (poly->length + 1), 2 * poly->exps_alloc);
            poly->exps = (ulong *)flint_realloc(poly->exps,
                new_alloc * sizeof(ulong));
            poly->exps_alloc = new_alloc;
        }
        poly->coeffs[poly->length] = c;
        mpoly_monomial_set(poly->exps + N*poly->length, packed_exp, N);
        poly->length++;
    }
    
    flint_free(packed_exp);
}

/* Macro for struct-based coefficient setting */
#define DEFINE_SET_COEFF_STRUCT(FIELD, COEFF_T) \
void FIELD##_mpoly_set_coeff_ui_ui(FIELD##_mpoly_t poly, const COEFF_T *c, \
    const ulong *exp, const FIELD##_mpoly_ctx_t ctx) \
{ \
    if (poly->bits == 0) poly->bits = MPOLY_MIN_BITS; \
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo); \
    ulong *packed_exp = (ulong *)flint_malloc(N * sizeof(ulong)); \
    mpoly_set_monomial_ui(packed_exp, exp, poly->bits, ctx->minfo); \
    \
    for (slong pos = 0; pos < poly->length; pos++) { \
        if (mpoly_monomial_equal(poly->exps + N*pos, packed_exp, N)) { \
            if (FIELD##_is_zero(c)) { \
                for (slong i = pos; i < poly->length - 1; i++) { \
                    poly->coeffs[i] = poly->coeffs[i + 1]; \
                    mpoly_monomial_set(poly->exps + N*i, poly->exps + N*(i+1), N); \
                } \
                poly->length--; \
            } else { \
                poly->coeffs[pos] = *c; \
            } \
            flint_free(packed_exp); \
            return; \
        } \
    } \
    \
    if (!FIELD##_is_zero(c)) { \
        if (poly->length >= poly->coeffs_alloc) { \
            slong new_alloc = FLINT_MAX(poly->length + 1, 2 * poly->coeffs_alloc); \
            poly->coeffs = (COEFF_T *)flint_realloc(poly->coeffs, \
                new_alloc * sizeof(COEFF_T)); \
            poly->coeffs_alloc = new_alloc; \
        } \
        if (N * poly->length >= poly->exps_alloc) { \
            slong new_alloc = FLINT_MAX(N * (poly->length + 1), 2 * poly->exps_alloc); \
            poly->exps = (ulong *)flint_realloc(poly->exps, \
                new_alloc * sizeof(ulong)); \
            poly->exps_alloc = new_alloc; \
        } \
        poly->coeffs[poly->length] = *c; \
        mpoly_monomial_set(poly->exps + N*poly->length, packed_exp, N); \
        poly->length++; \
    } \
    \
    flint_free(packed_exp); \
}

DEFINE_SET_COEFF_STRUCT(gf232, gf232_t)
DEFINE_SET_COEFF_STRUCT(gf264, gf264_t)
DEFINE_SET_COEFF_STRUCT(gf2128, gf2128_t)

/* ============================================================================
   PRINTING - Field-specific implementations
   ============================================================================ */

void gf28_mpoly_print(const gf28_mpoly_t poly, const char **vars,
    const gf28_mpoly_ctx_t ctx)
{
    if (poly->length == 0) { printf("0"); return; }
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
    
    for (slong i = 0; i < poly->length; i++) {
        if (i > 0) printf(" + ");
        printf("0x%02x", poly->coeffs[i]);
        mpoly_get_monomial_ui(exp, poly->exps + N*i, poly->bits, ctx->minfo);
        for (slong j = 0; j < nvars; j++) {
            if (exp[j] > 0) {
                printf("*%s", vars[j]);
                if (exp[j] > 1) printf("^%lu", exp[j]);
            }
        }
    }
    flint_free(exp);
}

void gf216_mpoly_print(const gf216_mpoly_t poly, const char **vars,
    const gf216_mpoly_ctx_t ctx)
{
    if (poly->length == 0) { printf("0"); return; }
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
    
    for (slong i = 0; i < poly->length; i++) {
        if (i > 0) printf(" + ");
        printf("0x%04x", poly->coeffs[i]);
        mpoly_get_monomial_ui(exp, poly->exps + N*i, poly->bits, ctx->minfo);
        for (slong j = 0; j < nvars; j++) {
            if (exp[j] > 0) {
                printf("*%s", vars[j]);
                if (exp[j] > 1) printf("^%lu", exp[j]);
            }
        }
    }
    flint_free(exp);
}

#define DEFINE_PRINT_STRUCT(FIELD, PRINT_FUNC) \
void FIELD##_mpoly_print(const FIELD##_mpoly_t poly, const char **vars, \
    const FIELD##_mpoly_ctx_t ctx) \
{ \
    if (poly->length == 0) { printf("0"); return; } \
    slong N = mpoly_words_per_exp(poly->bits, ctx->minfo); \
    slong nvars = ctx->minfo->nvars; \
    ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong)); \
    \
    for (slong i = 0; i < poly->length; i++) { \
        if (i > 0) printf(" + "); \
        PRINT_FUNC(&poly->coeffs[i]); \
        mpoly_get_monomial_ui(exp, poly->exps + N*i, poly->bits, ctx->minfo); \
        for (slong j = 0; j < nvars; j++) { \
            if (exp[j] > 0) { \
                printf("*%s", vars[j]); \
                if (exp[j] > 1) printf("^%lu", exp[j]); \
            } \
        } \
    } \
    flint_free(exp); \
}

DEFINE_PRINT_STRUCT(gf232, gf232_print)
DEFINE_PRINT_STRUCT(gf264, gf264_print)
DEFINE_PRINT_STRUCT(gf2128, gf2128_print)

/* ============================================================================
   RANDOM TESTING
   ============================================================================ */

void gf28_mpoly_randtest(gf28_mpoly_t poly, flint_rand_t state,
    slong length, slong exp_bound, const gf28_mpoly_ctx_t ctx)
{
    slong nvars = ctx->minfo->nvars;
    ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
    gf28_mpoly_zero(poly, ctx);
    
    for (slong i = 0; i < length; i++) {
        for (slong j = 0; j < nvars; j++) {
            exp[j] = n_randint(state, exp_bound);
        }
        uint8_t c = n_randint(state, 255) + 1;
        gf28_mpoly_set_coeff_ui_ui(poly, c, exp, ctx);
    }
    flint_free(exp);
}

void gf216_mpoly_randtest(gf216_mpoly_t poly, flint_rand_t state,
    slong length, slong exp_bound, const gf216_mpoly_ctx_t ctx)
{
    slong nvars = ctx->minfo->nvars;
    ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
    gf216_mpoly_zero(poly, ctx);
    
    for (slong i = 0; i < length; i++) {
        for (slong j = 0; j < nvars; j++) {
            exp[j] = n_randint(state, exp_bound);
        }
        uint16_t c = n_randint(state, 65535) + 1;
        gf216_mpoly_set_coeff_ui_ui(poly, c, exp, ctx);
    }
    flint_free(exp);
}

#define DEFINE_RANDTEST_STRUCT(FIELD, COEFF_T, CREATE_RANDOM) \
void FIELD##_mpoly_randtest(FIELD##_mpoly_t poly, flint_rand_t state, \
    slong length, slong exp_bound, const FIELD##_mpoly_ctx_t ctx) \
{ \
    slong nvars = ctx->minfo->nvars; \
    ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong)); \
    FIELD##_mpoly_zero(poly, ctx); \
    \
    for (slong i = 0; i < length; i++) { \
        for (slong j = 0; j < nvars; j++) { \
            exp[j] = n_randint(state, exp_bound); \
        } \
        COEFF_T c = CREATE_RANDOM; \
        if (!FIELD##_is_zero(&c)) { \
            FIELD##_mpoly_set_coeff_ui_ui(poly, &c, exp, ctx); \
        } \
    } \
    flint_free(exp); \
}

DEFINE_RANDTEST_STRUCT(gf232, gf232_t, gf232_create(n_randtest(state)))
DEFINE_RANDTEST_STRUCT(gf264, gf264_t, gf264_create(n_randtest(state)))

void gf2128_mpoly_randtest(gf2128_mpoly_t poly, flint_rand_t state,
    slong length, slong exp_bound, const gf2128_mpoly_ctx_t ctx)
{
    slong nvars = ctx->minfo->nvars;
    ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
    gf2128_mpoly_zero(poly, ctx);
    
    for (slong i = 0; i < length; i++) {
        for (slong j = 0; j < nvars; j++) {
            exp[j] = n_randint(state, exp_bound);
        }
        gf2128_t c;
        c.low = n_randtest(state);
        c.high = n_randtest(state);
        if (!gf2128_is_zero(&c)) {
            gf2128_mpoly_set_coeff_ui_ui(poly, &c, exp, ctx);
        }
    }
    flint_free(exp);
}

/* ============================================================================
   FIT LENGTH AND RESET BITS MACRO 
   ============================================================================ */

#define DEFINE_FIT_LENGTH_RESET_BITS(FIELD, COEFF_T) \
void FIELD##_mpoly_fit_length_reset_bits(FIELD##_mpoly_t poly, slong len, \
                                        flint_bitcnt_t bits, const FIELD##_mpoly_ctx_t ctx) \
{ \
    slong N = mpoly_words_per_exp(bits, ctx->minfo); \
    \
    /* Reallocate coefficients if needed */ \
    if (len > poly->coeffs_alloc) { \
        slong new_alloc = FLINT_MAX(len, 2 * poly->coeffs_alloc); \
        poly->coeffs = (COEFF_T *)flint_realloc(poly->coeffs, \
            new_alloc * sizeof(COEFF_T)); \
        poly->coeffs_alloc = new_alloc; \
    } \
    \
    /* Reallocate exponents if needed */ \
    if (N * len > poly->exps_alloc) { \
        slong new_alloc = FLINT_MAX(N * len, 2 * poly->exps_alloc); \
        poly->exps = (ulong *)flint_realloc(poly->exps, \
            new_alloc * sizeof(ulong)); \
        poly->exps_alloc = new_alloc; \
    } \
    \
    poly->bits = bits; \
}

/* Instantiate for all fields */
DEFINE_FIT_LENGTH_RESET_BITS(gf28, uint8_t)
DEFINE_FIT_LENGTH_RESET_BITS(gf216, uint16_t)
DEFINE_FIT_LENGTH_RESET_BITS(gf232, gf232_t)
DEFINE_FIT_LENGTH_RESET_BITS(gf264, gf264_t)
DEFINE_FIT_LENGTH_RESET_BITS(gf2128, gf2128_t)

/* ============================================================================
   CONVERSION FUNCTIONS IMPLEMENTATION
   ============================================================================ */

/* ============================================================================
   GF(2^8) CONVERSIONS (Already implemented in gf28_mpoly.c, shown for completeness)
   ============================================================================ */

void fq_nmod_mpoly_to_gf28_mpoly(gf28_mpoly_t res, const fq_nmod_mpoly_t poly,
    const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx) 
{
    gf28_mpoly_ctx_t ctx;
    gf28_mpoly_ctx_init(ctx, fq_mpoly_ctx->minfo->nvars, fq_mpoly_ctx->minfo->ord);
    
    gf28_mpoly_zero(res, ctx);
    
    slong len = fq_nmod_mpoly_length(poly, fq_mpoly_ctx);
    if (len == 0) {
        gf28_mpoly_ctx_clear(ctx);
        return;
    }
    
    /* Initialize conversion tables if needed */
    if (!g_gf28_conversion || !g_gf28_conversion->initialized) {
        init_gf28_conversion(fqctx);
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
            ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
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

void gf28_mpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t res, const gf28_mpoly_t poly,
    const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx) 
{
    /* Initialize conversion tables if needed */
    if (!g_gf28_conversion || !g_gf28_conversion->initialized) {
        init_gf28_conversion(fqctx);
    }
    
    fq_nmod_mpoly_zero(res, fq_mpoly_ctx);
    
    if (poly->length == 0) return;
    
    slong N = mpoly_words_per_exp(poly->bits, fq_mpoly_ctx->minfo);
    slong nvars = fq_mpoly_ctx->minfo->nvars;
    
    for (slong i = 0; i < poly->length; i++) {
        if (poly->coeffs[i] != 0) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, fqctx);
            
            gf28_elem_to_fq_nmod(coeff, poly->coeffs[i], fqctx);
            
            ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
            mpoly_get_monomial_ui(exp, poly->exps + N*i, poly->bits, fq_mpoly_ctx->minfo);
            
            fq_nmod_mpoly_set_coeff_fq_nmod_ui(res, coeff, exp, fq_mpoly_ctx);
            
            flint_free(exp);
            fq_nmod_clear(coeff, fqctx);
        }
    }
}

/* ============================================================================
   GF(2^16) CONVERSIONS
   ============================================================================ */

void fq_nmod_mpoly_to_gf216_mpoly(gf216_mpoly_t res, const fq_nmod_mpoly_t poly,
    const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx) 
{
    gf216_mpoly_ctx_t ctx;
    gf216_mpoly_ctx_init(ctx, fq_mpoly_ctx->minfo->nvars, fq_mpoly_ctx->minfo->ord);
    
    gf216_mpoly_zero(res, ctx);
    
    slong len = fq_nmod_mpoly_length(poly, fq_mpoly_ctx);
    if (len == 0) {
        gf216_mpoly_ctx_clear(ctx);
        return;
    }
    
    /* Initialize conversion tables if needed */
    if (!g_gf216_conversion || !g_gf216_conversion->initialized) {
        init_gf216_conversion(fqctx);
    }
    
    flint_bitcnt_t bits = FLINT_MAX(poly->bits, MPOLY_MIN_BITS);
    gf216_mpoly_fit_length_reset_bits(res, len, bits, ctx);
    res->length = 0;
    
    slong nvars = fq_mpoly_ctx->minfo->nvars;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    
    slong actual_terms = 0;
    for (slong i = 0; i < len; i++) {
        fq_nmod_t coeff;
        fq_nmod_init(coeff, fqctx);
        
        fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, poly, i, fq_mpoly_ctx);
        uint16_t gf216_coeff = fq_nmod_to_gf216_elem(coeff, fqctx);
        
        if (gf216_coeff != 0) {
            ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
            fq_nmod_mpoly_get_term_exp_ui(exp, poly, i, fq_mpoly_ctx);
            
            mpoly_set_monomial_ui(res->exps + N*actual_terms, exp, bits, ctx->minfo);
            res->coeffs[actual_terms] = gf216_coeff;
            actual_terms++;
            
            flint_free(exp);
        }
        
        fq_nmod_clear(coeff, fqctx);
    }
    
    res->length = actual_terms;
    gf216_mpoly_ctx_clear(ctx);
}

void gf216_mpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t res, const gf216_mpoly_t poly,
    const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx) 
{
    /* Initialize conversion tables if needed */
    if (!g_gf216_conversion || !g_gf216_conversion->initialized) {
        init_gf216_conversion(fqctx);
    }
    
    fq_nmod_mpoly_zero(res, fq_mpoly_ctx);
    
    if (poly->length == 0) return;
    
    slong N = mpoly_words_per_exp(poly->bits, fq_mpoly_ctx->minfo);
    slong nvars = fq_mpoly_ctx->minfo->nvars;
    
    for (slong i = 0; i < poly->length; i++) {
        if (poly->coeffs[i] != 0) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, fqctx);
            
            gf216_elem_to_fq_nmod(coeff, poly->coeffs[i], fqctx);
            
            ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
            mpoly_get_monomial_ui(exp, poly->exps + N*i, poly->bits, fq_mpoly_ctx->minfo);
            
            fq_nmod_mpoly_set_coeff_fq_nmod_ui(res, coeff, exp, fq_mpoly_ctx);
            
            flint_free(exp);
            fq_nmod_clear(coeff, fqctx);
        }
    }
}

/* ============================================================================
   GF(2^32) CONVERSIONS
   ============================================================================ */

void fq_nmod_mpoly_to_gf232_mpoly(gf232_mpoly_t res, const fq_nmod_mpoly_t poly,
    const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx) 
{
    gf232_mpoly_ctx_t ctx;
    gf232_mpoly_ctx_init(ctx, fq_mpoly_ctx->minfo->nvars, fq_mpoly_ctx->minfo->ord);
    
    gf232_mpoly_zero(res, ctx);
    
    slong len = fq_nmod_mpoly_length(poly, fq_mpoly_ctx);
    if (len == 0) {
        gf232_mpoly_ctx_clear(ctx);
        return;
    }
    
    /* Initialize conversion tables if needed */
    if (!g_gf232_conversion || !g_gf232_conversion->initialized) {
        init_gf232_conversion(fqctx);
    }
    
    flint_bitcnt_t bits = FLINT_MAX(poly->bits, MPOLY_MIN_BITS);
    gf232_mpoly_fit_length_reset_bits(res, len, bits, ctx);
    res->length = 0;
    
    slong nvars = fq_mpoly_ctx->minfo->nvars;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    
    slong actual_terms = 0;
    for (slong i = 0; i < len; i++) {
        fq_nmod_t coeff;
        fq_nmod_init(coeff, fqctx);
        
        fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, poly, i, fq_mpoly_ctx);
        gf232_t gf232_coeff = fq_nmod_to_gf232(coeff, fqctx);
        
        if (!gf232_is_zero(&gf232_coeff)) {
            ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
            fq_nmod_mpoly_get_term_exp_ui(exp, poly, i, fq_mpoly_ctx);
            
            mpoly_set_monomial_ui(res->exps + N*actual_terms, exp, bits, ctx->minfo);
            res->coeffs[actual_terms] = gf232_coeff;
            actual_terms++;
            
            flint_free(exp);
        }
        
        fq_nmod_clear(coeff, fqctx);
    }
    
    res->length = actual_terms;
    gf232_mpoly_ctx_clear(ctx);
}

void gf232_mpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t res, const gf232_mpoly_t poly,
    const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx) 
{
    /* Initialize conversion tables if needed */
    if (!g_gf232_conversion || !g_gf232_conversion->initialized) {
        init_gf232_conversion(fqctx);
    }
    
    fq_nmod_mpoly_zero(res, fq_mpoly_ctx);
    
    if (poly->length == 0) return;
    
    slong N = mpoly_words_per_exp(poly->bits, fq_mpoly_ctx->minfo);
    slong nvars = fq_mpoly_ctx->minfo->nvars;
    
    for (slong i = 0; i < poly->length; i++) {
        if (!gf232_is_zero(&poly->coeffs[i])) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, fqctx);
            
            gf232_to_fq_nmod(coeff, &poly->coeffs[i], fqctx);
            
            ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
            mpoly_get_monomial_ui(exp, poly->exps + N*i, poly->bits, fq_mpoly_ctx->minfo);
            
            fq_nmod_mpoly_set_coeff_fq_nmod_ui(res, coeff, exp, fq_mpoly_ctx);
            
            flint_free(exp);
            fq_nmod_clear(coeff, fqctx);
        }
    }
}

/* ============================================================================
   GF(2^64) CONVERSIONS
   ============================================================================ */

void fq_nmod_mpoly_to_gf264_mpoly(gf264_mpoly_t res, const fq_nmod_mpoly_t poly,
    const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx) 
{
    gf264_mpoly_ctx_t ctx;
    gf264_mpoly_ctx_init(ctx, fq_mpoly_ctx->minfo->nvars, fq_mpoly_ctx->minfo->ord);
    
    gf264_mpoly_zero(res, ctx);
    
    slong len = fq_nmod_mpoly_length(poly, fq_mpoly_ctx);
    if (len == 0) {
        gf264_mpoly_ctx_clear(ctx);
        return;
    }
    
    /* Initialize conversion tables if needed */
    if (!g_gf264_conversion || !g_gf264_conversion->initialized) {
        init_gf264_conversion(fqctx);
    }
    
    flint_bitcnt_t bits = FLINT_MAX(poly->bits, MPOLY_MIN_BITS);
    gf264_mpoly_fit_length_reset_bits(res, len, bits, ctx);
    res->length = 0;
    
    slong nvars = fq_mpoly_ctx->minfo->nvars;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    
    slong actual_terms = 0;
    for (slong i = 0; i < len; i++) {
        fq_nmod_t coeff;
        fq_nmod_init(coeff, fqctx);
        
        fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, poly, i, fq_mpoly_ctx);
        gf264_t gf264_coeff = fq_nmod_to_gf264(coeff, fqctx);
        
        if (!gf264_is_zero(&gf264_coeff)) {
            ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
            fq_nmod_mpoly_get_term_exp_ui(exp, poly, i, fq_mpoly_ctx);
            
            mpoly_set_monomial_ui(res->exps + N*actual_terms, exp, bits, ctx->minfo);
            res->coeffs[actual_terms] = gf264_coeff;
            actual_terms++;
            
            flint_free(exp);
        }
        
        fq_nmod_clear(coeff, fqctx);
    }
    
    res->length = actual_terms;
    gf264_mpoly_ctx_clear(ctx);
}

void gf264_mpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t res, const gf264_mpoly_t poly,
    const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx) 
{
    /* Initialize conversion tables if needed */
    if (!g_gf264_conversion || !g_gf264_conversion->initialized) {
        init_gf264_conversion(fqctx);
    }
    
    fq_nmod_mpoly_zero(res, fq_mpoly_ctx);
    
    if (poly->length == 0) return;
    
    slong N = mpoly_words_per_exp(poly->bits, fq_mpoly_ctx->minfo);
    slong nvars = fq_mpoly_ctx->minfo->nvars;
    
    for (slong i = 0; i < poly->length; i++) {
        if (!gf264_is_zero(&poly->coeffs[i])) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, fqctx);
            
            gf264_to_fq_nmod(coeff, &poly->coeffs[i], fqctx);
            
            ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
            mpoly_get_monomial_ui(exp, poly->exps + N*i, poly->bits, fq_mpoly_ctx->minfo);
            
            fq_nmod_mpoly_set_coeff_fq_nmod_ui(res, coeff, exp, fq_mpoly_ctx);
            
            flint_free(exp);
            fq_nmod_clear(coeff, fqctx);
        }
    }
}

/* ============================================================================
   GF(2^128) CONVERSIONS (Already implemented in gf2128_mpoly.c, shown for completeness)
   ============================================================================ */

void fq_nmod_mpoly_to_gf2128_mpoly(gf2128_mpoly_t res, const fq_nmod_mpoly_t poly,
    const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx) 
{
    gf2128_mpoly_ctx_t ctx;
    gf2128_mpoly_ctx_init(ctx, fq_mpoly_ctx->minfo->nvars, fq_mpoly_ctx->minfo->ord);
    
    gf2128_mpoly_zero(res, ctx);
    
    slong len = fq_nmod_mpoly_length(poly, fq_mpoly_ctx);
    if (len == 0) {
        gf2128_mpoly_ctx_clear(ctx);
        return;
    }
    
    /* Initialize conversion if needed */
    if (!g_gf2128_conversion || !g_gf2128_conversion->initialized) {
        init_gf2128_conversion(fqctx);
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
        gf2128_t gf2128_coeff = fq_nmod_to_gf2128(coeff, fqctx);
        
        if (!gf2128_is_zero(&gf2128_coeff)) {
            ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
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

void gf2128_mpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t res, const gf2128_mpoly_t poly,
    const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx) 
{
    /* Initialize conversion if needed */
    if (!g_gf2128_conversion || !g_gf2128_conversion->initialized) {
        init_gf2128_conversion(fqctx);
    }
    
    fq_nmod_mpoly_zero(res, fq_mpoly_ctx);
    
    if (poly->length == 0) return;
    
    slong N = mpoly_words_per_exp(poly->bits, fq_mpoly_ctx->minfo);
    slong nvars = fq_mpoly_ctx->minfo->nvars;
    
    for (slong i = 0; i < poly->length; i++) {
        if (!gf2128_is_zero(&poly->coeffs[i])) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, fqctx);
            
            gf2128_to_fq_nmod(coeff, &poly->coeffs[i], fqctx);
            
            ulong *exp = (ulong *)flint_malloc(nvars * sizeof(ulong));
            mpoly_get_monomial_ui(exp, poly->exps + N*i, poly->bits, fq_mpoly_ctx->minfo);
            
            fq_nmod_mpoly_set_coeff_fq_nmod_ui(res, coeff, exp, fq_mpoly_ctx);
            
            flint_free(exp);
            fq_nmod_clear(coeff, fqctx);
        }
    }
}