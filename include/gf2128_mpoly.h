/* gf2128_mpoly.h - GF(2^128) Multivariate Polynomial Header */
#ifndef GF2128_MPOLY_H
#define GF2128_MPOLY_H

/*
 * GF(2^128) Multivariate Polynomial Array Multiplication and Division
 * 
 * This implementation provides fast multiplication and division of multivariate 
 * polynomials over GF(2^128) using array-based methods with dynamic memory allocation.
 * 
 * Key features:
 * 1. Uses native PCLMUL instructions for GF(2^128) arithmetic when available
 * 2. Dynamic memory allocation to support arbitrarily large polynomials
 * 3. Optimized memory usage with chunking for very large polynomials
 * 4. Implements both multiplication and exact division
 * 5. Compatible with FLINT polynomial structures
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
#include "gf2n_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Block size for array operations (smaller for GF(2^128) due to larger element size) */
#define BLOCK 256

/* ============================================================================
   GF(2^128) MULTIVARIATE POLYNOMIAL STRUCTURES
   ============================================================================ */

typedef struct {
    gf2128_t *coeffs;
    ulong *exps;
    slong length;
    slong coeffs_alloc;
    slong exps_alloc;
    flint_bitcnt_t bits;
} gf2128_mpoly_struct;

typedef gf2128_mpoly_struct gf2128_mpoly_t[1];

typedef struct {
    mpoly_ctx_t minfo;
} gf2128_mpoly_ctx_struct;

typedef gf2128_mpoly_ctx_struct gf2128_mpoly_ctx_t[1];

/* ============================================================================
   DYNAMIC MEMORY MANAGEMENT FOR LARGE ARRAYS
   ============================================================================ */

typedef struct {
    gf2128_t *data;
    slong size;
    slong chunk_size;
    int use_chunks;
} gf2128_dynamic_array_t;

/* ============================================================================
   DENSE POLYNOMIAL STRUCTURES FOR DIVISION
   ============================================================================ */

typedef struct {
    gf2128_t *coeffs;
    slong alloc;
    slong *deg_bounds;
    slong nvars;
} gf2128_mpolyd_struct;

typedef gf2128_mpolyd_struct gf2128_mpolyd_t[1];

/* ============================================================================
   BASIC INITIALIZATION AND MEMORY MANAGEMENT
   ============================================================================ */

void gf2128_mpoly_init(gf2128_mpoly_t poly, const gf2128_mpoly_ctx_t ctx);
void gf2128_mpoly_clear(gf2128_mpoly_t poly, const gf2128_mpoly_ctx_t ctx);
void gf2128_mpoly_ctx_init(gf2128_mpoly_ctx_t ctx, slong nvars, const ordering_t ord);
void gf2128_mpoly_ctx_clear(gf2128_mpoly_ctx_t ctx);
void gf2128_mpoly_zero(gf2128_mpoly_t poly, const gf2128_mpoly_ctx_t ctx);

void _gf2128_mpoly_fit_length(gf2128_t **coeffs, slong *coeffs_alloc,
                             ulong **exps, slong *exps_alloc, slong N, slong length);

void gf2128_mpoly_fit_length_reset_bits(gf2128_mpoly_t poly, slong len, 
                                       flint_bitcnt_t bits, const gf2128_mpoly_ctx_t ctx);

void gf2128_mpoly_init3(gf2128_mpoly_t poly, slong alloc, flint_bitcnt_t bits,
                       const gf2128_mpoly_ctx_t ctx);

void gf2128_mpoly_swap(gf2128_mpoly_t poly1, gf2128_mpoly_t poly2, const gf2128_mpoly_ctx_t ctx);

void _gf2128_mpoly_set_length(gf2128_mpoly_t poly, slong len, const gf2128_mpoly_ctx_t ctx);

/* ============================================================================
   DYNAMIC MEMORY MANAGEMENT
   ============================================================================ */

void gf2128_dynamic_array_init(gf2128_dynamic_array_t *arr, slong size);
void gf2128_dynamic_array_clear(gf2128_dynamic_array_t *arr);

/* ============================================================================
   ARRAY MULTIPLICATION OPERATIONS
   ============================================================================ */

void _gf2128_mpoly_addmul_array1_safe(gf2128_t *poly1, slong array_size,
                                      const gf2128_t *poly2, const ulong *exp2, slong len2,
                                      const gf2128_t *poly3, const ulong *exp3, slong len3);

slong gf2128_mpoly_append_array_LEX_safe(gf2128_mpoly_t P, slong Plen, gf2128_t *coeff_array,
                                         const ulong *mults, slong num, slong array_size, 
                                         slong top, const gf2128_mpoly_ctx_t ctx);

slong gf2128_mpoly_append_array_DEGLEX_safe(gf2128_mpoly_t P, slong Plen, gf2128_t *coeff_array,
                                            slong top, slong nvars, slong degb,
                                            const gf2128_mpoly_ctx_t ctx);

slong gf2128_mpoly_append_array_DEGREVLEX_safe(gf2128_mpoly_t P, slong Plen, gf2128_t *coeff_array,
                                               slong top, slong nvars, slong degb,
                                               const gf2128_mpoly_ctx_t ctx);

void _gf2128_mpoly_mul_array_chunked_LEX(gf2128_mpoly_t P,
                                        const gf2128_mpoly_t A,
                                        const gf2128_mpoly_t B,
                                        const ulong *mults,
                                        const gf2128_mpoly_ctx_t ctx);

void _gf2128_mpoly_mul_array_chunked_DEG(gf2128_mpoly_t P,
                                        const gf2128_mpoly_t A,
                                        const gf2128_mpoly_t B,
                                        ulong degb,
                                        const gf2128_mpoly_ctx_t ctx);

int _gf2128_mpoly_mul_array_LEX(gf2128_mpoly_t A,
                               const gf2128_mpoly_t B,
                               fmpz *maxBfields,
                               const gf2128_mpoly_t C,
                               fmpz *maxCfields,
                               const gf2128_mpoly_ctx_t ctx);

int _gf2128_mpoly_mul_array_DEG(gf2128_mpoly_t A,
                               const gf2128_mpoly_t B,
                               fmpz *maxBfields,
                               const gf2128_mpoly_t C,
                               fmpz *maxCfields,
                               const gf2128_mpoly_ctx_t ctx);

int gf2128_mpoly_mul_array(gf2128_mpoly_t A, const gf2128_mpoly_t B,
                          const gf2128_mpoly_t C, const gf2128_mpoly_ctx_t ctx);

/* ============================================================================
   HELPER FUNCTIONS
   ============================================================================ */

void gf2128_mpoly_set(gf2128_mpoly_t res, const gf2128_mpoly_t poly, const gf2128_mpoly_ctx_t ctx);

void gf2128_mpoly_set_coeff_ui_ui(gf2128_mpoly_t poly, const gf2128_t *c, 
                                  const ulong *exp, const gf2128_mpoly_ctx_t ctx);

void gf2128_mpoly_print(const gf2128_mpoly_t poly, const char **vars, 
                       const gf2128_mpoly_ctx_t ctx);

/* ============================================================================
   DIVISION SUPPORT FUNCTIONS
   ============================================================================ */

void gf2128_mpolyd_init(gf2128_mpolyd_t A, slong nvars);
void gf2128_mpolyd_clear(gf2128_mpolyd_t A);
slong gf2128_mpolyd_offset(const gf2128_mpolyd_t A, const ulong *exp);

void gf2128_mpolyd_divrem_univar(gf2128_mpolyd_t Q, gf2128_mpolyd_t R,
                                 const gf2128_mpolyd_t A, const gf2128_mpolyd_t B);

int gf2128_mpolyd_is_zero(const gf2128_mpolyd_t A);
int gf2128_mpoly_is_monomial(const gf2128_mpoly_t poly);
gf2128_t gf2128_mpoly_get_monomial_coeff(const gf2128_mpoly_t poly);

void gf2128_mpoly_get_monomial_exp(ulong *exp, const gf2128_mpoly_t poly, 
                                   const gf2128_mpoly_ctx_t ctx);

void gf2128_mpoly_to_mpolyd(gf2128_mpolyd_t A, const gf2128_mpoly_t B, 
                            const gf2128_mpoly_ctx_t ctx);

void gf2128_mpolyd_to_mpoly(gf2128_mpoly_t A, const gf2128_mpolyd_t B,
                            const gf2128_mpoly_ctx_t ctx);

void gf2128_mpolyd_divrem_multivar(gf2128_mpolyd_t Q, gf2128_mpolyd_t R,
                                   const gf2128_mpolyd_t A, const gf2128_mpolyd_t B,
                                   const gf2128_mpoly_ctx_t ctx);

int gf2128_mpolyd_set_degbounds(gf2128_mpolyd_t A, const slong *bounds);

int gf2128_mpoly_divides_monomial(gf2128_mpoly_t Q, const gf2128_mpoly_t A, 
                                  const gf2128_mpoly_t B, const gf2128_mpoly_ctx_t ctx);

int gf2128_mpoly_divides_dense(gf2128_mpoly_t Q, const gf2128_mpoly_t A, 
                               const gf2128_mpoly_t B, const gf2128_mpoly_ctx_t ctx);

int gf2128_mpoly_divides(gf2128_mpoly_t Q, const gf2128_mpoly_t A, 
                        const gf2128_mpoly_t B, const gf2128_mpoly_ctx_t ctx);

/* ============================================================================
   ADDITIONAL HELPER FUNCTIONS
   ============================================================================ */

void gf2128_mpoly_mul_simple(gf2128_mpoly_t res, const gf2128_mpoly_t a, 
                             const gf2128_mpoly_t b, const gf2128_mpoly_ctx_t ctx);

int gf2128_mpoly_equal(const gf2128_mpoly_t A, const gf2128_mpoly_t B, 
                       const gf2128_mpoly_ctx_t ctx);

void gf2128_mpoly_randtest(gf2128_mpoly_t poly, flint_rand_t state,
                           slong length, slong exp_bound, 
                           const gf2128_mpoly_ctx_t ctx);

/* ============================================================================
   CONVERSION FUNCTIONS FOR FLINT COMPATIBILITY
   ============================================================================ */

void gf2128_mpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t res, const gf2128_mpoly_t poly,
                                   const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx);

void fq_nmod_mpoly_to_gf2128_mpoly(gf2128_mpoly_t res, const fq_nmod_mpoly_t poly,
                                   const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx);

/* ============================================================================
   TIMING UTILITIES
   ============================================================================ */

double get_wall_time_128(void);

#ifdef __cplusplus
}
#endif

#endif /* GF2128_MPOLY_H */