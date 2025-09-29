/* gf28_mpoly.h - GF(2^8) Multivariate Polynomial Header */
#ifndef GF28_MPOLY_H
#define GF28_MPOLY_H

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

/* Debug flag */
#define DEBUG_DIVISION 0

/* Block size for array operations */
#define BLOCK 256

/* ============================================================================
   GF(2^8) MULTIVARIATE POLYNOMIAL STRUCTURES
   ============================================================================ */

typedef struct {
    uint8_t *coeffs;
    ulong *exps;
    slong length;
    slong coeffs_alloc;
    slong exps_alloc;
    flint_bitcnt_t bits;
} gf28_mpoly_struct;

typedef gf28_mpoly_struct gf28_mpoly_t[1];

typedef struct {
    mpoly_ctx_t minfo;
} gf28_mpoly_ctx_struct;

typedef gf28_mpoly_ctx_struct gf28_mpoly_ctx_t[1];

/* ============================================================================
   DENSE POLYNOMIAL STRUCTURES FOR DIVISION
   ============================================================================ */

typedef struct {
    uint8_t *coeffs;
    slong alloc;
    slong *deg_bounds;
    slong nvars;
} gf28_mpolyd_struct;

typedef gf28_mpolyd_struct gf28_mpolyd_t[1];

/* ============================================================================
   DYNAMIC ARRAY STRUCTURE
   ============================================================================ */

typedef struct {
    uint8_t *data;
    slong size;
} gf28_dynamic_array_t;

/* ============================================================================
   BASIC INITIALIZATION AND MEMORY MANAGEMENT
   ============================================================================ */

void gf28_mpoly_init(gf28_mpoly_t poly, const gf28_mpoly_ctx_t ctx);
void gf28_mpoly_clear(gf28_mpoly_t poly, const gf28_mpoly_ctx_t ctx);
void gf28_mpoly_ctx_init(gf28_mpoly_ctx_t ctx, slong nvars, const ordering_t ord);
void gf28_mpoly_ctx_clear(gf28_mpoly_ctx_t ctx);
void gf28_mpoly_zero(gf28_mpoly_t poly, const gf28_mpoly_ctx_t ctx);

void _gf28_mpoly_fit_length(uint8_t **coeffs, slong *coeffs_alloc,
                           ulong **exps, slong *exps_alloc, slong N, slong length);

void gf28_mpoly_fit_length_reset_bits(gf28_mpoly_t poly, slong len, 
                                     flint_bitcnt_t bits, const gf28_mpoly_ctx_t ctx);

void gf28_mpoly_init3(gf28_mpoly_t poly, slong alloc, flint_bitcnt_t bits,
                     const gf28_mpoly_ctx_t ctx);

void gf28_mpoly_swap(gf28_mpoly_t poly1, gf28_mpoly_t poly2, const gf28_mpoly_ctx_t ctx);

void _gf28_mpoly_set_length(gf28_mpoly_t poly, slong len, const gf28_mpoly_ctx_t ctx);

/* ============================================================================
   HELPER FUNCTIONS
   ============================================================================ */

void gf28_mpoly_set(gf28_mpoly_t res, const gf28_mpoly_t poly, const gf28_mpoly_ctx_t ctx);

void gf28_mpoly_set_coeff_ui_ui(gf28_mpoly_t poly, uint8_t c, 
                                const ulong *exp, const gf28_mpoly_ctx_t ctx);

void gf28_mpoly_print(const gf28_mpoly_t poly, const char **vars, 
                     const gf28_mpoly_ctx_t ctx);

/* ============================================================================
   ARRAY MULTIPLICATION FUNCTIONS
   ============================================================================ */

void _gf28_mpoly_addmul_array1_safe(uint8_t *poly1, slong array_size,
                                   const uint8_t *poly2, const ulong *exp2, slong len2,
                                   const uint8_t *poly3, const ulong *exp3, slong len3);

slong gf28_mpoly_append_array_LEX_safe(gf28_mpoly_t P, slong Plen, uint8_t *coeff_array,
                                      const ulong *mults, slong num, slong array_size, 
                                      slong top, const gf28_mpoly_ctx_t ctx);

void gf28_dynamic_array_init(gf28_dynamic_array_t *arr, slong size);
void gf28_dynamic_array_clear(gf28_dynamic_array_t *arr);

void _gf28_mpoly_mul_array_chunked_LEX(gf28_mpoly_t P,
                                      const gf28_mpoly_t A,
                                      const gf28_mpoly_t B,
                                      const ulong *mults,
                                      const gf28_mpoly_ctx_t ctx);

int _gf28_mpoly_mul_array_LEX(gf28_mpoly_t A,
                             const gf28_mpoly_t B,
                             fmpz *maxBfields,
                             const gf28_mpoly_t C,
                             fmpz *maxCfields,
                             const gf28_mpoly_ctx_t ctx);

int gf28_mpoly_mul_array(gf28_mpoly_t A, const gf28_mpoly_t B,
                        const gf28_mpoly_t C, const gf28_mpoly_ctx_t ctx);

/* ============================================================================
   CONVERSION FUNCTIONS
   ============================================================================ */

void fq_nmod_mpoly_to_gf28_mpoly(gf28_mpoly_t res, const fq_nmod_mpoly_t poly,
                                 const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx);

void gf28_mpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t res, const gf28_mpoly_t poly,
                                const fq_nmod_ctx_t fqctx, const fq_nmod_mpoly_ctx_t fq_mpoly_ctx);

/* ============================================================================
   INTELLIGENT MULTIPLICATION AND UNIVARIATE OPTIMIZATION
   ============================================================================ */

int gf28_mpoly_is_univariate(const gf28_mpoly_t poly, const gf28_mpoly_ctx_t ctx, 
                            slong *main_var);

slong gf28_mpoly_univariate_degree(const gf28_mpoly_t poly, const gf28_mpoly_ctx_t ctx,
                                  slong main_var);

void gf28_mpoly_to_gf28_poly_univar(gf28_poly_t res, const gf28_mpoly_t poly,
                                   const gf28_mpoly_ctx_t ctx, slong main_var);

void gf28_poly_to_gf28_mpoly_univar(gf28_mpoly_t res, const gf28_poly_t poly,
                                   const gf28_mpoly_ctx_t ctx, slong main_var);

void gf28_mpoly_mul_flint_univar(gf28_mpoly_t res, const gf28_mpoly_t a, 
                                const gf28_mpoly_t b, const gf28_mpoly_ctx_t ctx,
                                slong main_var);

int gf28_mpoly_mul(gf28_mpoly_t res, const gf28_mpoly_t a, const gf28_mpoly_t b, 
                   const gf28_mpoly_ctx_t ctx);

/* ============================================================================
   DIVISION IMPLEMENTATION
   ============================================================================ */

void gf28_mpolyd_init(gf28_mpolyd_t A, slong nvars);
void gf28_mpolyd_clear(gf28_mpolyd_t A);
slong gf28_mpolyd_offset(const gf28_mpolyd_t A, const ulong *exp);

void batch_xor_sse2(uint8_t* dst, const uint8_t* src, slong len);

void gf28_mpolyd_divrem_univar(gf28_mpolyd_t Q, gf28_mpolyd_t R,
                              const gf28_mpolyd_t A, const gf28_mpolyd_t B);

int gf28_mpolyd_is_zero(const gf28_mpolyd_t A);
int gf28_mpoly_is_monomial(const gf28_mpoly_t poly);
uint8_t gf28_mpoly_get_monomial_coeff(const gf28_mpoly_t poly);

void gf28_mpoly_get_monomial_exp(ulong *exp, const gf28_mpoly_t poly, 
                                const gf28_mpoly_ctx_t ctx);

void gf28_mpoly_to_mpolyd(gf28_mpolyd_t A, const gf28_mpoly_t B, 
                         const gf28_mpoly_ctx_t ctx);

void gf28_mpolyd_to_mpoly_fast(gf28_mpoly_t A, const gf28_mpolyd_t B,
                              const gf28_mpoly_ctx_t ctx);

void gf28_mpolyd_to_mpoly_univariate(gf28_mpoly_t A, const gf28_mpolyd_t B,
                                    const gf28_mpoly_ctx_t ctx);

void gf28_mpolyd_to_mpoly(gf28_mpoly_t A, const gf28_mpolyd_t B,
                         const gf28_mpoly_ctx_t ctx);

void gf28_mpolyd_divrem_multivar(gf28_mpolyd_t Q, gf28_mpolyd_t R,
                                const gf28_mpolyd_t A, const gf28_mpolyd_t B,
                                const gf28_mpoly_ctx_t ctx);

int gf28_mpolyd_set_degbounds(gf28_mpolyd_t A, const slong *bounds);

int gf28_mpoly_divides_monomial(gf28_mpoly_t Q, const gf28_mpoly_t A, 
                               const gf28_mpoly_t B, const gf28_mpoly_ctx_t ctx);

int gf28_mpoly_divides_dense(gf28_mpoly_t Q, const gf28_mpoly_t A, 
                            const gf28_mpoly_t B, const gf28_mpoly_ctx_t ctx);

int gf28_mpoly_divides_flint_univar(gf28_mpoly_t Q, const gf28_mpoly_t A, 
                                   const gf28_mpoly_t B, const gf28_mpoly_ctx_t ctx,
                                   slong main_var);

int gf28_mpoly_divides(gf28_mpoly_t Q, const gf28_mpoly_t A, 
                      const gf28_mpoly_t B, const gf28_mpoly_ctx_t ctx);

/* ============================================================================
   ADDITIONAL UTILITY FUNCTIONS
   ============================================================================ */

void gf28_mpoly_mul_simple(gf28_mpoly_t res, const gf28_mpoly_t a, 
                          const gf28_mpoly_t b, const gf28_mpoly_ctx_t ctx);

int gf28_mpoly_equal(const gf28_mpoly_t A, const gf28_mpoly_t B, 
                    const gf28_mpoly_ctx_t ctx);

void gf28_mpoly_randtest(gf28_mpoly_t poly, flint_rand_t state,
                        slong length, slong exp_bound, 
                        const gf28_mpoly_ctx_t ctx);

/* ============================================================================
   TIMING UTILITIES
   ============================================================================ */

double get_wall_time_8(void);

#ifdef __cplusplus
}
#endif

#endif /* GF28_MPOLY_H */