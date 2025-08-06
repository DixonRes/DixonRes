/* unified_mpoly_interface.h - Complete single-file implementation */

#ifndef UNIFIED_MPOLY_INTERFACE_H
#define UNIFIED_MPOLY_INTERFACE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <flint/flint.h>
#include <flint/nmod_mpoly.h>
#include <flint/fq_nmod_mpoly.h>
#include <flint/fmpz.h>
#include <flint/mpoly.h>
#include <flint/longlong.h>
#include "fq_unified_interface.h"
#include "gf28_mpoly.h"
#include "gf2128_mpoly.h"

/* ============================================================================
   UNIFIED MULTIVARIATE POLYNOMIAL TYPES AND STRUCTURES
   ============================================================================ */

/* Context for multivariate polynomials */
typedef struct {
    field_ctx_t *field_ctx;
    slong nvars;
    ordering_t ord;
    
    /* Different contexts based on field type */
    union {
        nmod_mpoly_ctx_struct nmod_ctx;
        fq_nmod_mpoly_ctx_struct fq_ctx;
    } ctx;
    
    /* Method selection parameters */
    slong array_size_limit;
    slong dense_threshold;
} unified_mpoly_ctx_struct;

typedef unified_mpoly_ctx_struct *unified_mpoly_ctx_t;

/* Unified multivariate polynomial type */
typedef struct {
    /* Field type identifier */
    field_id_t field_id;
    
    /* Storage for different polynomial types */
    union {
        nmod_mpoly_struct nmod_poly;
        fq_nmod_mpoly_struct fq_poly;
    } data;
    
    /* Context pointer */
    unified_mpoly_ctx_t ctx_ptr;
} unified_mpoly_struct;

typedef unified_mpoly_struct *unified_mpoly_t;

/* ============================================================================
   HELPER MACROS
   ============================================================================ */

#define GET_NMOD_POLY(poly) (&(poly)->data.nmod_poly)
#define GET_FQ_POLY(poly) (&(poly)->data.fq_poly)

#define GET_NMOD_CTX(ctx) (&(ctx)->ctx.nmod_ctx)
#define GET_FQ_CTX(ctx) (&(ctx)->ctx.fq_ctx)

/* ============================================================================
   CONTEXT OPERATIONS
   ============================================================================ */

unified_mpoly_ctx_t unified_mpoly_ctx_init(slong nvars, const ordering_t ord, field_ctx_t *field_ctx) {
    unified_mpoly_ctx_t ctx = (unified_mpoly_ctx_t)malloc(sizeof(unified_mpoly_ctx_struct));
    if (!ctx) return NULL;
    
    ctx->field_ctx = field_ctx;
    ctx->nvars = nvars;
    ctx->ord = ord;
    ctx->array_size_limit = 1L << 26;  /* 64MB default */
    ctx->dense_threshold = 10;         /* 10% density threshold */
    
    switch (field_ctx->field_id) {
        case FIELD_ID_NMOD:
            {
                /* For prime fields, we need the modulus */
                ulong modulus = field_ctx->ctx.nmod_ctx.n;
                //printf("DEBUG: Initializing nmod context with modulus=%lu, nvars=%ld\n", modulus, nvars);
                nmod_mpoly_ctx_init(GET_NMOD_CTX(ctx), nvars, ord, modulus);
                
                // Verify the context
                nmod_mpoly_ctx_struct *nctx = GET_NMOD_CTX(ctx);
                //printf("DEBUG: nmod_ctx initialized: nvars=%ld, mod.n=%lu\n", nmod_mpoly_ctx_nvars(nctx), nctx->mod.n);
            }
            break;
            
        case FIELD_ID_GF28:
        case FIELD_ID_GF216:
        case FIELD_ID_GF232:
        case FIELD_ID_GF264:
        case FIELD_ID_GF2128:
            /* For GF(2^n) fields, use fq_nmod context */
            if (field_ctx->ctx.fq_ctx) {
                fq_nmod_mpoly_ctx_init(GET_FQ_CTX(ctx), nvars, ord, field_ctx->ctx.fq_ctx);
            }
            break;
            
        default:
            /* General finite fields */
            fq_nmod_mpoly_ctx_init(GET_FQ_CTX(ctx), nvars, ord, field_ctx->ctx.fq_ctx);
            break;
    }
    
    return ctx;
}

void unified_mpoly_ctx_clear(unified_mpoly_ctx_t ctx) {
    if (!ctx) return;
    
    switch (ctx->field_ctx->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_ctx_clear(GET_NMOD_CTX(ctx));
            break;
            
        default:
            fq_nmod_mpoly_ctx_clear(GET_FQ_CTX(ctx));
            break;
    }
    
    free(ctx);
}

/* ============================================================================
   POLYNOMIAL MEMORY MANAGEMENT
   ============================================================================ */

unified_mpoly_t unified_mpoly_init(unified_mpoly_ctx_t ctx) {
    unified_mpoly_t poly = (unified_mpoly_t)malloc(sizeof(unified_mpoly_struct));
    if (!poly) return NULL;
    
    poly->field_id = ctx->field_ctx->field_id;
    poly->ctx_ptr = ctx;
    
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_init(GET_NMOD_POLY(poly), GET_NMOD_CTX(ctx));
            break;
            
        default:
            fq_nmod_mpoly_init(GET_FQ_POLY(poly), GET_FQ_CTX(ctx));
            break;
    }
    
    return poly;
}

void unified_mpoly_clear(unified_mpoly_t poly) {
    if (!poly) return;
    
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_clear(GET_NMOD_POLY(poly), GET_NMOD_CTX(ctx));
            break;
            
        default:
            fq_nmod_mpoly_clear(GET_FQ_POLY(poly), GET_FQ_CTX(ctx));
            break;
    }
    
    free(poly);
}

/* ============================================================================
   BASIC OPERATIONS
   ============================================================================ */

void unified_mpoly_zero(unified_mpoly_t poly) {
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_zero(GET_NMOD_POLY(poly), GET_NMOD_CTX(ctx));
            break;
            
        default:
            fq_nmod_mpoly_zero(GET_FQ_POLY(poly), GET_FQ_CTX(ctx));
            break;
    }
}

void unified_mpoly_one(unified_mpoly_t poly) {
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_one(GET_NMOD_POLY(poly), GET_NMOD_CTX(ctx));
            break;
            
        default:
            fq_nmod_mpoly_one(GET_FQ_POLY(poly), GET_FQ_CTX(ctx));
            break;
    }
}

int unified_mpoly_is_zero(const unified_mpoly_t poly) {
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            return nmod_mpoly_is_zero(GET_NMOD_POLY(poly), GET_NMOD_CTX(ctx));
            
        default:
            return fq_nmod_mpoly_is_zero(GET_FQ_POLY(poly), GET_FQ_CTX(ctx));
    }
}

int unified_mpoly_is_one(const unified_mpoly_t poly) {
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            return nmod_mpoly_is_one(GET_NMOD_POLY(poly), GET_NMOD_CTX(ctx));
            
        default:
            return fq_nmod_mpoly_is_one(GET_FQ_POLY(poly), GET_FQ_CTX(ctx));
    }
}

slong unified_mpoly_length(const unified_mpoly_t poly) {
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            return nmod_mpoly_length(GET_NMOD_POLY(poly), GET_NMOD_CTX(ctx));
            
        default:
            return fq_nmod_mpoly_length(GET_FQ_POLY(poly), GET_FQ_CTX(ctx));
    }
}

void unified_mpoly_set(unified_mpoly_t poly1, const unified_mpoly_t poly2) {
    if (poly1 == poly2) return;
    
    unified_mpoly_ctx_t ctx = poly1->ctx_ptr;
    
    switch (poly1->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_set(GET_NMOD_POLY(poly1), GET_NMOD_POLY(poly2), GET_NMOD_CTX(ctx));
            break;
            
        default:
            fq_nmod_mpoly_set(GET_FQ_POLY(poly1), GET_FQ_POLY(poly2), GET_FQ_CTX(ctx));
            break;
    }
}

/* ============================================================================
   ARITHMETIC OPERATIONS
   ============================================================================ */

void unified_mpoly_add(unified_mpoly_t poly1, const unified_mpoly_t poly2,
                      const unified_mpoly_t poly3) {
    unified_mpoly_ctx_t ctx = poly1->ctx_ptr;
    
    switch (poly1->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_add(GET_NMOD_POLY(poly1), GET_NMOD_POLY(poly2),
                          GET_NMOD_POLY(poly3), GET_NMOD_CTX(ctx));
            break;
            
        default:
            fq_nmod_mpoly_add(GET_FQ_POLY(poly1), GET_FQ_POLY(poly2),
                             GET_FQ_POLY(poly3), GET_FQ_CTX(ctx));
            break;
    }
}

void unified_mpoly_sub(unified_mpoly_t poly1, const unified_mpoly_t poly2,
                      const unified_mpoly_t poly3) {
    unified_mpoly_ctx_t ctx = poly1->ctx_ptr;
    
    switch (poly1->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_sub(GET_NMOD_POLY(poly1), GET_NMOD_POLY(poly2),
                          GET_NMOD_POLY(poly3), GET_NMOD_CTX(ctx));
            break;
            
        case FIELD_ID_GF28:
        case FIELD_ID_GF216:
        case FIELD_ID_GF232:
        case FIELD_ID_GF264:
        case FIELD_ID_GF2128:
            /* In GF(2^n), subtraction is the same as addition */
            unified_mpoly_add(poly1, poly2, poly3);
            break;
            
        default:
            fq_nmod_mpoly_sub(GET_FQ_POLY(poly1), GET_FQ_POLY(poly2),
                             GET_FQ_POLY(poly3), GET_FQ_CTX(ctx));
            break;
    }
}

void unified_mpoly_neg(unified_mpoly_t poly1, const unified_mpoly_t poly2) {
    unified_mpoly_ctx_t ctx = poly1->ctx_ptr;
    
    switch (poly1->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_neg(GET_NMOD_POLY(poly1), GET_NMOD_POLY(poly2), GET_NMOD_CTX(ctx));
            break;
            
        case FIELD_ID_GF28:
        case FIELD_ID_GF216:
        case FIELD_ID_GF232:
        case FIELD_ID_GF264:
        case FIELD_ID_GF2128:
            /* In GF(2^n), -a = a */
            unified_mpoly_set(poly1, poly2);
            break;
            
        default:
            fq_nmod_mpoly_neg(GET_FQ_POLY(poly1), GET_FQ_POLY(poly2), GET_FQ_CTX(ctx));
            break;
    }
}

/* ============================================================================
   MULTIPLICATION - WITH OPTIMIZATION HOOKS
   ============================================================================ */

/* Global flags to enable/disable optimizations */
static int use_gf28_array_mul = 0;
static int use_gf2128_array_mul = 0;

/* Function to enable optimizations */
void unified_mpoly_enable_optimizations(field_id_t field_id, int enable) {
    switch (field_id) {
        case FIELD_ID_GF28:
            use_gf28_array_mul = enable;
            //printf("GF(2^8) array multiplication: %s\n", enable ? "enabled" : "disabled");
            break;
        case FIELD_ID_GF2128:
            use_gf2128_array_mul = enable;
            //printf("GF(2^128) array multiplication: %s\n", enable ? "enabled" : "disabled");
            break;
        default:
            break;
    }
}

int unified_mpoly_mul(unified_mpoly_t poly1, const unified_mpoly_t poly2,
                     const unified_mpoly_t poly3) {
    unified_mpoly_ctx_t ctx = poly1->ctx_ptr;

    switch (poly1->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_mul(GET_NMOD_POLY(poly1), GET_NMOD_POLY(poly2),
                          GET_NMOD_POLY(poly3), GET_NMOD_CTX(ctx));
            return 1;
            
        case FIELD_ID_GF28:
            if (1 || use_gf28_array_mul) {
                /* Use optimized GF(2^8) array multiplication */
                gf28_mpoly_t A, B, C;
                gf28_mpoly_ctx_t native_ctx;
                
                /* Initialize native context */
                gf28_mpoly_ctx_init(native_ctx, ctx->nvars, ctx->ord);
                
                /* Initialize native polynomials */
                gf28_mpoly_init(A, native_ctx);
                gf28_mpoly_init(B, native_ctx);
                gf28_mpoly_init(C, native_ctx);
                
                /* Convert from fq_nmod to native format */
                fq_nmod_mpoly_to_gf28_mpoly(A, GET_FQ_POLY(poly2), 
                                            ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                fq_nmod_mpoly_to_gf28_mpoly(B, GET_FQ_POLY(poly3), 
                                            ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                
                /* Try array multiplication */
                int success = gf28_mpoly_mul_array(C, A, B, native_ctx);
                
                if (success) {
                    /* Convert result back to fq_nmod format */
                    gf28_mpoly_to_fq_nmod_mpoly(GET_FQ_POLY(poly1), C, 
                                                ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
#if 0                    
                    /* CORRECTNESS TEST: Compare with standard multiplication */
                    fq_nmod_mpoly_t test_result;
                    fq_nmod_mpoly_init(test_result, GET_FQ_CTX(ctx));



                    fq_nmod_mpoly_mul(test_result, GET_FQ_POLY(poly2), GET_FQ_POLY(poly3), GET_FQ_CTX(ctx));


                    
                    if (!fq_nmod_mpoly_equal(GET_FQ_POLY(poly1), test_result, GET_FQ_CTX(ctx))) {
                        printf("\n=== GF(2^8) MULTIPLICATION ERROR ===\n");
                        printf("poly1 addr: %p\n", GET_FQ_POLY(poly1));
                        printf("poly2 addr: %p\n", GET_FQ_POLY(poly2));
                        printf("poly3 addr: %p\n", GET_FQ_POLY(poly3));
                        printf("Input A length: %ld, Input B length: %ld\n", 
                               fq_nmod_mpoly_length(GET_FQ_POLY(poly2), GET_FQ_CTX(ctx)),
                               fq_nmod_mpoly_length(GET_FQ_POLY(poly3), GET_FQ_CTX(ctx)));
                        printf("Array multiplication result length: %ld\n", 
                               fq_nmod_mpoly_length(GET_FQ_POLY(poly1), GET_FQ_CTX(ctx)));
                        printf("Standard multiplication result length: %ld\n", 
                               fq_nmod_mpoly_length(test_result, GET_FQ_CTX(ctx)));
                        
                        /* Print first few terms of each result */
                        printf("\nArray multiplication result (first 5 terms):\n");
                        const char *vars[] = {"x", "y", "z", "t", "u", "v", "w"};
                        fq_nmod_mpoly_t temp;
                        fq_nmod_mpoly_init(temp, GET_FQ_CTX(ctx));
                        for (slong i = 0; i < FLINT_MIN(5, fq_nmod_mpoly_length(GET_FQ_POLY(poly1), GET_FQ_CTX(ctx))); i++) {
                            fq_nmod_t coeff;
                            fq_nmod_init(coeff, ctx->field_ctx->ctx.fq_ctx);
                            fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, GET_FQ_POLY(poly1), i, GET_FQ_CTX(ctx));
                            
                            ulong *exp = (ulong *)malloc(ctx->nvars * sizeof(ulong));
                            fq_nmod_mpoly_get_term_exp_ui(exp, GET_FQ_POLY(poly1), i, GET_FQ_CTX(ctx));
                            
                            printf("  Term %ld: ", i);
                            fq_nmod_print_pretty(coeff, ctx->field_ctx->ctx.fq_ctx);
                            printf(" * ");
                            for (slong j = 0; j < ctx->nvars; j++) {
                                if (exp[j] > 0) {
                                    printf("%s^%lu ", vars[j], exp[j]);
                                }
                            }
                            printf("\n");
                            
                            free(exp);
                            fq_nmod_clear(coeff, ctx->field_ctx->ctx.fq_ctx);
                        }
                        
                        printf("\nStandard multiplication result (first 5 terms):\n");
                        for (slong i = 0; i < FLINT_MIN(5, fq_nmod_mpoly_length(test_result, GET_FQ_CTX(ctx))); i++) {
                            fq_nmod_t coeff;
                            fq_nmod_init(coeff, ctx->field_ctx->ctx.fq_ctx);
                            fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, test_result, i, GET_FQ_CTX(ctx));
                            
                            ulong *exp = (ulong *)malloc(ctx->nvars * sizeof(ulong));
                            fq_nmod_mpoly_get_term_exp_ui(exp, test_result, i, GET_FQ_CTX(ctx));
                            
                            printf("  Term %ld: ", i);
                            fq_nmod_print_pretty(coeff, ctx->field_ctx->ctx.fq_ctx);
                            printf(" * ");
                            for (slong j = 0; j < ctx->nvars; j++) {
                                if (exp[j] > 0) {
                                    printf("%s^%lu ", vars[j], exp[j]);
                                }
                            }
                            printf("\n");
                            
                            free(exp);
                            fq_nmod_clear(coeff, ctx->field_ctx->ctx.fq_ctx);
                        }
                        fq_nmod_mpoly_clear(temp, GET_FQ_CTX(ctx));
                        printf("===================================\n\n");
                    }
                    
                    fq_nmod_mpoly_clear(test_result, GET_FQ_CTX(ctx));
#endif                    
                    /* Clean up */
                    gf28_mpoly_clear(A, native_ctx);
                    gf28_mpoly_clear(B, native_ctx);
                    gf28_mpoly_clear(C, native_ctx);
                    gf28_mpoly_ctx_clear(native_ctx);
                    
                    return 1;
                } else {
                    /* Clean up and fall back to standard multiplication */
                    gf28_mpoly_clear(A, native_ctx);
                    gf28_mpoly_clear(B, native_ctx);
                    gf28_mpoly_clear(C, native_ctx);
                    gf28_mpoly_ctx_clear(native_ctx);
                    
                    printf("GF(2^8) array multiplication failed (array too large), using standard method\n");
                }
            }
            /* Fall through to standard multiplication */
            fq_nmod_mpoly_mul(GET_FQ_POLY(poly1), GET_FQ_POLY(poly2),
                             GET_FQ_POLY(poly3), GET_FQ_CTX(ctx));
            return 1;
            
        case FIELD_ID_GF2128:
            if (1 || use_gf2128_array_mul) {
                /* Use optimized GF(2^128) array multiplication */
                gf2128_mpoly_t A, B, C;
                gf2128_mpoly_ctx_t native_ctx;
                
                /* Initialize native context */
                gf2128_mpoly_ctx_init(native_ctx, ctx->nvars, ctx->ord);
                
                /* Initialize native polynomials */
                gf2128_mpoly_init(A, native_ctx);
                gf2128_mpoly_init(B, native_ctx);
                gf2128_mpoly_init(C, native_ctx);
                
                /* Convert from fq_nmod to native format */
                fq_nmod_mpoly_to_gf2128_mpoly(A, GET_FQ_POLY(poly2), 
                                              ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                fq_nmod_mpoly_to_gf2128_mpoly(B, GET_FQ_POLY(poly3), 
                                              ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                
                /* Try array multiplication */
                int success = gf2128_mpoly_mul_array(C, A, B, native_ctx);
                
                if (success) {
                    /* Convert result back to fq_nmod format */
                    gf2128_mpoly_to_fq_nmod_mpoly(GET_FQ_POLY(poly1), C, 
                                                  ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
#if 0                    
                    /* CORRECTNESS TEST: Compare with standard multiplication */
                    fq_nmod_mpoly_t test_result;
                    fq_nmod_mpoly_init(test_result, GET_FQ_CTX(ctx));
                    fq_nmod_mpoly_mul(test_result, GET_FQ_POLY(poly2), 
                                     GET_FQ_POLY(poly3), GET_FQ_CTX(ctx));
                    
                    if (!fq_nmod_mpoly_equal(GET_FQ_POLY(poly1), test_result, GET_FQ_CTX(ctx))) {
                        printf("\n=== GF(2^128) MULTIPLICATION ERROR ===\n");
                        printf("poly1 addr: %p\n", GET_FQ_POLY(poly1));
                        printf("poly2 addr: %p\n", GET_FQ_POLY(poly2));
                        printf("poly3 addr: %p\n", GET_FQ_POLY(poly3));
                        printf("Input A length: %ld, Input B length: %ld\n", 
                               fq_nmod_mpoly_length(GET_FQ_POLY(poly2), GET_FQ_CTX(ctx)),
                               fq_nmod_mpoly_length(GET_FQ_POLY(poly3), GET_FQ_CTX(ctx)));
                        printf("Array multiplication result length: %ld\n", 
                               fq_nmod_mpoly_length(GET_FQ_POLY(poly1), GET_FQ_CTX(ctx)));
                        printf("Standard multiplication result length: %ld\n", 
                               fq_nmod_mpoly_length(test_result, GET_FQ_CTX(ctx)));
                        
                        /* Print first few terms of each result */
                        printf("\nArray multiplication result (first 5 terms):\n");
                        const char *vars[] = {"x", "y", "z", "t", "u", "v", "w"};
                        for (slong i = 0; i < FLINT_MIN(5, fq_nmod_mpoly_length(GET_FQ_POLY(poly1), GET_FQ_CTX(ctx))); i++) {
                            fq_nmod_t coeff;
                            fq_nmod_init(coeff, ctx->field_ctx->ctx.fq_ctx);
                            fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, GET_FQ_POLY(poly1), i, GET_FQ_CTX(ctx));
                            
                            ulong *exp = (ulong *)malloc(ctx->nvars * sizeof(ulong));
                            fq_nmod_mpoly_get_term_exp_ui(exp, GET_FQ_POLY(poly1), i, GET_FQ_CTX(ctx));
                            
                            printf("  Term %ld: ", i);
                            fq_nmod_print_pretty(coeff, ctx->field_ctx->ctx.fq_ctx);
                            printf(" * ");
                            for (slong j = 0; j < ctx->nvars; j++) {
                                if (exp[j] > 0) {
                                    printf("%s^%lu ", vars[j], exp[j]);
                                }
                            }
                            printf("\n");
                            
                            free(exp);
                            fq_nmod_clear(coeff, ctx->field_ctx->ctx.fq_ctx);
                        }
                        
                        printf("\nStandard multiplication result (first 5 terms):\n");
                        for (slong i = 0; i < FLINT_MIN(5, fq_nmod_mpoly_length(test_result, GET_FQ_CTX(ctx))); i++) {
                            fq_nmod_t coeff;
                            fq_nmod_init(coeff, ctx->field_ctx->ctx.fq_ctx);
                            fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, test_result, i, GET_FQ_CTX(ctx));
                            
                            ulong *exp = (ulong *)malloc(ctx->nvars * sizeof(ulong));
                            fq_nmod_mpoly_get_term_exp_ui(exp, test_result, i, GET_FQ_CTX(ctx));
                            
                            printf("  Term %ld: ", i);
                            fq_nmod_print_pretty(coeff, ctx->field_ctx->ctx.fq_ctx);
                            printf(" * ");
                            for (slong j = 0; j < ctx->nvars; j++) {
                                if (exp[j] > 0) {
                                    printf("%s^%lu ", vars[j], exp[j]);
                                }
                            }
                            printf("\n");
                            
                            free(exp);
                            fq_nmod_clear(coeff, ctx->field_ctx->ctx.fq_ctx);
                        }
                        printf("=====================================\n\n");
                    }
                    
                    fq_nmod_mpoly_clear(test_result, GET_FQ_CTX(ctx));
#endif                    
                    /* Clean up */
                    gf2128_mpoly_clear(A, native_ctx);
                    gf2128_mpoly_clear(B, native_ctx);
                    gf2128_mpoly_clear(C, native_ctx);
                    gf2128_mpoly_ctx_clear(native_ctx);
                    
                    return 1;
                } else {
                    /* Clean up and fall back to standard multiplication */
                    gf2128_mpoly_clear(A, native_ctx);
                    gf2128_mpoly_clear(B, native_ctx);
                    gf2128_mpoly_clear(C, native_ctx);
                    gf2128_mpoly_ctx_clear(native_ctx);
                    
                    printf("GF(2^128) array multiplication failed (array too large), using standard method\n");
                }
            }
            /* Fall through to standard multiplication */
            fq_nmod_mpoly_mul(GET_FQ_POLY(poly1), GET_FQ_POLY(poly2),
                             GET_FQ_POLY(poly3), GET_FQ_CTX(ctx));
            return 1;
            
        default:
            fq_nmod_mpoly_mul(GET_FQ_POLY(poly1), GET_FQ_POLY(poly2),
                             GET_FQ_POLY(poly3), GET_FQ_CTX(ctx));
            return 1;
    }
}

int unified_mpoly_mul_old(unified_mpoly_t poly1, const unified_mpoly_t poly2,
                     const unified_mpoly_t poly3) {
    unified_mpoly_ctx_t ctx = poly1->ctx_ptr;
    
    switch (poly1->field_id) {
        case FIELD_ID_NMOD:
            //printf("DEBUG: nmod multiplication, input lengths: %ld * %ld\n", nmod_mpoly_length(GET_NMOD_POLY(poly2), GET_NMOD_CTX(ctx)), nmod_mpoly_length(GET_NMOD_POLY(poly3), GET_NMOD_CTX(ctx)));
            
            nmod_mpoly_mul(GET_NMOD_POLY(poly1), GET_NMOD_POLY(poly2),
                          GET_NMOD_POLY(poly3), GET_NMOD_CTX(ctx));
            
            //printf("DEBUG: nmod multiplication result length: %ld\n", nmod_mpoly_length(GET_NMOD_POLY(poly1), GET_NMOD_CTX(ctx)));
            return 1;
            
        case FIELD_ID_GF28:
            if (1 || use_gf28_array_mul) {
                /* Use optimized GF(2^8) array multiplication */
                gf28_mpoly_t A, B, C;
                gf28_mpoly_ctx_t native_ctx;
                
                /* Initialize native context */
                gf28_mpoly_ctx_init(native_ctx, ctx->nvars, ctx->ord);
                
                /* Initialize native polynomials */
                gf28_mpoly_init(A, native_ctx);
                gf28_mpoly_init(B, native_ctx);
                gf28_mpoly_init(C, native_ctx);
                
                /* Convert from fq_nmod to native format */
                fq_nmod_mpoly_to_gf28_mpoly(A, GET_FQ_POLY(poly2), 
                                            ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                fq_nmod_mpoly_to_gf28_mpoly(B, GET_FQ_POLY(poly3), 
                                            ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                
                /* Try array multiplication */
                int success = gf28_mpoly_mul_array(C, A, B, native_ctx);
                
                if (success) {
                    /* Convert result back to fq_nmod format */
                    gf28_mpoly_to_fq_nmod_mpoly(GET_FQ_POLY(poly1), C, 
                                                ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                    
                    /* Clean up */
                    gf28_mpoly_clear(A, native_ctx);
                    gf28_mpoly_clear(B, native_ctx);
                    gf28_mpoly_clear(C, native_ctx);
                    gf28_mpoly_ctx_clear(native_ctx);
                    
                    return 1;
                } else {
                    /* Clean up and fall back to standard multiplication */
                    gf28_mpoly_clear(A, native_ctx);
                    gf28_mpoly_clear(B, native_ctx);
                    gf28_mpoly_clear(C, native_ctx);
                    gf28_mpoly_ctx_clear(native_ctx);
                    
                    printf("GF(2^8) array multiplication failed (array too large), using standard method\n");
                }
            }
            /* Fall through to standard multiplication */
            fq_nmod_mpoly_mul(GET_FQ_POLY(poly1), GET_FQ_POLY(poly2),
                             GET_FQ_POLY(poly3), GET_FQ_CTX(ctx));
            return 1;
            
        case FIELD_ID_GF2128:
            if (1 || use_gf2128_array_mul) {
                /* Use optimized GF(2^128) array multiplication */
                gf2128_mpoly_t A, B, C;
                gf2128_mpoly_ctx_t native_ctx;
                
                /* Initialize native context */
                gf2128_mpoly_ctx_init(native_ctx, ctx->nvars, ctx->ord);
                
                /* Initialize native polynomials */
                gf2128_mpoly_init(A, native_ctx);
                gf2128_mpoly_init(B, native_ctx);
                gf2128_mpoly_init(C, native_ctx);
                
                /* Convert from fq_nmod to native format */
                fq_nmod_mpoly_to_gf2128_mpoly(A, GET_FQ_POLY(poly2), 
                                              ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                fq_nmod_mpoly_to_gf2128_mpoly(B, GET_FQ_POLY(poly3), 
                                              ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                
                /* Try array multiplication */
                int success = gf2128_mpoly_mul_array(C, A, B, native_ctx);
                
                if (success) {
                    /* Convert result back to fq_nmod format */
                    gf2128_mpoly_to_fq_nmod_mpoly(GET_FQ_POLY(poly1), C, 
                                                  ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                    
                    /* Clean up */
                    gf2128_mpoly_clear(A, native_ctx);
                    gf2128_mpoly_clear(B, native_ctx);
                    gf2128_mpoly_clear(C, native_ctx);
                    gf2128_mpoly_ctx_clear(native_ctx);
                    
                    return 1;
                } else {
                    /* Clean up and fall back to standard multiplication */
                    gf2128_mpoly_clear(A, native_ctx);
                    gf2128_mpoly_clear(B, native_ctx);
                    gf2128_mpoly_clear(C, native_ctx);
                    gf2128_mpoly_ctx_clear(native_ctx);
                    
                    printf("GF(2^128) array multiplication failed (array too large), using standard method\n");
                }
            }
            /* Fall through to standard multiplication */
            fq_nmod_mpoly_mul(GET_FQ_POLY(poly1), GET_FQ_POLY(poly2),
                             GET_FQ_POLY(poly3), GET_FQ_CTX(ctx));
            return 1;
            
        default:
            fq_nmod_mpoly_mul(GET_FQ_POLY(poly1), GET_FQ_POLY(poly2),
                             GET_FQ_POLY(poly3), GET_FQ_CTX(ctx));
            return 1;
    }
}

/* ============================================================================
   COEFFICIENT ACCESS
   ============================================================================ */

void unified_mpoly_set_coeff_ui(unified_mpoly_t poly, const field_elem_u *c,
                               const ulong *exp) {
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_set_coeff_ui_ui(GET_NMOD_POLY(poly), c->nmod, exp,
                                      GET_NMOD_CTX(ctx));
            break;
            
        default:
            {
                fq_nmod_t temp;
                fq_nmod_init(temp, ctx->field_ctx->ctx.fq_ctx);
                field_elem_to_fq_nmod(temp, c, ctx->field_ctx);
                fq_nmod_mpoly_set_coeff_fq_nmod_ui(GET_FQ_POLY(poly), temp, exp,
                                                   GET_FQ_CTX(ctx));
                fq_nmod_clear(temp, ctx->field_ctx->ctx.fq_ctx);
            }
            break;
    }
}

void unified_mpoly_get_coeff_ui(field_elem_u *c, const unified_mpoly_t poly,
                               const ulong *exp) {
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            c->nmod = nmod_mpoly_get_coeff_ui_ui(GET_NMOD_POLY(poly), exp,
                                                GET_NMOD_CTX(ctx));
            break;
            
        default:
            {
                fq_nmod_t temp;
                fq_nmod_init(temp, ctx->field_ctx->ctx.fq_ctx);
                fq_nmod_mpoly_get_coeff_fq_nmod_ui(temp, GET_FQ_POLY(poly), exp,
                                                   GET_FQ_CTX(ctx));
                fq_nmod_to_field_elem(c, temp, ctx->field_ctx);
                fq_nmod_clear(temp, ctx->field_ctx->ctx.fq_ctx);
            }
            break;
    }
}

/* ============================================================================
   DEGREE OPERATIONS
   ============================================================================ */

void unified_mpoly_degrees_si(slong *degs, const unified_mpoly_t poly) {
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_degrees_si(degs, GET_NMOD_POLY(poly), GET_NMOD_CTX(ctx));
            break;
            
        default:
            fq_nmod_mpoly_degrees_si(degs, GET_FQ_POLY(poly), GET_FQ_CTX(ctx));
            break;
    }
}

slong unified_mpoly_total_degree_si(const unified_mpoly_t poly) {
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            return nmod_mpoly_total_degree_si(GET_NMOD_POLY(poly), GET_NMOD_CTX(ctx));
            
        default:
            return fq_nmod_mpoly_total_degree_si(GET_FQ_POLY(poly), GET_FQ_CTX(ctx));
    }
}

/* ============================================================================
   PRINTING
   ============================================================================ */

void unified_mpoly_print_pretty(const unified_mpoly_t poly, const char **vars) {
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_print_pretty(GET_NMOD_POLY(poly), vars, GET_NMOD_CTX(ctx));
            break;
            
        default:
            fq_nmod_mpoly_print_pretty(GET_FQ_POLY(poly), vars, GET_FQ_CTX(ctx));
            break;
    }
}

/* ============================================================================
   SCALAR MULTIPLICATION
   ============================================================================ */

void unified_mpoly_scalar_mul_ui(unified_mpoly_t poly1, const unified_mpoly_t poly2,
                                ulong c) {
    unified_mpoly_ctx_t ctx = poly1->ctx_ptr;
    
    switch (poly1->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_scalar_mul_ui(GET_NMOD_POLY(poly1), GET_NMOD_POLY(poly2),
                                    c, GET_NMOD_CTX(ctx));
            break;
            
        default:
            {
                /* Convert scalar to field element */
                fq_nmod_t scalar;
                fq_nmod_init(scalar, ctx->field_ctx->ctx.fq_ctx);
                fq_nmod_set_ui(scalar, c, ctx->field_ctx->ctx.fq_ctx);
                
                /* Multiply */
                fq_nmod_mpoly_scalar_mul_fq_nmod(GET_FQ_POLY(poly1), GET_FQ_POLY(poly2),
                                               scalar, GET_FQ_CTX(ctx));
                
                fq_nmod_clear(scalar, ctx->field_ctx->ctx.fq_ctx);
            }
            break;
    }
}


/* ============================================================================
   MISSING RING OPERATIONS FOR UNIFIED MPOLY
   ============================================================================ */

void unified_mpoly_swap(unified_mpoly_t poly1, unified_mpoly_t poly2) {
    if (poly1 == poly2) return;
    
    unified_mpoly_struct temp = *poly1;
    *poly1 = *poly2;
    *poly2 = temp;
}

void unified_mpoly_set_fmpz(unified_mpoly_t poly, const fmpz_t c) {
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    unified_mpoly_zero(poly);
    
    if (fmpz_is_zero(c)) return;
    
    /* Set polynomial to constant c */
    field_elem_u coeff;
    ulong *exp = (ulong *)calloc(ctx->nvars, sizeof(ulong));
    
    switch (poly->field_id) {
        case FIELD_ID_NMOD:
            coeff.nmod = fmpz_get_ui(c) % GET_NMOD_CTX(ctx)->mod.n;
            break;
            
        default:
            {
                fq_nmod_t temp;
                fq_nmod_init(temp, ctx->field_ctx->ctx.fq_ctx);
                fq_nmod_set_fmpz(temp, c, ctx->field_ctx->ctx.fq_ctx);
                fq_nmod_to_field_elem(&coeff, temp, ctx->field_ctx);
                fq_nmod_clear(temp, ctx->field_ctx->ctx.fq_ctx);
            }
            break;
    }
    
    unified_mpoly_set_coeff_ui(poly, &coeff, exp);
    free(exp);
}

void unified_mpoly_mul_fmpz(unified_mpoly_t poly1, const unified_mpoly_t poly2,
                           const fmpz_t c) {
    if (fmpz_is_zero(c)) {
        unified_mpoly_zero(poly1);
        return;
    }
    
    if (fmpz_is_one(c)) {
        unified_mpoly_set(poly1, poly2);
        return;
    }
    
    /* Convert fmpz to ulong for scalar multiplication */
    unified_mpoly_ctx_t ctx = poly1->ctx_ptr;
    
    switch (poly1->field_id) {
        case FIELD_ID_NMOD:
            {
                ulong c_mod = fmpz_get_ui(c) % GET_NMOD_CTX(ctx)->mod.n;
                unified_mpoly_scalar_mul_ui(poly1, poly2, c_mod);
            }
            break;
            
        default:
            {
                /* For other fields, convert to field element and multiply */
                unified_mpoly_t temp;
                temp = unified_mpoly_init(ctx);
                unified_mpoly_set_fmpz(temp, c);
                unified_mpoly_mul(poly1, poly2, temp);
                unified_mpoly_clear(temp);
            }
            break;
    }
}

int unified_mpoly_pow_fmpz(unified_mpoly_t poly1, const unified_mpoly_t poly2,
                          const fmpz_t exp) {
    if (fmpz_sgn(exp) < 0) {
        fprintf(stderr, "unified_mpoly_pow_fmpz: negative exponent not supported\n");
        return 0;
    }
    
    if (fmpz_is_zero(exp)) {
        unified_mpoly_one(poly1);
        return 1;
    }
    
    if (fmpz_is_one(exp)) {
        unified_mpoly_set(poly1, poly2);
        return 1;
    }
    
    /* Binary exponentiation */
    unified_mpoly_ctx_t ctx = poly1->ctx_ptr;
    unified_mpoly_t temp = unified_mpoly_init(ctx);
    unified_mpoly_t result = unified_mpoly_init(ctx);
    fmpz_t e;
    
    fmpz_init_set(e, exp);
    unified_mpoly_set(temp, poly2);
    unified_mpoly_one(result);
    
    while (!fmpz_is_zero(e)) {
        if (fmpz_is_odd(e)) {
            unified_mpoly_mul(result, result, temp);
        }
        unified_mpoly_mul(temp, temp, temp);
        fmpz_fdiv_q_2exp(e, e, 1);
    }
    
    unified_mpoly_set(poly1, result);
    
    fmpz_clear(e);
    unified_mpoly_clear(temp);
    unified_mpoly_clear(result);
    
    return 1;
}

slong unified_mpoly_length_wrapper(const void *a, const void *ctx) {
    return unified_mpoly_length((const unified_mpoly_t)a);
}


/* ============================================================================
   DIVISION
   ============================================================================ */

void unified_mpoly_divrem(unified_mpoly_t Q, unified_mpoly_t R,
                        const unified_mpoly_t A, const unified_mpoly_t B) {
    unified_mpoly_ctx_t ctx = A->ctx_ptr;
    
    switch (A->field_id) {
        case FIELD_ID_NMOD:
            nmod_mpoly_divrem(GET_NMOD_POLY(Q), GET_NMOD_POLY(R),
                            GET_NMOD_POLY(A), GET_NMOD_POLY(B),
                            GET_NMOD_CTX(ctx));
            break;
            
        default:
            fq_nmod_mpoly_divrem(GET_FQ_POLY(Q), GET_FQ_POLY(R),
                               GET_FQ_POLY(A), GET_FQ_POLY(B),
                               GET_FQ_CTX(ctx));
            break;
    }
}


/* ============================================================================
   ENHANCED DIVISION OPERATIONS WITH OPTIMIZATIONS
   ============================================================================ */

/* Global flags for division optimizations */
static int use_gf28_div_opt = 0;
static int use_gf2128_div_opt = 0;

/* Enable/disable division optimizations */
void unified_mpoly_enable_div_optimizations(field_id_t field_id, int enable) {
    switch (field_id) {
        case FIELD_ID_GF28:
            use_gf28_div_opt = enable;
            printf("GF(2^8) optimized division: %s\n", enable ? "enabled" : "disabled");
            break;
        case FIELD_ID_GF2128:
            use_gf2128_div_opt = enable;
            printf("GF(2^128) optimized division: %s\n", enable ? "enabled" : "disabled");
            break;
        default:
            break;
    }
}

int unified_mpoly_divides(unified_mpoly_t Q, const unified_mpoly_t A,
                             const unified_mpoly_t B) {
    unified_mpoly_ctx_t ctx = Q->ctx_ptr;
    
    switch (Q->field_id) {
        case FIELD_ID_NMOD:
            return nmod_mpoly_divides(GET_NMOD_POLY(Q), GET_NMOD_POLY(A),
                                     GET_NMOD_POLY(B), GET_NMOD_CTX(ctx));
            
        case FIELD_ID_GF28:
            if (1 || use_gf28_div_opt) {
                /* Use optimized GF(2^8) division */
                gf28_mpoly_t A_native, B_native, Q_native;
                gf28_mpoly_ctx_t native_ctx;
                
                /* Initialize native context */
                gf28_mpoly_ctx_init(native_ctx, ctx->nvars, ctx->ord);
                
                /* Initialize native polynomials */
                gf28_mpoly_init(A_native, native_ctx);
                gf28_mpoly_init(B_native, native_ctx);
                gf28_mpoly_init(Q_native, native_ctx);
                
                /* Convert from fq_nmod to native format */
                fq_nmod_mpoly_to_gf28_mpoly(A_native, GET_FQ_POLY(A), 
                                            ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                fq_nmod_mpoly_to_gf28_mpoly(B_native, GET_FQ_POLY(B), 
                                            ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                
                /* Try optimized division */
                int success = gf28_mpoly_divides(Q_native, A_native, B_native, native_ctx);
                
                if (success) {
                    /* Convert result back to fq_nmod format */
                    gf28_mpoly_to_fq_nmod_mpoly(GET_FQ_POLY(Q), Q_native, 
                                                ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
#if 0                    
                    /* CORRECTNESS TEST: Compare with standard division */
                    fq_nmod_mpoly_t test_Q;
                    fq_nmod_mpoly_init(test_Q, GET_FQ_CTX(ctx));
                    int test_success = fq_nmod_mpoly_divides(test_Q, GET_FQ_POLY(A),
                                                             GET_FQ_POLY(B), GET_FQ_CTX(ctx));
                    
                    if (test_success != success || 
                        (success && !fq_nmod_mpoly_equal(GET_FQ_POLY(Q), test_Q, GET_FQ_CTX(ctx)))) {
                        printf("\n=== GF(2^8) DIVISION ERROR ===\n");
                        printf("Q addr: %p\n", GET_FQ_POLY(Q));
                        printf("A addr: %p\n", GET_FQ_POLY(A));
                        printf("B addr: %p\n", GET_FQ_POLY(B));
                        printf("A length: %ld, B length: %ld\n", 
                               fq_nmod_mpoly_length(GET_FQ_POLY(A), GET_FQ_CTX(ctx)),
                               fq_nmod_mpoly_length(GET_FQ_POLY(B), GET_FQ_CTX(ctx)));
                        printf("Array division success: %d, Standard division success: %d\n", 
                               success, test_success);
                        
                        if (success && test_success) {
                            printf("Array division result length: %ld\n", 
                                   fq_nmod_mpoly_length(GET_FQ_POLY(Q), GET_FQ_CTX(ctx)));
                            printf("Standard division result length: %ld\n", 
                                   fq_nmod_mpoly_length(test_Q, GET_FQ_CTX(ctx)));
                            
                            /* Print first few terms */
                            printf("\nArray division result (first 5 terms):\n");
                            const char *vars[] = {"x", "y", "z", "t", "u", "v", "w"};
                            for (slong i = 0; i < FLINT_MIN(5, fq_nmod_mpoly_length(GET_FQ_POLY(Q), GET_FQ_CTX(ctx))); i++) {
                                fq_nmod_t coeff;
                                fq_nmod_init(coeff, ctx->field_ctx->ctx.fq_ctx);
                                fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, GET_FQ_POLY(Q), i, GET_FQ_CTX(ctx));
                                
                                ulong *exp = (ulong *)malloc(ctx->nvars * sizeof(ulong));
                                fq_nmod_mpoly_get_term_exp_ui(exp, GET_FQ_POLY(Q), i, GET_FQ_CTX(ctx));
                                
                                printf("  Term %ld: ", i);
                                fq_nmod_print_pretty(coeff, ctx->field_ctx->ctx.fq_ctx);
                                printf(" * ");
                                for (slong j = 0; j < ctx->nvars; j++) {
                                    if (exp[j] > 0) {
                                        printf("%s^%lu ", vars[j], exp[j]);
                                    }
                                }
                                printf("\n");
                                
                                free(exp);
                                fq_nmod_clear(coeff, ctx->field_ctx->ctx.fq_ctx);
                            }
                            
                            printf("\nStandard division result (first 5 terms):\n");
                            for (slong i = 0; i < FLINT_MIN(5, fq_nmod_mpoly_length(test_Q, GET_FQ_CTX(ctx))); i++) {
                                fq_nmod_t coeff;
                                fq_nmod_init(coeff, ctx->field_ctx->ctx.fq_ctx);
                                fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, test_Q, i, GET_FQ_CTX(ctx));
                                
                                ulong *exp = (ulong *)malloc(ctx->nvars * sizeof(ulong));
                                fq_nmod_mpoly_get_term_exp_ui(exp, test_Q, i, GET_FQ_CTX(ctx));
                                
                                printf("  Term %ld: ", i);
                                fq_nmod_print_pretty(coeff, ctx->field_ctx->ctx.fq_ctx);
                                printf(" * ");
                                for (slong j = 0; j < ctx->nvars; j++) {
                                    if (exp[j] > 0) {
                                        printf("%s^%lu ", vars[j], exp[j]);
                                    }
                                }
                                printf("\n");
                                
                                free(exp);
                                fq_nmod_clear(coeff, ctx->field_ctx->ctx.fq_ctx);
                            }
                        }
                        printf("==============================\n\n");
                    }
                    
                    fq_nmod_mpoly_clear(test_Q, GET_FQ_CTX(ctx));
#endif
                }
                
                /* Clean up */
                gf28_mpoly_clear(A_native, native_ctx);
                gf28_mpoly_clear(B_native, native_ctx);
                gf28_mpoly_clear(Q_native, native_ctx);
                gf28_mpoly_ctx_clear(native_ctx);
                
                return success;
            }
            /* Fall through to standard division */
            printf("GF(2^8) dense divides failed (array too large), using standard method\n");
            return fq_nmod_mpoly_divides(GET_FQ_POLY(Q), GET_FQ_POLY(A),
                                        GET_FQ_POLY(B), GET_FQ_CTX(ctx));
            
        case FIELD_ID_GF2128:
            if (1 || use_gf2128_div_opt) {
                /* Use optimized GF(2^128) division */
                gf2128_mpoly_t A_native, B_native, Q_native;
                gf2128_mpoly_ctx_t native_ctx;
                
                /* Initialize native context */
                gf2128_mpoly_ctx_init(native_ctx, ctx->nvars, ctx->ord);
                
                /* Initialize native polynomials */
                gf2128_mpoly_init(A_native, native_ctx);
                gf2128_mpoly_init(B_native, native_ctx);
                gf2128_mpoly_init(Q_native, native_ctx);
                
                /* Convert from fq_nmod to native format */
                fq_nmod_mpoly_to_gf2128_mpoly(A_native, GET_FQ_POLY(A), 
                                              ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                fq_nmod_mpoly_to_gf2128_mpoly(B_native, GET_FQ_POLY(B), 
                                              ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                
                /* Try optimized division */
                int success = gf2128_mpoly_divides(Q_native, A_native, B_native, native_ctx);
                
                if (success) {
                    /* Convert result back to fq_nmod format */
                    gf2128_mpoly_to_fq_nmod_mpoly(GET_FQ_POLY(Q), Q_native, 
                                                  ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
#if 0                    
                    /* CORRECTNESS TEST: Compare with standard division */
                    fq_nmod_mpoly_t test_Q;
                    fq_nmod_mpoly_init(test_Q, GET_FQ_CTX(ctx));
                    int test_success = fq_nmod_mpoly_divides(test_Q, GET_FQ_POLY(A),
                                                             GET_FQ_POLY(B), GET_FQ_CTX(ctx));
                    
                    if (test_success != success || 
                        (success && !fq_nmod_mpoly_equal(GET_FQ_POLY(Q), test_Q, GET_FQ_CTX(ctx)))) {
                        printf("\n=== GF(2^128) DIVISION ERROR ===\n");
                        printf("Q addr: %p\n", GET_FQ_POLY(Q));
                        printf("A addr: %p\n", GET_FQ_POLY(A));
                        printf("B addr: %p\n", GET_FQ_POLY(B));
                        printf("A length: %ld, B length: %ld\n", 
                               fq_nmod_mpoly_length(GET_FQ_POLY(A), GET_FQ_CTX(ctx)),
                               fq_nmod_mpoly_length(GET_FQ_POLY(B), GET_FQ_CTX(ctx)));
                        printf("Array division success: %d, Standard division success: %d\n", 
                               success, test_success);
                        
                        if (success && test_success) {
                            printf("Array division result length: %ld\n", 
                                   fq_nmod_mpoly_length(GET_FQ_POLY(Q), GET_FQ_CTX(ctx)));
                            printf("Standard division result length: %ld\n", 
                                   fq_nmod_mpoly_length(test_Q, GET_FQ_CTX(ctx)));
                            
                            /* Print first few terms */
                            printf("\nArray division result (first 5 terms):\n");
                            const char *vars[] = {"x", "y", "z", "t", "u", "v", "w"};
                            for (slong i = 0; i < FLINT_MIN(5, fq_nmod_mpoly_length(GET_FQ_POLY(Q), GET_FQ_CTX(ctx))); i++) {
                                fq_nmod_t coeff;
                                fq_nmod_init(coeff, ctx->field_ctx->ctx.fq_ctx);
                                fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, GET_FQ_POLY(Q), i, GET_FQ_CTX(ctx));
                                
                                ulong *exp = (ulong *)malloc(ctx->nvars * sizeof(ulong));
                                fq_nmod_mpoly_get_term_exp_ui(exp, GET_FQ_POLY(Q), i, GET_FQ_CTX(ctx));
                                
                                printf("  Term %ld: ", i);
                                fq_nmod_print_pretty(coeff, ctx->field_ctx->ctx.fq_ctx);
                                printf(" * ");
                                for (slong j = 0; j < ctx->nvars; j++) {
                                    if (exp[j] > 0) {
                                        printf("%s^%lu ", vars[j], exp[j]);
                                    }
                                }
                                printf("\n");
                                
                                free(exp);
                                fq_nmod_clear(coeff, ctx->field_ctx->ctx.fq_ctx);
                            }
                            
                            printf("\nStandard division result (first 5 terms):\n");
                            for (slong i = 0; i < FLINT_MIN(5, fq_nmod_mpoly_length(test_Q, GET_FQ_CTX(ctx))); i++) {
                                fq_nmod_t coeff;
                                fq_nmod_init(coeff, ctx->field_ctx->ctx.fq_ctx);
                                fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, test_Q, i, GET_FQ_CTX(ctx));
                                
                                ulong *exp = (ulong *)malloc(ctx->nvars * sizeof(ulong));
                                fq_nmod_mpoly_get_term_exp_ui(exp, test_Q, i, GET_FQ_CTX(ctx));
                                
                                printf("  Term %ld: ", i);
                                fq_nmod_print_pretty(coeff, ctx->field_ctx->ctx.fq_ctx);
                                printf(" * ");
                                for (slong j = 0; j < ctx->nvars; j++) {
                                    if (exp[j] > 0) {
                                        printf("%s^%lu ", vars[j], exp[j]);
                                    }
                                }
                                printf("\n");
                                
                                free(exp);
                                fq_nmod_clear(coeff, ctx->field_ctx->ctx.fq_ctx);
                            }
                        }
                        printf("=================================\n\n");
                    }
                    
                    fq_nmod_mpoly_clear(test_Q, GET_FQ_CTX(ctx));
#endif                    
                }
                
                /* Clean up */
                gf2128_mpoly_clear(A_native, native_ctx);
                gf2128_mpoly_clear(B_native, native_ctx);
                gf2128_mpoly_clear(Q_native, native_ctx);
                gf2128_mpoly_ctx_clear(native_ctx);
                
                return success;
            }
            /* Fall through to standard division */
            printf("GF(2^128) dense divides failed (array too large), using standard method\n");
            return fq_nmod_mpoly_divides(GET_FQ_POLY(Q), GET_FQ_POLY(A),
                                        GET_FQ_POLY(B), GET_FQ_CTX(ctx));
            
        default:
            return fq_nmod_mpoly_divides(GET_FQ_POLY(Q), GET_FQ_POLY(A),
                                        GET_FQ_POLY(B), GET_FQ_CTX(ctx));
    }
}

/* Enhanced divides function with optimizations */
int unified_mpoly_divides_old(unified_mpoly_t Q, const unified_mpoly_t A,
                             const unified_mpoly_t B) {
    unified_mpoly_ctx_t ctx = Q->ctx_ptr;
    
    switch (Q->field_id) {
        case FIELD_ID_NMOD:
            return nmod_mpoly_divides(GET_NMOD_POLY(Q), GET_NMOD_POLY(A),
                                     GET_NMOD_POLY(B), GET_NMOD_CTX(ctx));
            
        case FIELD_ID_GF28:
            if (1 || use_gf28_div_opt) {
                /* Use optimized GF(2^8) division */
                gf28_mpoly_t A_native, B_native, Q_native;
                gf28_mpoly_ctx_t native_ctx;
                
                /* Initialize native context */
                gf28_mpoly_ctx_init(native_ctx, ctx->nvars, ctx->ord);
                
                /* Initialize native polynomials */
                gf28_mpoly_init(A_native, native_ctx);
                gf28_mpoly_init(B_native, native_ctx);
                gf28_mpoly_init(Q_native, native_ctx);
                
                /* Convert from fq_nmod to native format */
                fq_nmod_mpoly_to_gf28_mpoly(A_native, GET_FQ_POLY(A), 
                                            ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                fq_nmod_mpoly_to_gf28_mpoly(B_native, GET_FQ_POLY(B), 
                                            ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                
                /* Try optimized division */
                int success = gf28_mpoly_divides(Q_native, A_native, B_native, native_ctx);
                
                if (success) {
                    /* Convert result back to fq_nmod format */
                    gf28_mpoly_to_fq_nmod_mpoly(GET_FQ_POLY(Q), Q_native, 
                                                ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                }
                
                /* Clean up */
                gf28_mpoly_clear(A_native, native_ctx);
                gf28_mpoly_clear(B_native, native_ctx);
                gf28_mpoly_clear(Q_native, native_ctx);
                gf28_mpoly_ctx_clear(native_ctx);
                
                return success;
            }
            /* Fall through to standard division */
            printf("GF(2^8) dense divides failed (array too large), using standard method\n");
            return fq_nmod_mpoly_divides(GET_FQ_POLY(Q), GET_FQ_POLY(A),
                                        GET_FQ_POLY(B), GET_FQ_CTX(ctx));
            
        case FIELD_ID_GF2128:
            if (1 || use_gf2128_div_opt) {
                /* Use optimized GF(2^128) division */
                gf2128_mpoly_t A_native, B_native, Q_native;
                gf2128_mpoly_ctx_t native_ctx;
                
                /* Initialize native context */
                gf2128_mpoly_ctx_init(native_ctx, ctx->nvars, ctx->ord);
                
                /* Initialize native polynomials */
                gf2128_mpoly_init(A_native, native_ctx);
                gf2128_mpoly_init(B_native, native_ctx);
                gf2128_mpoly_init(Q_native, native_ctx);
                
                /* Convert from fq_nmod to native format */
                fq_nmod_mpoly_to_gf2128_mpoly(A_native, GET_FQ_POLY(A), 
                                              ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                fq_nmod_mpoly_to_gf2128_mpoly(B_native, GET_FQ_POLY(B), 
                                              ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                
                /* Try optimized division */
                int success = gf2128_mpoly_divides(Q_native, A_native, B_native, native_ctx);
                
                if (success) {
                    /* Convert result back to fq_nmod format */
                    gf2128_mpoly_to_fq_nmod_mpoly(GET_FQ_POLY(Q), Q_native, 
                                                  ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
                }
                
                /* Clean up */
                gf2128_mpoly_clear(A_native, native_ctx);
                gf2128_mpoly_clear(B_native, native_ctx);
                gf2128_mpoly_clear(Q_native, native_ctx);
                gf2128_mpoly_ctx_clear(native_ctx);
                
                return success;
            }
            /* Fall through to standard division */
            printf("GF(2^128) dense divides failed (array too large), using standard method\n");
            return fq_nmod_mpoly_divides(GET_FQ_POLY(Q), GET_FQ_POLY(A),
                                        GET_FQ_POLY(B), GET_FQ_CTX(ctx));
            
        default:
            return fq_nmod_mpoly_divides(GET_FQ_POLY(Q), GET_FQ_POLY(A),
                                        GET_FQ_POLY(B), GET_FQ_CTX(ctx));
    }
}

/* Enhanced divexact function with optimizations */
int unified_mpoly_divexact(unified_mpoly_t Q, const unified_mpoly_t A,
                              const unified_mpoly_t B) {
    /* First try optimized divides */
    int divides = unified_mpoly_divides(Q, A, B);
    
    if (!divides) {
        fprintf(stderr, "unified_mpoly_divexact: division is not exact\n");
        unified_mpoly_zero(Q);
    }
    
    return divides;
}



/* ============================================================================
   VOID RING INTERFACE FUNCTIONS
   ============================================================================ */

/* These functions adapt unified_mpoly operations to the void* interface */

static void unified_mpoly_void_init(void *a, const void *ctx) {
    unified_mpoly_ctx_t mpoly_ctx = (unified_mpoly_ctx_t)ctx;
    unified_mpoly_t *poly_ptr = (unified_mpoly_t *)a;
    *poly_ptr = unified_mpoly_init(mpoly_ctx);
}

static void unified_mpoly_void_clear(void *a, const void *ctx) {
    unified_mpoly_t *poly_ptr = (unified_mpoly_t *)a;
    unified_mpoly_clear(*poly_ptr);
}

static int unified_mpoly_void_is_zero(const void *a, const void *ctx) {
    const unified_mpoly_t *poly_ptr = (const unified_mpoly_t *)a;
    return unified_mpoly_is_zero(*poly_ptr);
}

static void unified_mpoly_void_zero(void *a, const void *ctx) {
    unified_mpoly_t *poly_ptr = (unified_mpoly_t *)a;
    unified_mpoly_zero(*poly_ptr);
}

static void unified_mpoly_void_one(void *a, const void *ctx) {
    unified_mpoly_t *poly_ptr = (unified_mpoly_t *)a;
    unified_mpoly_one(*poly_ptr);
}

static void unified_mpoly_void_set(void *a, const void *b, const void *ctx) {
    unified_mpoly_t *poly_a = (unified_mpoly_t *)a;
    const unified_mpoly_t *poly_b = (const unified_mpoly_t *)b;
    unified_mpoly_set(*poly_a, *poly_b);
}

static void unified_mpoly_void_set_fmpz(void *a, const fmpz_t b, const void *ctx) {
    unified_mpoly_t *poly_ptr = (unified_mpoly_t *)a;
    unified_mpoly_set_fmpz(*poly_ptr, b);
}

static void unified_mpoly_void_swap(void *a, void *b, const void *ctx) {
    unified_mpoly_t *poly_a = (unified_mpoly_t *)a;
    unified_mpoly_t *poly_b = (unified_mpoly_t *)b;
    unified_mpoly_swap(*poly_a, *poly_b);
}

static void unified_mpoly_void_neg(void *a, const void *b, const void *ctx) {
    unified_mpoly_t *poly_a = (unified_mpoly_t *)a;
    const unified_mpoly_t *poly_b = (const unified_mpoly_t *)b;
    unified_mpoly_neg(*poly_a, *poly_b);
}

static void unified_mpoly_void_add(void *a, const void *b, const void *c,
                                  const void *ctx) {
    unified_mpoly_t *poly_a = (unified_mpoly_t *)a;
    const unified_mpoly_t *poly_b = (const unified_mpoly_t *)b;
    const unified_mpoly_t *poly_c = (const unified_mpoly_t *)c;
    unified_mpoly_add(*poly_a, *poly_b, *poly_c);
}

static void unified_mpoly_void_sub(void *a, const void *b, const void *c,
                                  const void *ctx) {
    unified_mpoly_t *poly_a = (unified_mpoly_t *)a;
    const unified_mpoly_t *poly_b = (const unified_mpoly_t *)b;
    const unified_mpoly_t *poly_c = (const unified_mpoly_t *)c;
    unified_mpoly_sub(*poly_a, *poly_b, *poly_c);
}

static void unified_mpoly_void_mul(void *a, const void *b, const void *c,
                                  const void *ctx) {
    unified_mpoly_t *poly_a = (unified_mpoly_t *)a;
    const unified_mpoly_t *poly_b = (const unified_mpoly_t *)b;
    const unified_mpoly_t *poly_c = (const unified_mpoly_t *)c;
    unified_mpoly_mul(*poly_a, *poly_b, *poly_c);
}

static void unified_mpoly_void_mul_fmpz(void *a, const void *b,
                                       const fmpz_t c, const void *ctx) {
    unified_mpoly_t *poly_a = (unified_mpoly_t *)a;
    const unified_mpoly_t *poly_b = (const unified_mpoly_t *)b;
    unified_mpoly_mul_fmpz(*poly_a, *poly_b, c);
}

static void unified_mpoly_void_divexact(void *a, const void *b,
                                       const void *c, const void *ctx) {
    unified_mpoly_t *poly_a = (unified_mpoly_t *)a;
    const unified_mpoly_t *poly_b = (const unified_mpoly_t *)b;
    const unified_mpoly_t *poly_c = (const unified_mpoly_t *)c;
    if (!unified_mpoly_divexact(*poly_a, *poly_b, *poly_c)) {
        flint_throw(FLINT_ERROR, "unified_mpoly_void_divexact: nonexact");
    }
}

static int unified_mpoly_void_divides(void *a, const void *b,
                                     const void *c, const void *ctx) {
    unified_mpoly_t *poly_a = (unified_mpoly_t *)a;
    const unified_mpoly_t *poly_b = (const unified_mpoly_t *)b;
    const unified_mpoly_t *poly_c = (const unified_mpoly_t *)c;
    return unified_mpoly_divides(*poly_a, *poly_b, *poly_c);
}

static int unified_mpoly_void_pow_fmpz(void *a, const void *b,
                                      const fmpz_t c, const void *ctx) {
    unified_mpoly_t *poly_a = (unified_mpoly_t *)a;
    const unified_mpoly_t *poly_b = (const unified_mpoly_t *)b;
    return unified_mpoly_pow_fmpz(*poly_a, *poly_b, c);
}

static slong unified_mpoly_void_length(const void *a, const void *ctx) {
    const unified_mpoly_t *poly_ptr = (const unified_mpoly_t *)a;
    return unified_mpoly_length(*poly_ptr);
}

/* ============================================================================
   RING STRUCTURE DEFINITION
   ============================================================================ */

typedef struct {
    mpoly_void_ring_t ring;
    unified_mpoly_ctx_t ctx;
} unified_mpoly_ring_struct;

typedef unified_mpoly_ring_struct *unified_mpoly_ring_t;

/* Initialize the ring structure for unified multivariate polynomials */
void unified_mpoly_ring_init(mpoly_void_ring_t R, unified_mpoly_ctx_t ctx) {
    R->elem_size = sizeof(unified_mpoly_t);
    R->ctx = ctx;
    R->init = unified_mpoly_void_init;
    R->clear = unified_mpoly_void_clear;
    R->is_zero = unified_mpoly_void_is_zero;
    R->zero = unified_mpoly_void_zero;
    R->one = unified_mpoly_void_one;
    R->set = unified_mpoly_void_set;
    R->set_fmpz = unified_mpoly_void_set_fmpz;
    R->swap = unified_mpoly_void_swap;
    R->neg = unified_mpoly_void_neg;
    R->add = unified_mpoly_void_add;
    R->sub = unified_mpoly_void_sub;
    R->mul = unified_mpoly_void_mul;
    R->mul_fmpz = unified_mpoly_void_mul_fmpz;
    R->divexact = unified_mpoly_void_divexact;
    R->divides = unified_mpoly_void_divides;
    R->pow_fmpz = unified_mpoly_void_pow_fmpz;
    R->length = unified_mpoly_void_length;
}


void unified_mpoly_ring_clear(unified_mpoly_ring_t R) {
    if (!R) return;
    
    unified_mpoly_ctx_clear(R->ctx);
    free(R);
}

/* ============================================================================
   CONVENIENCE FUNCTIONS
   ============================================================================ */

/* Enable all optimizations for a specific field */
void unified_mpoly_enable_all_optimizations(field_id_t field_id) {
    unified_mpoly_enable_optimizations(field_id, 1);
    unified_mpoly_enable_div_optimizations(field_id, 1);
}

/* Disable all optimizations for a specific field */
void unified_mpoly_disable_all_optimizations(field_id_t field_id) {
    unified_mpoly_enable_optimizations(field_id, 0);
    unified_mpoly_enable_div_optimizations(field_id, 0);
}

/* ============================================================================
   TESTING UTILITIES
   ============================================================================ */

/* Generate random polynomial for testing */
void unified_mpoly_randtest(unified_mpoly_t poly, flint_rand_t state,
                           slong length, slong exp_bound) {
    unified_mpoly_ctx_t ctx = poly->ctx_ptr;
    slong nvars = ctx->nvars;
    ulong *exp = (ulong *)calloc(nvars, sizeof(ulong));
    field_elem_u coeff;
    
    unified_mpoly_zero(poly);
    
    for (slong i = 0; i < length; i++) {
        /* Generate random exponents */
        for (slong j = 0; j < nvars; j++) {
            exp[j] = n_randint(state, exp_bound);
        }
        
        /* Generate random coefficient */
        switch (poly->field_id) {
            case FIELD_ID_NMOD:
                coeff.nmod = n_randint(state, GET_NMOD_CTX(ctx)->mod.n);
                if (coeff.nmod == 0) coeff.nmod = 1;
                break;
                
            case FIELD_ID_GF28:
                coeff.gf28 = n_randint(state, 256);
                if (coeff.gf28 == 0) coeff.gf28 = 1;
                break;
                
            case FIELD_ID_GF2128:
                coeff.gf2128.low = n_randtest(state);
                coeff.gf2128.high = n_randtest(state);
                if (gf2128_is_zero(&coeff.gf2128)) {
                    coeff.gf2128 = gf2128_one();
                }
                break;
                
            default:
                /* For other fields, use generic method */
                {
                    fq_nmod_t temp;
                    fq_nmod_init(temp, ctx->field_ctx->ctx.fq_ctx);
                    fq_nmod_randtest_not_zero(temp, state, ctx->field_ctx->ctx.fq_ctx);
                    fq_nmod_to_field_elem(&coeff, temp, ctx->field_ctx);
                    fq_nmod_clear(temp, ctx->field_ctx->ctx.fq_ctx);
                }
                break;
        }
        
        unified_mpoly_set_coeff_ui(poly, &coeff, exp);
    }
    
    free(exp);
}

/* Check if two polynomials are equal */
int unified_mpoly_equal(const unified_mpoly_t A, const unified_mpoly_t B) {
    if (A->field_id != B->field_id) return 0;
    
    unified_mpoly_ctx_t ctx = A->ctx_ptr;
    
    switch (A->field_id) {
        case FIELD_ID_NMOD:
            return nmod_mpoly_equal(GET_NMOD_POLY(A), GET_NMOD_POLY(B),
                                   GET_NMOD_CTX(ctx));
            
        default:
            return fq_nmod_mpoly_equal(GET_FQ_POLY(A), GET_FQ_POLY(B),
                                      GET_FQ_CTX(ctx));
    }
}

#endif