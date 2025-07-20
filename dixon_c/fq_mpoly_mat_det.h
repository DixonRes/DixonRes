/*
 * Optimized polynomial matrix determinant for small matrices with dense polynomials
 * Focus on optimizing the polynomial operations themselves
 * Enhanced with prime field optimization using nmod_mpoly
 */

#ifndef FQ_MPOLY_MAT_DET_OPTIMIZED_H
#define FQ_MPOLY_MAT_DET_OPTIMIZED_H

#include <omp.h>
#include <flint/fq_nmod_mpoly.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_mat.h>
#include <flint/fq_nmod_poly.h>
#include <flint/nmod_mpoly.h>
#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include "fq_poly_mat_det.h"

// Configuration
#define PARALLEL_THRESHOLD 5
#define MAX_PARALLEL_DEPTH 2
#define UNIVARIATE_THRESHOLD 3

// Debug control
#define DEBUG_FQ_DET 0

#if DEBUG_FQ_DET
#define DET_PRINT(fmt, ...) printf("[FQ_DET] " fmt, ##__VA_ARGS__)
#else
#define DET_PRINT(fmt, ...)
#endif

// ============= Timing Utilities =============

typedef struct {
    double wall_time;
    double cpu_time;
} timing_info_t;

static double get_wall_time() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

static double get_cpu_time() {
    return ((double)clock()) / CLOCKS_PER_SEC;
}

static timing_info_t start_timing() {
    timing_info_t t;
    t.wall_time = get_wall_time();
    t.cpu_time = get_cpu_time();
    return t;
}

static timing_info_t end_timing(timing_info_t start) {
    timing_info_t elapsed;
    elapsed.wall_time = get_wall_time() - start.wall_time;
    elapsed.cpu_time = get_cpu_time() - start.cpu_time;
    return elapsed;
}

static void print_timing(const char* label, timing_info_t elapsed) {
    DET_PRINT("%s: Wall time: %.6f s, CPU time: %.6f s", 
              label, elapsed.wall_time, elapsed.cpu_time);
    if (elapsed.wall_time > 0) {
        DET_PRINT(" (CPU efficiency: %.1f%%)\n", 
                  (elapsed.cpu_time / elapsed.wall_time) * 100.0);
    } else {
        DET_PRINT("\n");
    }
}

// ============= Prime Field Detection =============

static inline int is_prime_field(const fq_nmod_ctx_t ctx) {
    return fq_nmod_ctx_degree(ctx) == 1;
}

// ============= Polynomial Operation Optimizations =============

// Custom multiplication for dense polynomials
static inline void poly_mul_dense_optimized(fq_nmod_mpoly_t c, 
                                           const fq_nmod_mpoly_t a,
                                           const fq_nmod_mpoly_t b,
                                           const fq_nmod_mpoly_ctx_t ctx) {
    slong alen = fq_nmod_mpoly_length(a, ctx);
    slong blen = fq_nmod_mpoly_length(b, ctx);
    
    // For very dense polynomials, try different multiplication algorithms
    if (alen > 100 && blen > 100) {
        // Try Johnson's multiplication for dense polynomials
        fq_nmod_mpoly_mul_johnson(c, a, b, ctx);
    } else {
        // Default multiplication
        fq_nmod_mpoly_mul(c, a, b, ctx);
    }
}

// ============= Prime Field Conversion Functions =============

void fq_mvpoly_to_nmod_mpoly(nmod_mpoly_t mpoly, const fq_mvpoly_t *poly, 
                             nmod_mpoly_ctx_t mpoly_ctx) {
    nmod_mpoly_zero(mpoly, mpoly_ctx);
    
    if (poly->nterms == 0) return;
    
    slong total_vars = poly->nvars + poly->npars;
    
    // Pre-allocate space for better performance
    nmod_mpoly_fit_length(mpoly, poly->nterms, mpoly_ctx);
    
    for (slong i = 0; i < poly->nterms; i++) {
        ulong *exps = (ulong*) flint_calloc(total_vars, sizeof(ulong));
        
        if (poly->terms[i].var_exp && poly->nvars > 0) {
            for (slong j = 0; j < poly->nvars; j++) {
                exps[j] = (ulong)poly->terms[i].var_exp[j];
            }
        }
        
        if (poly->terms[i].par_exp && poly->npars > 0) {
            for (slong j = 0; j < poly->npars; j++) {
                exps[poly->nvars + j] = (ulong)poly->terms[i].par_exp[j];
            }
        }
        
        // For prime fields, extract the coefficient as ulong
        ulong coeff_val = nmod_poly_get_coeff_ui(poly->terms[i].coeff, 0);
        nmod_mpoly_push_term_ui_ui(mpoly, coeff_val, exps, mpoly_ctx);
        flint_free(exps);
    }
    
    nmod_mpoly_sort_terms(mpoly, mpoly_ctx);
    nmod_mpoly_combine_like_terms(mpoly, mpoly_ctx);
}

void nmod_mpoly_to_fq_mvpoly(fq_mvpoly_t *poly, const nmod_mpoly_t mpoly,
                             slong nvars, slong npars, 
                             nmod_mpoly_ctx_t mpoly_ctx, const fq_nmod_ctx_t ctx) {
    fq_mvpoly_init(poly, nvars, npars, ctx);
    
    slong nterms = nmod_mpoly_length(mpoly, mpoly_ctx);
    if (nterms == 0) return;
    
    slong total_vars = nvars + npars;
    
    poly->alloc = nterms;
    poly->terms = (fq_monomial_t*) flint_realloc(poly->terms, poly->alloc * sizeof(fq_monomial_t));
    poly->nterms = nterms;
    
    ulong *exp_buffer = (ulong*) flint_malloc(total_vars * sizeof(ulong));
    
    for (slong i = 0; i < nterms; i++) {
        fq_nmod_init(poly->terms[i].coeff, ctx);
        
        // Get coefficient and convert to fq_nmod
        ulong coeff_val = nmod_mpoly_get_term_coeff_ui(mpoly, i, mpoly_ctx);
        fq_nmod_set_ui(poly->terms[i].coeff, coeff_val, ctx);
        
        nmod_mpoly_get_term_exp_ui(exp_buffer, mpoly, i, mpoly_ctx);
        
        if (nvars > 0) {
            poly->terms[i].var_exp = (slong*) flint_calloc(nvars, sizeof(slong));
            for (slong j = 0; j < nvars; j++) {
                poly->terms[i].var_exp[j] = (slong)exp_buffer[j];
            }
        } else {
            poly->terms[i].var_exp = NULL;
        }
        
        if (npars > 0) {
            poly->terms[i].par_exp = (slong*) flint_calloc(npars, sizeof(slong));
            for (slong j = 0; j < npars; j++) {
                poly->terms[i].par_exp[j] = (slong)exp_buffer[nvars + j];
            }
        } else {
            poly->terms[i].par_exp = NULL;
        }
    }
    
    flint_free(exp_buffer);
}

void fq_matrix_mvpoly_to_nmod_mpoly(nmod_mpoly_t **mpoly_matrix, 
                                    fq_mvpoly_t **mvpoly_matrix, 
                                    slong size, 
                                    nmod_mpoly_ctx_t mpoly_ctx) {
    DET_PRINT("Converting %ld x %ld matrix to nmod_mpoly\n", size, size);
    
    timing_info_t start = start_timing();
    
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            nmod_mpoly_init(mpoly_matrix[i][j], mpoly_ctx);
            fq_mvpoly_to_nmod_mpoly(mpoly_matrix[i][j], &mvpoly_matrix[i][j], mpoly_ctx);
        }
    }
    
    timing_info_t elapsed = end_timing(start);
    print_timing("Matrix conversion to nmod_mpoly", elapsed);
}

// ============= Prime Field Determinant Computation =============

// Hand-optimized 3x3 determinant for nmod_mpoly
static void compute_det_3x3_nmod_optimized(nmod_mpoly_t det, 
                                          nmod_mpoly_t **m,
                                          nmod_mpoly_ctx_t ctx) {
    nmod_mpoly_t t1, t2, t3, t4, t5, t6, sum;
    
    // Initialize temporaries
    nmod_mpoly_init(t1, ctx);
    nmod_mpoly_init(t2, ctx);
    nmod_mpoly_init(t3, ctx);
    nmod_mpoly_init(t4, ctx);
    nmod_mpoly_init(t5, ctx);
    nmod_mpoly_init(t6, ctx);
    nmod_mpoly_init(sum, ctx);
    
    // Compute 6 products in parallel if beneficial
    #pragma omp parallel sections if(omp_get_max_threads() > 2)
    {
        #pragma omp section
        {
            nmod_mpoly_mul(t1, m[1][1], m[2][2], ctx);
            nmod_mpoly_mul(t1, m[0][0], t1, ctx);
        }
        #pragma omp section
        {
            nmod_mpoly_mul(t2, m[1][2], m[2][0], ctx);
            nmod_mpoly_mul(t2, m[0][1], t2, ctx);
        }
        #pragma omp section
        {
            nmod_mpoly_mul(t3, m[1][0], m[2][1], ctx);
            nmod_mpoly_mul(t3, m[0][2], t3, ctx);
        }
        #pragma omp section
        {
            nmod_mpoly_mul(t4, m[1][0], m[2][2], ctx);
            nmod_mpoly_mul(t4, m[0][1], t4, ctx);
        }
        #pragma omp section
        {
            nmod_mpoly_mul(t5, m[1][1], m[2][0], ctx);
            nmod_mpoly_mul(t5, m[0][2], t5, ctx);
        }
        #pragma omp section
        {
            nmod_mpoly_mul(t6, m[1][2], m[2][1], ctx);
            nmod_mpoly_mul(t6, m[0][0], t6, ctx);
        }
    }
    
    // Sum with signs
    nmod_mpoly_add(sum, t1, t2, ctx);
    nmod_mpoly_add(sum, sum, t3, ctx);
    nmod_mpoly_sub(sum, sum, t4, ctx);
    nmod_mpoly_sub(sum, sum, t5, ctx);
    nmod_mpoly_sub(det, sum, t6, ctx);
    
    // Cleanup
    nmod_mpoly_clear(t1, ctx);
    nmod_mpoly_clear(t2, ctx);
    nmod_mpoly_clear(t3, ctx);
    nmod_mpoly_clear(t4, ctx);
    nmod_mpoly_clear(t5, ctx);
    nmod_mpoly_clear(t6, ctx);
    nmod_mpoly_clear(sum, ctx);
}

// Recursive determinant for nmod_mpoly
void compute_nmod_mpoly_det_recursive(nmod_mpoly_t det_result, 
                                     nmod_mpoly_t **mpoly_matrix, 
                                     slong size, 
                                     nmod_mpoly_ctx_t mpoly_ctx) {
    if (size <= 0) {
        nmod_mpoly_one(det_result, mpoly_ctx);
        return;
    }
    
    if (size == 1) {
        nmod_mpoly_set(det_result, mpoly_matrix[0][0], mpoly_ctx);
        return;
    }
    
    if (size == 2) {
        nmod_mpoly_t ad, bc;
        nmod_mpoly_init(ad, mpoly_ctx);
        nmod_mpoly_init(bc, mpoly_ctx);
        
        nmod_mpoly_mul(ad, mpoly_matrix[0][0], mpoly_matrix[1][1], mpoly_ctx);
        nmod_mpoly_mul(bc, mpoly_matrix[0][1], mpoly_matrix[1][0], mpoly_ctx);
        nmod_mpoly_sub(det_result, ad, bc, mpoly_ctx);
        
        nmod_mpoly_clear(ad, mpoly_ctx);
        nmod_mpoly_clear(bc, mpoly_ctx);
        return;
    }
    
    if (size == 3) {
        compute_det_3x3_nmod_optimized(det_result, mpoly_matrix, mpoly_ctx);
        return;
    }
    
    // General case: Laplace expansion
    nmod_mpoly_zero(det_result, mpoly_ctx);
    
    nmod_mpoly_t temp_result, cofactor, subdet;
    nmod_mpoly_init(temp_result, mpoly_ctx);
    nmod_mpoly_init(cofactor, mpoly_ctx);
    nmod_mpoly_init(subdet, mpoly_ctx);
    
    for (slong col = 0; col < size; col++) {
        if (nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
            continue;
        }
        
        // Create submatrix
        nmod_mpoly_t **submatrix = (nmod_mpoly_t**) flint_malloc((size-1) * sizeof(nmod_mpoly_t*));
        for (slong i = 0; i < size-1; i++) {
            submatrix[i] = (nmod_mpoly_t*) flint_malloc((size-1) * sizeof(nmod_mpoly_t));
            for (slong j = 0; j < size-1; j++) {
                nmod_mpoly_init(submatrix[i][j], mpoly_ctx);
            }
        }
        
        // Fill submatrix
        for (slong i = 1; i < size; i++) {
            slong sub_j = 0;
            for (slong j = 0; j < size; j++) {
                if (j != col) {
                    nmod_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j], mpoly_ctx);
                    sub_j++;
                }
            }
        }
        
        // Recursive computation
        compute_nmod_mpoly_det_recursive(subdet, submatrix, size-1, mpoly_ctx);
        
        // Compute cofactor
        nmod_mpoly_mul(cofactor, mpoly_matrix[0][col], subdet, mpoly_ctx);
        
        // Add/subtract to result
        if (col % 2 == 0) {
            nmod_mpoly_add(temp_result, det_result, cofactor, mpoly_ctx);
        } else {
            nmod_mpoly_sub(temp_result, det_result, cofactor, mpoly_ctx);
        }
        nmod_mpoly_set(det_result, temp_result, mpoly_ctx);
        
        // Cleanup submatrix
        for (slong i = 0; i < size-1; i++) {
            for (slong j = 0; j < size-1; j++) {
                nmod_mpoly_clear(submatrix[i][j], mpoly_ctx);
            }
            flint_free(submatrix[i]);
        }
        flint_free(submatrix);
    }
    
    nmod_mpoly_clear(temp_result, mpoly_ctx);
    nmod_mpoly_clear(cofactor, mpoly_ctx);
    nmod_mpoly_clear(subdet, mpoly_ctx);
}

// Parallel determinant computation for nmod_mpoly
void compute_nmod_mpoly_det_parallel_optimized(nmod_mpoly_t det_result, 
                                              nmod_mpoly_t **mpoly_matrix, 
                                              slong size, 
                                              nmod_mpoly_ctx_t mpoly_ctx,
                                              slong depth) {
    if (size < PARALLEL_THRESHOLD || depth >= MAX_PARALLEL_DEPTH) {
        compute_nmod_mpoly_det_recursive(det_result, mpoly_matrix, size, mpoly_ctx);
        return;
    }
    
    DET_PRINT("Parallel nmod computation for %ld x %ld matrix (depth %ld)\n", size, size, depth);
    
    if (size <= 3) {
        compute_nmod_mpoly_det_recursive(det_result, mpoly_matrix, size, mpoly_ctx);
        return;
    }
    
    nmod_mpoly_zero(det_result, mpoly_ctx);
    
    // Count non-zero entries in first row
    slong nonzero_count = 0;
    for (slong col = 0; col < size; col++) {
        if (!nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
            nonzero_count++;
        }
    }
    
    if (nonzero_count < 2) {
        compute_nmod_mpoly_det_recursive(det_result, mpoly_matrix, size, mpoly_ctx);
        return;
    }
    
    // Allocate space for partial results
    nmod_mpoly_t *partial_results = (nmod_mpoly_t*) flint_malloc(size * sizeof(nmod_mpoly_t));
    for (slong i = 0; i < size; i++) {
        nmod_mpoly_init(partial_results[i], mpoly_ctx);
        nmod_mpoly_zero(partial_results[i], mpoly_ctx);
    }
    
    // Parallel computation of cofactors
    #pragma omp parallel for schedule(dynamic, 1) if(nonzero_count >= 3)
    for (slong col = 0; col < size; col++) {
        if (nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
            continue;
        }
        
        nmod_mpoly_t cofactor, subdet;
        nmod_mpoly_init(cofactor, mpoly_ctx);
        nmod_mpoly_init(subdet, mpoly_ctx);
        
        // Create submatrix
        nmod_mpoly_t **submatrix = (nmod_mpoly_t**) flint_malloc((size-1) * sizeof(nmod_mpoly_t*));
        for (slong i = 0; i < size-1; i++) {
            submatrix[i] = (nmod_mpoly_t*) flint_malloc((size-1) * sizeof(nmod_mpoly_t));
            for (slong j = 0; j < size-1; j++) {
                nmod_mpoly_init(submatrix[i][j], mpoly_ctx);
            }
        }
        
        // Fill submatrix
        for (slong i = 1; i < size; i++) {
            slong sub_j = 0;
            for (slong j = 0; j < size; j++) {
                if (j != col) {
                    nmod_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j], mpoly_ctx);
                    sub_j++;
                }
            }
        }
        
        // Recursive call
        if (size-1 >= PARALLEL_THRESHOLD && depth < MAX_PARALLEL_DEPTH) {
            compute_nmod_mpoly_det_parallel_optimized(subdet, submatrix, size-1, mpoly_ctx, depth+1);
        } else {
            compute_nmod_mpoly_det_recursive(subdet, submatrix, size-1, mpoly_ctx);
        }
        
        // Compute cofactor

        if (!nmod_mpoly_mul_array(cofactor, mpoly_matrix[0][col], subdet, mpoly_ctx)) {
            //printf("Back to normal method");
            nmod_mpoly_mul(cofactor, mpoly_matrix[0][col], subdet, mpoly_ctx);
        }
        // Store with sign
        if (col % 2 == 0) {
            nmod_mpoly_set(partial_results[col], cofactor, mpoly_ctx);
        } else {
            nmod_mpoly_neg(partial_results[col], cofactor, mpoly_ctx);
        }
        
        // Cleanup
        for (slong i = 0; i < size-1; i++) {
            for (slong j = 0; j < size-1; j++) {
                nmod_mpoly_clear(submatrix[i][j], mpoly_ctx);
            }
            flint_free(submatrix[i]);
        }
        flint_free(submatrix);
        
        nmod_mpoly_clear(cofactor, mpoly_ctx);
        nmod_mpoly_clear(subdet, mpoly_ctx);
    }
    
    // Sum results
    nmod_mpoly_t temp_sum;
    nmod_mpoly_init(temp_sum, mpoly_ctx);
    
    for (slong col = 0; col < size; col++) {
        if (!nmod_mpoly_is_zero(partial_results[col], mpoly_ctx)) {
            nmod_mpoly_add(temp_sum, det_result, partial_results[col], mpoly_ctx);
            nmod_mpoly_set(det_result, temp_sum, mpoly_ctx);
        }
        nmod_mpoly_clear(partial_results[col], mpoly_ctx);
    }
    
    nmod_mpoly_clear(temp_sum, mpoly_ctx);
    flint_free(partial_results);
}

// ============= Univariate Optimization =============

int is_univariate_matrix(fq_mvpoly_t **matrix, slong size) {
    if (size == 0) return 0;
    slong nvars = matrix[0][0].nvars;
    slong npars = matrix[0][0].npars;
    return (nvars == 1 && npars == 0);
}

void compute_fq_det_univariate_optimized(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size) {
    if (size <= 0) {
        fq_mvpoly_init(result, matrix[0][0].nvars, matrix[0][0].npars, matrix[0][0].ctx);
        return;
    }
    
    DET_PRINT("Using univariate polynomial matrix optimization for %ldx%ld matrix\n", size, size);
    
    const fq_nmod_ctx_struct *ctx = matrix[0][0].ctx;
    fq_mvpoly_init(result, 1, 0, ctx);
    
    timing_info_t start = start_timing();
    
    fq_nmod_poly_mat_t poly_mat;
    fq_nmod_poly_mat_init(poly_mat, size, size, ctx);
    
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            fq_nmod_poly_struct *entry = fq_nmod_poly_mat_entry(poly_mat, i, j);
            fq_nmod_poly_zero(entry, ctx);
            
            for (slong k = 0; k < matrix[i][j].nterms; k++) {
                fq_monomial_t *term = &matrix[i][j].terms[k];
                slong degree = 0;
                if (term->var_exp && matrix[i][j].nvars > 0) {
                    degree = term->var_exp[0];
                }
                fq_nmod_poly_set_coeff(entry, degree, term->coeff, ctx);
            }
        }
    }
    
    fq_nmod_poly_t det_poly;
    fq_nmod_poly_init(det_poly, ctx);
    
    fq_nmod_poly_mat_det_iter(det_poly, poly_mat, ctx);
    
    timing_info_t conv_elapsed = end_timing(start);
    print_timing("Univariate matrix determinant", conv_elapsed);
    
    slong degree = fq_nmod_poly_degree(det_poly, ctx);
    if (degree >= 0) {
        for (slong d = 0; d <= degree; d++) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ctx);
            fq_nmod_poly_get_coeff(coeff, det_poly, d, ctx);
            
            if (!fq_nmod_is_zero(coeff, ctx)) {
                slong *var_exp = (slong*) flint_calloc(1, sizeof(slong));
                var_exp[0] = d;
                fq_mvpoly_add_term(result, var_exp, NULL, coeff);
                flint_free(var_exp);
            }
            
            fq_nmod_clear(coeff, ctx);
        }
    }
    
    fq_nmod_poly_clear(det_poly, ctx);
    fq_nmod_poly_mat_clear(poly_mat, ctx);
}

// ============= Conversion Functions =============

void fq_mvpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t mpoly, const fq_mvpoly_t *poly, 
                               fq_nmod_mpoly_ctx_t mpoly_ctx) {
    fq_nmod_mpoly_zero(mpoly, mpoly_ctx);
    
    if (poly->nterms == 0) return;
    
    slong total_vars = poly->nvars + poly->npars;
    
    // Pre-allocate space for better performance
    fq_nmod_mpoly_fit_length(mpoly, poly->nterms, mpoly_ctx);
    
    for (slong i = 0; i < poly->nterms; i++) {
        ulong *exps = (ulong*) flint_calloc(total_vars, sizeof(ulong));
        
        if (poly->terms[i].var_exp && poly->nvars > 0) {
            for (slong j = 0; j < poly->nvars; j++) {
                exps[j] = (ulong)poly->terms[i].var_exp[j];
            }
        }
        
        if (poly->terms[i].par_exp && poly->npars > 0) {
            for (slong j = 0; j < poly->npars; j++) {
                exps[poly->nvars + j] = (ulong)poly->terms[i].par_exp[j];
            }
        }
        
        fq_nmod_mpoly_push_term_fq_nmod_ui(mpoly, poly->terms[i].coeff, exps, mpoly_ctx);
        flint_free(exps);
    }
    
    fq_nmod_mpoly_sort_terms(mpoly, mpoly_ctx);
    fq_nmod_mpoly_combine_like_terms(mpoly, mpoly_ctx);
}

void fq_nmod_mpoly_to_fq_mvpoly(fq_mvpoly_t *poly, const fq_nmod_mpoly_t mpoly,
                                slong nvars, slong npars, 
                                fq_nmod_mpoly_ctx_t mpoly_ctx, const fq_nmod_ctx_t ctx) {
    fq_mvpoly_init(poly, nvars, npars, ctx);
    
    slong nterms = fq_nmod_mpoly_length(mpoly, mpoly_ctx);
    if (nterms == 0) return;
    
    slong total_vars = nvars + npars;
    
    poly->alloc = nterms;
    poly->terms = (fq_monomial_t*) flint_realloc(poly->terms, poly->alloc * sizeof(fq_monomial_t));
    poly->nterms = nterms;
    
    ulong *exp_buffer = (ulong*) flint_malloc(total_vars * sizeof(ulong));
    
    for (slong i = 0; i < nterms; i++) {
        fq_nmod_init(poly->terms[i].coeff, ctx);
        
        fq_nmod_mpoly_get_term_coeff_fq_nmod(poly->terms[i].coeff, mpoly, i, mpoly_ctx);
        fq_nmod_mpoly_get_term_exp_ui(exp_buffer, mpoly, i, mpoly_ctx);
        
        if (nvars > 0) {
            poly->terms[i].var_exp = (slong*) flint_calloc(nvars, sizeof(slong));
            for (slong j = 0; j < nvars; j++) {
                poly->terms[i].var_exp[j] = (slong)exp_buffer[j];
            }
        } else {
            poly->terms[i].var_exp = NULL;
        }
        
        if (npars > 0) {
            poly->terms[i].par_exp = (slong*) flint_calloc(npars, sizeof(slong));
            for (slong j = 0; j < npars; j++) {
                poly->terms[i].par_exp[j] = (slong)exp_buffer[nvars + j];
            }
        } else {
            poly->terms[i].par_exp = NULL;
        }
    }
    
    flint_free(exp_buffer);
}

void fq_matrix_mvpoly_to_mpoly(fq_nmod_mpoly_t **mpoly_matrix, 
                              fq_mvpoly_t **mvpoly_matrix, 
                              slong size, 
                              fq_nmod_mpoly_ctx_t mpoly_ctx) {
    DET_PRINT("Converting %ld x %ld matrix\n", size, size);
    
    timing_info_t start = start_timing();
    
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            fq_nmod_mpoly_init(mpoly_matrix[i][j], mpoly_ctx);
            fq_mvpoly_to_fq_nmod_mpoly(mpoly_matrix[i][j], &mvpoly_matrix[i][j], mpoly_ctx);
        }
    }
    
    timing_info_t elapsed = end_timing(start);
    print_timing("Matrix conversion", elapsed);
}

// ============= Optimized Determinant Computation =============

// Hand-optimized 3x3 determinant
static void compute_det_3x3_optimized(fq_nmod_mpoly_t det, 
                                     fq_nmod_mpoly_t **m,
                                     fq_nmod_mpoly_ctx_t ctx) {
    fq_nmod_mpoly_t t1, t2, t3, t4, t5, t6, sum;
    
    // Initialize temporaries
    fq_nmod_mpoly_init(t1, ctx);
    fq_nmod_mpoly_init(t2, ctx);
    fq_nmod_mpoly_init(t3, ctx);
    fq_nmod_mpoly_init(t4, ctx);
    fq_nmod_mpoly_init(t5, ctx);
    fq_nmod_mpoly_init(t6, ctx);
    fq_nmod_mpoly_init(sum, ctx);
    
    // Compute 6 products in parallel if beneficial
    #pragma omp parallel sections if(omp_get_max_threads() > 2)
    {
        #pragma omp section
        {
            fq_nmod_mpoly_mul(t1, m[1][1], m[2][2], ctx);
            poly_mul_dense_optimized(t1, m[0][0], t1, ctx);
        }
        #pragma omp section
        {
            fq_nmod_mpoly_mul(t2, m[1][2], m[2][0], ctx);
            poly_mul_dense_optimized(t2, m[0][1], t2, ctx);
        }
        #pragma omp section
        {
            fq_nmod_mpoly_mul(t3, m[1][0], m[2][1], ctx);
            poly_mul_dense_optimized(t3, m[0][2], t3, ctx);
        }
        #pragma omp section
        {
            fq_nmod_mpoly_mul(t4, m[1][0], m[2][2], ctx);
            poly_mul_dense_optimized(t4, m[0][1], t4, ctx);
        }
        #pragma omp section
        {
            fq_nmod_mpoly_mul(t5, m[1][1], m[2][0], ctx);
            poly_mul_dense_optimized(t5, m[0][2], t5, ctx);
        }
        #pragma omp section
        {
            fq_nmod_mpoly_mul(t6, m[1][2], m[2][1], ctx);
            poly_mul_dense_optimized(t6, m[0][0], t6, ctx);
        }
    }
    
    // Sum with signs
    fq_nmod_mpoly_add(sum, t1, t2, ctx);
    fq_nmod_mpoly_add(sum, sum, t3, ctx);
    fq_nmod_mpoly_sub(sum, sum, t4, ctx);
    fq_nmod_mpoly_sub(sum, sum, t5, ctx);
    fq_nmod_mpoly_sub(det, sum, t6, ctx);
    
    // Cleanup
    fq_nmod_mpoly_clear(t1, ctx);
    fq_nmod_mpoly_clear(t2, ctx);
    fq_nmod_mpoly_clear(t3, ctx);
    fq_nmod_mpoly_clear(t4, ctx);
    fq_nmod_mpoly_clear(t5, ctx);
    fq_nmod_mpoly_clear(t6, ctx);
    fq_nmod_mpoly_clear(sum, ctx);
}

// Recursive determinant with optimizations
void compute_fq_nmod_mpoly_det_recursive(fq_nmod_mpoly_t det_result, 
                                        fq_nmod_mpoly_t **mpoly_matrix, 
                                        slong size, 
                                        fq_nmod_mpoly_ctx_t mpoly_ctx) {
    if (size <= 0) {
        fq_nmod_mpoly_one(det_result, mpoly_ctx);
        return;
    }
    
    if (size == 1) {
        fq_nmod_mpoly_set(det_result, mpoly_matrix[0][0], mpoly_ctx);
        return;
    }
    
    if (size == 2) {
        fq_nmod_mpoly_t ad, bc;
        fq_nmod_mpoly_init(ad, mpoly_ctx);
        fq_nmod_mpoly_init(bc, mpoly_ctx);
        
        poly_mul_dense_optimized(ad, mpoly_matrix[0][0], mpoly_matrix[1][1], mpoly_ctx);
        poly_mul_dense_optimized(bc, mpoly_matrix[0][1], mpoly_matrix[1][0], mpoly_ctx);
        fq_nmod_mpoly_sub(det_result, ad, bc, mpoly_ctx);
        
        fq_nmod_mpoly_clear(ad, mpoly_ctx);
        fq_nmod_mpoly_clear(bc, mpoly_ctx);
        return;
    }
    
    if (size == 3) {
        compute_det_3x3_optimized(det_result, mpoly_matrix, mpoly_ctx);
        return;
    }
    
    // General case: Laplace expansion
    fq_nmod_mpoly_zero(det_result, mpoly_ctx);
    
    fq_nmod_mpoly_t temp_result, cofactor, subdet;
    fq_nmod_mpoly_init(temp_result, mpoly_ctx);
    fq_nmod_mpoly_init(cofactor, mpoly_ctx);
    fq_nmod_mpoly_init(subdet, mpoly_ctx);
    
    for (slong col = 0; col < size; col++) {
        if (fq_nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
            continue;
        }
        
        // Create submatrix
        fq_nmod_mpoly_t **submatrix = (fq_nmod_mpoly_t**) flint_malloc((size-1) * sizeof(fq_nmod_mpoly_t*));
        for (slong i = 0; i < size-1; i++) {
            submatrix[i] = (fq_nmod_mpoly_t*) flint_malloc((size-1) * sizeof(fq_nmod_mpoly_t));
            for (slong j = 0; j < size-1; j++) {
                fq_nmod_mpoly_init(submatrix[i][j], mpoly_ctx);
            }
        }
        
        // Fill submatrix
        for (slong i = 1; i < size; i++) {
            slong sub_j = 0;
            for (slong j = 0; j < size; j++) {
                if (j != col) {
                    fq_nmod_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j], mpoly_ctx);
                    sub_j++;
                }
            }
        }
        
        // Recursive computation
        compute_fq_nmod_mpoly_det_recursive(subdet, submatrix, size-1, mpoly_ctx);
        
        // Compute cofactor
        poly_mul_dense_optimized(cofactor, mpoly_matrix[0][col], subdet, mpoly_ctx);
        
        // Add/subtract to result
        if (col % 2 == 0) {
            fq_nmod_mpoly_add(temp_result, det_result, cofactor, mpoly_ctx);
        } else {
            fq_nmod_mpoly_sub(temp_result, det_result, cofactor, mpoly_ctx);
        }
        fq_nmod_mpoly_set(det_result, temp_result, mpoly_ctx);
        
        // Cleanup submatrix
        for (slong i = 0; i < size-1; i++) {
            for (slong j = 0; j < size-1; j++) {
                fq_nmod_mpoly_clear(submatrix[i][j], mpoly_ctx);
            }
            flint_free(submatrix[i]);
        }
        flint_free(submatrix);
    }
    
    fq_nmod_mpoly_clear(temp_result, mpoly_ctx);
    fq_nmod_mpoly_clear(cofactor, mpoly_ctx);
    fq_nmod_mpoly_clear(subdet, mpoly_ctx);
}

// Parallel determinant computation
void compute_fq_nmod_mpoly_det_parallel_optimized(fq_nmod_mpoly_t det_result, 
                                                  fq_nmod_mpoly_t **mpoly_matrix, 
                                                  slong size, 
                                                  fq_nmod_mpoly_ctx_t mpoly_ctx,
                                                  slong depth) {
    if (size < PARALLEL_THRESHOLD || depth >= MAX_PARALLEL_DEPTH) {
        compute_fq_nmod_mpoly_det_recursive(det_result, mpoly_matrix, size, mpoly_ctx);
        return;
    }
    
    DET_PRINT("Parallel computation for %ld x %ld matrix (depth %ld)\n", size, size, depth);
    
    if (size <= 3) {
        compute_fq_nmod_mpoly_det_recursive(det_result, mpoly_matrix, size, mpoly_ctx);
        return;
    }
    
    fq_nmod_mpoly_zero(det_result, mpoly_ctx);
    
    // Count non-zero entries in first row
    slong nonzero_count = 0;
    for (slong col = 0; col < size; col++) {
        if (!fq_nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
            nonzero_count++;
        }
    }
    
    if (nonzero_count < 2) {
        compute_fq_nmod_mpoly_det_recursive(det_result, mpoly_matrix, size, mpoly_ctx);
        return;
    }
    
    // Allocate space for partial results
    fq_nmod_mpoly_t *partial_results = (fq_nmod_mpoly_t*) flint_malloc(size * sizeof(fq_nmod_mpoly_t));
    for (slong i = 0; i < size; i++) {
        fq_nmod_mpoly_init(partial_results[i], mpoly_ctx);
        fq_nmod_mpoly_zero(partial_results[i], mpoly_ctx);
    }
    
    // Parallel computation of cofactors
    #pragma omp parallel for schedule(dynamic, 1) if(nonzero_count >= 3)
    for (slong col = 0; col < size; col++) {
        if (fq_nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
            continue;
        }
        
        fq_nmod_mpoly_t cofactor, subdet;
        fq_nmod_mpoly_init(cofactor, mpoly_ctx);
        fq_nmod_mpoly_init(subdet, mpoly_ctx);
        
        // Create submatrix
        fq_nmod_mpoly_t **submatrix = (fq_nmod_mpoly_t**) flint_malloc((size-1) * sizeof(fq_nmod_mpoly_t*));
        for (slong i = 0; i < size-1; i++) {
            submatrix[i] = (fq_nmod_mpoly_t*) flint_malloc((size-1) * sizeof(fq_nmod_mpoly_t));
            for (slong j = 0; j < size-1; j++) {
                fq_nmod_mpoly_init(submatrix[i][j], mpoly_ctx);
            }
        }
        
        // Fill submatrix
        for (slong i = 1; i < size; i++) {
            slong sub_j = 0;
            for (slong j = 0; j < size; j++) {
                if (j != col) {
                    fq_nmod_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j], mpoly_ctx);
                    sub_j++;
                }
            }
        }
        
        // Recursive call
        if (size-1 >= PARALLEL_THRESHOLD && depth < MAX_PARALLEL_DEPTH) {
            compute_fq_nmod_mpoly_det_parallel_optimized(subdet, submatrix, size-1, mpoly_ctx, depth+1);
        } else {
            compute_fq_nmod_mpoly_det_recursive(subdet, submatrix, size-1, mpoly_ctx);
        }
        
        // Compute cofactor
        poly_mul_dense_optimized(cofactor, mpoly_matrix[0][col], subdet, mpoly_ctx);
        
        // Store with sign
        if (col % 2 == 0) {
            fq_nmod_mpoly_set(partial_results[col], cofactor, mpoly_ctx);
        } else {
            fq_nmod_mpoly_neg(partial_results[col], cofactor, mpoly_ctx);
        }
        
        // Cleanup
        for (slong i = 0; i < size-1; i++) {
            for (slong j = 0; j < size-1; j++) {
                fq_nmod_mpoly_clear(submatrix[i][j], mpoly_ctx);
            }
            flint_free(submatrix[i]);
        }
        flint_free(submatrix);
        
        fq_nmod_mpoly_clear(cofactor, mpoly_ctx);
        fq_nmod_mpoly_clear(subdet, mpoly_ctx);
    }
    
    // Sum results
    fq_nmod_mpoly_t temp_sum;
    fq_nmod_mpoly_init(temp_sum, mpoly_ctx);
    
    for (slong col = 0; col < size; col++) {
        if (!fq_nmod_mpoly_is_zero(partial_results[col], mpoly_ctx)) {
            fq_nmod_mpoly_add(temp_sum, det_result, partial_results[col], mpoly_ctx);
            fq_nmod_mpoly_set(det_result, temp_sum, mpoly_ctx);
        }
        fq_nmod_mpoly_clear(partial_results[col], mpoly_ctx);
    }
    
    fq_nmod_mpoly_clear(temp_sum, mpoly_ctx);
    flint_free(partial_results);
}

// ============= Main Interface with Prime Field Detection =============

void compute_fq_det_recursive_flint(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size) {
    if (size <= 0) {
        fq_mvpoly_init(result, matrix[0][0].nvars, matrix[0][0].npars, matrix[0][0].ctx);
        return;
    }
    
    timing_info_t total_start = start_timing();
    
    slong max_threads = omp_get_max_threads();
    DET_PRINT("Computing %ldx%ld determinant (OpenMP: %ld threads available)\n", 
              size, size, max_threads);
    
    slong nvars = matrix[0][0].nvars;
    slong npars = matrix[0][0].npars;
    slong total_vars = nvars + npars;
    const fq_nmod_ctx_struct *ctx = matrix[0][0].ctx;
    
    fq_mvpoly_init(result, nvars, npars, ctx);
    
    // Check for univariate optimization
    if (is_univariate_matrix(matrix, size) && size >= UNIVARIATE_THRESHOLD) {
        DET_PRINT("Detected univariate matrix, using specialized optimization\n");
        compute_fq_det_univariate_optimized(result, matrix, size);
        
        timing_info_t total_elapsed = end_timing(total_start);
        print_timing("Total univariate computation", total_elapsed);
        return;
    }
    
    // Check if we can use prime field optimization
    if (is_prime_field(ctx)) {
        DET_PRINT("Detected prime field, using nmod_mpoly optimization\n");
        
        // Get the modulus
        mp_limb_t modulus = fq_nmod_ctx_modulus(ctx)->mod.n;
        
        // Create nmod_mpoly context
        nmod_mpoly_ctx_t nmod_ctx;
        nmod_mpoly_ctx_init(nmod_ctx, total_vars, ORD_LEX, modulus);
        
        // Allocate nmod_mpoly matrix
        nmod_mpoly_t **nmod_matrix = (nmod_mpoly_t**) flint_malloc(size * sizeof(nmod_mpoly_t*));
        for (slong i = 0; i < size; i++) {
            nmod_matrix[i] = (nmod_mpoly_t*) flint_malloc(size * sizeof(nmod_mpoly_t));
        }
        
        // Convert matrix
        fq_matrix_mvpoly_to_nmod_mpoly(nmod_matrix, matrix, size, nmod_ctx);
        
        // Compute determinant
        nmod_mpoly_t det_nmod;
        nmod_mpoly_init(det_nmod, nmod_ctx);
        
        timing_info_t det_start = start_timing();
        if (size >= PARALLEL_THRESHOLD && max_threads > 1) {
            DET_PRINT("Using parallel nmod determinant computation\n");
            compute_nmod_mpoly_det_parallel_optimized(det_nmod, nmod_matrix, size, nmod_ctx, 0);
        } else {
            DET_PRINT("Using serial nmod determinant computation\n");
            compute_nmod_mpoly_det_recursive(det_nmod, nmod_matrix, size, nmod_ctx);
        }
        timing_info_t det_elapsed = end_timing(det_start);
        print_timing("Prime field determinant computation", det_elapsed);
        
        // Convert result back
        timing_info_t result_start = start_timing();
        nmod_mpoly_to_fq_mvpoly(result, det_nmod, nvars, npars, nmod_ctx, ctx);
        timing_info_t result_elapsed = end_timing(result_start);
        print_timing("Result conversion from nmod", result_elapsed);
        
        DET_PRINT("Final result: %ld terms\n", result->nterms);
        
        // Cleanup
        for (slong i = 0; i < size; i++) {
            for (slong j = 0; j < size; j++) {
                nmod_mpoly_clear(nmod_matrix[i][j], nmod_ctx);
            }
            flint_free(nmod_matrix[i]);
        }
        flint_free(nmod_matrix);
        
        nmod_mpoly_clear(det_nmod, nmod_ctx);
        nmod_mpoly_ctx_clear(nmod_ctx);
        
        timing_info_t total_elapsed = end_timing(total_start);
        printf("Total computation (prime field): Wall time: %.6f s, CPU time: %.6f s", 
               total_elapsed.wall_time, total_elapsed.cpu_time);
        if (total_elapsed.wall_time > 0) {
            printf(" (CPU efficiency: %.1f%%)\n", 
                   (total_elapsed.cpu_time / total_elapsed.wall_time) * 100.0);
        } else {
            printf("\n");
        }
        return;
    }
    
    // General case: Use fq_nmod_mpoly
    DET_PRINT("Using general multivariate polynomial matrix method\n");
    
    // Create mpoly context
    fq_nmod_mpoly_ctx_t mpoly_ctx;
    fq_nmod_mpoly_ctx_init(mpoly_ctx, total_vars, ORD_LEX, ctx);
    
    // Allocate mpoly matrix
    fq_nmod_mpoly_t **mpoly_matrix = (fq_nmod_mpoly_t**) flint_malloc(size * sizeof(fq_nmod_mpoly_t*));
    for (slong i = 0; i < size; i++) {
        mpoly_matrix[i] = (fq_nmod_mpoly_t*) flint_malloc(size * sizeof(fq_nmod_mpoly_t));
    }
    
    // Convert matrix
    fq_matrix_mvpoly_to_mpoly(mpoly_matrix, matrix, size, mpoly_ctx);
    
    // Compute determinant
    fq_nmod_mpoly_t det_mpoly;
    fq_nmod_mpoly_init(det_mpoly, mpoly_ctx);
    
    timing_info_t det_start = start_timing();
    if (size >= PARALLEL_THRESHOLD && max_threads > 1) {
        DET_PRINT("Using parallel determinant computation\n");
        compute_fq_nmod_mpoly_det_parallel_optimized(det_mpoly, mpoly_matrix, size, mpoly_ctx, 0);
    } else {
        DET_PRINT("Using serial determinant computation\n");
        compute_fq_nmod_mpoly_det_recursive(det_mpoly, mpoly_matrix, size, mpoly_ctx);
    }
    timing_info_t det_elapsed = end_timing(det_start);
    print_timing("Determinant computation", det_elapsed);
    
    // Convert result back
    timing_info_t result_start = start_timing();
    fq_nmod_mpoly_to_fq_mvpoly(result, det_mpoly, nvars, npars, mpoly_ctx, ctx);
    timing_info_t result_elapsed = end_timing(result_start);
    print_timing("Result conversion", result_elapsed);
    
    DET_PRINT("Final result: %ld terms\n", result->nterms);
    
    // Cleanup
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            fq_nmod_mpoly_clear(mpoly_matrix[i][j], mpoly_ctx);
        }
        flint_free(mpoly_matrix[i]);
    }
    flint_free(mpoly_matrix);
    
    fq_nmod_mpoly_clear(det_mpoly, mpoly_ctx);
    fq_nmod_mpoly_ctx_clear(mpoly_ctx);
    
    timing_info_t total_elapsed = end_timing(total_start);
    printf("Total computation: Wall time: %.6f s, CPU time: %.6f s", 
           total_elapsed.wall_time, total_elapsed.cpu_time);
    if (total_elapsed.wall_time > 0) {
        printf(" (CPU efficiency: %.1f%%)\n", 
               (total_elapsed.cpu_time / total_elapsed.wall_time) * 100.0);
    } else {
        printf("\n");
    }
}

// Compatibility interfaces
void compute_fq_det_recursive(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size) {
    compute_fq_det_recursive_flint(result, matrix, size);
}

void compute_fq_det_polynomial_matrix_simple(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size) {
    if (is_univariate_matrix(matrix, size)) {
        DET_PRINT("Using simple polynomial matrix method for univariate case\n");
        compute_fq_det_univariate_optimized(result, matrix, size);
    } else {
        DET_PRINT("Matrix is not univariate, falling back to general method\n");
        compute_fq_det_recursive_flint(result, matrix, size);
    }
}

#endif // FQ_MPOLY_MAT_DET_OPTIMIZED_H
