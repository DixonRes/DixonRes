/*
 * Optimized polynomial matrix determinant for small matrices with dense polynomials
 * Focus on optimizing the polynomial operations themselves
 * Enhanced with prime field optimization using nmod_mpoly
 * Added multiple algorithm options including polynomial recursive
 */

#ifndef FQ_MPOLY_MAT_DET_H
#define FQ_MPOLY_MAT_DET_H

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
#include "fq_multivariate_interpolation.h"
// Algorithm selection options
#define DET_ALGORITHM_RECURSIVE       0  // Original recursive expansion algorithm
#define DET_ALGORITHM_INTERPOLATION   1  // Multivariate interpolation algorithm
#define DET_ALGORITHM_KRONECKER       2  // Kronecker substitution to univariate
#define DET_ALGORITHM_POLY_RECURSIVE  3  // Convert to fq_nmod_poly and use recursive

// Default algorithm selection
#ifndef DET_ALGORITHM
#define DET_ALGORITHM DET_ALGORITHM_RECURSIVE
#endif

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

// Enable nested parallelism at program start
void init_nested_parallelism(void) {
    #ifdef _OPENMP
    omp_set_nested(1);  // Enable nested parallelism
    omp_set_max_active_levels(2);  // Allow 2 levels of parallelism
    printf("Nested parallelism enabled with max 2 levels\n");
    printf("Max threads available: %d\n", omp_get_max_threads());
    #endif
}

static void compute_kronecker_bounds(slong *var_bounds, fq_mvpoly_t **matrix, 
                                    slong size, slong nvars, slong npars);
static void mvpoly_to_univariate_kronecker(fq_nmod_poly_t uni_poly,
                                          const fq_mvpoly_t *mv_poly,
                                          const slong *substitution_powers,
                                          const fq_nmod_ctx_t ctx);
// Convert univariate polynomial back to multivariate
static void univariate_to_mvpoly_kronecker(fq_mvpoly_t *mv_poly,
                                          const fq_nmod_poly_t uni_poly,
                                          const slong *substitution_powers,
                                          const slong *var_bounds,
                                          slong nvars, slong npars,
                                          const fq_nmod_ctx_t ctx);
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

// ============= Conversion Functions for Polynomial Recursive =============

// Check if all polynomials in matrix use only the first variable
static int is_essentially_univariate(fq_mvpoly_t **matrix, slong size, slong *active_var) {
    if (size == 0) return 0;
    
    *active_var = -1;
    
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            fq_mvpoly_t *poly = &matrix[i][j];
            
            for (slong t = 0; t < poly->nterms; t++) {
                // Check variables
                if (poly->terms[t].var_exp) {
                    for (slong v = 0; v < poly->nvars; v++) {
                        if (poly->terms[t].var_exp[v] > 0) {
                            if (*active_var == -1) {
                                *active_var = v;
                            } else if (*active_var != v) {
                                return 0; // Multiple variables used
                            }
                        }
                    }
                }
                
                // Check parameters - treat them as additional variables
                if (poly->terms[t].par_exp) {
                    for (slong p = 0; p < poly->npars; p++) {
                        if (poly->terms[t].par_exp[p] > 0) {
                            slong var_idx = poly->nvars + p;
                            if (*active_var == -1) {
                                *active_var = var_idx;
                            } else if (*active_var != var_idx) {
                                return 0; // Multiple variables used
                            }
                        }
                    }
                }
            }
        }
    }
    
    return 1; // All entries use at most one variable
}

// Convert fq_mvpoly to fq_nmod_poly for a specific variable
static void mvpoly_to_fq_nmod_poly(fq_nmod_poly_t poly, const fq_mvpoly_t *mvpoly, 
                                   slong var_index, const fq_nmod_ctx_t ctx) {
    fq_nmod_poly_zero(poly, ctx);
    
    slong nvars = mvpoly->nvars;
    
    for (slong t = 0; t < mvpoly->nterms; t++) {
        slong degree = 0;
        int other_vars_zero = 1;
        
        // Check if this term has non-zero exponents in other variables
        if (mvpoly->terms[t].var_exp) {
            for (slong v = 0; v < nvars; v++) {
                if (v == var_index && var_index < nvars) {
                    degree = mvpoly->terms[t].var_exp[v];
                } else if (mvpoly->terms[t].var_exp[v] > 0) {
                    other_vars_zero = 0;
                    break;
                }
            }
        }
        
        if (mvpoly->terms[t].par_exp && other_vars_zero) {
            for (slong p = 0; p < mvpoly->npars; p++) {
                slong idx = nvars + p;
                if (idx == var_index) {
                    degree = mvpoly->terms[t].par_exp[p];
                } else if (mvpoly->terms[t].par_exp[p] > 0) {
                    other_vars_zero = 0;
                    break;
                }
            }
        }
        
        // Only include terms where all other variables have zero exponent
        if (other_vars_zero) {
            fq_nmod_t existing;
            fq_nmod_init(existing, ctx);
            fq_nmod_poly_get_coeff(existing, poly, degree, ctx);
            fq_nmod_add(existing, existing, mvpoly->terms[t].coeff, ctx);
            fq_nmod_poly_set_coeff(poly, degree, existing, ctx);
            fq_nmod_clear(existing, ctx);
        }
    }
}

// Convert fq_nmod_poly back to fq_mvpoly
static void fq_nmod_poly_to_mvpoly(fq_mvpoly_t *mvpoly, const fq_nmod_poly_t poly,
                                   slong var_index, slong nvars, slong npars,
                                   const fq_nmod_ctx_t ctx) {
    fq_mvpoly_init(mvpoly, nvars, npars, ctx);
    
    slong degree = fq_nmod_poly_degree(poly, ctx);
    if (degree < 0) return;
    
    for (slong d = 0; d <= degree; d++) {
        fq_nmod_t coeff;
        fq_nmod_init(coeff, ctx);
        fq_nmod_poly_get_coeff(coeff, poly, d, ctx);
        
        if (!fq_nmod_is_zero(coeff, ctx)) {
            slong *var_exp = NULL;
            slong *par_exp = NULL;
            
            if (var_index < nvars && nvars > 0) {
                var_exp = (slong*) flint_calloc(nvars, sizeof(slong));
                var_exp[var_index] = d;
            } else if (var_index >= nvars && npars > 0) {
                par_exp = (slong*) flint_calloc(npars, sizeof(slong));
                par_exp[var_index - nvars] = d;
            }
            
            fq_mvpoly_add_term(mvpoly, var_exp, par_exp, coeff);
            
            if (var_exp) flint_free(var_exp);
            if (par_exp) flint_free(par_exp);
        }
        
        fq_nmod_clear(coeff, ctx);
    }
}

// ============= Polynomial Recursive Determinant =============

// Recursive determinant computation using fq_nmod_poly operations
static void compute_det_poly_recursive_helper(fq_nmod_poly_t det, 
                                             fq_nmod_poly_t **matrix,
                                             slong size, 
                                             const fq_nmod_ctx_t ctx) {
    if (size == 0) {
        fq_nmod_poly_one(det, ctx);
        return;
    }
    
    if (size == 1) {
        fq_nmod_poly_set(det, matrix[0][0], ctx);
        return;
    }
    
    if (size == 2) {
        fq_nmod_poly_t ad, bc;
        fq_nmod_poly_init(ad, ctx);
        fq_nmod_poly_init(bc, ctx);
        
        fq_nmod_poly_mul(ad, matrix[0][0], matrix[1][1], ctx);
        fq_nmod_poly_mul(bc, matrix[0][1], matrix[1][0], ctx);
        fq_nmod_poly_sub(det, ad, bc, ctx);
        
        fq_nmod_poly_clear(ad, ctx);
        fq_nmod_poly_clear(bc, ctx);
        return;
    }
    
    // General case: Laplace expansion
    fq_nmod_poly_zero(det, ctx);
    
    fq_nmod_poly_t cofactor, subdet;
    fq_nmod_poly_init(cofactor, ctx);
    fq_nmod_poly_init(subdet, ctx);
    
    // Allocate submatrix
    fq_nmod_poly_t **submatrix = (fq_nmod_poly_t**) malloc((size-1) * sizeof(fq_nmod_poly_t*));
    for (slong i = 0; i < size-1; i++) {
        submatrix[i] = (fq_nmod_poly_t*) malloc((size-1) * sizeof(fq_nmod_poly_t));
        for (slong j = 0; j < size-1; j++) {
            fq_nmod_poly_init(submatrix[i][j], ctx);
        }
    }
    
    for (slong col = 0; col < size; col++) {
        if (fq_nmod_poly_is_zero(matrix[0][col], ctx)) {
            continue;
        }
        
        // Build submatrix
        for (slong i = 1; i < size; i++) {
            slong sub_j = 0;
            for (slong j = 0; j < size; j++) {
                if (j != col) {
                    fq_nmod_poly_set(submatrix[i-1][sub_j], matrix[i][j], ctx);
                    sub_j++;
                }
            }
        }
        
        // Recursive call
        compute_det_poly_recursive_helper(subdet, submatrix, size-1, ctx);
        
        // Multiply by matrix element
        fq_nmod_poly_mul(cofactor, matrix[0][col], subdet, ctx);
        
        // Add or subtract based on sign
        if (col % 2 == 0) {
            fq_nmod_poly_add(det, det, cofactor, ctx);
        } else {
            fq_nmod_poly_sub(det, det, cofactor, ctx);
        }
    }
    
    // Cleanup
    for (slong i = 0; i < size-1; i++) {
        for (slong j = 0; j < size-1; j++) {
            fq_nmod_poly_clear(submatrix[i][j], ctx);
        }
        free(submatrix[i]);
    }
    free(submatrix);
    
    fq_nmod_poly_clear(cofactor, ctx);
    fq_nmod_poly_clear(subdet, ctx);
}

// Main function for polynomial recursive algorithm (FIXED)
void compute_fq_det_poly_recursive(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size) {
    if (size <= 0) {
        fq_mvpoly_init(result, matrix[0][0].nvars, matrix[0][0].npars, matrix[0][0].ctx);
        return;
    }
    
    timing_info_t total_start = start_timing();
    
    const fq_nmod_ctx_struct *ctx = matrix[0][0].ctx;
    slong nvars = matrix[0][0].nvars;
    slong npars = matrix[0][0].npars;
    slong total_vars = nvars + npars;
    
    DET_PRINT("Computing %ldx%ld determinant via polynomial recursive with Kronecker\n", size, size);
    DET_PRINT("Variables: %ld, Parameters: %ld\n", nvars, npars);
    
    // Special case: if already univariate, no Kronecker needed
    if (total_vars == 1) {
        DET_PRINT("Already univariate, using direct polynomial recursive\n");
        
        // Convert to fq_nmod_poly format
        fq_nmod_poly_t **poly_matrix = (fq_nmod_poly_t**) malloc(size * sizeof(fq_nmod_poly_t*));
        for (slong i = 0; i < size; i++) {
            poly_matrix[i] = (fq_nmod_poly_t*) malloc(size * sizeof(fq_nmod_poly_t));
            for (slong j = 0; j < size; j++) {
                fq_nmod_poly_init(poly_matrix[i][j], ctx);
                
                // Convert mvpoly to poly (direct conversion for univariate)
                for (slong t = 0; t < matrix[i][j].nterms; t++) {
                    slong degree = 0;
                    if (matrix[i][j].terms[t].var_exp && nvars > 0) {
                        degree = matrix[i][j].terms[t].var_exp[0];
                    } else if (matrix[i][j].terms[t].par_exp && npars > 0) {
                        degree = matrix[i][j].terms[t].par_exp[0];
                    }
                    fq_nmod_poly_set_coeff(poly_matrix[i][j], degree, 
                                          matrix[i][j].terms[t].coeff, ctx);
                }
            }
        }
        
        // Compute determinant
        fq_nmod_poly_t det_poly;
        fq_nmod_poly_init(det_poly, ctx);
        compute_det_poly_recursive_helper(det_poly, poly_matrix, size, ctx);
        
        // Convert back to mvpoly
        fq_mvpoly_init(result, nvars, npars, ctx);
        slong degree = fq_nmod_poly_degree(det_poly, ctx);
        for (slong d = 0; d <= degree; d++) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ctx);
            fq_nmod_poly_get_coeff(coeff, det_poly, d, ctx);
            
            if (!fq_nmod_is_zero(coeff, ctx)) {
                if (nvars > 0) {
                    slong *var_exp = (slong*) flint_calloc(nvars, sizeof(slong));
                    var_exp[0] = d;
                    fq_mvpoly_add_term(result, var_exp, NULL, coeff);
                    flint_free(var_exp);
                } else {
                    slong *par_exp = (slong*) flint_calloc(npars, sizeof(slong));
                    par_exp[0] = d;
                    fq_mvpoly_add_term(result, NULL, par_exp, coeff);
                    flint_free(par_exp);
                }
            }
            
            fq_nmod_clear(coeff, ctx);
        }
        
        // Cleanup
        for (slong i = 0; i < size; i++) {
            for (slong j = 0; j < size; j++) {
                fq_nmod_poly_clear(poly_matrix[i][j], ctx);
            }
            free(poly_matrix[i]);
        }
        free(poly_matrix);
        fq_nmod_poly_clear(det_poly, ctx);
        
        timing_info_t total_elapsed = end_timing(total_start);
        print_timing("Total poly recursive (univariate)", total_elapsed);
        return;
    }
    
    // Multivariate case: use Kronecker substitution
    
    // Step 1: Compute variable bounds (same as in Kronecker algorithm)
    timing_info_t bounds_start = start_timing();
    slong *var_bounds = (slong*) malloc(total_vars * sizeof(slong));
    compute_kronecker_bounds(var_bounds, matrix, size, nvars, npars);
    timing_info_t bounds_elapsed = end_timing(bounds_start);
    
    // Step 2: Compute substitution powers
    slong *substitution_powers = (slong*) malloc(total_vars * sizeof(slong));
    substitution_powers[0] = 1;
    for (slong v = 1; v < total_vars; v++) {
        substitution_powers[v] = substitution_powers[v-1] * var_bounds[v-1];
    }
    
    DET_PRINT("Substitution powers: ");
    for (slong v = 0; v < total_vars; v++) {
        DET_PRINT("%ld ", substitution_powers[v]);
    }
    DET_PRINT("\n");
    
    // Step 3: Convert matrix to univariate using Kronecker
    timing_info_t convert_start = start_timing();
    fq_nmod_poly_t **poly_matrix = (fq_nmod_poly_t**) malloc(size * sizeof(fq_nmod_poly_t*));
    for (slong i = 0; i < size; i++) {
        poly_matrix[i] = (fq_nmod_poly_t*) malloc(size * sizeof(fq_nmod_poly_t));
        for (slong j = 0; j < size; j++) {
            fq_nmod_poly_init(poly_matrix[i][j], ctx);
            mvpoly_to_univariate_kronecker(poly_matrix[i][j], &matrix[i][j], 
                                            substitution_powers, ctx);
        }
    }
    timing_info_t convert_elapsed = end_timing(convert_start);
    
    // Step 4: Compute determinant using recursive algorithm
    timing_info_t det_start = start_timing();
    fq_nmod_poly_t det_poly;
    fq_nmod_poly_init(det_poly, ctx);
    
    compute_det_poly_recursive_helper(det_poly, poly_matrix, size, ctx);
    
    timing_info_t det_elapsed = end_timing(det_start);
    DET_PRINT("Univariate determinant degree: %ld\n", fq_nmod_poly_degree(det_poly, ctx));
    
    // Step 5: Convert back to multivariate
    timing_info_t back_start = start_timing();
    univariate_to_mvpoly_kronecker(result, det_poly, substitution_powers, 
                                  var_bounds, nvars, npars, ctx);
    timing_info_t back_elapsed = end_timing(back_start);
    
    // Cleanup
    free(var_bounds);
    free(substitution_powers);
    
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            fq_nmod_poly_clear(poly_matrix[i][j], ctx);
        }
        free(poly_matrix[i]);
    }
    free(poly_matrix);
    fq_nmod_poly_clear(det_poly, ctx);
    
    timing_info_t total_elapsed = end_timing(total_start);
    
    printf("\n=== Polynomial Recursive Time Statistics ===\n");
    print_timing("Compute bounds", bounds_elapsed);
    print_timing("Convert to univariate", convert_elapsed);
    print_timing("Recursive determinant", det_elapsed);
    print_timing("Convert back", back_elapsed);
    print_timing("Total poly recursive", total_elapsed);
    printf("Final result: %ld terms\n", result->nterms);
    printf("============================================\n");
}

// ============= Kronecker Substitution Implementation =============

// Compute the Kronecker bound for a multivariate polynomial matrix
static void compute_kronecker_bounds(slong *var_bounds, fq_mvpoly_t **matrix, 
                                    slong size, slong nvars, slong npars) {
    slong total_vars = nvars + npars;
    
    // Initialize bounds
    for (slong v = 0; v < total_vars; v++) {
        var_bounds[v] = 0;
    }
    
    // Find maximum degree in each variable across all matrix entries
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            fq_mvpoly_t *poly = &matrix[i][j];
            
            for (slong t = 0; t < poly->nterms; t++) {
                // Check variable degrees
                if (poly->terms[t].var_exp && nvars > 0) {
                    for (slong v = 0; v < nvars && v < poly->nvars; v++) {
                        if (poly->terms[t].var_exp[v] > var_bounds[v]) {
                            var_bounds[v] = poly->terms[t].var_exp[v];
                        }
                    }
                }
                
                // Check parameter degrees
                if (poly->terms[t].par_exp && npars > 0) {
                    for (slong p = 0; p < npars && p < poly->npars; p++) {
                        if (poly->terms[t].par_exp[p] > var_bounds[nvars + p]) {
                            var_bounds[nvars + p] = poly->terms[t].par_exp[p];
                        }
                    }
                }
            }
        }
    }
    
    // Compute degree bound for determinant (sum of row maximums)
    for (slong v = 0; v < total_vars; v++) {
        slong det_bound = 0;
        
        for (slong row = 0; row < size; row++) {
            slong row_max = 0;
            
            for (slong col = 0; col < size; col++) {
                fq_mvpoly_t *poly = &matrix[row][col];
                
                for (slong t = 0; t < poly->nterms; t++) {
                    slong deg = 0;
                    
                    if (v < nvars && poly->terms[t].var_exp && v < poly->nvars) {
                        deg = poly->terms[t].var_exp[v];
                    } else if (v >= nvars && poly->terms[t].par_exp && 
                              v - nvars < poly->npars) {
                        deg = poly->terms[t].par_exp[v - nvars];
                    }
                    
                    if (deg > row_max) row_max = deg;
                }
            }
            det_bound += row_max;
        }
        
        var_bounds[v] = det_bound + 1;  // Add 1 for safety
    }
    for (slong v = 0; v < nvars/2; v++) {
        var_bounds[v] = var_bounds[v + nvars/2];
    }
}

// Convert multivariate polynomial to univariate using Kronecker substitution
static void mvpoly_to_univariate_kronecker(fq_nmod_poly_t uni_poly,
                                          const fq_mvpoly_t *mv_poly,
                                          const slong *substitution_powers,
                                          const fq_nmod_ctx_t ctx) {
    fq_nmod_poly_zero(uni_poly, ctx);
    
    if (mv_poly->nterms == 0) return;
    
    slong total_vars = mv_poly->nvars + mv_poly->npars;
    
    for (slong t = 0; t < mv_poly->nterms; t++) {
        slong uni_exp = 0;
        
        // Compute univariate exponent: sum of var_exp[i] * substitution_powers[i]
        if (mv_poly->terms[t].var_exp) {
            for (slong v = 0; v < mv_poly->nvars; v++) {
                uni_exp += mv_poly->terms[t].var_exp[v] * substitution_powers[v];
            }
        }
        
        if (mv_poly->terms[t].par_exp) {
            for (slong p = 0; p < mv_poly->npars; p++) {
                uni_exp += mv_poly->terms[t].par_exp[p] * 
                          substitution_powers[mv_poly->nvars + p];
            }
        }
        
        // Add coefficient at computed degree
        fq_nmod_t existing;
        fq_nmod_init(existing, ctx);
        fq_nmod_poly_get_coeff(existing, uni_poly, uni_exp, ctx);
        fq_nmod_add(existing, existing, mv_poly->terms[t].coeff, ctx);
        fq_nmod_poly_set_coeff(uni_poly, uni_exp, existing, ctx);
        fq_nmod_clear(existing, ctx);
    }
}

// Convert univariate polynomial back to multivariate
static void univariate_to_mvpoly_kronecker(fq_mvpoly_t *mv_poly,
                                          const fq_nmod_poly_t uni_poly,
                                          const slong *substitution_powers,
                                          const slong *var_bounds,
                                          slong nvars, slong npars,
                                          const fq_nmod_ctx_t ctx) {
    fq_mvpoly_init(mv_poly, nvars, npars, ctx);
    
    slong degree = fq_nmod_poly_degree(uni_poly, ctx);
    if (degree < 0) return;
    
    slong total_vars = nvars + npars;
    
    for (slong d = 0; d <= degree; d++) {
        fq_nmod_t coeff;
        fq_nmod_init(coeff, ctx);
        fq_nmod_poly_get_coeff(coeff, uni_poly, d, ctx);
        
        if (!fq_nmod_is_zero(coeff, ctx)) {
            // Decompose d into multivariate exponents
            slong *var_exp = NULL;
            slong *par_exp = NULL;
            
            if (nvars > 0) {
                var_exp = (slong*) flint_calloc(nvars, sizeof(slong));
            }
            if (npars > 0) {
                par_exp = (slong*) flint_calloc(npars, sizeof(slong));
            }
            
            slong remaining = d;
            
            // Extract exponents in reverse order (largest substitution power first)
            for (slong v = total_vars - 1; v >= 0; v--) {
                slong exp = remaining / substitution_powers[v];
                remaining = remaining % substitution_powers[v];
                
                if (v < nvars && var_exp) {
                    var_exp[v] = exp;
                } else if (v >= nvars && par_exp) {
                    par_exp[v - nvars] = exp;
                }
            }
            
            fq_mvpoly_add_term(mv_poly, var_exp, par_exp, coeff);
            
            if (var_exp) flint_free(var_exp);
            if (par_exp) flint_free(par_exp);
        }
        
        fq_nmod_clear(coeff, ctx);
    }
}

// Compute determinant using Kronecker substitution
void compute_fq_det_kronecker(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size) {
    if (size <= 0) {
        fq_mvpoly_init(result, matrix[0][0].nvars, matrix[0][0].npars, matrix[0][0].ctx);
        return;
    }
    
    timing_info_t total_start = start_timing();
    
    const fq_nmod_ctx_struct *ctx = matrix[0][0].ctx;
    slong nvars = matrix[0][0].nvars;
    slong npars = matrix[0][0].npars;
    slong total_vars = nvars + npars;
    
    DET_PRINT("Computing %ldx%ld determinant via Kronecker substitution\n", size, size);
    DET_PRINT("Variables: %ld, Parameters: %ld\n", nvars, npars);
    
    // Step 1: Compute variable bounds
    timing_info_t bounds_start = start_timing();
    slong *var_bounds = (slong*) malloc(total_vars * sizeof(slong));
    compute_kronecker_bounds(var_bounds, matrix, size, nvars, npars);
    timing_info_t bounds_elapsed = end_timing(bounds_start);
    
    DET_PRINT("Variable bounds: ");
    for (slong v = 0; v < total_vars; v++) {
        DET_PRINT("%ld ", var_bounds[v]);
    }
    DET_PRINT("\n");
    
    // Step 2: Compute substitution powers
    slong *substitution_powers = (slong*) malloc(total_vars * sizeof(slong));
    substitution_powers[0] = 1;
    for (slong v = 1; v < total_vars; v++) {
        substitution_powers[v] = substitution_powers[v-1] * var_bounds[v-1];
    }
    
    DET_PRINT("Substitution powers: ");
    for (slong v = 0; v < total_vars; v++) {
        DET_PRINT("%ld ", substitution_powers[v]);
    }
    DET_PRINT("\n");
    
    // Step 3: Convert matrix to univariate
    timing_info_t convert_start = start_timing();
    fq_nmod_poly_mat_t uni_mat;
    fq_nmod_poly_mat_init(uni_mat, size, size, ctx);
    
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            mvpoly_to_univariate_kronecker(fq_nmod_poly_mat_entry(uni_mat, i, j),
                                          &matrix[i][j], substitution_powers, ctx);
        }
    }
    timing_info_t convert_elapsed = end_timing(convert_start);
    
    // Check maximum degree
    slong max_degree = 0;
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            slong deg = fq_nmod_poly_degree(fq_nmod_poly_mat_entry(uni_mat, i, j), ctx);
            if (deg > max_degree) max_degree = deg;
        }
    }
    DET_PRINT("Maximum univariate degree after conversion: %ld\n", max_degree);
    
    // Step 4: Compute univariate determinant
    timing_info_t det_start = start_timing();
    fq_nmod_poly_t det_poly;
    fq_nmod_poly_init(det_poly, ctx);
    
    // Use the optimized univariate determinant function
    fq_nmod_poly_mat_det_iter(det_poly, uni_mat, ctx);
    
    timing_info_t det_elapsed = end_timing(det_start);
    DET_PRINT("Univariate determinant degree: %ld\n", fq_nmod_poly_degree(det_poly, ctx));
    
    // Step 5: Convert back to multivariate
    timing_info_t back_start = start_timing();
    univariate_to_mvpoly_kronecker(result, det_poly, substitution_powers, 
                                  var_bounds, nvars, npars, ctx);
    timing_info_t back_elapsed = end_timing(back_start);
    
    // Cleanup
    free(var_bounds);
    free(substitution_powers);
    fq_nmod_poly_clear(det_poly, ctx);
    fq_nmod_poly_mat_clear(uni_mat, ctx);
    
    timing_info_t total_elapsed = end_timing(total_start);
    
    printf("\n=== Kronecker Substitution Time Statistics ===\n");
    print_timing("Compute bounds", bounds_elapsed);
    print_timing("Convert to univariate", convert_elapsed);
    print_timing("Univariate determinant", det_elapsed);
    print_timing("Convert back", back_elapsed);
    print_timing("Total Kronecker", total_elapsed);
    printf("Final result: %ld terms\n", result->nterms);
    printf("==============================================\n");
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

void nmod_mpoly_to_fq_mvpoly_old(fq_mvpoly_t *poly, const nmod_mpoly_t mpoly,
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

void nmod_mpoly_to_fq_mvpoly_old2(fq_mvpoly_t *result, const nmod_mpoly_t poly,
                             slong nvars, slong npars,
                             const nmod_mpoly_ctx_t mpoly_ctx,
                             const fq_nmod_ctx_t field_ctx) {
    fq_mvpoly_init(result, nvars, npars, field_ctx);
    
    slong nterms = nmod_mpoly_length(poly, mpoly_ctx);
    slong total_vars = nmod_mpoly_ctx_nvars(mpoly_ctx);
    
    for (slong i = 0; i < nterms; i++) {
        /* Get coefficient */
        mp_limb_t coeff_limb = nmod_mpoly_get_term_coeff_ui(poly, i, mpoly_ctx);
        fq_nmod_t coeff;
        fq_nmod_init(coeff, field_ctx);
        fq_nmod_set_ui(coeff, coeff_limb, field_ctx);
        
        /* Get exponents */
        ulong *exp = (ulong*) flint_malloc(total_vars * sizeof(ulong));
        nmod_mpoly_get_term_exp_ui(exp, poly, i, mpoly_ctx);
        
        /* Split exponents into variables and parameters */
        slong *var_exp = NULL;
        slong *par_exp = NULL;
        
        if (nvars > 0) {
            var_exp = (slong*) flint_calloc(nvars, sizeof(slong));
            for (slong j = 0; j < nvars && j < total_vars; j++) {
                var_exp[j] = exp[j];
            }
        }
        
        if (npars > 0 && total_vars > nvars) {
            par_exp = (slong*) flint_calloc(npars, sizeof(slong));
            for (slong j = 0; j < npars && (nvars + j) < total_vars; j++) {
                par_exp[j] = exp[nvars + j];
            }
        }
        
        /* Add term to result */
        fq_mvpoly_add_term_fast(result, var_exp, par_exp, coeff);
        
        /* Cleanup */
        fq_nmod_clear(coeff, field_ctx);
        flint_free(exp);
        if (var_exp) flint_free(var_exp);
        if (par_exp) flint_free(par_exp);
    }
}

void nmod_mpoly_to_fq_mvpoly(fq_mvpoly_t *result, const nmod_mpoly_t poly,
                                       slong nvars, slong npars,
                                       const nmod_mpoly_ctx_t mpoly_ctx,
                                       const fq_nmod_ctx_t field_ctx) {
    fq_mvpoly_init(result, nvars, npars, field_ctx);
    
    slong nterms = nmod_mpoly_length(poly, mpoly_ctx);
    if (nterms == 0) return;
    
    slong total_vars = nmod_mpoly_ctx_nvars(mpoly_ctx);
    
    if (result->alloc < nterms) {
        result->alloc = nterms;
        result->terms = (fq_monomial_t*) flint_realloc(result->terms, 
                                                        result->alloc * sizeof(fq_monomial_t));
    }
    
    // 批量分配指数缓冲区
    ulong *exp_buffer = (ulong*) flint_malloc(total_vars * sizeof(ulong));
    
    // 批量转换，但使用正确的初始化
    for (slong i = 0; i < nterms; i++) {
        // 获取系数
        mp_limb_t coeff_limb = nmod_mpoly_get_term_coeff_ui(poly, i, mpoly_ctx);
        
        // 初始化系数（重要！）
        fq_nmod_init(result->terms[i].coeff, field_ctx);
        fq_nmod_set_ui(result->terms[i].coeff, coeff_limb, field_ctx);
        
        // 获取指数
        nmod_mpoly_get_term_exp_ui(exp_buffer, poly, i, mpoly_ctx);
        
        // 分配并设置变量指数
        if (nvars > 0) {
            result->terms[i].var_exp = (slong*) flint_calloc(nvars, sizeof(slong));
            for (slong j = 0; j < nvars && j < total_vars; j++) {
                result->terms[i].var_exp[j] = (slong)exp_buffer[j];
            }
        } else {
            result->terms[i].var_exp = NULL;
        }
        
        // 分配并设置参数指数
        if (npars > 0 && total_vars > nvars) {
            result->terms[i].par_exp = (slong*) flint_calloc(npars, sizeof(slong));
            for (slong j = 0; j < npars && (nvars + j) < total_vars; j++) {
                result->terms[i].par_exp[j] = (slong)exp_buffer[nvars + j];
            }
        } else {
            result->terms[i].par_exp = NULL;
        }
    }
    
    // 设置项数
    result->nterms = nterms;
    
    // 清理
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

// Parallel determinant computation for nmod_mpoly with proper nested parallelism
void compute_nmod_mpoly_det_parallel_optimized(nmod_mpoly_t det_result, 
                                              nmod_mpoly_t **mpoly_matrix, 
                                              slong size, 
                                              nmod_mpoly_ctx_t mpoly_ctx,
                                              slong depth) {
    // For deep recursion or small matrices, use sequential
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
    
    // Determine parallelism strategy based on depth
    if (depth == 0) {
        // First level: use parallel for with nested parallelism enabled
        #pragma omp parallel for schedule(static) num_threads(FLINT_MIN(nonzero_count, omp_get_max_threads()))
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
            
            // Recursive call - this will use nested parallelism at depth 1
            compute_nmod_mpoly_det_parallel_optimized(subdet, submatrix, size-1, mpoly_ctx, depth+1);
            
            // Compute cofactor
            nmod_mpoly_mul(cofactor, mpoly_matrix[0][col], subdet, mpoly_ctx);
            
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
    } else if (depth == 1 && size >= PARALLEL_THRESHOLD) {
        // Second level: also use parallel for, but with fewer threads
        slong max_threads_level2 = FLINT_MAX(1, omp_get_max_threads() / size);
        
        #pragma omp parallel for schedule(static) num_threads(FLINT_MIN(nonzero_count, max_threads_level2)) if(nonzero_count >= 3)
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
            
            // At depth > 1, use sequential computation
            compute_nmod_mpoly_det_recursive(subdet, submatrix, size-1, mpoly_ctx);
            
            // Compute cofactor
            nmod_mpoly_mul(cofactor, mpoly_matrix[0][col], subdet, mpoly_ctx);
            
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
    } else {
        // Sequential fallback for deeper levels
        for (slong col = 0; col < size; col++) {
            if (nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
                continue;
            }
            
            nmod_mpoly_t cofactor, subdet;
            nmod_mpoly_init(cofactor, mpoly_ctx);
            nmod_mpoly_init(subdet, mpoly_ctx);
            
            // Create and fill submatrix
            nmod_mpoly_t **submatrix = (nmod_mpoly_t**) flint_malloc((size-1) * sizeof(nmod_mpoly_t*));
            for (slong i = 0; i < size-1; i++) {
                submatrix[i] = (nmod_mpoly_t*) flint_malloc((size-1) * sizeof(nmod_mpoly_t));
                for (slong j = 0; j < size-1; j++) {
                    nmod_mpoly_init(submatrix[i][j], mpoly_ctx);
                }
            }
            
            for (slong i = 1; i < size; i++) {
                slong sub_j = 0;
                for (slong j = 0; j < size; j++) {
                    if (j != col) {
                        nmod_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j], mpoly_ctx);
                        sub_j++;
                    }
                }
            }
            
            // Sequential computation
            compute_nmod_mpoly_det_recursive(subdet, submatrix, size-1, mpoly_ctx);
            
            // Compute cofactor and store
            nmod_mpoly_mul(cofactor, mpoly_matrix[0][col], subdet, mpoly_ctx);
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
    }
    
    // Sum results (sequential to avoid race conditions)
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

void fq_nmod_mpoly_to_fq_mvpoly_old(fq_mvpoly_t *poly, const fq_nmod_mpoly_t mpoly,
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

void fq_nmod_mpoly_to_fq_mvpoly(fq_mvpoly_t *poly, const fq_nmod_mpoly_t mpoly,
                                slong nvars, slong npars, 
                                fq_nmod_mpoly_ctx_t mpoly_ctx, const fq_nmod_ctx_t ctx) {
    fq_mvpoly_init(poly, nvars, npars, ctx);
    
    slong nterms = fq_nmod_mpoly_length(mpoly, mpoly_ctx);
    if (nterms == 0) return;
    
    slong total_vars = fq_nmod_mpoly_ctx_nvars(mpoly_ctx);
    
    // Pre-allocate the terms array
    if (poly->alloc < nterms) {
        poly->alloc = nterms;
        poly->terms = (fq_monomial_t*) flint_realloc(poly->terms, 
                                                      poly->alloc * sizeof(fq_monomial_t));
    }
    
    // Allocate a single exponent buffer for reading
    ulong *exp_buffer = (ulong*) flint_malloc(total_vars * sizeof(ulong));
    
    // Process all terms - but allocate individually for compatibility
    for (slong i = 0; i < nterms; i++) {
        // Initialize coefficient
        fq_nmod_init(poly->terms[i].coeff, ctx);
        fq_nmod_mpoly_get_term_coeff_fq_nmod(poly->terms[i].coeff, mpoly, i, mpoly_ctx);
        
        // Get exponents for this term
        fq_nmod_mpoly_get_term_exp_ui(exp_buffer, mpoly, i, mpoly_ctx);
        
        // Allocate and set variable exponents
        if (nvars > 0) {
            poly->terms[i].var_exp = (slong*) flint_calloc(nvars, sizeof(slong));
            for (slong j = 0; j < nvars && j < total_vars; j++) {
                poly->terms[i].var_exp[j] = (slong)exp_buffer[j];
            }
        } else {
            poly->terms[i].var_exp = NULL;
        }
        
        // Allocate and set parameter exponents
        if (npars > 0 && total_vars > nvars) {
            poly->terms[i].par_exp = (slong*) flint_calloc(npars, sizeof(slong));
            for (slong j = 0; j < npars && (nvars + j) < total_vars; j++) {
                poly->terms[i].par_exp[j] = (slong)exp_buffer[nvars + j];
            }
        } else {
            poly->terms[i].par_exp = NULL;
        }
    }
    
    // Set the number of terms
    poly->nterms = nterms;
    
    // Cleanup
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
    
    // Determine parallelism strategy based on depth
    if (depth == 0) {
        // First level: use parallel for with nested parallelism enabled
        #pragma omp parallel for schedule(static) num_threads(FLINT_MIN(nonzero_count, omp_get_max_threads()))
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
            
            // Recursive call - this will use nested parallelism at depth 1
            compute_fq_nmod_mpoly_det_parallel_optimized(subdet, submatrix, size-1, mpoly_ctx, depth+1);
            
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
    } else if (depth == 1 && size >= PARALLEL_THRESHOLD) {
        // Second level: also use parallel for, but with fewer threads
        slong max_threads_level2 = FLINT_MAX(1, omp_get_max_threads() / size);
        
        #pragma omp parallel for schedule(static) num_threads(FLINT_MIN(nonzero_count, max_threads_level2)) if(nonzero_count >= 3)
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
            
            // At depth > 1, use sequential computation
            compute_fq_nmod_mpoly_det_recursive(subdet, submatrix, size-1, mpoly_ctx);
            
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
    } else {
        // Sequential fallback for deeper levels or small matrices
        for (slong col = 0; col < size; col++) {
            if (fq_nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
                continue;
            }
            
            fq_nmod_mpoly_t cofactor, subdet;
            fq_nmod_mpoly_init(cofactor, mpoly_ctx);
            fq_nmod_mpoly_init(subdet, mpoly_ctx);
            
            // Create and fill submatrix
            fq_nmod_mpoly_t **submatrix = (fq_nmod_mpoly_t**) flint_malloc((size-1) * sizeof(fq_nmod_mpoly_t*));
            for (slong i = 0; i < size-1; i++) {
                submatrix[i] = (fq_nmod_mpoly_t*) flint_malloc((size-1) * sizeof(fq_nmod_mpoly_t));
                for (slong j = 0; j < size-1; j++) {
                    fq_nmod_mpoly_init(submatrix[i][j], mpoly_ctx);
                }
            }
            
            for (slong i = 1; i < size; i++) {
                slong sub_j = 0;
                for (slong j = 0; j < size; j++) {
                    if (j != col) {
                        fq_nmod_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j], mpoly_ctx);
                        sub_j++;
                    }
                }
            }
            
            // Sequential computation
            compute_fq_nmod_mpoly_det_recursive(subdet, submatrix, size-1, mpoly_ctx);
            
            // Compute cofactor and store
            poly_mul_dense_optimized(cofactor, mpoly_matrix[0][col], subdet, mpoly_ctx);
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
    }
    
    // Sum results (sequential to avoid race conditions)
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


// ============= Main Interface with Algorithm Selection =============

void compute_fq_det_recursive_flint(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size) {
    if (size <= 0) {
        fq_mvpoly_init(result, matrix[0][0].nvars, matrix[0][0].npars, matrix[0][0].ctx);
        return;
    }
    
    timing_info_t total_start = start_timing();
    
    slong nvars = matrix[0][0].nvars;
    slong npars = matrix[0][0].npars;
    const fq_nmod_ctx_struct *ctx = matrix[0][0].ctx;
    
    // Choose algorithm based on configuration
    #if DET_ALGORITHM == DET_ALGORITHM_INTERPOLATION
    {
        printf("Using multivariate interpolation algorithm\n");
        
        // Include the interpolation header if not already included
        #ifndef FQ_NMOD_INTERPOLATION_OPTIMIZED_H
        #include "fq_multivariate_interpolation.h"
        #endif
        slong total_vars = nvars + npars;
        slong *var_bounds = (slong*) malloc(total_vars * sizeof(slong));
        compute_kronecker_bounds(var_bounds, matrix, size, nvars, npars);
        // Use interpolation algorithm
        fq_compute_det_by_interpolation_optimized(result, matrix, size, 
                                                 nvars, npars, ctx, var_bounds);
        return;
    }
    #elif DET_ALGORITHM == DET_ALGORITHM_KRONECKER
    {
        printf("Using Kronecker substitution algorithm\n");
        compute_fq_det_kronecker(result, matrix, size);
        return;
    }
    #elif DET_ALGORITHM == DET_ALGORITHM_POLY_RECURSIVE
    {
        printf("Using polynomial recursive algorithm\n");
        compute_fq_det_poly_recursive(result, matrix, size);
        return;
    }
    #else
    {
        // Original recursive algorithm (default)
        slong max_threads = omp_get_max_threads();
        DET_PRINT("Computing %ldx%ld determinant (OpenMP: %ld threads available)\n", 
                  size, size, max_threads);
        
        slong total_vars = nvars + npars;
        
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
                static int nested_init = 0;
                if (!nested_init) {
                    init_nested_parallelism();
                    nested_init = 1;
                }  
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
            fq_mvpoly_clear(result);  // Clear the initialization from the beginning
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
            // Initialize nested parallelism if not already done
            static int nested_init = 0;
            if (!nested_init) {
                init_nested_parallelism();
                nested_init = 1;
            }  
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
        fq_mvpoly_clear(result);
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
    #endif
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