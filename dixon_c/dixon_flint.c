
/*
 * Complete Dixon Resultant Implementation for Finite Extension Fields
 * 
 * Fixed compilation issues:
 * - Context pointer handling
 * - Function declaration order
 * - Assignment to array types
 * 
 * Compile with:
 * gcc -o dixon_fq_nmod dixon_main.c -lflint -lmpfr -lgmp -lpthread
     --param max-inline-insns-single=1000 \
    --param max-inline-insns-auto=2000 \
    --param inline-unit-growth=50 \
    --param large-function-growth=200 \
    gcc -O3 -march=native --param max-inline-insns-single=1000 --param max-inline-insns-auto=2000 --param inline-unit-growth=50 --param large-function-growth=200 -Winline -fdump-ipa-inline -o dixon_flint dixon_flint.c -lflint -lmpfr -lgmp -lpthread -L/home/suohaohai02/mylinks -lflint -lstdc++ -lpml2
 //gcc -O3 -march=native -o dixon_flint dixon_flint.c -lflint -lmpfr -lgmp -lpthread -L/home/suohaohai02/mylinks -lflint -lstdc++ -lpml2
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include <unistd.h>
#include <ctype.h>

#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_poly.h>
#include <flint/fq_nmod_mat.h>
#include <flint/fq_nmod_mpoly.h>
#include <flint/nmod_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/profiler.h>

// Method enumeration for determinant computation
typedef enum {
    DET_METHOD_KRONECKER = 0,
    DET_METHOD_MBOT = 1,
    DET_METHOD_INTERPOLATION = 2
} det_method_t;

// Monomial structure for finite extension fields
typedef struct {
    slong *var_exp;      // exponent vector for variables (including duals)
    slong *par_exp;      // exponent vector for parameters
    fq_nmod_t coeff;     // coefficient in F_{p^d}
} fq_monomial_t;

// Multivariate polynomial with parameters over F_{p^d}
typedef struct {
    slong nvars;         // number of variables (not including duals)
    slong npars;         // number of parameters
    slong nterms;        // number of terms
    slong alloc;         // allocated space
    fq_monomial_t *terms; // array of terms
    const fq_nmod_ctx_struct *ctx;   // finite field context (stored as pointer to struct)
} fq_mvpoly_t;

// Comparison function for sorting degrees
typedef struct {
    slong index;
    slong degree;
} fq_index_degree_pair;


// Forward declarations
void fq_mvpoly_init(fq_mvpoly_t *p, slong nvars, slong npars, const fq_nmod_ctx_t ctx);
void fq_mvpoly_clear(fq_mvpoly_t *p);
void fq_mvpoly_copy(fq_mvpoly_t *dest, const fq_mvpoly_t *src);
void fq_mvpoly_add_term(fq_mvpoly_t *p, const slong *var_exp, const slong *par_exp, const fq_nmod_t coeff);
void fq_mvpoly_add_term_fast(fq_mvpoly_t *p, const slong *var_exp, const slong *par_exp, const fq_nmod_t coeff);
void fq_mvpoly_print(const fq_mvpoly_t *p, const char *name);
void fq_mvpoly_print_expanded(const fq_mvpoly_t *p, const char *name, int use_dual);
void fq_mvpoly_mul(fq_mvpoly_t *result, const fq_mvpoly_t *a, const fq_mvpoly_t *b);
void fq_mvpoly_pow(fq_mvpoly_t *result, const fq_mvpoly_t *base, slong power);
void fq_mvpoly_scalar_mul(fq_mvpoly_t *result, const fq_mvpoly_t *p, const fq_nmod_t scalar);
void fq_mvpoly_add(fq_mvpoly_t *result, const fq_mvpoly_t *a, const fq_mvpoly_t *b);
void evaluate_fq_mvpoly_at_params(fq_nmod_t result, const fq_mvpoly_t *poly, const fq_nmod_t *param_vals);
void compute_fq_cancel_matrix_det(fq_mvpoly_t *result, fq_mvpoly_t **modified_M_mvpoly,
                                 slong nvars, slong npars, det_method_t method);
void divide_by_fq_linear_factor_improved(fq_mvpoly_t *quotient, const fq_mvpoly_t *dividend,
                                        slong var_idx, slong nvars, slong npars);
void fq_dixon_resultant(fq_mvpoly_t *result, fq_mvpoly_t *polys, 
                       slong nvars, slong npars);

#include "fq_unified_interface.h"
#include "fq_mpoly_mat_det.h"
#include "fq_multivariate_interpolation.h"
#include "dixon_interface_flint.h"



// ============ Basic fq_mvpoly operations ============

void fq_mvpoly_init(fq_mvpoly_t *p, slong nvars, slong npars, const fq_nmod_ctx_t ctx) {
    p->nvars = nvars;
    p->npars = npars;
    p->nterms = 0;
    p->alloc = 16;
    p->terms = (fq_monomial_t*) flint_malloc(p->alloc * sizeof(fq_monomial_t));
    p->ctx = ctx;  // Store pointer to the context struct
}

void fq_mvpoly_clear(fq_mvpoly_t *p) {
    if (p->nterms > 0 && p->terms) {
        // Clear all coefficient elements
        for (slong i = 0; i < p->nterms; i++) {
            fq_nmod_clear(p->terms[i].coeff, p->ctx);
            if (p->terms[i].var_exp) flint_free(p->terms[i].var_exp);
            if (p->terms[i].par_exp) flint_free(p->terms[i].par_exp);
        }
        flint_free(p->terms);
    }
    // Don't clear context as it's not owned by this polynomial
    p->terms = NULL;
    p->nterms = 0;
    p->alloc = 0;
}

void fq_nmod_mat_transpose(fq_nmod_mat_t dest, const fq_nmod_mat_t src, const fq_nmod_ctx_t ctx) {
    slong nrows = fq_nmod_mat_nrows(src, ctx);
    slong ncols = fq_nmod_mat_ncols(src, ctx);
    
    // Ensure dest has correct dimensions (ncols x nrows)
    if (fq_nmod_mat_nrows(dest, ctx) != ncols || fq_nmod_mat_ncols(dest, ctx) != nrows) {
        // If dimensions don't match, we can't safely transpose
        // This should be handled by proper initialization before calling
        printf("Error: Matrix dimensions don't match for transpose\n");
        return;
    }
    
    for (slong i = 0; i < nrows; i++) {
        for (slong j = 0; j < ncols; j++) {
            fq_nmod_set(fq_nmod_mat_entry(dest, j, i), 
                        fq_nmod_mat_entry(src, i, j), ctx);
        }
    }
}

static int compare_fq_degrees(const void *a, const void *b) {
    fq_index_degree_pair *pa = (fq_index_degree_pair *)a;
    fq_index_degree_pair *pb = (fq_index_degree_pair *)b;
    
    // Sort by degree (ascending)
    if (pa->degree < pb->degree) return -1;
    if (pa->degree > pb->degree) return 1;
    
    // Same degree, sort by index for stability
    if (pa->index < pb->index) return -1;
    if (pa->index > pb->index) return 1;
    return 0;
}


void fq_mvpoly_add_term_fast(fq_mvpoly_t *p, const slong *var_exp, const slong *par_exp, const fq_nmod_t coeff) {
    if (fq_nmod_is_zero(coeff, p->ctx)) return;
    
    // Add new term without checking for duplicates
    if (p->nterms >= p->alloc) {
        p->alloc = p->alloc ? p->alloc * 2 : 16;
        p->terms = (fq_monomial_t*) flint_realloc(p->terms, p->alloc * sizeof(fq_monomial_t));
    }
    
    // Initialize coefficient
    fq_nmod_init(p->terms[p->nterms].coeff, p->ctx);
    fq_nmod_set(p->terms[p->nterms].coeff, coeff, p->ctx);
    
    // Allocate and copy exponents
    if (p->nvars > 0) {
        p->terms[p->nterms].var_exp = (slong*) flint_calloc(p->nvars, sizeof(slong));
        if (var_exp) {
            memcpy(p->terms[p->nterms].var_exp, var_exp, p->nvars * sizeof(slong));
        }
    } else {
        p->terms[p->nterms].var_exp = NULL;
    }
    
    if (p->npars > 0) {
        p->terms[p->nterms].par_exp = (slong*) flint_calloc(p->npars, sizeof(slong));
        if (par_exp) {
            memcpy(p->terms[p->nterms].par_exp, par_exp, p->npars * sizeof(slong));
        }
    } else {
        p->terms[p->nterms].par_exp = NULL;
    }
    
    p->nterms++;
}

void fq_mvpoly_copy(fq_mvpoly_t *dest, const fq_mvpoly_t *src) {
    fq_mvpoly_init(dest, src->nvars, src->npars, src->ctx);
    
    for (slong i = 0; i < src->nterms; i++) {
        fq_mvpoly_add_term_fast(dest, src->terms[i].var_exp, src->terms[i].par_exp, src->terms[i].coeff);
    }
}

void fq_mvpoly_add_term(fq_mvpoly_t *p, const slong *var_exp, const slong *par_exp, const fq_nmod_t coeff) {
    if (fq_nmod_is_zero(coeff, p->ctx)) return;
    
    // Check if term already exists
    for (slong i = 0; i < p->nterms; i++) {
        int same = 1;
        
        // Check variable exponents
        if (p->nvars > 0) {
            for (slong j = 0; j < p->nvars; j++) {
                slong e1 = var_exp ? var_exp[j] : 0;
                slong e2 = p->terms[i].var_exp ? p->terms[i].var_exp[j] : 0;
                if (e1 != e2) {
                    same = 0;
                    break;
                }
            }
        }
        
        // Check parameter exponents
        if (same && p->npars > 0) {
            for (slong j = 0; j < p->npars; j++) {
                slong e1 = par_exp ? par_exp[j] : 0;
                slong e2 = p->terms[i].par_exp ? p->terms[i].par_exp[j] : 0;
                if (e1 != e2) {
                    same = 0;
                    break;
                }
            }
        }
        
        if (same) {
            fq_nmod_add(p->terms[i].coeff, p->terms[i].coeff, coeff, p->ctx);
            if (fq_nmod_is_zero(p->terms[i].coeff, p->ctx)) {
                // Remove zero term
                fq_nmod_clear(p->terms[i].coeff, p->ctx);
                if (p->terms[i].var_exp) flint_free(p->terms[i].var_exp);
                if (p->terms[i].par_exp) flint_free(p->terms[i].par_exp);
                for (slong j = i; j < p->nterms - 1; j++) {
                    p->terms[j] = p->terms[j + 1];
                }
                p->nterms--;
            }
            return;
        }
    }
    
    // Add new term
    fq_mvpoly_add_term_fast(p, var_exp, par_exp, coeff);
}

void fq_mvpoly_print(const fq_mvpoly_t *p, const char *name) {
    fq_mvpoly_print_expanded(p, name, 0);
}

void fq_mvpoly_print_expanded(const fq_mvpoly_t *p, const char *name, int use_dual) {
    printf("%s = ", name);
    if (p->nterms == 0) {
        printf("0\n");
        return;
    }
    
    char var_names[] = {'x', 'y', 'z', 'w', 'v', 'u'};
    char par_names[] = {'a', 'b', 'c', 'd'};
    
    for (slong i = 0; i < p->nterms; i++) {
        if (i > 0) printf(" + ");
        
        // Handle coefficient
        fq_nmod_t neg_one, one;
        fq_nmod_init(neg_one, p->ctx);
        fq_nmod_init(one, p->ctx);
        fq_nmod_one(one, p->ctx);
        fq_nmod_neg(neg_one, one, p->ctx);
        
        // Check if we have any variables or parameters for this term
        int has_vars = 0;
        if (p->terms[i].var_exp) {
            
            for (slong j = 0; j < p->nvars; j++) {
                if (p->terms[i].var_exp[j] > 0) {
                    has_vars = 1;
                    break;
                }
            }
        }
        if (!has_vars && p->terms[i].par_exp) {
            for (slong j = 0; j < p->npars; j++) {
                if (p->terms[i].par_exp[j] > 0) {
                    has_vars = 1;
                    break;
                }
            }
        }
        
        // Print coefficient correctly
        if (fq_nmod_equal(p->terms[i].coeff, neg_one, p->ctx)) {
            printf("-");
            if (!has_vars) printf("1");
        } else if (fq_nmod_is_one(p->terms[i].coeff, p->ctx)) {
            if (!has_vars) printf("1");
            // For coefficient 1 with variables, don't print anything
        } else {
            // For all other coefficients, show using fq_nmod_print_pretty
            fq_nmod_print_pretty(p->terms[i].coeff, p->ctx);
            // Only add * if we have variables/parameters to follow
            if (has_vars) printf("*");
        }
        
        fq_nmod_clear(neg_one, p->ctx);
        fq_nmod_clear(one, p->ctx);
        int var_printed = 0;
        // Print variables - simplified version
        if (use_dual) {
            // For expanded polynomials with dual variables
            slong actual_nvars = p->nvars / 2;
            
            // First print regular variables
            for (slong j = 0; j < actual_nvars; j++) {
                if (p->terms[i].var_exp && p->terms[i].var_exp[j] > 0) {
                    if (var_printed) printf("*");
                    if (j < 6) printf("%c", var_names[j]);
                    else printf("x_%ld", j);
                    
                    if (p->terms[i].var_exp[j] > 1) {
                        printf("^%ld", p->terms[i].var_exp[j]);
                    }
                    var_printed = 1;
                }
            }
            
            // Then print dual variables with tilde
            for (slong j = actual_nvars; j < p->nvars; j++) {
                if (p->terms[i].var_exp && p->terms[i].var_exp[j] > 0) {
                    if (var_printed) printf("*");
                    slong orig_idx = j - actual_nvars;
                    if (orig_idx < 6) printf("~%c", var_names[orig_idx]);
                    else printf("~x_%ld", orig_idx);
                    
                    if (p->terms[i].var_exp[j] > 1) {
                        printf("^%ld", p->terms[i].var_exp[j]);
                    }
                    var_printed = 1;
                }
            }
        } else {
            // Normal variables
            for (slong j = 0; j < p->nvars; j++) {
                if (p->terms[i].var_exp && p->terms[i].var_exp[j] > 0) {
                    if (var_printed) printf("*");
                    if (j < 6) printf("%c", var_names[j]);
                    else printf("x_%ld", j);
                    
                    if (p->terms[i].var_exp[j] > 1) {
                        printf("^%ld", p->terms[i].var_exp[j]);
                    }
                    var_printed = 1;
                }
            }
        }
        
        // Print parameters
        for (slong j = 0; j < p->npars; j++) {
            if (p->terms[i].par_exp && p->terms[i].par_exp[j] > 0) {
                if (var_printed) printf("*");
                if (j < 4) printf("%c", par_names[j]);
                else printf("p_%ld", j);
                
                if (p->terms[i].par_exp[j] > 1) {
                    printf("^%ld", p->terms[i].par_exp[j]);
                }
                var_printed = 1;
            }
        }
    }
    printf("\n");
}

// Multiply two multivariate polynomials
void fq_mvpoly_mul(fq_mvpoly_t *result, const fq_mvpoly_t *a, const fq_mvpoly_t *b) {
    fq_mvpoly_init(result, a->nvars, a->npars, a->ctx);
    
    for (slong i = 0; i < a->nterms; i++) {
        for (slong j = 0; j < b->nterms; j++) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, a->ctx);
            fq_nmod_mul(coeff, a->terms[i].coeff, b->terms[j].coeff, a->ctx);
            
            slong *var_exp = NULL;
            slong *par_exp = NULL;
            
            if (a->nvars > 0) {
                var_exp = (slong*) flint_calloc(a->nvars, sizeof(slong));
                for (slong k = 0; k < a->nvars; k++) {
                    var_exp[k] = (a->terms[i].var_exp ? a->terms[i].var_exp[k] : 0) +
                                 (b->terms[j].var_exp ? b->terms[j].var_exp[k] : 0);
                }
            }
            
            if (a->npars > 0) {
                par_exp = (slong*) flint_calloc(a->npars, sizeof(slong));
                for (slong k = 0; k < a->npars; k++) {
                    par_exp[k] = (a->terms[i].par_exp ? a->terms[i].par_exp[k] : 0) +
                                 (b->terms[j].par_exp ? b->terms[j].par_exp[k] : 0);
                }
            }
            
            fq_mvpoly_add_term(result, var_exp, par_exp, coeff);
            
            fq_nmod_clear(coeff, a->ctx);
            if (var_exp) flint_free(var_exp);
            if (par_exp) flint_free(par_exp);
        }
    }
}

// Compute polynomial to a power
void fq_mvpoly_pow(fq_mvpoly_t *result, const fq_mvpoly_t *base, slong power) {
    fq_mvpoly_init(result, base->nvars, base->npars, base->ctx);
    
    fq_nmod_t one;
    fq_nmod_init(one, base->ctx);
    fq_nmod_one(one, base->ctx);
    fq_mvpoly_add_term(result, NULL, NULL, one); // result = 1
    fq_nmod_clear(one, base->ctx);
    
    if (power == 0) return;
    
    fq_mvpoly_t base_copy, temp;
    fq_mvpoly_copy(&base_copy, base);
    
    while (power > 0) {
        if (power & 1) {
            fq_mvpoly_mul(&temp, result, &base_copy);
            fq_mvpoly_clear(result);
            fq_mvpoly_copy(result, &temp);
            fq_mvpoly_clear(&temp);
        }
        
        if (power > 1) {
            fq_mvpoly_mul(&temp, &base_copy, &base_copy);
            fq_mvpoly_clear(&base_copy);
            fq_mvpoly_copy(&base_copy, &temp);
            fq_mvpoly_clear(&temp);
        }
        
        power >>= 1;
    }
    
    fq_mvpoly_clear(&base_copy);
}

// Multiply polynomial by scalar
void fq_mvpoly_scalar_mul(fq_mvpoly_t *result, const fq_mvpoly_t *p, const fq_nmod_t scalar) {
    fq_mvpoly_init(result, p->nvars, p->npars, p->ctx);
    
    for (slong i = 0; i < p->nterms; i++) {
        fq_nmod_t new_coeff;
        fq_nmod_init(new_coeff, p->ctx);
        fq_nmod_mul(new_coeff, p->terms[i].coeff, scalar, p->ctx);
        fq_mvpoly_add_term(result, p->terms[i].var_exp, p->terms[i].par_exp, new_coeff);
        fq_nmod_clear(new_coeff, p->ctx);
    }
}

// Add polynomial b to polynomial a
void fq_mvpoly_add(fq_mvpoly_t *result, const fq_mvpoly_t *a, const fq_mvpoly_t *b) {
    fq_mvpoly_init(result, a->nvars, a->npars, a->ctx);
    
    // Add all terms from a
    for (slong i = 0; i < a->nterms; i++) {
        fq_mvpoly_add_term(result, a->terms[i].var_exp, a->terms[i].par_exp, a->terms[i].coeff);
    }
    
    // Add all terms from b
    for (slong i = 0; i < b->nterms; i++) {
        fq_mvpoly_add_term(result, b->terms[i].var_exp, b->terms[i].par_exp, b->terms[i].coeff);
    }
}

void print_fq_matrix_mvpoly(fq_mvpoly_t **matrix, slong nrows, slong ncols, 
                            const char *matrix_name, int show_details) {
    printf("\n=== %s Matrix (%ld x %ld) ===\n", matrix_name, nrows, ncols);
    
    if (show_details) {
        for (slong i = 0; i < nrows; i++) {
            for (slong j = 0; j < ncols; j++) {
                printf("M[%ld][%ld]: ", i, j);
                if (matrix[i][j].nterms == 0) {
                    printf("0\n");
                } else if (matrix[i][j].nterms <= 5) {
                    // 如果项数不多，显示完整表达式
                    fq_mvpoly_print_expanded(&matrix[i][j], "", 1);
                } else {
                    // 如果项数很多，只显示项数和最大次数
                    printf("%ld terms", matrix[i][j].nterms);
                    
                    // 计算最大次数
                    slong max_var_deg = 0, max_par_deg = 0;
                    for (slong t = 0; t < matrix[i][j].nterms; t++) {
                        if (matrix[i][j].terms[t].var_exp) {
                            for (slong k = 0; k < matrix[i][j].nvars; k++) {
                                if (matrix[i][j].terms[t].var_exp[k] > max_var_deg) {
                                    max_var_deg = matrix[i][j].terms[t].var_exp[k];
                                }
                            }
                        }
                        if (matrix[i][j].terms[t].par_exp && matrix[i][j].npars > 0) {
                            for (slong k = 0; k < matrix[i][j].npars; k++) {
                                if (matrix[i][j].terms[t].par_exp[k] > max_par_deg) {
                                    max_par_deg = matrix[i][j].terms[t].par_exp[k];
                                }
                            }
                        }
                    }
                    printf(" (max var deg: %ld, max par deg: %ld)\n", max_var_deg, max_par_deg);
                }
            }
        }
    } else {
        // 简化显示：只显示每个元素的项数
        printf("Matrix term counts:\n");
        for (slong i = 0; i < nrows; i++) {
            printf("  Row %ld: ", i);
            for (slong j = 0; j < ncols; j++) {
                printf("%3ld", matrix[i][j].nterms);
                if (j < ncols - 1) printf(" ");
            }
            printf("\n");
        }
    }
    printf("=== End %s Matrix ===\n\n", matrix_name);
}

// 分析矩阵的函数
void analyze_fq_matrix_mvpoly(fq_mvpoly_t **matrix, slong nrows, slong ncols, 
                              const char *matrix_name) {
    printf("\n--- Analysis of %s Matrix ---\n", matrix_name);
    
    slong total_terms = 0;
    slong zero_entries = 0;
    slong max_terms = 0;
    slong max_var_deg = 0;
    slong max_par_deg = 0;
    
    for (slong i = 0; i < nrows; i++) {
        for (slong j = 0; j < ncols; j++) {
            total_terms += matrix[i][j].nterms;
            if (matrix[i][j].nterms == 0) {
                zero_entries++;
            }
            if (matrix[i][j].nterms > max_terms) {
                max_terms = matrix[i][j].nterms;
            }
            
            // 分析次数
            for (slong t = 0; t < matrix[i][j].nterms; t++) {
                if (matrix[i][j].terms[t].var_exp) {
                    for (slong k = 0; k < matrix[i][j].nvars; k++) {
                        if (matrix[i][j].terms[t].var_exp[k] > max_var_deg) {
                            max_var_deg = matrix[i][j].terms[t].var_exp[k];
                        }
                    }
                }
                if (matrix[i][j].terms[t].par_exp && matrix[i][j].npars > 0) {
                    for (slong k = 0; k < matrix[i][j].npars; k++) {
                        if (matrix[i][j].terms[t].par_exp[k] > max_par_deg) {
                            max_par_deg = matrix[i][j].terms[t].par_exp[k];
                        }
                    }
                }
            }
        }
    }
    
    printf("  Total entries: %ld\n", nrows * ncols);
    printf("  Zero entries: %ld (%.1f%%)\n", zero_entries, 
           100.0 * zero_entries / (nrows * ncols));
    printf("  Non-zero entries: %ld\n", nrows * ncols - zero_entries);
    printf("  Total terms: %ld\n", total_terms);
    printf("  Average terms per entry: %.2f\n", (double)total_terms / (nrows * ncols));
    printf("  Maximum terms in single entry: %ld\n", max_terms);
    printf("  Maximum variable degree: %ld\n", max_var_deg);
    printf("  Maximum parameter degree: %ld\n", max_par_deg);
    printf("--- End Analysis ---\n\n");
}
// ============ Kronecker substitution helpers ============

slong exp_to_kronecker_index(const slong *exp, const slong *degs, slong n) {
    slong index = 0;
    slong stride = 1;
    
    for (slong i = 0; i < n; i++) {
        index += exp[i] * stride;
        stride *= (degs[i] + 1);
    }
    
    return index;
}

void kronecker_index_to_exp(slong index, slong *exp, const slong *degs, slong n) {
    for (slong i = 0; i < n; i++) {
        exp[i] = index % (degs[i] + 1);
        index /= (degs[i] + 1);
    }
}

// Convert fq_mvpoly to univariate via Kronecker
void fq_mvpoly_to_kronecker_full(fq_nmod_poly_t out, const fq_mvpoly_t *p, 
                                const slong *var_degs, const slong *par_degs) {
    // Calculate total degree
    slong total_vars = p->nvars + p->npars;
    slong *all_degs = (slong*) flint_malloc(total_vars * sizeof(slong));
    
    // Copy degree bounds
    for (slong i = 0; i < p->nvars; i++) {
        all_degs[i] = var_degs[i];
    }
    for (slong i = 0; i < p->npars; i++) {
        all_degs[p->nvars + i] = par_degs[i];
    }
    
    // Calculate maximum degree
    slong max_deg = 1;
    for (slong i = 0; i < total_vars; i++) {
        max_deg *= (all_degs[i] + 1);
    }
    max_deg--;
    
    fq_nmod_poly_fit_length(out, max_deg + 1, p->ctx);
    fq_nmod_poly_zero(out, p->ctx);
    
    // Convert each term
    slong *combined_exp = (slong*) flint_calloc(total_vars, sizeof(slong));
    
    for (slong i = 0; i < p->nterms; i++) {
        // Combine variable and parameter exponents
        for (slong j = 0; j < p->nvars; j++) {
            combined_exp[j] = p->terms[i].var_exp ? p->terms[i].var_exp[j] : 0;
        }
        for (slong j = 0; j < p->npars; j++) {
            combined_exp[p->nvars + j] = p->terms[i].par_exp ? p->terms[i].par_exp[j] : 0;
        }
        
        slong idx = exp_to_kronecker_index(combined_exp, all_degs, total_vars);
        
        fq_nmod_t old_coeff, new_coeff;
        fq_nmod_init(old_coeff, p->ctx);
        fq_nmod_init(new_coeff, p->ctx);
        
        fq_nmod_poly_get_coeff(old_coeff, out, idx, p->ctx);
        fq_nmod_add(new_coeff, old_coeff, p->terms[i].coeff, p->ctx);
        fq_nmod_poly_set_coeff(out, idx, new_coeff, p->ctx);
        
        fq_nmod_clear(old_coeff, p->ctx);
        fq_nmod_clear(new_coeff, p->ctx);
    }
    
    _fq_nmod_poly_normalise(out, p->ctx);
    
    flint_free(combined_exp);
    flint_free(all_degs);
}

void kronecker_to_fq_mvpoly_full(fq_mvpoly_t *out, const fq_nmod_poly_t in, 
                                const slong *var_degs, slong nvars,
                                const slong *par_degs, slong npars, const fq_nmod_ctx_t ctx) {
    fq_mvpoly_init(out, nvars, npars, ctx);
    
    slong deg = fq_nmod_poly_degree(in, ctx);
    if (deg < 0) return;
    
    // Pre-calculate non-zero terms
    slong nterms = 0;
    fq_nmod_t coeff;
    fq_nmod_init(coeff, ctx);
    
    for (slong i = 0; i <= deg; i++) {
        fq_nmod_poly_get_coeff(coeff, in, i, ctx);
        if (!fq_nmod_is_zero(coeff, ctx)) {
            nterms++;
        }
    }
    
    if (nterms == 0) {
        fq_nmod_clear(coeff, ctx);
        return;
    }
    
    // Pre-allocate space
    if (out->alloc < nterms) {
        out->alloc = nterms;
        out->terms = (fq_monomial_t*) flint_realloc(out->terms, out->alloc * sizeof(fq_monomial_t));
    }
    
    // Prepare work arrays
    slong total_vars = nvars + npars;
    slong *all_degs = (slong*) flint_malloc(total_vars * sizeof(slong));
    
    for (slong i = 0; i < nvars; i++) {
        all_degs[i] = var_degs[i];
    }
    for (slong i = 0; i < npars; i++) {
        all_degs[nvars + i] = par_degs[i];
    }
    
    slong *combined_exp = (slong*) flint_calloc(total_vars, sizeof(slong));
    
    // Fill terms
    slong term_idx = 0;
    for (slong i = 0; i <= deg; i++) {
        fq_nmod_poly_get_coeff(coeff, in, i, ctx);
        if (!fq_nmod_is_zero(coeff, ctx)) {
            // Decode exponents
            slong idx = i;
            for (slong j = 0; j < total_vars; j++) {
                combined_exp[j] = idx % (all_degs[j] + 1);
                idx /= (all_degs[j] + 1);
            }
            
            // Initialize and set coefficient
            fq_nmod_init(out->terms[term_idx].coeff, ctx);
            fq_nmod_set(out->terms[term_idx].coeff, coeff, ctx);
            
            if (nvars > 0) {
                out->terms[term_idx].var_exp = (slong*) flint_calloc(nvars, sizeof(slong));
                for (slong j = 0; j < nvars; j++) {
                    out->terms[term_idx].var_exp[j] = combined_exp[j];
                }
            } else {
                out->terms[term_idx].var_exp = NULL;
            }
            
            if (npars > 0) {
                out->terms[term_idx].par_exp = (slong*) flint_calloc(npars, sizeof(slong));
                for (slong j = 0; j < npars; j++) {
                    out->terms[term_idx].par_exp[j] = combined_exp[nvars + j];
                }
            } else {
                out->terms[term_idx].par_exp = NULL;
            }
            
            term_idx++;
        }
    }
    
    out->nterms = term_idx;
    
    fq_nmod_clear(coeff, ctx);
    flint_free(combined_exp);
    flint_free(all_degs);
}


// Division using FLINT's built-in functions (alternative approach)
void divide_by_fq_linear_factor_flint(fq_mvpoly_t *quotient, const fq_mvpoly_t *dividend,
                                     slong var_idx, slong nvars, slong npars) {
    // Initialize the quotient
    fq_mvpoly_init(quotient, nvars, npars, dividend->ctx);
    
    if (dividend->nterms == 0) return;
    
    // Create mpoly context for the combined variable and parameter space
    slong total_vars = nvars + npars;
    fq_nmod_mpoly_ctx_t mpoly_ctx;
    fq_nmod_mpoly_ctx_init(mpoly_ctx, total_vars, ORD_LEX, dividend->ctx);
    
    // Convert dividend to FLINT mpoly format
    fq_nmod_mpoly_t dividend_mpoly;
    fq_nmod_mpoly_init(dividend_mpoly, mpoly_ctx);
    fq_mvpoly_to_fq_nmod_mpoly(dividend_mpoly, dividend, mpoly_ctx);
    
    // Create the linear factor (x_i - ~x_i) as an mpoly
    fq_nmod_mpoly_t divisor;
    fq_nmod_mpoly_init(divisor, mpoly_ctx);
    
    slong actual_nvars = nvars / 2;  // Since we have dual variables
    
    // Create coefficient constants
    fq_nmod_t one, neg_one;
    fq_nmod_init(one, dividend->ctx);
    fq_nmod_init(neg_one, dividend->ctx);
    fq_nmod_one(one, dividend->ctx);
    fq_nmod_neg(neg_one, one, dividend->ctx);
    
    // Allocate exponent array
    ulong *exps = (ulong*) flint_calloc(total_vars, sizeof(ulong));
    
    // Add x_i term: coefficient = 1, exponent vector has 1 at position var_idx
    memset(exps, 0, total_vars * sizeof(ulong));
    exps[var_idx] = 1;  // x_i^1
    fq_nmod_mpoly_push_term_fq_nmod_ui(divisor, one, exps, mpoly_ctx);
    
    // Add -~x_i term: coefficient = -1, exponent vector has 1 at dual position
    memset(exps, 0, total_vars * sizeof(ulong));
    exps[var_idx + actual_nvars] = 1;  // ~x_i^1
    fq_nmod_mpoly_push_term_fq_nmod_ui(divisor, neg_one, exps, mpoly_ctx);
    
    // Finalize the divisor polynomial
    fq_nmod_mpoly_sort_terms(divisor, mpoly_ctx);
    fq_nmod_mpoly_combine_like_terms(divisor, mpoly_ctx);
    
    // Perform the division using FLINT's exact division
    fq_nmod_mpoly_t quotient_mpoly;
    fq_nmod_mpoly_init(quotient_mpoly, mpoly_ctx);
    
    // Try exact division first
    int divides = fq_nmod_mpoly_divides(quotient_mpoly, dividend_mpoly, divisor, mpoly_ctx);
    
    if (divides) {
        // Successful exact division - convert result back to fq_mvpoly_t
        fq_nmod_mpoly_to_fq_mvpoly(quotient, quotient_mpoly, nvars, npars, mpoly_ctx, dividend->ctx);
        
        //printf("    FLINT exact division successful\n");
    } else {
        // Exact division failed - try general division and check remainder
        fq_nmod_mpoly_t remainder;
        fq_nmod_mpoly_init(remainder, mpoly_ctx);
        
        fq_nmod_mpoly_divrem(quotient_mpoly, remainder, dividend_mpoly, divisor, mpoly_ctx);
        
        if (fq_nmod_mpoly_is_zero(remainder, mpoly_ctx)) {
            // Division is exact despite divides() returning false
            fq_nmod_mpoly_to_fq_mvpoly(quotient, quotient_mpoly, nvars, npars, mpoly_ctx, dividend->ctx);
            //printf("    FLINT general division successful (zero remainder)\n");
        } else {
            // Division is not exact - fall back to direct method
            printf("    Warning: FLINT division not exact, falling back to direct method\n");
            fq_mvpoly_clear(quotient);
            divide_by_fq_linear_factor_improved(quotient, dividend, var_idx, nvars, npars);
        }
        
        fq_nmod_mpoly_clear(remainder, mpoly_ctx);
    }
    
    // Cleanup
    flint_free(exps);
    fq_nmod_clear(one, dividend->ctx);
    fq_nmod_clear(neg_one, dividend->ctx);
    fq_nmod_mpoly_clear(dividend_mpoly, mpoly_ctx);
    fq_nmod_mpoly_clear(divisor, mpoly_ctx);
    fq_nmod_mpoly_clear(quotient_mpoly, mpoly_ctx);
    fq_nmod_mpoly_ctx_clear(mpoly_ctx);
}

// Improved division by linear factor (x_i - ~x_i) - direct method
void divide_by_fq_linear_factor_improved(fq_mvpoly_t *quotient, const fq_mvpoly_t *dividend,
                                        slong var_idx, slong nvars, slong npars) {
    fq_mvpoly_init(quotient, nvars, npars, dividend->ctx);
    
    if (dividend->nterms == 0) return;
    
    slong actual_nvars = nvars / 2;  // Since we have dual variables
    
    // Process each term in the dividend
    for (slong i = 0; i < dividend->nterms; i++) {
        slong deg_x = dividend->terms[i].var_exp ? dividend->terms[i].var_exp[var_idx] : 0;
        slong deg_dual = dividend->terms[i].var_exp ? dividend->terms[i].var_exp[var_idx + actual_nvars] : 0;
        
        if (deg_x == deg_dual) {
            // This term is divisible by (x_i - ~x_i), the result has the same exponents
            fq_mvpoly_add_term(quotient, dividend->terms[i].var_exp, dividend->terms[i].par_exp, 
                              dividend->terms[i].coeff);
        } else if (deg_x > deg_dual) {
            // Reduce x degree by 1
            slong *new_exp = (slong*) flint_calloc(nvars, sizeof(slong));
            if (dividend->terms[i].var_exp) {
                memcpy(new_exp, dividend->terms[i].var_exp, nvars * sizeof(slong));
            }
            new_exp[var_idx]--;
            
            fq_mvpoly_add_term(quotient, new_exp, dividend->terms[i].par_exp, 
                              dividend->terms[i].coeff);
            flint_free(new_exp);
        } else { // deg_dual > deg_x
            // Reduce ~x degree by 1, with negative coefficient
            slong *new_exp = (slong*) flint_calloc(nvars, sizeof(slong));
            if (dividend->terms[i].var_exp) {
                memcpy(new_exp, dividend->terms[i].var_exp, nvars * sizeof(slong));
            }
            new_exp[var_idx + actual_nvars]--;
            
            fq_nmod_t neg_coeff;
            fq_nmod_init(neg_coeff, dividend->ctx);
            fq_nmod_neg(neg_coeff, dividend->terms[i].coeff, dividend->ctx);
            
            fq_mvpoly_add_term(quotient, new_exp, dividend->terms[i].par_exp, neg_coeff);
            
            fq_nmod_clear(neg_coeff, dividend->ctx);
            flint_free(new_exp);
        }
    }
}


// ============ Evaluation functions ============

void evaluate_fq_mvpoly_at_params(fq_nmod_t result, const fq_mvpoly_t *poly, const fq_nmod_t *param_vals) {
    fq_nmod_zero(result, poly->ctx);
    
    for (slong i = 0; i < poly->nterms; i++) {
        fq_nmod_t term_val;
        fq_nmod_init(term_val, poly->ctx);
        fq_nmod_set(term_val, poly->terms[i].coeff, poly->ctx);
        
        // Multiply by parameter powers
        if (poly->npars > 0 && poly->terms[i].par_exp) {
            for (slong j = 0; j < poly->npars; j++) {
                for (slong k = 0; k < poly->terms[i].par_exp[j]; k++) {
                    fq_nmod_mul(term_val, term_val, param_vals[j], poly->ctx);
                }
            }
        }
        
        fq_nmod_add(result, result, term_val, poly->ctx);
        fq_nmod_clear(term_val, poly->ctx);
    }
}


// ============ Matrix operations ============

// Build cancellation matrix in multivariate form
void build_fq_cancellation_matrix_mvpoly(fq_mvpoly_t ***M, fq_mvpoly_t *polys, 
                                        slong nvars, slong npars) {
    slong n = nvars + 1;
    
    // Allocate matrix space
    *M = (fq_mvpoly_t**) flint_malloc(n * sizeof(fq_mvpoly_t*));
    for (slong i = 0; i < n; i++) {
        (*M)[i] = (fq_mvpoly_t*) flint_malloc(n * sizeof(fq_mvpoly_t));
    }
    
    printf("\nBuilding %ld x %ld cancellation matrix (multivariate form)\n", n, n);
    
    // Build matrix entries
    for (slong i = 0; i < n; i++) {
        for (slong j = 0; j < n; j++) {
            fq_mvpoly_init(&(*M)[i][j], 2 * nvars, npars, polys[0].ctx);
            
            // Substitute variables according to row index
            for (slong t = 0; t < polys[j].nterms; t++) {
                slong *new_var_exp = (slong*) flint_calloc(2 * nvars, sizeof(slong));
                
                for (slong k = 0; k < nvars; k++) {
                    slong orig_exp = polys[j].terms[t].var_exp ? polys[j].terms[t].var_exp[k] : 0;
                    
                    if (k < i) {
                        // Use dual variable ~x_k
                        new_var_exp[nvars + k] = orig_exp;
                    } else {
                        // Use original variable x_k
                        new_var_exp[k] = orig_exp;
                    }
                }
                
                fq_mvpoly_add_term(&(*M)[i][j], new_var_exp, polys[j].terms[t].par_exp, 
                                  polys[j].terms[t].coeff);
                flint_free(new_var_exp);
            }
        }
    }
}

void perform_fq_matrix_row_operations_mvpoly(fq_mvpoly_t ***new_matrix, fq_mvpoly_t ***original_matrix,
                                            slong nvars, slong npars) {
    slong n = nvars + 1;
    
    // Allocate new matrix space
    *new_matrix = (fq_mvpoly_t**) flint_malloc(n * sizeof(fq_mvpoly_t*));
    for (slong i = 0; i < n; i++) {
        (*new_matrix)[i] = (fq_mvpoly_t*) flint_malloc(n * sizeof(fq_mvpoly_t));
    }
    
    printf("\nPerforming row operations on matrix:\n");
    
    // First row remains unchanged
    printf("Row 0: unchanged\n");
    for (slong j = 0; j < n; j++) {
        fq_mvpoly_copy(&(*new_matrix)[0][j], &(*original_matrix)[0][j]);
    }
    
    // Process remaining rows
    for (slong i = 1; i < n; i++) {
        printf("Row %ld: (row[%ld] - row[%ld]) / (x%ld - ~x%ld)\n", i, i, i-1, i-1, i-1);
        
        for (slong j = 0; j < n; j++) {
            // Calculate difference
            fq_mvpoly_t diff;
            fq_mvpoly_init(&diff, 2 * nvars, npars, (*original_matrix)[0][0].ctx);
            
            // diff = M[i][j] - M[i-1][j]
            fq_mvpoly_copy(&diff, &(*original_matrix)[i][j]);
            for (slong t = 0; t < (*original_matrix)[i-1][j].nterms; t++) {
                fq_nmod_t neg_coeff;
                fq_nmod_init(neg_coeff, diff.ctx);
                fq_nmod_neg(neg_coeff, (*original_matrix)[i-1][j].terms[t].coeff, diff.ctx);
                
                fq_mvpoly_add_term(&diff, 
                                  (*original_matrix)[i-1][j].terms[t].var_exp,
                                  (*original_matrix)[i-1][j].terms[t].par_exp,
                                  neg_coeff);
                                  
                fq_nmod_clear(neg_coeff, diff.ctx);
            }
            
            // Perform division using improved method
            fq_mvpoly_init(&(*new_matrix)[i][j], 2 * nvars, npars, diff.ctx);
            
            if (diff.nterms > 0) {
                // Use the improved division method (could also use divide_by_fq_linear_factor_flint)
                divide_by_fq_linear_factor_flint(&(*new_matrix)[i][j], &diff, 
                                                   i-1, 2*nvars, npars);
            }
            
            fq_mvpoly_clear(&diff);
        }
    }
}




// Compute Dixon resultant degree bound
slong compute_fq_dixon_resultant_degree_bound(fq_mvpoly_t *polys, slong npolys, slong nvars, slong npars) {
    slong degree_product = 1;
    
    for (slong i = 0; i < npolys; i++) {
        slong max_total_deg = 0;
        
        // Find maximum total degree of polynomial i
        for (slong t = 0; t < polys[i].nterms; t++) {
            slong total_deg = 0;
            
            // Sum variable degrees
            if (polys[i].terms[t].var_exp) {
                for (slong j = 0; j < nvars; j++) {
                    total_deg += polys[i].terms[t].var_exp[j];
                }
            }
            
            // Sum parameter degrees
            if (polys[i].terms[t].par_exp && npars > 0) {
                for (slong j = 0; j < npars; j++) {
                    total_deg += polys[i].terms[t].par_exp[j];
                }
            }
            
            if (total_deg > max_total_deg) {
                max_total_deg = total_deg;
            }
        }
        
        degree_product *= max_total_deg;
    }
    
    return degree_product + 1;
}

// 修复的系数矩阵行列式计算函数
// 修复的系数矩阵行列式计算函数 - 使用Mulders-Storjohann算法处理单变量情况
void compute_fq_coefficient_matrix_det(fq_mvpoly_t *result, fq_mvpoly_t **coeff_matrix,
                                       slong size, slong npars, const fq_nmod_ctx_t ctx,
                                       det_method_t method, slong res_deg_bound) {
    if (size == 0) {
        fq_mvpoly_init(result, 0, npars, ctx);
        return;
    }
    
    fq_mvpoly_init(result, 0, npars, ctx);
    
    if (npars == 0) {
        // No parameters - scalar entries
        fq_nmod_mat_t scalar_mat;
        fq_nmod_mat_init(scalar_mat, size, size, ctx);
        
        for (slong i = 0; i < size; i++) {
            for (slong j = 0; j < size; j++) {
                if (coeff_matrix[i][j].nterms > 0) {
                    fq_nmod_set(fq_nmod_mat_entry(scalar_mat, i, j), 
                                coeff_matrix[i][j].terms[0].coeff, ctx);
                } else {
                    fq_nmod_zero(fq_nmod_mat_entry(scalar_mat, i, j), ctx);
                }
            }
        }
        
        printf("\nComputing Resultant\n");
        clock_t start = clock();
        
        fq_nmod_t det;
        fq_nmod_init(det, ctx);
        fq_nmod_mat_det(det, scalar_mat, ctx);
        
        clock_t end = clock();
        double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
        
        printf("End (%.3f seconds)\n", elapsed);
        
        if (!fq_nmod_is_zero(det, ctx)) {
            fq_mvpoly_add_term(result, NULL, NULL, det);
        }
        
        fq_nmod_clear(det, ctx);
        fq_nmod_mat_clear(scalar_mat, ctx);
        
    } else if (npars == 1) {
        // One parameter - use Mulders-Storjohann algorithm
        printf("\nComputing Resultant using Mulders-Storjohann algorithm (univariate case)\n");
        clock_t start = clock();
        
        // Create polynomial matrix using our fq_nmod_poly_mat_t structure
        fq_nmod_poly_mat_t poly_mat;
        fq_nmod_poly_mat_init(poly_mat, size, size, ctx);
        
        // Convert mvpoly coefficients to polynomial matrix
        for (slong i = 0; i < size; i++) {
            for (slong j = 0; j < size; j++) {
                fq_nmod_poly_struct *entry = fq_nmod_poly_mat_entry(poly_mat, i, j);
                fq_nmod_poly_zero(entry, ctx);
                
                // Convert mvpoly to fq_nmod_poly (using parameter as polynomial variable)
                for (slong t = 0; t < coeff_matrix[i][j].nterms; t++) {
                    slong deg = coeff_matrix[i][j].terms[t].par_exp ? 
                               coeff_matrix[i][j].terms[t].par_exp[0] : 0;
                    fq_nmod_poly_set_coeff(entry, deg, 
                                          coeff_matrix[i][j].terms[t].coeff, ctx);
                }
            }
        }
        
        // Compute determinant using Mulders-Storjohann algorithm
        fq_nmod_poly_t det_poly;
        fq_nmod_poly_init(det_poly, ctx);
        
        printf("  Matrix size: %ld x %ld\n", size, size);
        printf("  Using weak Popov form method...\n");
        
        fq_nmod_poly_mat_det_iter(det_poly, poly_mat, ctx);
        
        // Convert result back to fq_mvpoly
        slong det_deg = fq_nmod_poly_degree(det_poly, ctx);
        if (det_deg >= 0) {
            for (slong i = 0; i <= det_deg; i++) {
                fq_nmod_t coeff;
                fq_nmod_init(coeff, ctx);
                fq_nmod_poly_get_coeff(coeff, det_poly, i, ctx);
                if (!fq_nmod_is_zero(coeff, ctx)) {
                    slong par_exp[1] = {i};
                    fq_mvpoly_add_term(result, NULL, par_exp, coeff);
                }
                fq_nmod_clear(coeff, ctx);
            }
        }
        
        printf("  Determinant degree: %ld\n", det_deg);
        printf("  Result terms: %ld\n", result->nterms);
        
        // Cleanup
        fq_nmod_poly_clear(det_poly, ctx);
        fq_nmod_poly_mat_clear(poly_mat, ctx);
        
        clock_t end = clock();
        double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
        printf("End (%.3f seconds)\n", elapsed);
        
    } else {
        // Multiple parameters (npars >= 2) - use interpolation method
        printf("\nMultiple parameters detected. Using fq_nmod interpolation method.\n");
        printf("  Parameters: %ld\n", npars);
        printf("  Matrix size: %ld x %ld\n", size, size);
        
        clock_t start = clock();
        // Use the fixed fq interpolation
        fq_compute_det_by_interpolation(result, coeff_matrix, size,
                                       0, npars, ctx, res_deg_bound);
        
        clock_t end = clock();
        double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
        printf("End (%.3f seconds)\n", elapsed);
    }
}



// ============ Submatrix finding functions ============

// Check if a row is linearly dependent
static int is_fq_linearly_dependent_row(fq_nmod_mat_t mat, slong new_row_idx, 
                                       slong *selected_rows, slong num_selected_rows,
                                       slong ncols, const fq_nmod_ctx_t ctx) {
    if (num_selected_rows == 0) return 0;
    
    // Build matrix with selected rows and new row
    fq_nmod_mat_t test_mat;
    fq_nmod_mat_init(test_mat, num_selected_rows + 1, ncols, ctx);
    
    // Copy selected rows
    for (slong i = 0; i < num_selected_rows; i++) {
        for (slong j = 0; j < ncols; j++) {
            fq_nmod_set(fq_nmod_mat_entry(test_mat, i, j), 
                       fq_nmod_mat_entry(mat, selected_rows[i], j), ctx);
        }
    }
    
    // Add new row
    for (slong j = 0; j < ncols; j++) {
        fq_nmod_set(fq_nmod_mat_entry(test_mat, num_selected_rows, j),
                   fq_nmod_mat_entry(mat, new_row_idx, j), ctx);
    }
    
    // Compute rank
    slong rank = fq_nmod_mat_rank(test_mat, ctx);
    fq_nmod_mat_clear(test_mat, ctx);
    
    // If rank didn't increase, new row is linearly dependent
    return (rank == num_selected_rows);
}

// Check if a column is linearly dependent
static int is_fq_linearly_dependent_col(fq_nmod_mat_t mat, slong new_col_idx,
                                       slong *selected_cols, slong num_selected_cols,
                                       slong *selected_rows, slong num_selected_rows,
                                       const fq_nmod_ctx_t ctx) {
    if (num_selected_cols == 0 || num_selected_rows == 0) return 0;
    
    // Build submatrix with selected rows
    fq_nmod_mat_t test_mat;
    fq_nmod_mat_init(test_mat, num_selected_rows, num_selected_cols + 1, ctx);
    
    // Copy selected columns
    for (slong i = 0; i < num_selected_rows; i++) {
        for (slong j = 0; j < num_selected_cols; j++) {
            fq_nmod_set(fq_nmod_mat_entry(test_mat, i, j),
                       fq_nmod_mat_entry(mat, selected_rows[i], selected_cols[j]), ctx);
        }
        // Add new column
        fq_nmod_set(fq_nmod_mat_entry(test_mat, i, num_selected_cols),
                   fq_nmod_mat_entry(mat, selected_rows[i], new_col_idx), ctx);
    }
    
    // Compute rank
    slong rank = fq_nmod_mat_rank(test_mat, ctx);
    fq_nmod_mat_clear(test_mat, ctx);
    
    // If rank didn't increase, new column is linearly dependent
    return (rank == num_selected_cols);
}



// Helper function to compute total degree of a polynomial
static slong compute_fq_polynomial_total_degree(fq_mvpoly_t *poly, slong npars) {
    if (poly->nterms == 0) {
        return 0;
    }
    
    slong max_degree = -1;
    
    for (slong t = 0; t < poly->nterms; t++) {
        slong term_total_degree = 0;
        
        // Variable degrees
        if (poly->terms[t].var_exp) {
            for (slong v = 0; v < poly->nvars; v++) {
                term_total_degree += poly->terms[t].var_exp[v];
            }
        }
        
        // Parameter degrees
        if (poly->terms[t].par_exp && npars > 0) {
            for (slong p = 0; p < npars; p++) {
                term_total_degree += poly->terms[t].par_exp[p];
            }
        }
        
        if (max_degree == -1 || term_total_degree > max_degree) {
            max_degree = term_total_degree;
        }
    }
    
    return max_degree;
}

// Compute row maximum total degree
static slong compute_fq_row_max_total_degree(fq_mvpoly_t **matrix_row, slong ncols, slong npars) {
    slong max_degree = -1;
    
    for (slong j = 0; j < ncols; j++) {
        slong poly_deg = compute_fq_polynomial_total_degree(matrix_row[j], npars);
        if (poly_deg > max_degree) {
            max_degree = poly_deg;
        }
    }
    
    return max_degree;
}

// Compute column maximum total degree
static slong compute_fq_col_max_total_degree(fq_mvpoly_t ***matrix, slong col_idx, slong nrows, slong npars) {
    slong max_degree = -1;
    
    for (slong i = 0; i < nrows; i++) {
        slong poly_deg = compute_fq_polynomial_total_degree(matrix[i][col_idx], npars);
        if (poly_deg > max_degree) {
            max_degree = poly_deg;
        }
    }
    
    return max_degree;
}

// 修复后的函数
// 简化优化版本 - 直接替换原有函数
// 只需要在 dixon_flint.c 中替换 find_fq_optimal_maximal_rank_submatrix 函数

// 辅助结构：维护当前行空间的基
typedef struct {
    fq_nmod_mat_t reduced_rows;    // 已约化的行向量
    slong *pivot_cols;             // 每行的主元列位置  
    slong *selected_indices;       // 已选择的原始行索引
    slong current_rank;            // 当前秩
} row_basis_tracker_t;

// 初始化行基跟踪器
static void row_basis_tracker_init(row_basis_tracker_t *tracker, slong max_size, slong ncols, 
                                  const fq_nmod_ctx_struct *ctx) {
    fq_nmod_mat_init(tracker->reduced_rows, max_size, ncols, ctx);
    tracker->pivot_cols = (slong*) flint_malloc(max_size * sizeof(slong));
    tracker->selected_indices = (slong*) flint_malloc(max_size * sizeof(slong));
    tracker->current_rank = 0;
    
    for (slong i = 0; i < max_size; i++) {
        tracker->pivot_cols[i] = -1;
    }
}

// 清理行基跟踪器
static void row_basis_tracker_clear(row_basis_tracker_t *tracker, const fq_nmod_ctx_struct *ctx) {
    fq_nmod_mat_clear(tracker->reduced_rows, ctx);
    flint_free(tracker->pivot_cols);
    flint_free(tracker->selected_indices);
}

// 尝试添加新行到基中，返回1表示成功（线性无关），0表示失败（线性相关）
static int try_add_row_to_basis(row_basis_tracker_t *tracker, const fq_nmod_mat_t eval_mat, 
                               slong new_row_idx, slong ncols, const fq_nmod_ctx_struct *ctx) {
    // 复制新行到工作向量
    fq_nmod_t *work_row = (fq_nmod_t*) flint_malloc(ncols * sizeof(fq_nmod_t));
    for (slong j = 0; j < ncols; j++) {
        fq_nmod_init(work_row[j], ctx);
        fq_nmod_set(work_row[j], fq_nmod_mat_entry(eval_mat, new_row_idx, j), ctx);
    }
    
    // 用现有基对新行进行高斯消元
    for (slong i = 0; i < tracker->current_rank; i++) {
        slong pivot_col = tracker->pivot_cols[i];
        if (pivot_col >= 0 && !fq_nmod_is_zero(work_row[pivot_col], ctx)) {
            // 计算消元因子
            fq_nmod_t factor;
            fq_nmod_init(factor, ctx);
            fq_nmod_div(factor, work_row[pivot_col], 
                       fq_nmod_mat_entry(tracker->reduced_rows, i, pivot_col), ctx);
            
            // 消元: work_row = work_row - factor * reduced_rows[i]
            for (slong j = 0; j < ncols; j++) {
                fq_nmod_t temp;
                fq_nmod_init(temp, ctx);
                fq_nmod_mul(temp, factor, fq_nmod_mat_entry(tracker->reduced_rows, i, j), ctx);
                fq_nmod_sub(work_row[j], work_row[j], temp, ctx);
                fq_nmod_clear(temp, ctx);
            }
            fq_nmod_clear(factor, ctx);
        }
    }
    
    // 检查消元后是否为零向量
    slong first_nonzero = -1;
    for (slong j = 0; j < ncols; j++) {
        if (!fq_nmod_is_zero(work_row[j], ctx)) {
            first_nonzero = j;
            break;
        }
    }
    
    if (first_nonzero == -1) {
        // 零向量，线性相关
        for (slong j = 0; j < ncols; j++) {
            fq_nmod_clear(work_row[j], ctx);
        }
        flint_free(work_row);
        return 0;
    }
    
    // 非零向量，线性无关，标准化并添加到基中
    fq_nmod_t leading;
    fq_nmod_init(leading, ctx);
    fq_nmod_set(leading, work_row[first_nonzero], ctx);
    
    // 标准化：首个非零元素变为1
    for (slong j = 0; j < ncols; j++) {
        fq_nmod_div(work_row[j], work_row[j], leading, ctx);
        fq_nmod_set(fq_nmod_mat_entry(tracker->reduced_rows, tracker->current_rank, j), 
                   work_row[j], ctx);
    }
    
    // 更新基信息
    tracker->pivot_cols[tracker->current_rank] = first_nonzero;
    tracker->selected_indices[tracker->current_rank] = new_row_idx;
    tracker->current_rank++;
    
    // 清理
    fq_nmod_clear(leading, ctx);
    for (slong j = 0; j < ncols; j++) {
        fq_nmod_clear(work_row[j], ctx);
    }
    flint_free(work_row);
    
    return 1;
}

// 优化的主函数 - 直接替换原有版本
void find_fq_optimal_maximal_rank_submatrix(fq_mvpoly_t ***full_matrix, 
                                           slong nrows, slong ncols,
                                           slong **row_indices_out, 
                                           slong **col_indices_out,
                                           slong *num_rows, slong *num_cols,
                                           slong npars) {
    printf("Finding maximal rank submatrix using optimized incremental method...\n");
    
    // 获取上下文
    const fq_nmod_ctx_struct *ctx = full_matrix[0][0]->ctx;
    
    // 1. 计算度数并排序（复用原有函数）
    fq_index_degree_pair *row_degrees = (fq_index_degree_pair*) flint_malloc(nrows * sizeof(fq_index_degree_pair));
    fq_index_degree_pair *col_degrees = (fq_index_degree_pair*) flint_malloc(ncols * sizeof(fq_index_degree_pair));
    
    for (slong i = 0; i < nrows; i++) {
        row_degrees[i].index = i;
        row_degrees[i].degree = compute_fq_row_max_total_degree(full_matrix[i], ncols, npars);
    }
    
    for (slong j = 0; j < ncols; j++) {
        col_degrees[j].index = j;
        col_degrees[j].degree = compute_fq_col_max_total_degree(full_matrix, j, nrows, npars);
    }
    
    qsort(row_degrees, nrows, sizeof(fq_index_degree_pair), compare_fq_degrees);
    qsort(col_degrees, ncols, sizeof(fq_index_degree_pair), compare_fq_degrees);
    
    // 2. 评估矩阵
    printf("Evaluating matrix at parameter values...\n");
    fq_nmod_t *param_vals = (fq_nmod_t*) flint_malloc(npars * sizeof(fq_nmod_t));
    for (slong i = 0; i < npars; i++) {
        fq_nmod_init(param_vals[i], ctx);
        fq_nmod_set_si(param_vals[i], 2 + i, ctx);
    }
    
    fq_nmod_mat_t eval_mat;
    fq_nmod_mat_init(eval_mat, nrows, ncols, ctx);
    
    for (slong i = 0; i < nrows; i++) {
        for (slong j = 0; j < ncols; j++) {
            fq_nmod_t val;
            fq_nmod_init(val, ctx);
            evaluate_fq_mvpoly_at_params(val, full_matrix[i][j], param_vals);
            fq_nmod_set(fq_nmod_mat_entry(eval_mat, i, j), val, ctx);
            fq_nmod_clear(val, ctx);
        }
    }
    
    // 3. 使用优化方法选择行
    printf("Selecting rows using incremental Gaussian elimination...\n");
    row_basis_tracker_t row_tracker;
    row_basis_tracker_init(&row_tracker, FLINT_MIN(nrows, ncols), ncols, ctx);
    
    for (slong i = 0; i < nrows && row_tracker.current_rank < FLINT_MIN(nrows, ncols); i++) {
        slong row_idx = row_degrees[i].index;
        
        if (try_add_row_to_basis(&row_tracker, eval_mat, row_idx, ncols, ctx)) {
            /*
            if (row_tracker.current_rank % 10 == 0 || row_tracker.current_rank <= 10) {
                printf("  Selected row %ld (total: %ld)\n", row_idx, row_tracker.current_rank);
            }
            */
        }
    }
    
    printf("Selected %ld linearly independent rows\n", row_tracker.current_rank);
    
    // 4. 选择列（简化版本，基于所选行的子矩阵）
    printf("Selecting columns...\n");
    
    // 构建所选行的子矩阵
    fq_nmod_mat_t selected_submat;
    fq_nmod_mat_init(selected_submat, row_tracker.current_rank, ncols, ctx);
    
    for (slong i = 0; i < row_tracker.current_rank; i++) {
        for (slong j = 0; j < ncols; j++) {
            fq_nmod_set(fq_nmod_mat_entry(selected_submat, i, j),
                       fq_nmod_mat_entry(eval_mat, row_tracker.selected_indices[i], j), ctx);
        }
    }
    
    // 选择列基
    row_basis_tracker_t col_tracker;
    row_basis_tracker_init(&col_tracker, ncols, row_tracker.current_rank, ctx);
    
    // 转置思考：选择列等价于在转置矩阵上选择行
    fq_nmod_mat_t transposed;
    fq_nmod_mat_init(transposed, ncols, row_tracker.current_rank, ctx);
    fq_nmod_mat_transpose(transposed, selected_submat, ctx);
    
    for (slong j = 0; j < ncols && col_tracker.current_rank < row_tracker.current_rank; j++) {
        slong col_idx = col_degrees[j].index;
        
        if (try_add_row_to_basis(&col_tracker, transposed, col_idx, row_tracker.current_rank, ctx)) {
            /*
            if (col_tracker.current_rank % 10 == 0 || col_tracker.current_rank <= 10) {
                printf("  Selected col %ld (total: %ld)\n", col_idx, col_tracker.current_rank);
            }
            */
        }
    }
    
    printf("Selected %ld linearly independent columns\n", col_tracker.current_rank);
    
    // 5. 构建最终结果
    slong final_size = FLINT_MIN(row_tracker.current_rank, col_tracker.current_rank);
    
    *row_indices_out = (slong*) flint_malloc(final_size * sizeof(slong));
    *col_indices_out = (slong*) flint_malloc(final_size * sizeof(slong));
    
    for (slong i = 0; i < final_size; i++) {
        (*row_indices_out)[i] = row_tracker.selected_indices[i];
        (*col_indices_out)[i] = col_tracker.selected_indices[i];
    }
    
    *num_rows = final_size;
    *num_cols = final_size;
    
    // 6. 验证结果
    printf("Verifying final %ld x %ld submatrix...\n", final_size, final_size);
    fq_nmod_mat_t final_mat;
    fq_nmod_mat_init(final_mat, final_size, final_size, ctx);
    
    for (slong i = 0; i < final_size; i++) {
        for (slong j = 0; j < final_size; j++) {
            fq_nmod_set(fq_nmod_mat_entry(final_mat, i, j),
                       fq_nmod_mat_entry(eval_mat, (*row_indices_out)[i], (*col_indices_out)[j]), ctx);
        }
    }
    
    slong final_rank = fq_nmod_mat_rank(final_mat, ctx);
    printf("Final rank: %ld/%ld %s\n", final_rank, final_size, 
           (final_rank == final_size) ? "✓" : "✗");
    
    // 清理
    row_basis_tracker_clear(&row_tracker, ctx);
    row_basis_tracker_clear(&col_tracker, ctx);
    flint_free(row_degrees);
    flint_free(col_degrees);
    
    for (slong i = 0; i < npars; i++) {
        fq_nmod_clear(param_vals[i], ctx);
    }
    flint_free(param_vals);
    
    fq_nmod_mat_clear(eval_mat, ctx);
    fq_nmod_mat_clear(selected_submat, ctx);
    fq_nmod_mat_clear(transposed, ctx);
    fq_nmod_mat_clear(final_mat, ctx);
}
// 优化版本的 find_fq_optimal_maximal_rank_submatrix
// ============ Extract coefficient matrix ============

void extract_fq_coefficient_matrix_from_dixon(fq_mvpoly_t ***coeff_matrix,
                                              slong *row_indices, slong *col_indices,
                                              slong *matrix_size,
                                              const fq_mvpoly_t *dixon_poly,
                                              slong nvars, slong npars) {
    printf("\nExtracting coefficient matrix from Dixon polynomial\n");
    
    // First compute degree bounds
    slong *d0 = (slong*) flint_calloc(nvars, sizeof(slong));
    slong *d1 = (slong*) flint_calloc(nvars, sizeof(slong));
    
    // Find maximum degree for each variable
    for (slong i = 0; i < dixon_poly->nterms; i++) {
        if (dixon_poly->terms[i].var_exp) {
            // Check x variables (first nvars)
            for (slong j = 0; j < nvars; j++) {
                if (dixon_poly->terms[i].var_exp[j] > d0[j]) {
                    d0[j] = dixon_poly->terms[i].var_exp[j];
                }
            }
            // Check ~x variables (next nvars)
            for (slong j = 0; j < nvars; j++) {
                if (dixon_poly->terms[i].var_exp[nvars + j] > d1[j]) {
                    d1[j] = dixon_poly->terms[i].var_exp[nvars + j];
                }
            }
        }
    }
    
    // Add 1 to each degree bound
    for (slong i = 0; i < nvars; i++) {
        d0[i]++;
        d1[i]++;
    }
    
    printf("Degree bounds - x vars: ");
    for (slong i = 0; i < nvars; i++) printf("%ld ", d0[i]);
    printf("\nDegree bounds - ~x vars: ");
    for (slong i = 0; i < nvars; i++) printf("%ld ", d1[i]);
    printf("\n");
    
    // Calculate expected matrix dimensions
    slong expected_rows = 1;
    slong expected_cols = 1;
    for (slong i = 0; i < nvars; i++) {
        expected_rows *= d0[i];
        expected_cols *= d1[i];
    }
    printf("Expected matrix size: %ld x %ld\n", expected_rows, expected_cols);
    
    // Collect monomials
    typedef struct {
        slong *exp;
        slong idx;
    } monom_t;
    
    monom_t *x_monoms = (monom_t*) flint_malloc(dixon_poly->nterms * sizeof(monom_t));
    monom_t *dual_monoms = (monom_t*) flint_malloc(dixon_poly->nterms * sizeof(monom_t));
    slong nx_monoms = 0, ndual_monoms = 0;
    
    // Collect unique monomials that satisfy degree bounds
    for (slong i = 0; i < dixon_poly->nterms; i++) {
        // Check if this term satisfies degree bounds
        int valid = 1;
        if (dixon_poly->terms[i].var_exp) {
            for (slong k = 0; k < nvars; k++) {
                if (dixon_poly->terms[i].var_exp[k] >= d0[k] || 
                    dixon_poly->terms[i].var_exp[nvars + k] >= d1[k]) {
                    valid = 0;
                    break;
                }
            }
        }
        
        if (!valid) continue;
        
        // Check x-monomial
        int found = 0;
        for (slong j = 0; j < nx_monoms; j++) {
            int same = 1;
            for (slong k = 0; k < nvars; k++) {
                if (x_monoms[j].exp[k] != dixon_poly->terms[i].var_exp[k]) {
                    same = 0;
                    break;
                }
            }
            if (same) {
                found = 1;
                break;
            }
        }
        if (!found) {
            x_monoms[nx_monoms].exp = (slong*) flint_calloc(nvars, sizeof(slong));
            x_monoms[nx_monoms].idx = nx_monoms;
            for (slong k = 0; k < nvars; k++) {
                x_monoms[nx_monoms].exp[k] = dixon_poly->terms[i].var_exp[k];
            }
            nx_monoms++;
        }
        
        // Check ~x-monomial
        found = 0;
        for (slong j = 0; j < ndual_monoms; j++) {
            int same = 1;
            for (slong k = 0; k < nvars; k++) {
                if (dual_monoms[j].exp[k] != dixon_poly->terms[i].var_exp[nvars + k]) {
                    same = 0;
                    break;
                }
            }
            if (same) {
                found = 1;
                break;
            }
        }
        if (!found) {
            dual_monoms[ndual_monoms].exp = (slong*) flint_calloc(nvars, sizeof(slong));
            dual_monoms[ndual_monoms].idx = ndual_monoms;
            for (slong k = 0; k < nvars; k++) {
                dual_monoms[ndual_monoms].exp[k] = dixon_poly->terms[i].var_exp[nvars + k];
            }
            ndual_monoms++;
        }
    }
    
    printf("Found %ld x-monomials and %ld ~x-monomials (after degree filtering)\n", 
           nx_monoms, ndual_monoms);
    
    if (nx_monoms == 0 || ndual_monoms == 0) {
        printf("Warning: Empty coefficient matrix\n");
        *matrix_size = 0;
        flint_free(d0);
        flint_free(d1);
        flint_free(x_monoms);
        flint_free(dual_monoms);
        return;
    }
    
    // Build full coefficient matrix (entries are polynomials in parameters)
    fq_mvpoly_t ***full_matrix = (fq_mvpoly_t***) flint_malloc(nx_monoms * sizeof(fq_mvpoly_t**));
    for (slong i = 0; i < nx_monoms; i++) {
        full_matrix[i] = (fq_mvpoly_t**) flint_malloc(ndual_monoms * sizeof(fq_mvpoly_t*));
        for (slong j = 0; j < ndual_monoms; j++) {
            full_matrix[i][j] = (fq_mvpoly_t*) flint_malloc(sizeof(fq_mvpoly_t));
            fq_mvpoly_init(full_matrix[i][j], 0, npars, dixon_poly->ctx);
        }
    }
    
    // Fill the coefficient matrix
    for (slong t = 0; t < dixon_poly->nterms; t++) {
        // Check degree bounds again
        int valid = 1;
        if (dixon_poly->terms[t].var_exp) {
            for (slong k = 0; k < nvars; k++) {
                if (dixon_poly->terms[t].var_exp[k] >= d0[k] || 
                    dixon_poly->terms[t].var_exp[nvars + k] >= d1[k]) {
                    valid = 0;
                    break;
                }
            }
        }
        
        if (!valid) continue;
        
        // Find row index (x-monomial)
        slong row = -1;
        for (slong i = 0; i < nx_monoms; i++) {
            int same = 1;
            for (slong k = 0; k < nvars; k++) {
                if (x_monoms[i].exp[k] != dixon_poly->terms[t].var_exp[k]) {
                    same = 0;
                    break;
                }
            }
            if (same) {
                row = i;
                break;
            }
        }
        
        // Find column index (~x-monomial)
        slong col = -1;
        for (slong j = 0; j < ndual_monoms; j++) {
            int same = 1;
            for (slong k = 0; k < nvars; k++) {
                if (dual_monoms[j].exp[k] != dixon_poly->terms[t].var_exp[nvars + k]) {
                    same = 0;
                    break;
                }
            }
            if (same) {
                col = j;
                break;
            }
        }
        
        // Add the parameter polynomial to the matrix entry
        if (row >= 0 && col >= 0) {
            fq_mvpoly_add_term(full_matrix[row][col], NULL, 
                              dixon_poly->terms[t].par_exp, dixon_poly->terms[t].coeff);
        }
    }
    
    // Find maximal rank submatrix
    slong *row_idx_array = NULL;
    slong *col_idx_array = NULL;
    slong num_rows, num_cols;
    
    if (npars == 0) {
        // Special case: no parameters, directly find rank
        fq_nmod_mat_t eval_mat;
        fq_nmod_mat_init(eval_mat, nx_monoms, ndual_monoms, dixon_poly->ctx);
        
        for (slong i = 0; i < nx_monoms; i++) {
            for (slong j = 0; j < ndual_monoms; j++) {
                if (full_matrix[i][j]->nterms > 0) {
                    fq_nmod_set(fq_nmod_mat_entry(eval_mat, i, j), 
                               full_matrix[i][j]->terms[0].coeff, dixon_poly->ctx);
                } else {
                    fq_nmod_zero(fq_nmod_mat_entry(eval_mat, i, j), dixon_poly->ctx);
                }
            }
        }
        
        slong rank = fq_nmod_mat_rank(eval_mat, dixon_poly->ctx);
        
        // Create identity selection for square case
        slong min_size = FLINT_MIN(nx_monoms, ndual_monoms);
        slong actual_size = FLINT_MIN(rank, min_size);
        
        row_idx_array = (slong*) flint_malloc(actual_size * sizeof(slong));
        col_idx_array = (slong*) flint_malloc(actual_size * sizeof(slong));
        
        for (slong i = 0; i < actual_size; i++) {
            row_idx_array[i] = i;
            col_idx_array[i] = i;
        }
        
        num_rows = actual_size;
        num_cols = actual_size;
        
        fq_nmod_mat_clear(eval_mat, dixon_poly->ctx);
    } else {
        // Call submatrix finding function
        find_fq_optimal_maximal_rank_submatrix(full_matrix, nx_monoms, ndual_monoms,
                                              &row_idx_array, &col_idx_array, 
                                              &num_rows, &num_cols,
                                              npars);
    }
    
    // Take square part if needed
    slong submat_rank = FLINT_MIN(num_rows, num_cols);
    
    if (submat_rank == 0) {
        printf("Warning: Matrix has rank 0\n");
        *matrix_size = 0;
        
        // Cleanup
        for (slong i = 0; i < nx_monoms; i++) {
            for (slong j = 0; j < ndual_monoms; j++) {
                fq_mvpoly_clear(full_matrix[i][j]);
                flint_free(full_matrix[i][j]);
            }
            flint_free(full_matrix[i]);
        }
        flint_free(full_matrix);
        
        if (row_idx_array) flint_free(row_idx_array);
        if (col_idx_array) flint_free(col_idx_array);
        
        for (slong i = 0; i < nx_monoms; i++) {
            flint_free(x_monoms[i].exp);
        }
        for (slong j = 0; j < ndual_monoms; j++) {
            flint_free(dual_monoms[j].exp);
        }
        flint_free(x_monoms);
        flint_free(dual_monoms);
        flint_free(d0);
        flint_free(d1);
        return;
    }
    
    printf("\nExtracted submatrix of size %ld x %ld\n", submat_rank, submat_rank);
    
    // Build the output submatrix
    *coeff_matrix = (fq_mvpoly_t**) flint_malloc(submat_rank * sizeof(fq_mvpoly_t*));
    for (slong i = 0; i < submat_rank; i++) {
        (*coeff_matrix)[i] = (fq_mvpoly_t*) flint_malloc(submat_rank * sizeof(fq_mvpoly_t));
        for (slong j = 0; j < submat_rank; j++) {
            fq_mvpoly_copy(&(*coeff_matrix)[i][j], full_matrix[row_idx_array[i]][col_idx_array[j]]);
        }
    }
    
    // Copy indices
    for (slong i = 0; i < submat_rank; i++) {
        row_indices[i] = row_idx_array[i];
        col_indices[i] = col_idx_array[i];
    }
    *matrix_size = submat_rank;
    
    // Cleanup
    for (slong i = 0; i < nx_monoms; i++) {
        for (slong j = 0; j < ndual_monoms; j++) {
            fq_mvpoly_clear(full_matrix[i][j]);
            flint_free(full_matrix[i][j]);
        }
        flint_free(full_matrix[i]);
    }
    flint_free(full_matrix);
    
    if (row_idx_array) flint_free(row_idx_array);
    if (col_idx_array) flint_free(col_idx_array);
    
    for (slong i = 0; i < nx_monoms; i++) {
        flint_free(x_monoms[i].exp);
    }
    for (slong j = 0; j < ndual_monoms; j++) {
        flint_free(dual_monoms[j].exp);
    }
    flint_free(x_monoms);
    flint_free(dual_monoms);
    flint_free(d0);
    flint_free(d1);
}

// Compute determinant of cancellation matrix
void compute_fq_cancel_matrix_det(fq_mvpoly_t *result, fq_mvpoly_t **modified_M_mvpoly,
                                 slong nvars, slong npars, det_method_t method) {
    printf("Computing cancellation matrix determinant using recursive expansion\n");
    
    clock_t start = clock();
    compute_fq_det_recursive(result, modified_M_mvpoly, nvars + 1);
    clock_t end = clock();
    
    printf("End (%.3f seconds)\n", (double)(end - start) / CLOCKS_PER_SEC);
}

// ============ Compute determinant of coefficient matrix ============

// ============ Main Dixon resultant function ============
void fq_dixon_resultant(fq_mvpoly_t *result, fq_mvpoly_t *polys, 
                       slong nvars, slong npars) {
    printf("\n=== Dixon Resultant over F_{p^d} ===\n");
    
    // Step 1: Build cancellation matrix
    printf("\nStep 1: Build Cancellation Matrix\n");
    fq_mvpoly_t **M_mvpoly;
    build_fq_cancellation_matrix_mvpoly(&M_mvpoly, polys, nvars, npars);
    
    // 显示原始矩阵的分析
    analyze_fq_matrix_mvpoly(M_mvpoly, nvars + 1, nvars + 1, "Original Cancellation");
    
    // Step 2: Perform row operations
    printf("\nStep 2: Perform Matrix Row Operations\n");
    fq_mvpoly_t **modified_M_mvpoly;
    perform_fq_matrix_row_operations_mvpoly(&modified_M_mvpoly, &M_mvpoly, nvars, npars);
    
    // **新增：显示除法后的矩阵**
    printf("\n=== Matrix After Row Operations ===\n");
    analyze_fq_matrix_mvpoly(modified_M_mvpoly, nvars + 1, nvars + 1, "After Row Operations");
    
    // 如果矩阵不太大，显示详细内容
    if (nvars <= 3) {
        print_fq_matrix_mvpoly(modified_M_mvpoly, nvars + 1, nvars + 1, 
                               "After Row Operations (Detailed)", 1);
    } else {
        print_fq_matrix_mvpoly(modified_M_mvpoly, nvars + 1, nvars + 1, 
                               "After Row Operations (Summary)", 0);
    }
    
    // 比较原始矩阵和修改后矩阵的复杂度
    printf("\n=== Complexity Comparison ===\n");
    
    // 计算原始矩阵的总项数
    slong orig_total_terms = 0;
    for (slong i = 0; i <= nvars; i++) {
        for (slong j = 0; j <= nvars; j++) {
            orig_total_terms += M_mvpoly[i][j].nterms;
        }
    }
    
    // 计算修改后矩阵的总项数
    slong modified_total_terms = 0;
    for (slong i = 0; i <= nvars; i++) {
        for (slong j = 0; j <= nvars; j++) {
            modified_total_terms += modified_M_mvpoly[i][j].nterms;
        }
    }
    
    printf("Original matrix total terms: %ld\n", orig_total_terms);
    printf("Modified matrix total terms: %ld\n", modified_total_terms);
    printf("Complexity ratio: %.2f\n", (double)modified_total_terms / orig_total_terms);
    
    // Step 3: Compute determinant of modified matrix
    printf("\nStep 3: Compute determinant of modified matrix\n");
    fq_mvpoly_t d_poly;
    compute_fq_cancel_matrix_det(&d_poly, modified_M_mvpoly, nvars, npars, DET_METHOD_MBOT);
    
    printf("Dixon polynomial has %ld terms\n", d_poly.nterms);
    
    if (d_poly.nterms <= 100) {
        fq_mvpoly_print_expanded(&d_poly, "Dixon", 1);
    } else {
        printf("Dixon polynomial too large to display (%ld terms)\n", d_poly.nterms);
        
        // 显示部分项作为示例
        printf("First 5 terms of Dixon polynomial:\n");
        for (slong i = 0; i < FLINT_MIN(5, d_poly.nterms); i++) {
            printf("  Term %ld: ", i);
            fq_nmod_print_pretty(d_poly.terms[i].coeff, d_poly.ctx);
            
            // 显示变量指数
            if (d_poly.terms[i].var_exp) {
                for (slong j = 0; j < d_poly.nvars; j++) {
                    if (d_poly.terms[i].var_exp[j] > 0) {
                        if (j < d_poly.nvars / 2) {
                            printf(" * x%ld^%ld", j, d_poly.terms[i].var_exp[j]);
                        } else {
                            printf(" * ~x%ld^%ld", j - d_poly.nvars / 2, d_poly.terms[i].var_exp[j]);
                        }
                    }
                }
            }
            
            // 显示参数指数
            if (d_poly.terms[i].par_exp && d_poly.npars > 0) {
                for (slong j = 0; j < d_poly.npars; j++) {
                    if (d_poly.terms[i].par_exp[j] > 0) {
                        printf(" * p%ld^%ld", j, d_poly.terms[i].par_exp[j]);
                    }
                }
            }
            printf("\n");
        }
    }
    
    // Step 4: Extract coefficient matrix
    printf("\nStep 4: Extract coefficient matrix\n");
    
    fq_mvpoly_t **coeff_matrix;
    slong *row_indices = (slong*) flint_malloc(d_poly.nterms * sizeof(slong));
    slong *col_indices = (slong*) flint_malloc(d_poly.nterms * sizeof(slong));
    slong matrix_size;
    
    extract_fq_coefficient_matrix_from_dixon(&coeff_matrix, row_indices, col_indices,
                                            &matrix_size, &d_poly, nvars, npars);
    
    if (matrix_size > 0) {
        printf("\nStep 5: Compute determinant of coefficient matrix (fq_nmod)\n");
        
        slong res_deg_bound = compute_fq_dixon_resultant_degree_bound(polys, nvars+1, nvars, npars);
        printf("Resultant degree bound: %ld\n", res_deg_bound);
        
        det_method_t coeff_method = DET_METHOD_INTERPOLATION;
        char *method_env = getenv("DIXON_DET_METHOD");
        if (method_env) {
            coeff_method = atoi(method_env);
        }
        
        printf("Using coefficient matrix determinant method: %d\n", coeff_method);
        
        compute_fq_coefficient_matrix_det(result, coeff_matrix, matrix_size,
                                         npars, polys[0].ctx, coeff_method, res_deg_bound);
        
        printf("\nResultant polynomial has %ld terms\n", result->nterms);
        if (result->nterms <= 20) {
            fq_mvpoly_print(result, "Final Resultant");
        } else {
            printf("Final resultant too large to display (%ld terms)\n", result->nterms);
            
            // 显示次数信息
            slong max_par_deg = 0;
            for (slong i = 0; i < result->nterms; i++) {
                if (result->terms[i].par_exp && result->npars > 0) {
                    for (slong j = 0; j < result->npars; j++) {
                        if (result->terms[i].par_exp[j] > max_par_deg) {
                            max_par_deg = result->terms[i].par_exp[j];
                        }
                    }
                }
            }
            printf("Maximum parameter degree in resultant: %ld\n", max_par_deg);
        }
        
        // Cleanup coefficient matrix
        for (slong i = 0; i < matrix_size; i++) {
            for (slong j = 0; j < matrix_size; j++) {
                fq_mvpoly_clear(&coeff_matrix[i][j]);
            }
            flint_free(coeff_matrix[i]);
        }
        flint_free(coeff_matrix);
    } else {
        fq_mvpoly_init(result, 0, npars, polys[0].ctx);
        printf("Warning: Empty coefficient matrix, resultant is 0\n");
    }
    
    // Cleanup
    flint_free(row_indices);
    flint_free(col_indices);
    
    for (slong i = 0; i <= nvars; i++) {
        for (slong j = 0; j <= nvars; j++) {
            fq_mvpoly_clear(&M_mvpoly[i][j]);
            fq_mvpoly_clear(&modified_M_mvpoly[i][j]);
        }
        flint_free(M_mvpoly[i]);
        flint_free(modified_M_mvpoly[i]);
    }
    flint_free(M_mvpoly);
    flint_free(modified_M_mvpoly);
    fq_mvpoly_clear(&d_poly);
    
    printf("\n=== Dixon Resultant Computation Complete ===\n");
}



//#include "gf2n_flint_comparison.h"
//#include "test_gf2n_unified.h"
//#include "gf2n_extended_tests.h"
//#include "gf2n_fft_mul.h"
//#include "test_gf2n_fields.h" 
//#include "test_unified_interface.h" 
//#include "fq_test_performance_fixed.h"
//#include "test_fft_debug.h"
// Main function

int main() {
    printf("Dixon Resultant Computation over Finite Extension Fields\n");
    printf("========================================================\n");
    
    //int result = test_fq_dixon_with_analysis();
    //quick_extension_test();
        
//fq_test_performance_fixed();
    //test_gf232_fft();
    //test_gf2n_fields();
    //gf2128_fft_mul();
    omp_set_num_threads(64);
    //test_gf2128_fft_specific();
    //gf2n_tests();
    //gf2n_extended_tests();
    //test_gf2n_unified();
    int testdixon = test_fq_string_interface_enhanced();
    //test_gf2n_flint_comparison();
    // 清理
    //cleanup_gf28_tables();

    //test_unified_interface_in_main();

    //test_fq_nmod_poly_mat_det();
    //result = main_string_parser();
    return 0;
}

