#ifndef FQ_MVPOLY_H
#define FQ_MVPOLY_H

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

void fq_mvpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t mpoly, const fq_mvpoly_t *poly, 
                               fq_nmod_mpoly_ctx_t mpoly_ctx);
void fq_nmod_mpoly_to_fq_mvpoly(fq_mvpoly_t *poly, const fq_nmod_mpoly_t mpoly,
                                slong nvars, slong npars, 
                                fq_nmod_mpoly_ctx_t mpoly_ctx, const fq_nmod_ctx_t ctx);

// ============ Basic fq_mvpoly operations ============

void fq_mvpoly_init(fq_mvpoly_t *p, slong nvars, slong npars, const fq_nmod_ctx_t ctx) {
    p->nvars = nvars;
    p->npars = npars;
    p->nterms = 0;
    p->alloc = 16;
    p->terms = (fq_monomial_t*) flint_malloc(p->alloc * sizeof(fq_monomial_t));
    
    memset(p->terms, 0, p->alloc * sizeof(fq_monomial_t));
    
    p->ctx = ctx; 
}


void fq_mvpoly_clear(fq_mvpoly_t *p) {
    if (p->terms) {
        for (slong i = 0; i < p->nterms; i++) {
            fq_nmod_clear(p->terms[i].coeff, p->ctx);
            if (p->terms[i].var_exp) flint_free(p->terms[i].var_exp);
            if (p->terms[i].par_exp) flint_free(p->terms[i].par_exp);
        }

        flint_free(p->terms);
        p->terms = NULL;
    }
    p->nterms = 0;
    p->alloc = 0;
    p->nvars = 0;
    p->npars = 0;
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
    if (dest == src) {
        return;
    }
    
    fq_mvpoly_t temp;
    fq_mvpoly_init(&temp, src->nvars, src->npars, src->ctx);
    
    for (slong i = 0; i < src->nterms; i++) {
        fq_mvpoly_add_term_fast(&temp, src->terms[i].var_exp, src->terms[i].par_exp, src->terms[i].coeff);
    }
    
    *dest = temp;
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
    fq_mvpoly_t temp;
    fq_mvpoly_init(&temp, a->nvars, a->npars, a->ctx);
    
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
            
            fq_mvpoly_add_term(&temp, var_exp, par_exp, coeff);
            
            fq_nmod_clear(coeff, a->ctx);
            if (var_exp) flint_free(var_exp);
            if (par_exp) flint_free(par_exp);
        }
    }
    
    // 处理 result 可能已初始化的情况
    if (result == a || result == b) {
        fq_mvpoly_clear(result);
    }
    *result = temp;
}

// Compute polynomial to a power
void fq_mvpoly_pow(fq_mvpoly_t *result, const fq_mvpoly_t *base, slong power) {
    fq_mvpoly_t temp;
    fq_mvpoly_init(&temp, base->nvars, base->npars, base->ctx);
    
    fq_nmod_t one;
    fq_nmod_init(one, base->ctx);
    fq_nmod_one(one, base->ctx);
    fq_mvpoly_add_term(&temp, NULL, NULL, one); // temp = 1
    fq_nmod_clear(one, base->ctx);
    
    if (power == 0) {
        *result = temp;
        return;
    }
    
    fq_mvpoly_t base_copy, temp2;
    fq_mvpoly_copy(&base_copy, base);
    
    while (power > 0) {
        if (power & 1) {
            fq_mvpoly_mul(&temp2, &temp, &base_copy);
            fq_mvpoly_clear(&temp);
            temp = temp2;  // 结构体赋值
        }
        
        if (power > 1) {
            fq_mvpoly_mul(&temp2, &base_copy, &base_copy);
            fq_mvpoly_clear(&base_copy);
            base_copy = temp2;  // 结构体赋值
        }
        
        power >>= 1;
    }
    
    fq_mvpoly_clear(&base_copy);
    
    // 处理 result 可能已初始化的情况
    if (result == base) {
        fq_mvpoly_clear(result);
    }
    *result = temp;
}

// Multiply polynomial by scalar
void fq_mvpoly_scalar_mul(fq_mvpoly_t *result, const fq_mvpoly_t *p, const fq_nmod_t scalar) {
    fq_mvpoly_t temp;
    fq_mvpoly_init(&temp, p->nvars, p->npars, p->ctx);
    
    for (slong i = 0; i < p->nterms; i++) {
        fq_nmod_t new_coeff;
        fq_nmod_init(new_coeff, p->ctx);
        fq_nmod_mul(new_coeff, p->terms[i].coeff, scalar, p->ctx);
        fq_mvpoly_add_term(&temp, p->terms[i].var_exp, p->terms[i].par_exp, new_coeff);
        fq_nmod_clear(new_coeff, p->ctx);
    }
    
    // 处理 result 可能已初始化的情况
    if (result == p) {
        fq_mvpoly_clear(result);
    }
    *result = temp;
}

// Add polynomial b to polynomial a
void fq_mvpoly_add(fq_mvpoly_t *result, const fq_mvpoly_t *a, const fq_mvpoly_t *b) {
    // 检查 result 是否是输入参数之一
    if (result == a || result == b) {
        // 需要使用临时变量
        fq_mvpoly_t temp;
        fq_mvpoly_init(&temp, a->nvars, a->npars, a->ctx);
        
        // 添加所有项从 a
        for (slong i = 0; i < a->nterms; i++) {
            fq_mvpoly_add_term(&temp, a->terms[i].var_exp, a->terms[i].par_exp, a->terms[i].coeff);
        }
        
        // 添加所有项从 b
        for (slong i = 0; i < b->nterms; i++) {
            fq_mvpoly_add_term(&temp, b->terms[i].var_exp, b->terms[i].par_exp, b->terms[i].coeff);
        }
        
        // 清理旧的 result 并复制新内容
        fq_mvpoly_clear(result);
        fq_mvpoly_init(result, a->nvars, a->npars, a->ctx);
        
        // 复制 temp 到 result
        for (slong i = 0; i < temp.nterms; i++) {
            fq_mvpoly_add_term(result, temp.terms[i].var_exp, temp.terms[i].par_exp, temp.terms[i].coeff);
        }
        
        fq_mvpoly_clear(&temp);
    } else {
        // result 不是输入参数，可以直接操作
        
        // 如果 result 已经包含数据，先清理
        if (result->terms != NULL) {
            fq_mvpoly_clear(result);
        }
        
        // 初始化 result
        fq_mvpoly_init(result, a->nvars, a->npars, a->ctx);
        
        // 添加所有项从 a
        for (slong i = 0; i < a->nterms; i++) {
            fq_mvpoly_add_term(result, a->terms[i].var_exp, a->terms[i].par_exp, a->terms[i].coeff);
        }
        
        // 添加所有项从 b
        for (slong i = 0; i < b->nterms; i++) {
            fq_mvpoly_add_term(result, b->terms[i].var_exp, b->terms[i].par_exp, b->terms[i].coeff);
        }
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


// Improved division by linear factor (x_i - ~x_i) - direct method
void divide_by_fq_linear_factor_improved(fq_mvpoly_t *quotient, const fq_mvpoly_t *dividend,
                                        slong var_idx, slong nvars, slong npars) {
    fq_mvpoly_t temp;
    fq_mvpoly_init(&temp, nvars, npars, dividend->ctx);
    
    if (dividend->nterms == 0) {
        *quotient = temp;
        return;
    }
    
    slong actual_nvars = nvars / 2;  // Since we have dual variables
    
    // Process each term in the dividend
    for (slong i = 0; i < dividend->nterms; i++) {
        slong deg_x = dividend->terms[i].var_exp ? dividend->terms[i].var_exp[var_idx] : 0;
        slong deg_dual = dividend->terms[i].var_exp ? dividend->terms[i].var_exp[var_idx + actual_nvars] : 0;
        
        if (deg_x == deg_dual) {
            // This term is divisible by (x_i - ~x_i), the result has the same exponents
            fq_mvpoly_add_term(&temp, dividend->terms[i].var_exp, dividend->terms[i].par_exp, 
                              dividend->terms[i].coeff);
        } else if (deg_x > deg_dual) {
            // Reduce x degree by 1
            slong *new_exp = (slong*) flint_calloc(nvars, sizeof(slong));
            if (dividend->terms[i].var_exp) {
                memcpy(new_exp, dividend->terms[i].var_exp, nvars * sizeof(slong));
            }
            new_exp[var_idx]--;
            
            fq_mvpoly_add_term(&temp, new_exp, dividend->terms[i].par_exp, 
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
            
            fq_mvpoly_add_term(&temp, new_exp, dividend->terms[i].par_exp, neg_coeff);
            
            fq_nmod_clear(neg_coeff, dividend->ctx);
            flint_free(new_exp);
        }
    }
    
    // 处理 quotient 可能已初始化的情况
    //if (quotient == dividend) {
        //fq_mvpoly_clear(quotient);
    //}
    if (quotient->terms != NULL) {
        // 已初始化，需要先清理
        fq_mvpoly_clear(quotient);
    }
    *quotient = temp;
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

#endif