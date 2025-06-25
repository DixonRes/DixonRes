/*
 * Modified Dixon Resultant Implementation with Row Operations
 * 
 * This implementation modifies the Dixon resultant algorithm to:
 * - First compute row differences and divisions
 * - Then compute the determinant of the modified matrix
 * - Uses Kronecker substitution for univariate polynomial division
 * 
 * gcc -o dixon_modified dixon_resultant.c -lflint -lmpfr -lgmp -lpthread -L/home/suohaohai02/mylinks -lflint -lstdc++ -lpml2
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/nmod.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/nmod_mat.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/nmod_types.h>
#include <flint/profiler.h>
#include <nmod_poly_mat_utils.h>
#include <nmod_poly_mat_extra.h>

// Structure for multivariate monomials with parameters
typedef struct {
    slong *var_exp;      // exponent vector for variables (including duals)
    slong *par_exp;      // exponent vector for parameters
    mp_limb_t coeff;     // coefficient
} monomial_t;

// Structure for multivariate polynomial with parameters
typedef struct {
    slong nvars;         // number of variables (not including duals)
    slong npars;         // number of parameters
    slong nterms;        // number of terms
    slong alloc;         // allocated space
    monomial_t *terms;   // array of terms
    nmod_t mod;          // modulus
} mvpoly_t;

// Forward declarations
void mvpoly_init(mvpoly_t *p, slong nvars, slong npars, mp_limb_t mod);
void mvpoly_clear(mvpoly_t *p);
void mvpoly_copy(mvpoly_t *dest, const mvpoly_t *src);
void mvpoly_add_term(mvpoly_t *p, const slong *var_exp, const slong *par_exp, mp_limb_t coeff);
void mvpoly_print(const mvpoly_t *p, const char *name);
void mvpoly_print_expanded(const mvpoly_t *p, const char *name, int use_dual);

// Initialize multivariate polynomial with parameters
void mvpoly_init(mvpoly_t *p, slong nvars, slong npars, mp_limb_t mod)
{
    p->nvars = nvars;
    p->npars = npars;
    p->nterms = 0;
    p->alloc = 16;
    p->terms = (monomial_t*) flint_malloc(p->alloc * sizeof(monomial_t));
    nmod_init(&p->mod, mod);
}

void mvpoly_clear(mvpoly_t *p)
{
    for (slong i = 0; i < p->nterms; i++)
    {
        if (p->terms[i].var_exp) flint_free(p->terms[i].var_exp);
        if (p->terms[i].par_exp) flint_free(p->terms[i].par_exp);
    }
    if (p->terms) flint_free(p->terms);
}

// Copy polynomial
void mvpoly_copy(mvpoly_t *dest, const mvpoly_t *src)
{
    mvpoly_init(dest, src->nvars, src->npars, src->mod.n);
    
    for (slong i = 0; i < src->nterms; i++)
    {
        mvpoly_add_term(dest, src->terms[i].var_exp, src->terms[i].par_exp, src->terms[i].coeff);
    }
}


// Print multivariate polynomial (basic version)
void mvpoly_print(const mvpoly_t *p, const char *name)
{
    mvpoly_print_expanded(p, name, 0);
}

// Print multivariate polynomial with expanded variable names
void mvpoly_print_expanded(const mvpoly_t *p, const char *name, int use_dual)
{
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
        if (p->terms[i].coeff == p->mod.n - 1) {
            printf("-");
            // Check if we need to print 1
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
            if (!has_vars) printf("1");
        } else {
            printf("%lu", (unsigned long)p->terms[i].coeff);
        }
        
        // Print variables
        if (use_dual) {
            // For expanded polynomials with dual variables
            slong actual_nvars = p->nvars / 2;
            
            // First print regular variables
            for (slong j = 0; j < actual_nvars; j++) {
                if (p->terms[i].var_exp && p->terms[i].var_exp[j] > 0) {
                    printf("*");
                    if (j < 6) printf("%c", var_names[j]);
                    else printf("x_%ld", j);
                    
                    if (p->terms[i].var_exp[j] > 1) {
                        printf("^%ld", p->terms[i].var_exp[j]);
                    }
                }
            }
            
            // Then print dual variables with tilde
            for (slong j = actual_nvars; j < p->nvars; j++) {
                if (p->terms[i].var_exp && p->terms[i].var_exp[j] > 0) {
                    printf("*");
                    slong orig_idx = j - actual_nvars;
                    if (orig_idx < 6) printf("~%c", var_names[orig_idx]);
                    else printf("~x_%ld", orig_idx);
                    
                    if (p->terms[i].var_exp[j] > 1) {
                        printf("^%ld", p->terms[i].var_exp[j]);
                    }
                }
            }
        } else {
            // Normal variables
            for (slong j = 0; j < p->nvars; j++) {
                if (p->terms[i].var_exp && p->terms[i].var_exp[j] > 0) {
                    printf("*");
                    if (j < 6) printf("%c", var_names[j]);
                    else printf("x_%ld", j);
                    
                    if (p->terms[i].var_exp[j] > 1) {
                        printf("^%ld", p->terms[i].var_exp[j]);
                    }
                }
            }
        }
        
        // Print parameters
        for (slong j = 0; j < p->npars; j++) {
            if (p->terms[i].par_exp && p->terms[i].par_exp[j] > 0) {
                printf("*");
                if (j < 4) printf("%c", par_names[j]);
                else printf("p_%ld", j);
                
                if (p->terms[i].par_exp[j] > 1) {
                    printf("^%ld", p->terms[i].par_exp[j]);
                }
            }
        }
    }
    printf("\n");
}

// Add a term to multivariate polynomial
void mvpoly_add_term(mvpoly_t *p, const slong *var_exp, const slong *par_exp, mp_limb_t coeff)
{
    if (coeff == 0) return;
    
    // Check if term already exists
    for (slong i = 0; i < p->nterms; i++)
    {
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
        
        if (same)
        {
            p->terms[i].coeff = nmod_add(p->terms[i].coeff, coeff, p->mod);
            if (p->terms[i].coeff == 0)
            {
                // Remove zero term
                if (p->terms[i].var_exp) flint_free(p->terms[i].var_exp);
                if (p->terms[i].par_exp) flint_free(p->terms[i].par_exp);
                for (slong j = i; j < p->nterms - 1; j++)
                {
                    p->terms[j] = p->terms[j + 1];
                }
                p->nterms--;
            }
            return;
        }
    }
    
    // Add new term
    if (p->nterms >= p->alloc)
    {
        p->alloc *= 2;
        p->terms = (monomial_t*) flint_realloc(p->terms, p->alloc * sizeof(monomial_t));
    }
    
    // Allocate and copy exponents
    if (p->nvars > 0) {
        p->terms[p->nterms].var_exp = (slong*) flint_calloc(p->nvars, sizeof(slong));
        if (var_exp) {
            for (slong i = 0; i < p->nvars; i++) {
                p->terms[p->nterms].var_exp[i] = var_exp[i];
            }
        }
    } else {
        p->terms[p->nterms].var_exp = NULL;
    }
    
    if (p->npars > 0) {
        p->terms[p->nterms].par_exp = (slong*) flint_calloc(p->npars, sizeof(slong));
        if (par_exp) {
            for (slong i = 0; i < p->npars; i++) {
                p->terms[p->nterms].par_exp[i] = par_exp[i];
            }
        }
    } else {
        p->terms[p->nterms].par_exp = NULL;
    }
    
    p->terms[p->nterms].coeff = coeff;
    p->nterms++;
}


// Kronecker substitution helpers
slong exp_to_kronecker_index(const slong *exp, const slong *degs, slong n)
{
    slong index = 0;
    slong stride = 1;
    
    for (slong i = 0; i < n; i++)
    {
        index += exp[i] * stride;
        stride *= (degs[i] + 1);
    }
    
    return index;
}

void kronecker_index_to_exp(slong index, slong *exp, const slong *degs, slong n)
{
    for (slong i = 0; i < n; i++)
    {
        exp[i] = index % (degs[i] + 1);
        index /= (degs[i] + 1);
    }
}

// Convert mvpoly to univariate via Kronecker (for both variables and parameters)
void mvpoly_to_kronecker_full(nmod_poly_t out, const mvpoly_t *p, 
                             const slong *var_degs, const slong *par_degs)
{
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
    
    nmod_poly_fit_length(out, max_deg + 1);
    nmod_poly_zero(out);
    
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
        mp_limb_t old_coeff = nmod_poly_get_coeff_ui(out, idx);
        mp_limb_t new_coeff = nmod_add(old_coeff, p->terms[i].coeff, p->mod);
        nmod_poly_set_coeff_ui(out, idx, new_coeff);
    }
    
    _nmod_poly_normalise(out);
    
    flint_free(combined_exp);
    flint_free(all_degs);
}

void kronecker_to_mvpoly_full(mvpoly_t *out, const nmod_poly_t in, 
                                         const slong *var_degs, slong nvars,
                                         const slong *par_degs, slong npars, mp_limb_t mod)
{
    mvpoly_init(out, nvars, npars, mod);
    
    slong deg = nmod_poly_degree(in);
    if (deg < 0) return;
    
    // 预先计算非零项数
    slong nterms = 0;
    for (slong i = 0; i <= deg; i++) {
        if (nmod_poly_get_coeff_ui(in, i) != 0) {
            nterms++;
        }
    }
    
    if (nterms == 0) return;
    
    // 预分配足够的空间
    if (out->alloc < nterms) {
        out->alloc = nterms;
        out->terms = (monomial_t*) flint_realloc(out->terms, out->alloc * sizeof(monomial_t));
    }
    
    // 准备工作数组
    slong total_vars = nvars + npars;
    slong *all_degs = (slong*) flint_malloc(total_vars * sizeof(slong));
    
    for (slong i = 0; i < nvars; i++) {
        all_degs[i] = var_degs[i];
    }
    for (slong i = 0; i < npars; i++) {
        all_degs[nvars + i] = par_degs[i];
    }
    
    slong *combined_exp = (slong*) flint_calloc(total_vars, sizeof(slong));
    
    // 填充项 - 使用标准的内存分配方式以保持兼容性
    slong term_idx = 0;
    for (slong i = 0; i <= deg; i++) {
        mp_limb_t coeff = nmod_poly_get_coeff_ui(in, i);
        if (coeff != 0) {
            // 解码指数
            slong idx = i;
            for (slong j = 0; j < total_vars; j++) {
                combined_exp[j] = idx % (all_degs[j] + 1);
                idx /= (all_degs[j] + 1);
            }
            
            // 使用标准方式分配内存（与原始代码兼容）
            out->terms[term_idx].coeff = coeff;
            
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
    
    flint_free(combined_exp);
    flint_free(all_degs);
}

// Build cancellation matrix entry
void build_cancellation_entry(nmod_poly_t entry, const mvpoly_t *poly, 
                             slong row, slong nvars_orig, const slong *var_degs, 
                             const slong *par_degs, mp_limb_t prime)
{
    // Create polynomial with dual variables
    mvpoly_t expanded;
    mvpoly_init(&expanded, 2 * nvars_orig, poly->npars, prime);
    
    // Substitute variables according to row index
    for (slong t = 0; t < poly->nterms; t++) {
        slong *new_var_exp = (slong*) flint_calloc(2 * nvars_orig, sizeof(slong));
        
        for (slong k = 0; k < nvars_orig; k++) {
            slong orig_exp = poly->terms[t].var_exp ? poly->terms[t].var_exp[k] : 0;
            
            if (k < row) {
                // Use dual variable ~x_k
                new_var_exp[nvars_orig + k] = orig_exp;
            } else {
                // Use original variable x_k
                new_var_exp[k] = orig_exp;
            }
        }
        
        mvpoly_add_term(&expanded, new_var_exp, poly->terms[t].par_exp, poly->terms[t].coeff);
        flint_free(new_var_exp);
    }
    
    // Convert to univariate via Kronecker
    mvpoly_to_kronecker_full(entry, &expanded, var_degs, par_degs);
    
    mvpoly_clear(&expanded);
}

// Helper function to compute tight degree bounds for determinant
void compute_determinant_degree_bounds(slong *det_var_degs, slong *det_par_degs,
                                      mvpoly_t *polys, slong nvars, slong npars,
                                      slong matrix_size)
{
    // Initialize output degree bounds
    for (slong i = 0; i < 2 * nvars; i++) {
        det_var_degs[i] = 0;
    }
    for (slong i = 0; i < npars; i++) {
        det_par_degs[i] = 0;
    }
    
    // For each row of the cancellation matrix
    for (slong row = 0; row < matrix_size; row++) {
        // Track maximum degree in this row for each variable/parameter
        slong *row_var_max = (slong*) flint_calloc(2 * nvars, sizeof(slong));
        slong *row_par_max = (slong*) flint_calloc(npars, sizeof(slong));
        
        // For each column (polynomial) in this row
        for (slong col = 0; col < matrix_size; col++) {
            mvpoly_t *poly = &polys[col];
            
            // For each term in the polynomial
            for (slong t = 0; t < poly->nterms; t++) {
                // Calculate degrees after substitution for this row
                // Variables x_k with k < row become ~x_k (dual variables)
                
                for (slong k = 0; k < nvars; k++) {
                    slong orig_deg = poly->terms[t].var_exp ? poly->terms[t].var_exp[k] : 0;
                    
                    if (k < row) {
                        // This variable becomes dual
                        if (orig_deg > row_var_max[nvars + k]) {
                            row_var_max[nvars + k] = orig_deg;
                        }
                    } else {
                        // This variable stays original
                        if (orig_deg > row_var_max[k]) {
                            row_var_max[k] = orig_deg;
                        }
                    }
                }
                
                // Parameters don't change
                for (slong k = 0; k < npars; k++) {
                    slong deg = poly->terms[t].par_exp ? poly->terms[t].par_exp[k] : 0;
                    if (deg > row_par_max[k]) {
                        row_par_max[k] = deg;
                    }
                }
            }
        }
        
        // Add this row's maximum degrees to the total
        for (slong i = 0; i < 2 * nvars; i++) {
            det_var_degs[i] += row_var_max[i];
        }
        for (slong i = 0; i < npars; i++) {
            det_par_degs[i] += row_par_max[i];
        }
        
        flint_free(row_var_max);
        flint_free(row_par_max);
    }
    
    printf("Determinant degree bounds:\n");
    printf("  Variables: ");
    for (slong i = 0; i < nvars; i++) {
        if (det_var_degs[i] > 0) printf("x%ld:%ld ", i, det_var_degs[i]);
    }
    for (slong i = 0; i < nvars; i++) {
        if (det_var_degs[nvars + i] > 0) printf("~x%ld:%ld ", i, det_var_degs[nvars + i]);
    }
    printf("\n  Parameters: ");
    for (slong i = 0; i < npars; i++) {
        if (det_par_degs[i] > 0) printf("p%ld:%ld ", i, det_par_degs[i]);
    }
    printf("\n");
}

// Build cancellation matrix with degree bounds
void build_cancellation_matrix_with_bounds(nmod_poly_mat_t M, mvpoly_t *polys, 
                                          slong nvars, slong npars, mp_limb_t prime,
                                          slong *det_var_degs, slong *det_par_degs)
{
    slong n = nvars + 1;
    nmod_poly_mat_init(M, n, n, prime);
    
    printf("\nBuilding %ld x %ld cancellation matrix\n", n, n);
    
    // Compute tight degree bounds for the determinant
    compute_determinant_degree_bounds(det_var_degs, det_par_degs, polys, nvars, npars, n);
    
    // Build matrix entries using the computed degree bounds
    for (slong i = 0; i < n; i++) {
        for (slong j = 0; j < n; j++) {
            printf("M[%ld,%ld]: poly[%ld] with ", i, j, j);
            for (slong k = 0; k < nvars; k++) {
                if (k < i) printf("x%ld->~x%ld ", k, k);
                else printf("x%ld->x%ld ", k, k);
            }
            
            build_cancellation_entry(nmod_poly_mat_entry(M, i, j), &polys[j], 
                                   i, nvars, det_var_degs, det_par_degs, prime);
            
            printf("degree=%ld\n", nmod_poly_degree(nmod_poly_mat_entry(M, i, j)));
        }
    }
}

// NEW FUNCTION: Divide univariate polynomial by linear factor using Kronecker substitution
// Divides poly by (x_var - x_dual) where x_var and x_dual are represented in Kronecker form
void divide_by_linear_factor_kronecker(nmod_poly_t quotient, const nmod_poly_t dividend,
                                      slong var_idx, const slong *var_degs, slong nvars,
                                      const slong *par_degs, slong npars, mp_limb_t prime)
{
    printf("  Dividing by (x%ld - ~x%ld) using Kronecker division\n", var_idx, var_idx);
    
    // Convert back to multivariate to understand structure
    mvpoly_t dividend_mv;
    kronecker_to_mvpoly_full(&dividend_mv, dividend, var_degs, nvars, par_degs, npars, prime);
    
    printf("    Dividend has %ld terms\n", dividend_mv.nterms);
    
    // Create quotient polynomial
    mvpoly_t quotient_mv;
    mvpoly_init(&quotient_mv, nvars, npars, prime);
    
    // Perform symbolic division
    // We need to divide by (x_var - x_dual)
    // This is a multivariate polynomial division
    
    mvpoly_t work;
    mvpoly_copy(&work, &dividend_mv);
    
    int iterations = 0;
    while (1) {
        iterations++;
        if (iterations > 100000) {
            printf("    ERROR: Maximum iterations reached\n");
            break;
        }
        
        // Find a term where deg(x_i) != deg(~x_i)
        slong found_idx = -1;
        slong actual_nvars = nvars / 2;  // Since we have dual variables
        
        for (slong i = 0; i < work.nterms; i++) {
            slong deg_x = work.terms[i].var_exp ? work.terms[i].var_exp[var_idx] : 0;
            slong deg_dual = work.terms[i].var_exp ? work.terms[i].var_exp[var_idx + actual_nvars] : 0;
            
            if (deg_x != deg_dual) {
                found_idx = i;
                break;
            }
        }
        
        if (found_idx == -1) break;  // No more terms to process
        
        // Get the term to process
        slong deg_x = work.terms[found_idx].var_exp[var_idx];
        slong deg_dual = work.terms[found_idx].var_exp[var_idx + actual_nvars];
        mp_limb_t coeff = work.terms[found_idx].coeff;
        
        // Create exponent vectors for the operations
        slong *exp_orig = (slong*) flint_calloc(nvars, sizeof(slong));
        slong *exp_quot = (slong*) flint_calloc(nvars, sizeof(slong));
        slong *par_exp = (slong*) flint_calloc(npars, sizeof(slong));
        
        // Copy exponents
        for (slong j = 0; j < nvars; j++) {
            exp_orig[j] = work.terms[found_idx].var_exp[j];
            exp_quot[j] = work.terms[found_idx].var_exp[j];
        }
        if (npars > 0 && work.terms[found_idx].par_exp) {
            for (slong j = 0; j < npars; j++) {
                par_exp[j] = work.terms[found_idx].par_exp[j];
            }
        }
        
        if (deg_x > deg_dual) {
            // Quotient term: reduce x degree by 1
            exp_quot[var_idx]--;
            mvpoly_add_term(&quotient_mv, exp_quot, par_exp, coeff);
            
            // Update work: subtract (x_i - ~x_i) * quotient_term
            mvpoly_add_term(&work, exp_orig, par_exp, nmod_neg(coeff, work.mod));
            exp_quot[var_idx + actual_nvars]++;
            mvpoly_add_term(&work, exp_quot, par_exp, coeff);
            
        } else { // deg_dual > deg_x
            // Quotient term: reduce ~x degree by 1, with negative coefficient
            exp_quot[var_idx + actual_nvars]--;
            mvpoly_add_term(&quotient_mv, exp_quot, par_exp, nmod_neg(coeff, work.mod));
            
            // Update work
            mvpoly_add_term(&work, exp_orig, par_exp, nmod_neg(coeff, work.mod));
            exp_quot[var_idx]++;
            mvpoly_add_term(&work, exp_quot, par_exp, coeff);
        }
        
        flint_free(exp_orig);
        flint_free(exp_quot);
        flint_free(par_exp);
    }
    
    // Add remaining terms to quotient
    for (slong i = 0; i < work.nterms; i++) {
        mvpoly_add_term(&quotient_mv, work.terms[i].var_exp, work.terms[i].par_exp, work.terms[i].coeff);
    }
    
    printf("    Division result: %ld terms (iterations: %d)\n", quotient_mv.nterms, iterations);
    
    // Convert back to Kronecker form
    mvpoly_to_kronecker_full(quotient, &quotient_mv, var_degs, par_degs);
    
    mvpoly_clear(&dividend_mv);
    mvpoly_clear(&quotient_mv);
    mvpoly_clear(&work);
}

// NEW FUNCTION: Perform row operations on matrix - compute (row[i+1] - row[i]) / (x_i - ~x_i)
void perform_matrix_row_operations(nmod_poly_mat_t new_matrix, const nmod_poly_mat_t original_matrix,
                                  slong nvars, const slong *var_degs, const slong *par_degs, 
                                  slong npars, mp_limb_t prime)
{
    slong n = nvars + 1;
    nmod_poly_mat_init(new_matrix, n, n, prime);
    
    printf("\nPerforming row operations on matrix:\n");
    
    // First row remains unchanged
    printf("Row 0: unchanged\n");
    for (slong j = 0; j < n; j++) {
        nmod_poly_set(nmod_poly_mat_entry(new_matrix, 0, j), 
                     nmod_poly_mat_entry(original_matrix, 0, j));
    }
    
    // For rows 1 to nvars: new_row[i] = (old_row[i+1] - old_row[i]) / (x_{i-1} - ~x_{i-1})
    for (slong i = 1; i < n; i++) {
        printf("Row %ld: (row[%ld] - row[%ld]) / (x%ld - ~x%ld)\n", i, i, i-1, i-1, i-1);
        
        for (slong j = 0; j < n; j++) {
            // Compute row[i] - row[i-1]
            nmod_poly_t diff;
            nmod_poly_init(diff, prime);
            nmod_poly_sub(diff, nmod_poly_mat_entry(original_matrix, i, j),
                               nmod_poly_mat_entry(original_matrix, i-1, j));
            
            printf("  Column %ld: difference degree = %ld", j, nmod_poly_degree(diff));
            
            // Divide by (x_{i-1} - ~x_{i-1})
            nmod_poly_t quotient;
            nmod_poly_init(quotient, prime);
            
            if (nmod_poly_degree(diff) >= 0) {
                divide_by_linear_factor_kronecker(quotient, diff, i-1, var_degs, 2*nvars, 
                                                 par_degs, npars, prime);
            }
            
            printf(", quotient degree = %ld\n", nmod_poly_degree(quotient));
            
            nmod_poly_set(nmod_poly_mat_entry(new_matrix, i, j), quotient);
            
            nmod_poly_clear(diff);
            nmod_poly_clear(quotient);
        }
    }
}

// Extract coefficient matrix for final determinant
void extract_coefficient_matrix(nmod_mat_t M, const mvpoly_t *dixon_poly, 
                               slong *row_monoms, slong *col_monoms,
                               slong *nrows, slong *ncols, slong nvars)
{
    // Find all distinct monomials in x variables and ~x variables
    slong max_monoms = dixon_poly->nterms;
    slong *x_monoms = (slong*) flint_malloc(max_monoms * nvars * sizeof(slong));
    slong *dual_monoms = (slong*) flint_malloc(max_monoms * nvars * sizeof(slong));
    
    *nrows = 0;
    *ncols = 0;
    
    // Extract unique monomials
    for (slong i = 0; i < dixon_poly->nterms; i++) {
        // Check if x-monomial is new
        int found = 0;
        for (slong j = 0; j < *nrows; j++) {
            int same = 1;
            for (slong k = 0; k < nvars; k++) {
                if (x_monoms[j * nvars + k] != dixon_poly->terms[i].var_exp[k]) {
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
            for (slong k = 0; k < nvars; k++) {
                x_monoms[(*nrows) * nvars + k] = dixon_poly->terms[i].var_exp[k];
            }
            (*nrows)++;
        }
        
        // Check if ~x-monomial is new
        found = 0;
        for (slong j = 0; j < *ncols; j++) {
            int same = 1;
            for (slong k = 0; k < nvars; k++) {
                if (dual_monoms[j * nvars + k] != dixon_poly->terms[i].var_exp[nvars + k]) {
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
            for (slong k = 0; k < nvars; k++) {
                dual_monoms[(*ncols) * nvars + k] = dixon_poly->terms[i].var_exp[nvars + k];
            }
            (*ncols)++;
        }
    }
    
    printf("  Found %ld x-monomials and %ld ~x-monomials\n", *nrows, *ncols);
    
    // Initialize matrix
    nmod_mat_init(M, *nrows, *ncols, dixon_poly->mod.n);
    nmod_mat_zero(M);
    
    // Fill matrix
    for (slong t = 0; t < dixon_poly->nterms; t++) {
        // Find row index
        slong row = -1;
        for (slong i = 0; i < *nrows; i++) {
            int same = 1;
            for (slong k = 0; k < nvars; k++) {
                if (x_monoms[i * nvars + k] != dixon_poly->terms[t].var_exp[k]) {
                    same = 0;
                    break;
                }
            }
            if (same) {
                row = i;
                break;
            }
        }
        
        // Find column index
        slong col = -1;
        for (slong j = 0; j < *ncols; j++) {
            int same = 1;
            for (slong k = 0; k < nvars; k++) {
                if (dual_monoms[j * nvars + k] != dixon_poly->terms[t].var_exp[nvars + k]) {
                    same = 0;
                    break;
                }
            }
            if (same) {
                col = j;
                break;
            }
        }
        
        if (row >= 0 && col >= 0) {
            nmod_mat_entry(M, row, col) = dixon_poly->terms[t].coeff;
        }
    }
    
    // Copy monomials for output
    for (slong i = 0; i < (*nrows) * nvars; i++) {
        row_monoms[i] = x_monoms[i];
    }
    for (slong i = 0; i < (*ncols) * nvars; i++) {
        col_monoms[i] = dual_monoms[i];
    }
    
    flint_free(x_monoms);
    flint_free(dual_monoms);
}
// Helper function to evaluate a parametric polynomial at specific parameter values
mp_limb_t evaluate_mvpoly_at_params(const mvpoly_t *poly, const mp_limb_t *param_vals)
{
    mp_limb_t result = 0;
    
    for (slong i = 0; i < poly->nterms; i++) {
        mp_limb_t term_val = poly->terms[i].coeff;
        
        // Multiply by parameter powers
        if (poly->npars > 0 && poly->terms[i].par_exp) {
            for (slong j = 0; j < poly->npars; j++) {
                for (slong k = 0; k < poly->terms[i].par_exp[j]; k++) {
                    term_val = nmod_mul(term_val, param_vals[j], poly->mod);
                }
            }
        }
        
        result = nmod_add(result, term_val, poly->mod);
    }
    
    return result;
}

// Find pivot columns in row echelon form
slong find_pivot_columns(nmod_mat_t mat, slong *pivot_cols)
{
    slong nrows = nmod_mat_nrows(mat);
    slong ncols = nmod_mat_ncols(mat);
    slong rank = 0;
    
    // Convert to row echelon form (modifying the matrix)
    nmod_mat_t work;
    nmod_mat_init(work, nrows, ncols, mat->mod.n);
    nmod_mat_set(work, mat);
    
    slong current_row = 0;
    for (slong col = 0; col < ncols && current_row < nrows; col++) {
        // Find pivot in current column
        slong pivot_row = -1;
        for (slong row = current_row; row < nrows; row++) {
            if (nmod_mat_entry(work, row, col) != 0) {
                pivot_row = row;
                break;
            }
        }
        
        if (pivot_row == -1) continue; // No pivot in this column
        
        // Swap rows if necessary
        if (pivot_row != current_row) {
            for (slong j = 0; j < ncols; j++) {
                mp_limb_t temp = nmod_mat_entry(work, current_row, j);
                nmod_mat_entry(work, current_row, j) = nmod_mat_entry(work, pivot_row, j);
                nmod_mat_entry(work, pivot_row, j) = temp;
            }
        }
        
        // Record pivot column
        pivot_cols[rank] = col;
        rank++;
        
        // Eliminate below
        mp_limb_t pivot = nmod_mat_entry(work, current_row, col);
        mp_limb_t pivot_inv = nmod_inv(pivot, work->mod);
        
        for (slong row = current_row + 1; row < nrows; row++) {
            mp_limb_t val = nmod_mat_entry(work, row, col);
            if (val != 0) {
                mp_limb_t factor = nmod_mul(val, pivot_inv, work->mod);
                for (slong j = col; j < ncols; j++) {
                    mp_limb_t sub = nmod_mul(factor, nmod_mat_entry(work, current_row, j), work->mod);
                    nmod_mat_entry(work, row, j) = nmod_sub(nmod_mat_entry(work, row, j), sub, work->mod);
                }
            }
        }
        
        current_row++;
    }
    
    nmod_mat_clear(work);
    return rank;
}

// Find maximal rank submatrix by evaluating at primes
void find_maximal_rank_submatrix(mvpoly_t ***full_matrix, slong nrows, slong ncols,
                                 slong *row_indices, slong *col_indices, slong *rank,
                                 slong npars, mp_limb_t prime)
{
    printf("Finding maximal rank submatrix by evaluation...\n");
    
    // Primes to use for parameter evaluation
    mp_limb_t eval_primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};
    slong n_eval_primes = sizeof(eval_primes) / sizeof(eval_primes[0]);
    
    *rank = 0;
    
    // Try different parameter evaluations
    for (slong trial = 0; trial < 5 && *rank < FLINT_MIN(nrows, ncols); trial++) {
        // Set parameter values
        mp_limb_t *param_vals = (mp_limb_t*) flint_malloc(npars * sizeof(mp_limb_t));
        for (slong i = 0; i < npars; i++) {
            param_vals[i] = eval_primes[(trial * npars + i) % n_eval_primes];
        }
        
        printf("  Trial %ld: parameters = (", trial);
        for (slong i = 0; i < npars; i++) {
            if (i > 0) printf(", ");
            printf("%lu", (unsigned long)param_vals[i]);
        }
        printf(")\n");
        
        // Evaluate matrix at these parameter values
        nmod_mat_t eval_mat;
        nmod_mat_init(eval_mat, nrows, ncols, prime);
        
        for (slong i = 0; i < nrows; i++) {
            for (slong j = 0; j < ncols; j++) {
                mp_limb_t val = evaluate_mvpoly_at_params(full_matrix[i][j], param_vals);
                nmod_mat_entry(eval_mat, i, j) = val;
            }
        }
        
        // Find pivot columns (these give us linearly independent columns)
        slong *pivot_cols = (slong*) flint_malloc(ncols * sizeof(slong));
        slong col_rank = find_pivot_columns(eval_mat, pivot_cols);
        
        printf("    Column rank: %ld\n", col_rank);
        
        // Now work with transpose to find linearly independent rows
        nmod_mat_t eval_mat_t;
        nmod_mat_init(eval_mat_t, ncols, nrows, prime);
        nmod_mat_transpose(eval_mat_t, eval_mat);
        
        // Extract submatrix with pivot columns (which are rows in transpose)
        nmod_mat_t submat_t;
        nmod_mat_init(submat_t, col_rank, nrows, prime);
        
        for (slong i = 0; i < col_rank; i++) {
            for (slong j = 0; j < nrows; j++) {
                nmod_mat_entry(submat_t, i, j) = nmod_mat_entry(eval_mat_t, pivot_cols[i], j);
            }
        }
        
        // Find pivot columns in this submatrix (gives us linearly independent rows)
        slong *pivot_rows = (slong*) flint_malloc(nrows * sizeof(slong));
        slong row_rank = find_pivot_columns(submat_t, pivot_rows);
        
        printf("    Row rank in selected columns: %ld\n", row_rank);
        
        if (row_rank > *rank) {
            *rank = row_rank;
            
            // Store the selected rows and columns
            for (slong i = 0; i < row_rank; i++) {
                row_indices[i] = pivot_rows[i];
                col_indices[i] = pivot_cols[i];
            }
            
            printf("    Found submatrix of rank %ld\n", *rank);
        }
        
        flint_free(pivot_cols);
        flint_free(pivot_rows);
        flint_free(param_vals);
        nmod_mat_clear(eval_mat);
        nmod_mat_clear(eval_mat_t);
        nmod_mat_clear(submat_t);
        
        if (*rank == FLINT_MIN(nrows, ncols)) {
            printf("  Found full rank submatrix\n");
            break;
        }
    }
}

// Similarly for extract_parametric_resultant - compute tight bounds
void compute_coeff_matrix_det_bounds(slong *det_par_degs, 
                                    mvpoly_t **coeff_matrix,
                                    slong mat_size, slong npars)
{
    // Initialize
    for (slong p = 0; p < npars; p++) {
        det_par_degs[p] = 0;
    }
    
    // For each row, find maximum degree for each parameter
    for (slong row = 0; row < mat_size; row++) {
        slong *row_max = (slong*) flint_calloc(npars, sizeof(slong));
        
        // Check all entries in this row
        for (slong col = 0; col < mat_size; col++) {
            for (slong t = 0; t < coeff_matrix[row][col].nterms; t++) {
                for (slong p = 0; p < npars; p++) {
                    slong deg = coeff_matrix[row][col].terms[t].par_exp ? 
                               coeff_matrix[row][col].terms[t].par_exp[p] : 0;
                    if (deg > row_max[p]) {
                        row_max[p] = deg;
                    }
                }
            }
        }
        
        // Add to total
        for (slong p = 0; p < npars; p++) {
            det_par_degs[p] += row_max[p];
        }
        
        flint_free(row_max);
    }
    
    printf("Coefficient matrix determinant parameter bounds: ");
    for (slong p = 0; p < npars; p++) {
        printf("p%ld:%ld ", p, det_par_degs[p]);
    }
    printf("\n");
}

// Updated extract_parametric_resultant function
void extract_parametric_resultant(mvpoly_t *result, const mvpoly_t *dixon_poly, 
                                 slong nvars, slong npars, mp_limb_t prime)
{
    printf("\nExtracting parametric resultant from Dixon polynomial\n");
    
    // First, collect all distinct x-monomials and ~x-monomials
    typedef struct {
        slong *exp;
        slong idx;
    } monom_t;
    
    monom_t *x_monoms = (monom_t*) flint_malloc(dixon_poly->nterms * sizeof(monom_t));
    monom_t *dual_monoms = (monom_t*) flint_malloc(dixon_poly->nterms * sizeof(monom_t));
    slong nx_monoms = 0, ndual_monoms = 0;
    
    // Collect unique monomials
    for (slong i = 0; i < dixon_poly->nterms; i++) {
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
    
    printf("Found %ld x-monomials and %ld ~x-monomials\n", nx_monoms, ndual_monoms);
    
    if (nx_monoms == 0 || ndual_monoms == 0) {
        printf("Warning: Empty coefficient matrix\n");
        mvpoly_init(result, 0, npars, prime);
        goto cleanup;
    }
    
    // Build full coefficient matrix (entries are polynomials in parameters)
    mvpoly_t ***full_matrix = (mvpoly_t***) flint_malloc(nx_monoms * sizeof(mvpoly_t**));
    for (slong i = 0; i < nx_monoms; i++) {
        full_matrix[i] = (mvpoly_t**) flint_malloc(ndual_monoms * sizeof(mvpoly_t*));
        for (slong j = 0; j < ndual_monoms; j++) {
            full_matrix[i][j] = (mvpoly_t*) flint_malloc(sizeof(mvpoly_t));
            mvpoly_init(full_matrix[i][j], 0, npars, prime);
        }
    }
    
    // Fill the coefficient matrix
    for (slong t = 0; t < dixon_poly->nterms; t++) {
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
            mvpoly_add_term(full_matrix[row][col], NULL, 
                          dixon_poly->terms[t].par_exp, dixon_poly->terms[t].coeff);
        }
    }
    
    // Find maximal rank submatrix using evaluation
    slong *row_indices = (slong*) flint_malloc(FLINT_MIN(nx_monoms, ndual_monoms) * sizeof(slong));
    slong *col_indices = (slong*) flint_malloc(FLINT_MIN(nx_monoms, ndual_monoms) * sizeof(slong));
    slong submat_rank;
    
    if (npars == 0) {
        // Special case: no parameters, directly find rank
        nmod_mat_t eval_mat;
        nmod_mat_init(eval_mat, nx_monoms, ndual_monoms, prime);
        
        for (slong i = 0; i < nx_monoms; i++) {
            for (slong j = 0; j < ndual_monoms; j++) {
                mp_limb_t val = 0;
                if (full_matrix[i][j]->nterms > 0) {
                    val = full_matrix[i][j]->terms[0].coeff;
                }
                nmod_mat_entry(eval_mat, i, j) = val;
            }
        }
        
        slong *pivot_cols = (slong*) flint_malloc(ndual_monoms * sizeof(slong));
        submat_rank = find_pivot_columns(eval_mat, pivot_cols);
        
        // For square case, just use first rank indices
        for (slong i = 0; i < submat_rank; i++) {
            row_indices[i] = i;
            col_indices[i] = pivot_cols[i];
        }
        
        flint_free(pivot_cols);
        nmod_mat_clear(eval_mat);
    } else {
        find_maximal_rank_submatrix(full_matrix, nx_monoms, ndual_monoms,
                                    row_indices, col_indices, &submat_rank,
                                    npars, prime);
    }
    
    if (submat_rank == 0) {
        printf("Warning: Matrix has rank 0\n");
        mvpoly_init(result, 0, npars, prime);
        goto cleanup_matrix;
    }
    
    printf("\nExtracted submatrix of size %ld x %ld\n", submat_rank, submat_rank);
    
    // Build the submatrix
    mvpoly_t **coeff_matrix = (mvpoly_t**) flint_malloc(submat_rank * sizeof(mvpoly_t*));
    for (slong i = 0; i < submat_rank; i++) {
        coeff_matrix[i] = (mvpoly_t*) flint_malloc(submat_rank * sizeof(mvpoly_t));
        for (slong j = 0; j < submat_rank; j++) {
            mvpoly_copy(&coeff_matrix[i][j], full_matrix[row_indices[i]][col_indices[j]]);
        }
    }
    
    // Print the coefficient submatrix for debugging
    if (submat_rank <= 10) {
        printf("\nCoefficient submatrix (%ld x %ld):\n", submat_rank, submat_rank);
        for (slong i = 0; i < submat_rank; i++) {
            printf("[");
            for (slong j = 0; j < submat_rank; j++) {
                if (j > 0) printf(", ");
                if (coeff_matrix[i][j].nterms == 0) {
                    printf("0");
                } else {
                    mvpoly_print(&coeff_matrix[i][j], "");
                }
            }
            printf("]\n");
        }
    }
    
    // Now compute the determinant of the submatrix (rest of the code remains the same)
    mvpoly_init(result, 0, npars, prime);
    
    if (npars == 0) {
        // No parameters - scalar entries
        nmod_mat_t scalar_mat;
        nmod_mat_init(scalar_mat, submat_rank, submat_rank, prime);
        
        for (slong i = 0; i < submat_rank; i++) {
            for (slong j = 0; j < submat_rank; j++) {
                mp_limb_t val = 0;
                if (coeff_matrix[i][j].nterms > 0) {
                    val = coeff_matrix[i][j].terms[0].coeff;
                }
                nmod_mat_entry(scalar_mat, i, j) = val;
            }
        }

        printf("\nComputing Resultant\n");
        clock_t start = clock();  // 记录开始时间
        
        mp_limb_t det = nmod_mat_det(scalar_mat);
                
        clock_t end = clock();    // 记录结束时间
        double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
        
        printf("End (%.3f seconds)\n", elapsed);
        
        if (det != 0) {
            mvpoly_add_term(result, NULL, NULL, det);
        }
        
        nmod_mat_clear(scalar_mat);
        
    } else if (npars == 1) {
        // One parameter - use nmod_poly_mat_det directly
        nmod_poly_mat_t poly_mat;
        nmod_poly_mat_init(poly_mat, submat_rank, submat_rank, prime);
        
        // Convert mvpoly to nmod_poly
        for (slong i = 0; i < submat_rank; i++) {
            for (slong j = 0; j < submat_rank; j++) {
                nmod_poly_t entry;
                nmod_poly_init(entry, prime);
                
                for (slong k = 0; k < coeff_matrix[i][j].nterms; k++) {
                    slong deg = coeff_matrix[i][j].terms[k].par_exp ? 
                               coeff_matrix[i][j].terms[k].par_exp[0] : 0;
                    mp_limb_t coeff = coeff_matrix[i][j].terms[k].coeff;
                    nmod_poly_set_coeff_ui(entry, deg, coeff);
                }
                
                nmod_poly_set(nmod_poly_mat_entry(poly_mat, i, j), entry);
                nmod_poly_clear(entry);
            }
        }
        
        // Compute determinant
        nmod_poly_t det_poly;
        nmod_poly_init(det_poly, prime);

        printf("\nComputing Resultant\n");
        clock_t start = clock();  // 记录开始时间
        
        nmod_poly_mat_det_iter(det_poly, poly_mat);
                
        clock_t end = clock();    // 记录结束时间
        double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
        
        printf("End (%.3f seconds)\n", elapsed);
        // Convert back to mvpoly
        slong det_deg = nmod_poly_degree(det_poly);
        for (slong i = 0; i <= det_deg; i++) {
            mp_limb_t coeff = nmod_poly_get_coeff_ui(det_poly, i);
            if (coeff != 0) {
                slong par_exp[1] = {i};
                mvpoly_add_term(result, NULL, par_exp, coeff);
            }
        }
        
        nmod_poly_clear(det_poly);
        nmod_poly_mat_clear(poly_mat);
        
    } else {
        // Multiple parameters - use Kronecker substitution
        nmod_poly_mat_t poly_mat;
        nmod_poly_mat_init(poly_mat, submat_rank, submat_rank, prime);
        
        // Compute tight degree bounds for determinant
        slong *par_degs = (slong*) flint_calloc(npars, sizeof(slong));
        compute_coeff_matrix_det_bounds(par_degs, coeff_matrix, submat_rank, npars);
        
        // Convert each matrix entry via Kronecker
        slong dummy_var_degs[1] = {0}; // No variables, only parameters
        for (slong i = 0; i < submat_rank; i++) {
            for (slong j = 0; j < submat_rank; j++) {
                nmod_poly_t entry;
                nmod_poly_init(entry, prime);
                mvpoly_to_kronecker_full(entry, &coeff_matrix[i][j], 
                                       dummy_var_degs, par_degs);
                nmod_poly_set(nmod_poly_mat_entry(poly_mat, i, j), entry);
                nmod_poly_clear(entry);
            }
        }
            
        // Compute determinant
        nmod_poly_t det_poly;
        nmod_poly_init(det_poly, prime);

        printf("\nComputing Resultant\n");
        clock_t start = clock();  // 记录开始时间
        
        nmod_poly_mat_det(det_poly, poly_mat); //nmod_poly_mat_det_iter
                
        clock_t end = clock();    // 记录结束时间
        double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
        
        printf("End (%.3f seconds)\n", elapsed);
        
        // Convert back from Kronecker
        kronecker_to_mvpoly_full(result, det_poly, dummy_var_degs, 0, 
                               par_degs, npars, prime);
        
        flint_free(par_degs);
        nmod_poly_clear(det_poly);
        nmod_poly_mat_clear(poly_mat);
    }
    
    printf("\nResultant polynomial has %ld terms\n", result->nterms);
    
    // Cleanup submatrix
    for (slong i = 0; i < submat_rank; i++) {
        for (slong j = 0; j < submat_rank; j++) {
            mvpoly_clear(&coeff_matrix[i][j]);
        }
        flint_free(coeff_matrix[i]);
    }
    flint_free(coeff_matrix);
    
cleanup_matrix:
    // Cleanup full matrix
    for (slong i = 0; i < nx_monoms; i++) {
        for (slong j = 0; j < ndual_monoms; j++) {
            mvpoly_clear(full_matrix[i][j]);
            flint_free(full_matrix[i][j]);
        }
        flint_free(full_matrix[i]);
    }
    flint_free(full_matrix);
    flint_free(row_indices);
    flint_free(col_indices);
    
cleanup:
    for (slong i = 0; i < nx_monoms; i++) {
        flint_free(x_monoms[i].exp);
    }
    for (slong j = 0; j < ndual_monoms; j++) {
        flint_free(dual_monoms[j].exp);
    }
    flint_free(x_monoms);
    flint_free(dual_monoms);
}

// Modified dixon_resultant function using row operations
void dixon_resultant_with_row_ops(mvpoly_t *result, mvpoly_t *polys, slong nvars, slong npars, mp_limb_t prime)
{
    printf("\n=== Modified Dixon Resultant with Row Operations ===\n");
    printf("Variables: %ld, Parameters: %ld, Prime: %lu\n", nvars, npars, (unsigned long)prime);
    
    slong npolys = nvars + 1;
    
    // Allocate space for degree bounds
    slong *det_var_degs = (slong*) flint_calloc(2 * nvars, sizeof(slong));
    slong *det_par_degs = (slong*) flint_calloc(npars, sizeof(slong));
    
    // Step 1: Build cancellation matrix and get degree bounds
    printf("\nStep 1: Build Cancellation Matrix\n");
    nmod_poly_mat_t M;
    build_cancellation_matrix_with_bounds(M, polys, nvars, npars, prime, 
                                         det_var_degs, det_par_degs);
    
    // Print small matrices
    if (nvars <= 3) {
        printf("\nOriginal Cancellation Matrix:\n");
        for (slong i = 0; i < npolys; i++) {
            printf("[");
            for (slong j = 0; j < npolys; j++) {
                if (j > 0) printf(", ");
                printf("deg=%ld", nmod_poly_degree(nmod_poly_mat_entry(M, i, j)));
            }
            printf("]\n");
        }
    }
    
    // Step 2: Perform row operations on matrix (NEW APPROACH)
    printf("\nStep 2: Perform Matrix Row Operations\n");
    nmod_poly_mat_t modified_M;
    perform_matrix_row_operations(modified_M, M, nvars, det_var_degs, det_par_degs, 
                                 npars, prime);
    
    // Print the modified matrix
    if (nvars <= 3) {
        printf("\nModified Matrix after row operations:\n");
        for (slong i = 0; i < npolys; i++) {
            printf("[");
            for (slong j = 0; j < npolys; j++) {
                if (j > 0) printf(", ");
                printf("deg=%ld", nmod_poly_degree(nmod_poly_mat_entry(modified_M, i, j)));
            }
            printf("]\n");
        }
    }
    
    // Step 3: Compute determinant of modified matrix
    printf("\nStep 3: Compute Determinant of Modified Matrix\n");
    nmod_poly_t det_poly;
    nmod_poly_init(det_poly, prime);
    
    printf("\nComputing Dixon Polynomial\n");
    clock_t start = clock();  // 记录开始时间
    
    nmod_poly_mat_det(det_poly, modified_M); //nmod_poly_mat_det_iter
    
    clock_t end = clock();    // 记录结束时间
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    
    printf("End (%.3f seconds)\n", elapsed);
    
    printf("Modified matrix determinant degree: %ld\n", nmod_poly_degree(det_poly));
    
    if (nmod_poly_degree(det_poly) < 0) {
        printf("Warning: Determinant is zero!\n");
        mvpoly_init(result, 0, npars, prime);
        goto cleanup;
    }

    mvpoly_t d_poly;
    
    // Step 4: Convert result back to multivariate form
    printf("\nStep 4: Convert Result to Multivariate Form\n");
    kronecker_to_mvpoly_full(&d_poly, det_poly, det_var_degs, 2 * nvars, 
                            det_par_degs, npars, prime);
    
    printf("Dixon polynomial has %ld terms\n", d_poly.nterms);

    if (d_poly.nterms <= 100) {
        mvpoly_print_expanded(&d_poly, "Dixon", 1);
    }
    
    // Step 4: Extract coefficient matrix (if no parameters)
    if (npars == 0) {
        printf("\nStep 4: Extract Coefficient Matrix\n");
        
        slong nrows, ncols;
        slong *row_monoms = (slong*) flint_malloc(d_poly.nterms * nvars * sizeof(slong));
        slong *col_monoms = (slong*) flint_malloc(d_poly.nterms * nvars * sizeof(slong));
        
        nmod_mat_t coeff_mat;
        extract_coefficient_matrix(coeff_mat, &d_poly, row_monoms, col_monoms, 
                                  &nrows, &ncols, nvars);
        
        printf("Coefficient matrix: %ld x %ld\n", nrows, ncols);
        
        // Find maximal non-singular submatrix
        slong size = FLINT_MIN(nrows, ncols);
        if (size > 0) {
            nmod_mat_t submat;
            nmod_mat_init(submat, size, size, prime);
            
            for (slong i = 0; i < size; i++) {
                for (slong j = 0; j < size; j++) {
                    nmod_mat_entry(submat, i, j) = nmod_mat_entry(coeff_mat, i, j);
                }
            }
            
            mp_limb_t det = nmod_mat_det(submat);
            printf("Dixon resultant: %lu\n", (unsigned long)det);
            
            // Store result
            mvpoly_init(result, 0, 0, prime);
            if (det != 0) {
                mvpoly_add_term(result, NULL, NULL, det);
            }
            
            nmod_mat_clear(submat);
        }
        
        nmod_mat_clear(coeff_mat);
        flint_free(row_monoms);
        flint_free(col_monoms);
    } else {
        // With parameters, extract the parametric resultant
        extract_parametric_resultant(result, &d_poly, nvars, npars, prime);
    }
cleanup:
    // Cleanup
    mvpoly_clear(&d_poly);
    nmod_poly_clear(det_poly);
    nmod_poly_mat_clear(M);
    flint_free(det_var_degs);
    flint_free(det_par_degs);
    
}

// Test case 1: No parameters
void test_case_1()
{
    printf("\n========== Test Case 1: System with No Parameters ==========\n");
    
    mp_limb_t prime = 101;
    slong nvars = 2;
    slong npars = 0;
    slong npolys = nvars + 1;
    
    mvpoly_t *polys = (mvpoly_t*) flint_malloc(npolys * sizeof(mvpoly_t));
    
    for (slong i = 0; i < npolys; i++) {
        mvpoly_init(&polys[i], nvars, npars, prime);
    }
    
    slong var_exp[2];
    
    // p0 = x^2 + y^2 - 1
    var_exp[0] = 2; var_exp[1] = 0;
    mvpoly_add_term(&polys[0], var_exp, NULL, 1);
    var_exp[0] = 0; var_exp[1] = 2;
    mvpoly_add_term(&polys[0], var_exp, NULL, 1);
    var_exp[0] = 0; var_exp[1] = 0;
    mvpoly_add_term(&polys[0], var_exp, NULL, prime - 1);
    
    // p1 = x + y - 1
    var_exp[0] = 1; var_exp[1] = 0;
    mvpoly_add_term(&polys[1], var_exp, NULL, 1);
    var_exp[0] = 0; var_exp[1] = 1;
    mvpoly_add_term(&polys[1], var_exp, NULL, 1);
    var_exp[0] = 0; var_exp[1] = 0;
    mvpoly_add_term(&polys[1], var_exp, NULL, prime - 1);
    
    // p2 = 2x - y
    var_exp[0] = 1; var_exp[1] = 0;
    mvpoly_add_term(&polys[2], var_exp, NULL, 2);
    var_exp[0] = 0; var_exp[1] = 1;
    mvpoly_add_term(&polys[2], var_exp, NULL, prime - 1);
    
    printf("\nInput polynomials:\n");
    for (slong i = 0; i < npolys; i++) {
        char name[16];
        sprintf(name, "p%ld", i);
        mvpoly_print(&polys[i], name);
    }
    
    mvpoly_t result;
    dixon_resultant_with_row_ops(&result, polys, nvars, npars, prime);
    
    printf("\nFinal result: ");
    mvpoly_print(&result, "R");
    
    mvpoly_clear(&result);
    for (slong i = 0; i < npolys; i++) {
        mvpoly_clear(&polys[i]);
    }
    flint_free(polys);
}

// Test case 2: One parameter
void test_case_2()
{
    printf("\n========== Test Case 2: System with One Parameter ==========\n");
    
    mp_limb_t prime = 101;
    slong nvars = 2;
    slong npars = 1;
    slong npolys = nvars + 1;
    
    mvpoly_t *polys = (mvpoly_t*) flint_malloc(npolys * sizeof(mvpoly_t));
    
    for (slong i = 0; i < npolys; i++) {
        mvpoly_init(&polys[i], nvars, npars, prime);
    }
    
    slong var_exp[2], par_exp[1];
    
    // p0 = x^2 + y^2 - 1
    var_exp[0] = 2; var_exp[1] = 0; par_exp[0] = 0;
    mvpoly_add_term(&polys[0], var_exp, par_exp, 1);
    var_exp[0] = 0; var_exp[1] = 2;
    mvpoly_add_term(&polys[0], var_exp, par_exp, 1);
    var_exp[0] = 0; var_exp[1] = 0;
    mvpoly_add_term(&polys[0], var_exp, par_exp, prime - 1);
    
    // p1 = x + y - a
    var_exp[0] = 1; var_exp[1] = 0; par_exp[0] = 0;
    mvpoly_add_term(&polys[1], var_exp, par_exp, 1);
    var_exp[0] = 0; var_exp[1] = 1;
    mvpoly_add_term(&polys[1], var_exp, par_exp, 1);
    var_exp[0] = 0; var_exp[1] = 0; par_exp[0] = 1;
    mvpoly_add_term(&polys[1], var_exp, par_exp, prime - 1);
    
    // p2 = 2x - y
    var_exp[0] = 1; var_exp[1] = 0; par_exp[0] = 0;
    mvpoly_add_term(&polys[2], var_exp, par_exp, 2);
    var_exp[0] = 0; var_exp[1] = 1;
    mvpoly_add_term(&polys[2], var_exp, par_exp, prime - 1);
    
    printf("\nInput polynomials:\n");
    for (slong i = 0; i < npolys; i++) {
        char name[16];
        sprintf(name, "p%ld", i);
        mvpoly_print(&polys[i], name);
    }
    
    mvpoly_t result;
    dixon_resultant_with_row_ops(&result, polys, nvars, npars, prime);
    
    printf("\nFinal result: ");
    mvpoly_print(&result, "R(a)");
    
    mvpoly_clear(&result);
    for (slong i = 0; i < npolys; i++) {
        mvpoly_clear(&polys[i]);
    }
    flint_free(polys);
}


// ============= String Parser Interface =============

// Token types for the parser
typedef enum {
    TOK_NUMBER,
    TOK_VARIABLE,
    TOK_PLUS,
    TOK_MINUS,
    TOK_MULT,
    TOK_POWER,
    TOK_LPAREN,
    TOK_RPAREN,
    TOK_EOF
} token_type_t;

typedef struct {
    token_type_t type;
    char *str;
    mp_limb_t value;
} token_t;

typedef struct {
    const char *input;
    size_t pos;
    token_t current;
    
    // Variable and parameter tracking
    char **var_names;
    slong nvars;
    char **par_names;
    slong npars;
    slong max_pars;
    
    // Current parsing context
    mp_limb_t mod;
} parser_state_t;

// Find or add a parameter
static slong find_or_add_parameter(parser_state_t *state, const char *name)
{
    // First check if it's a variable
    for (slong i = 0; i < state->nvars; i++) {
        if (strcmp(state->var_names[i], name) == 0) {
            return -1; // It's a variable, not a parameter
        }
    }
    
    // Check existing parameters
    for (slong i = 0; i < state->npars; i++) {
        if (strcmp(state->par_names[i], name) == 0) {
            return i;
        }
    }
    
    // Add new parameter
    if (state->npars >= state->max_pars) {
        state->max_pars *= 2;
        state->par_names = (char**) realloc(state->par_names, 
                                           state->max_pars * sizeof(char*));
    }
    
    state->par_names[state->npars] = strdup(name);
    return state->npars++;
}

// Get variable index
static slong get_variable_index(parser_state_t *state, const char *name)
{
    for (slong i = 0; i < state->nvars; i++) {
        if (strcmp(state->var_names[i], name) == 0) {
            return i;
        }
    }
    return -1;
}

// Tokenizer
static void next_token(parser_state_t *state)
{
    // Skip whitespace
    while (state->input[state->pos] && isspace(state->input[state->pos])) {
        state->pos++;
    }
    
    // Free previous token string
    if (state->current.str) {
        free(state->current.str);
        state->current.str = NULL;
    }
    
    if (!state->input[state->pos]) {
        state->current.type = TOK_EOF;
        return;
    }
    
    char c = state->input[state->pos];
    
    // Single character tokens
    if (c == '+') {
        state->current.type = TOK_PLUS;
        state->pos++;
        return;
    } else if (c == '-') {
        state->current.type = TOK_MINUS;
        state->pos++;
        return;
    } else if (c == '*') {
        state->current.type = TOK_MULT;
        state->pos++;
        return;
    } else if (c == '^') {
        state->current.type = TOK_POWER;
        state->pos++;
        return;
    } else if (c == '(') {
        state->current.type = TOK_LPAREN;
        state->pos++;
        return;
    } else if (c == ')') {
        state->current.type = TOK_RPAREN;
        state->pos++;
        return;
    }
    
    // Numbers
    if (isdigit(c)) {
        size_t start = state->pos;
        while (isdigit(state->input[state->pos])) {
            state->pos++;
        }
        
        size_t len = state->pos - start;
        state->current.str = (char*) malloc(len + 1);
        strncpy(state->current.str, state->input + start, len);
        state->current.str[len] = '\0';
        state->current.value = atol(state->current.str);
        state->current.type = TOK_NUMBER;
        return;
    }
    
    // Variables/parameters (alphabetic)
    if (isalpha(c) || c == '_') {
        size_t start = state->pos;
        while (isalnum(state->input[state->pos]) || state->input[state->pos] == '_') {
            state->pos++;
        }
        
        size_t len = state->pos - start;
        state->current.str = (char*) malloc(len + 1);
        strncpy(state->current.str, state->input + start, len);
        state->current.str[len] = '\0';
        state->current.type = TOK_VARIABLE;
        return;
    }
    
    // Unknown character
    fprintf(stderr, "Unknown character: %c at position %zu\n", c, state->pos);
    state->pos++;
    next_token(state);
}

// Forward declarations
static void parse_polynomial(parser_state_t *state, mvpoly_t *poly);
static void parse_expression(parser_state_t *state, mvpoly_t *poly);
static void parse_term(parser_state_t *state, mvpoly_t *poly);
static void parse_factor(parser_state_t *state, mvpoly_t *poly);
static void parse_base(parser_state_t *state, mvpoly_t *poly);

// Parse base (number, variable, or parenthesized expression)
static void parse_base(parser_state_t *state, mvpoly_t *poly)
{
    if (state->current.type == TOK_NUMBER) {
        mp_limb_t coeff = state->current.value % state->mod;
        mvpoly_add_term(poly, NULL, NULL, coeff);
        next_token(state);
        
    } else if (state->current.type == TOK_VARIABLE) {
        char *name = strdup(state->current.str);
        next_token(state);
        
        slong var_idx = get_variable_index(state, name);
        if (var_idx >= 0) {
            // It's a variable
            slong *var_exp = (slong*) calloc(state->nvars, sizeof(slong));
            var_exp[var_idx] = 1;
            mvpoly_add_term(poly, var_exp, NULL, 1);
            free(var_exp);
        } else {
            // It's a parameter
            slong par_idx = find_or_add_parameter(state, name);
            if (par_idx >= 0) {
                slong *par_exp = (slong*) calloc(state->max_pars, sizeof(slong));
                par_exp[par_idx] = 1;
                mvpoly_add_term(poly, NULL, par_exp, 1);
                free(par_exp);
            }
        }
        free(name);
        
    } else if (state->current.type == TOK_LPAREN) {
        next_token(state);
        parse_expression(state, poly);
        if (state->current.type != TOK_RPAREN) {
            fprintf(stderr, "Expected ')'\n");
        } else {
            next_token(state);
        }
        
    } else if (state->current.type == TOK_MINUS) {
        next_token(state);
        mvpoly_t temp;
        mvpoly_init(&temp, state->nvars, state->max_pars, state->mod);
        parse_base(state, &temp);
        
        // Negate all terms
        for (slong i = 0; i < temp.nterms; i++) {
            mp_limb_t neg_coeff = nmod_neg(temp.terms[i].coeff, temp.mod);
            mvpoly_add_term(poly, temp.terms[i].var_exp, temp.terms[i].par_exp, neg_coeff);
        }
        mvpoly_clear(&temp);
        
    } else {
        fprintf(stderr, "Unexpected token in parse_base\n");
    }
}

// Parse factor (base with optional exponent)
static void parse_factor(parser_state_t *state, mvpoly_t *poly)
{
    mvpoly_t base;
    mvpoly_init(&base, state->nvars, state->max_pars, state->mod);
    parse_base(state, &base);
    
    if (state->current.type == TOK_POWER) {
        next_token(state);
        if (state->current.type != TOK_NUMBER) {
            fprintf(stderr, "Expected number after ^\n");
            mvpoly_clear(&base);
            return;
        }
        
        slong exp = state->current.value;
        next_token(state);
        
        // Compute base^exp
        mvpoly_t result;
        mvpoly_init(&result, state->nvars, state->max_pars, state->mod);
        mvpoly_add_term(&result, NULL, NULL, 1); // Start with 1
        
        for (slong i = 0; i < exp; i++) {
            mvpoly_t temp;
            mvpoly_init(&temp, state->nvars, state->max_pars, state->mod);
            
            // Multiply result by base
            for (slong j = 0; j < result.nterms; j++) {
                for (slong k = 0; k < base.nterms; k++) {
                    // Multiply terms
                    mp_limb_t coeff = nmod_mul(result.terms[j].coeff, base.terms[k].coeff, result.mod);
                    
                    slong *var_exp = NULL;
                    slong *par_exp = NULL;
                    
                    if (state->nvars > 0) {
                        var_exp = (slong*) calloc(state->nvars, sizeof(slong));
                        for (slong m = 0; m < state->nvars; m++) {
                            slong e1 = result.terms[j].var_exp ? result.terms[j].var_exp[m] : 0;
                            slong e2 = base.terms[k].var_exp ? base.terms[k].var_exp[m] : 0;
                            var_exp[m] = e1 + e2;
                        }
                    }
                    
                    if (state->max_pars > 0) {
                        par_exp = (slong*) calloc(state->max_pars, sizeof(slong));
                        for (slong m = 0; m < state->max_pars; m++) {
                            slong e1 = result.terms[j].par_exp ? result.terms[j].par_exp[m] : 0;
                            slong e2 = base.terms[k].par_exp ? base.terms[k].par_exp[m] : 0;
                            par_exp[m] = e1 + e2;
                        }
                    }
                    
                    mvpoly_add_term(&temp, var_exp, par_exp, coeff);
                    if (var_exp) free(var_exp);
                    if (par_exp) free(par_exp);
                }
            }
            
            mvpoly_clear(&result);
            result = temp;
        }
        
        // Add result to poly
        for (slong i = 0; i < result.nterms; i++) {
            mvpoly_add_term(poly, result.terms[i].var_exp, result.terms[i].par_exp, result.terms[i].coeff);
        }
        mvpoly_clear(&result);
        
    } else {
        // Just copy base to poly
        for (slong i = 0; i < base.nterms; i++) {
            mvpoly_add_term(poly, base.terms[i].var_exp, base.terms[i].par_exp, base.terms[i].coeff);
        }
    }
    
    mvpoly_clear(&base);
}

// Parse term (product of factors)
static void parse_term(parser_state_t *state, mvpoly_t *poly)
{
    mvpoly_t result;
    mvpoly_init(&result, state->nvars, state->max_pars, state->mod);
    
    parse_factor(state, &result);
    
    while (state->current.type == TOK_MULT) {
        next_token(state);
        
        mvpoly_t factor;
        mvpoly_init(&factor, state->nvars, state->max_pars, state->mod);
        parse_factor(state, &factor);
        
        // Multiply result by factor
        mvpoly_t temp;
        mvpoly_init(&temp, state->nvars, state->max_pars, state->mod);
        
        for (slong i = 0; i < result.nterms; i++) {
            for (slong j = 0; j < factor.nterms; j++) {
                mp_limb_t coeff = nmod_mul(result.terms[i].coeff, factor.terms[j].coeff, result.mod);
                
                slong *var_exp = NULL;
                slong *par_exp = NULL;
                
                if (state->nvars > 0) {
                    var_exp = (slong*) calloc(state->nvars, sizeof(slong));
                    for (slong k = 0; k < state->nvars; k++) {
                        slong e1 = result.terms[i].var_exp ? result.terms[i].var_exp[k] : 0;
                        slong e2 = factor.terms[j].var_exp ? factor.terms[j].var_exp[k] : 0;
                        var_exp[k] = e1 + e2;
                    }
                }
                
                if (state->max_pars > 0) {
                    par_exp = (slong*) calloc(state->max_pars, sizeof(slong));
                    for (slong k = 0; k < state->max_pars; k++) {
                        slong e1 = result.terms[i].par_exp ? result.terms[i].par_exp[k] : 0;
                        slong e2 = factor.terms[j].par_exp ? factor.terms[j].par_exp[k] : 0;
                        par_exp[k] = e1 + e2;
                    }
                }
                
                mvpoly_add_term(&temp, var_exp, par_exp, coeff);
                if (var_exp) free(var_exp);
                if (par_exp) free(par_exp);
            }
        }
        
        mvpoly_clear(&result);
        mvpoly_clear(&factor);
        result = temp;
    }
    
    // Add result to poly
    for (slong i = 0; i < result.nterms; i++) {
        mvpoly_add_term(poly, result.terms[i].var_exp, result.terms[i].par_exp, result.terms[i].coeff);
    }
    mvpoly_clear(&result);
}

// Parse expression (sum of terms)
static void parse_expression(parser_state_t *state, mvpoly_t *poly)
{
    // Handle leading sign
    int negate = 0;
    if (state->current.type == TOK_MINUS) {
        negate = 1;
        next_token(state);
    } else if (state->current.type == TOK_PLUS) {
        next_token(state);
    }
    
    mvpoly_t first_term;
    mvpoly_init(&first_term, state->nvars, state->max_pars, state->mod);
    parse_term(state, &first_term);
    
    if (negate) {
        for (slong i = 0; i < first_term.nterms; i++) {
            mp_limb_t neg_coeff = nmod_neg(first_term.terms[i].coeff, first_term.mod);
            mvpoly_add_term(poly, first_term.terms[i].var_exp, first_term.terms[i].par_exp, neg_coeff);
        }
    } else {
        for (slong i = 0; i < first_term.nterms; i++) {
            mvpoly_add_term(poly, first_term.terms[i].var_exp, first_term.terms[i].par_exp, 
                           first_term.terms[i].coeff);
        }
    }
    mvpoly_clear(&first_term);
    
    // Additional terms
    while (state->current.type == TOK_PLUS || state->current.type == TOK_MINUS) {
        int subtract = (state->current.type == TOK_MINUS);
        next_token(state);
        
        mvpoly_t term;
        mvpoly_init(&term, state->nvars, state->max_pars, state->mod);
        parse_term(state, &term);
        
        for (slong i = 0; i < term.nterms; i++) {
            mp_limb_t coeff = term.terms[i].coeff;
            if (subtract) {
                coeff = nmod_neg(coeff, term.mod);
            }
            mvpoly_add_term(poly, term.terms[i].var_exp, term.terms[i].par_exp, coeff);
        }
        mvpoly_clear(&term);
    }
}

// Main polynomial parser
static void parse_polynomial(parser_state_t *state, mvpoly_t *poly)
{
    parse_expression(state, poly);
}

// Main string-based interface
void compute_dixon_resultant_string(const char **poly_strings, slong npoly_strings,
                                   const char **var_names, slong nvars,
                                   mp_limb_t prime)
{
    printf("\n=== Dixon Resultant String Interface ===\n");
    printf("Prime: %lu\n", (unsigned long)prime);
    printf("Variables (%ld): ", nvars);
    for (slong i = 0; i < nvars; i++) {
        if (i > 0) printf(", ");
        printf("%s", var_names[i]);
    }
    printf("\n");
    
    if (npoly_strings != nvars + 1) {
        fprintf(stderr, "Error: Need exactly %ld polynomials for %ld variables\n",
                nvars + 1, nvars);
        return;
    }
    
    // Initialize parser state
    parser_state_t state;
    state.var_names = (char**) malloc(nvars * sizeof(char*));
    for (slong i = 0; i < nvars; i++) {
        state.var_names[i] = strdup(var_names[i]);
    }
    state.nvars = nvars;
    state.npars = 0;
    state.max_pars = 16;
    state.par_names = (char**) malloc(state.max_pars * sizeof(char*));
    state.mod = prime;
    state.current.str = NULL;
    
    // First pass: identify all parameters
    printf("\nFirst pass: identifying parameters...\n");
    for (slong i = 0; i < npoly_strings; i++) {
        mvpoly_t temp;
        mvpoly_init(&temp, nvars, state.max_pars, prime);
        
        state.input = poly_strings[i];
        state.pos = 0;
        next_token(&state);
        
        parse_polynomial(&state, &temp);
        mvpoly_clear(&temp);
        
        if (state.current.type != TOK_EOF) {
            fprintf(stderr, "Warning: Extra characters at end of polynomial %ld\n", i);
        }
    }
    
    // Report detected parameters
    printf("\nDetected parameters (%ld): ", state.npars);
    if (state.npars == 0) {
        printf("none");
    } else {
        for (slong i = 0; i < state.npars; i++) {
            if (i > 0) printf(", ");
            printf("%s", state.par_names[i]);
        }
    }
    printf("\n");
    
    // Second pass: parse polynomials with correct parameter count
    mvpoly_t *polys = (mvpoly_t*) malloc(npoly_strings * sizeof(mvpoly_t));
    printf("\nParsing polynomials:\n");
    
    for (slong i = 0; i < npoly_strings; i++) {
        mvpoly_init(&polys[i], nvars, state.npars, prime);
        
        state.input = poly_strings[i];
        state.pos = 0;
        if (state.current.str) {
            free(state.current.str);
            state.current.str = NULL;
        }
        next_token(&state);
        
        parse_polynomial(&state, &polys[i]);
        
        printf("  p%ld = %s => ", i, poly_strings[i]);
        mvpoly_print(&polys[i], "");
    }
    
    // Compute Dixon resultant
    printf("\nComputing Dixon resultant...\n");
    mvpoly_t result;
    dixon_resultant_with_row_ops(&result, polys, nvars, state.npars, prime);
    
    printf("\n=== Final Result ===\n");
    printf("Dixon Resultant: ");
    mvpoly_print(&result, "");
    printf("\n");
    
    // Cleanup
    mvpoly_clear(&result);
    for (slong i = 0; i < npoly_strings; i++) {
        mvpoly_clear(&polys[i]);
    }
    free(polys);
    
    for (slong i = 0; i < nvars; i++) {
        free(state.var_names[i]);
    }
    free(state.var_names);
    
    for (slong i = 0; i < state.npars; i++) {
        free(state.par_names[i]);
    }
    free(state.par_names);
    
    if (state.current.str) {
        free(state.current.str);
    }
}

// Example usage
int main()
{
    printf("Dixon Resultant Computation with String Parser\n");
    printf("=============================================\n");
    
    // Example 1: Your requested example
    {
        const char *polys[] = {
            "x^2 + y^2 + z^2 - 1",
            "x + y + z - a",
            "a*x + b*y - z",
            "x*y + y*z + z*x - b"
        };
        const char *vars[] = {"x", "y", "z"};
        
        compute_dixon_resultant_string(polys, 4, vars, 3, 101);
    }
    
    // Example 2: Simple 2D system
    {
        printf("\n\n=== Example 2: Simple 2D System ===\n");
        const char *polys[] = {
            "x^2 + y^2 - 1",
            "x + y - 1", 
            "2*x - y"
        };
        const char *vars[] = {"x", "y"};
        
        compute_dixon_resultant_string(polys, 3, vars, 2, 101);
    }
    
    // Example 3: With complex expressions
    {
        printf("\n\n=== Example 3: Complex Expressions ===\n");
        const char *polys[] = {
            "x^13 - 3*x*y^2 + 1",
            "3*x^2*y - y^13",
            "x^12 + y^2 - a^2"
        };
        const char *vars[] = {"x", "y"};
        
        compute_dixon_resultant_string(polys, 3, vars, 2, 101);
    }

    // Example 4: With more complex expressions
    {
        printf("\n\n=== Example 4: Complex Expressions ===\n");
        const char *polys[] = {"18284*x0^10 - 31925*x0^9*x1 + 10273*x0^9*x2 + 21254*x0^9 + 31834*x0^8*x1^2 - 7283*x0^8*x1*x2 - 2904*x0^8*x1 + 8729*x0^8*x2^2 + 9717*x0^8*x2 + 21720*x0^8 + 4348*x0^7*x1^3 - 10388*x0^7*x1^2*x2 + 5000*x0^7*x1^2 + 16643*x0^7*x1*x2^2 + 21708*x0^7*x1*x2 - 11969*x0^7*x1 + 4188*x0^7*x2^3 - 25424*x0^7*x2^2 + 19409*x0^7*x2 - 5865*x0^7 - 8367*x0^6*x1^4 - 5132*x0^6*x1^3*x2 - 29514*x0^6*x1^3 + 14710*x0^6*x1^2*x2^2 - 5285*x0^6*x1^2*x2 - 4907*x0^6*x1^2 - 3302*x0^6*x1*x2^3 - 9961*x0^6*x1*x2^2 + 560*x0^6*x1*x2 - 28327*x0^6*x1 + 12177*x0^6*x2^4 - 5119*x0^6*x2^3 + 31236*x0^6*x2^2 + 12387*x0^6*x2 - 22594*x0^6 - 8063*x0^5*x1^5 + 22487*x0^5*x1^4*x2 + 9438*x0^5*x1^4 - 2395*x0^5*x1^3*x2^2 - 144*x0^5*x1^3*x2 - 13472*x0^5*x1^3 - 24289*x0^5*x1^2*x2^3 - 28235*x0^5*x1^2*x2^2 + 767*x0^5*x1^2*x2 - 10200*x0^5*x1^2 - 11823*x0^5*x1*x2^4 - 7499*x0^5*x1*x2^3 + 28165*x0^5*x1*x2^2 + 22264*x0^5*x1*x2 + 6855*x0^5*x1 + 18137*x0^5*x2^5 + 4239*x0^5*x2^4 - 8457*x0^5*x2^3 - 18094*x0^5*x2^2 - 18209*x0^5*x2 + 23437*x0^5 + 40*x0^4*x1^6 + 30025*x0^4*x1^5*x2 - 12713*x0^4*x1^5 - 15712*x0^4*x1^4*x2^2 - 1744*x0^4*x1^4*x2 - 18056*x0^4*x1^4 - 19269*x0^4*x1^3*x2^3 - 11695*x0^4*x1^3*x2^2 + 16148*x0^4*x1^3*x2 + 22693*x0^4*x1^3 - 18584*x0^4*x1^2*x2^4 - 19429*x0^4*x1^2*x2^3 + 4245*x0^4*x1^2*x2^2 + 16105*x0^4*x1^2*x2 + 24164*x0^4*x1^2 - 14871*x0^4*x1*x2^5 - 1290*x0^4*x1*x2^4 - 17390*x0^4*x1*x2^3 + 12108*x0^4*x1*x2^2 + 11435*x0^4*x1*x2 + 21602*x0^4*x1 + 6930*x0^4*x2^6 - 7210*x0^4*x2^5 + 16728*x0^4*x2^4 - 619*x0^4*x2^3 - 1549*x0^4*x2^2 - 4789*x0^4*x2 - 6886*x0^4 - 31787*x0^3*x1^7 + 19033*x0^3*x1^6*x2 - 15981*x0^3*x1^6 - 3967*x0^3*x1^5*x2^2 + 27459*x0^3*x1^5*x2 - 11338*x0^3*x1^5 - 4226*x0^3*x1^4*x2^3 + 25431*x0^3*x1^4*x2^2 - 29484*x0^3*x1^4*x2 + 13674*x0^3*x1^4 + 10389*x0^3*x1^3*x2^4 + 13489*x0^3*x1^3*x2^3 + 511*x0^3*x1^3*x2^2 - 23225*x0^3*x1^3*x2 - 5616*x0^3*x1^3 + 23197*x0^3*x1^2*x2^5 - 31688*x0^3*x1^2*x2^4 - 30228*x0^3*x1^2*x2^3 + 101*x0^3*x1^2*x2^2 - 29757*x0^3*x1^2*x2 - 23620*x0^3*x1^2 + 28673*x0^3*x1*x2^6 + 5135*x0^3*x1*x2^5 - 26908*x0^3*x1*x2^4 + 26636*x0^3*x1*x2^3 - 28809*x0^3*x1*x2^2 - 22432*x0^3*x1*x2 - 19530*x0^3*x1 + 18328*x0^3*x2^7 - 10208*x0^3*x2^6 - 30284*x0^3*x2^5 - 10813*x0^3*x2^4 - 1051*x0^3*x2^3 - 26322*x0^3*x2^2 + 2971*x0^3*x2 - 4953*x0^3 - 7816*x0^2*x1^8 - 1902*x0^2*x1^7*x2 + 25248*x0^2*x1^7 - 11069*x0^2*x1^6*x2^2 + 32150*x0^2*x1^6*x2 + 25608*x0^2*x1^6 - 5281*x0^2*x1^5*x2^3 + 29789*x0^2*x1^5*x2^2 + 1238*x0^2*x1^5*x2 + 13900*x0^2*x1^5 + 18013*x0^2*x1^4*x2^4 - 32603*x0^2*x1^4*x2^3 - 10153*x0^2*x1^4*x2^2 + 24064*x0^2*x1^4*x2 - 32531*x0^2*x1^4 - 8418*x0^2*x1^3*x2^5 + 7952*x0^2*x1^3*x2^4 - 29769*x0^2*x1^3*x2^3 + 5467*x0^2*x1^3*x2^2 + 1112*x0^2*x1^3*x2 + 26836*x0^2*x1^3 + 26098*x0^2*x1^2*x2^6 + 21585*x0^2*x1^2*x2^5 - 16751*x0^2*x1^2*x2^4 - 30264*x0^2*x1^2*x2^3 - 6023*x0^2*x1^2*x2^2 - 3556*x0^2*x1^2*x2 - 15916*x0^2*x1^2 - 692*x0^2*x1*x2^7 - 22372*x0^2*x1*x2^6 + 31298*x0^2*x1*x2^5 - 14211*x0^2*x1*x2^4 - 11406*x0^2*x1*x2^3 - 25272*x0^2*x1*x2^2 - 7759*x0^2*x1*x2 - 2045*x0^2*x1 + 2344*x0^2*x2^8 - 32051*x0^2*x2^7 - 9410*x0^2*x2^6 + 14331*x0^2*x2^5 + 617*x0^2*x2^4 + 29671*x0^2*x2^3 + 29138*x0^2*x2^2 - 14421*x0^2*x2 - 19335*x0^2 + 3579*x0*x1^9 - 28011*x0*x1^8*x2 - 18868*x0*x1^8 + 26582*x0*x1^7*x2^2 + 5795*x0*x1^7*x2 + 8124*x0*x1^7 + 27627*x0*x1^6*x2^3 - 29607*x0*x1^6*x2^2 - 22058*x0*x1^6*x2 + 21807*x0*x1^6 + 19922*x0*x1^5*x2^4 + 20675*x0*x1^5*x2^3 + 27904*x0*x1^5*x2^2 + 9184*x0*x1^5*x2 - 28242*x0*x1^5 - 25689*x0*x1^4*x2^5 + 6364*x0*x1^4*x2^4 - 10599*x0*x1^4*x2^3 - 32615*x0*x1^4*x2^2 - 1093*x0*x1^4*x2 - 8567*x0*x1^4 + 30629*x0*x1^3*x2^6 - 4977*x0*x1^3*x2^5 + 10193*x0*x1^3*x2^4 + 12720*x0*x1^3*x2^3 - 677*x0*x1^3*x2^2 - 18111*x0*x1^3*x2 + 26753*x0*x1^3 + 30785*x0*x1^2*x2^7 + 14506*x0*x1^2*x2^6 + 11953*x0*x1^2*x2^5 + 20371*x0*x1^2*x2^4 + 2948*x0*x1^2*x2^3 - 19956*x0*x1^2*x2^2 + 26853*x0*x1^2*x2 + 26334*x0*x1^2 + 7823*x0*x1*x2^8 - 21766*x0*x1*x2^7 + 30645*x0*x1*x2^6 - 6661*x0*x1*x2^5 - 11887*x0*x1*x2^4 - 27768*x0*x1*x2^3 + 26192*x0*x1*x2^2 + 18411*x0*x1*x2 - 21656*x0*x1 + 18440*x0*x2^9 + 22528*x0*x2^8 - 26506*x0*x2^7 + 7822*x0*x2^6 + 29716*x0*x2^5 - 20664*x0*x2^4 - 10672*x0*x2^3 + 7638*x0*x2^2 - 7327*x0*x2 + 32525*x0 + 7459*x1^10 - 2836*x1^9*x2 + 30104*x1^9 - 5865*x1^8*x2^2 + 21775*x1^8*x2 + 1600*x1^8 + 18160*x1^7*x2^3 + 19968*x1^7*x2^2 + 16947*x1^7*x2 - 1802*x1^7 - 11767*x1^6*x2^4 + 22622*x1^6*x2^3 - 29767*x1^6*x2^2 + 32389*x1^6*x2 + 761*x1^6 + 5271*x1^5*x2^5 + 27422*x1^5*x2^4 + 22266*x1^5*x2^3 - 20423*x1^5*x2^2 - 9811*x1^5*x2 - 10437*x1^5 - 1166*x1^4*x2^6 + 16693*x1^4*x2^5 + 15110*x1^4*x2^4 + 16687*x1^4*x2^3 + 10868*x1^4*x2^2 + 4197*x1^4*x2 - 12100*x1^4 - 1691*x1^3*x2^7 - 19213*x1^3*x2^6 - 810*x1^3*x2^5 + 19577*x1^3*x2^4 - 2143*x1^3*x2^3 + 22081*x1^3*x2^2 + 21590*x1^3*x2 + 32363*x1^3 + 29734*x1^2*x2^8 + 3232*x1^2*x2^7 + 17203*x1^2*x2^6 + 20808*x1^2*x2^5 + 9097*x1^2*x2^4 + 12492*x1^2*x2^3 + 19191*x1^2*x2^2 - 9490*x1^2*x2 - 2848*x1^2 + 29970*x1*x2^9 - 4485*x1*x2^8 - 31058*x1*x2^7 + 8711*x1*x2^6 + 14343*x1*x2^5 - 28836*x1*x2^4 + 18185*x1*x2^3 - 17771*x1*x2^2 + 16063*x1*x2 + 24560*x1 + 10560*x2^10 - 25503*x2^9 + 2950*x2^8 - 19684*x2^7 + 21903*x2^6 + 31438*x2^5 - 29787*x2^4 - 30434*x2^3 - 1328*x2^2 + 30087*x2 + 13955",
 "-18123*x0^10 + 9236*x0^9*x1 - 16587*x0^9*x2 - 31326*x0^9 + 9588*x0^8*x1^2 + 16510*x0^8*x1*x2 + 18350*x0^8*x1 - 20417*x0^8*x2^2 + 27939*x0^8*x2 - 68*x0^8 + 31071*x0^7*x1^3 + 9736*x0^7*x1^2*x2 + 17414*x0^7*x1^2 + 28383*x0^7*x1*x2^2 - 29307*x0^7*x1*x2 + 15090*x0^7*x1 - 25095*x0^7*x2^3 - 25190*x0^7*x2^2 + 18465*x0^7*x2 + 32472*x0^7 + 24633*x0^6*x1^4 - 7959*x0^6*x1^3*x2 - 5077*x0^6*x1^3 - 10424*x0^6*x1^2*x2^2 - 23754*x0^6*x1^2*x2 - 23570*x0^6*x1^2 - 2959*x0^6*x1*x2^3 - 10981*x0^6*x1*x2^2 + 8892*x0^6*x1*x2 + 14173*x0^6*x1 + 24927*x0^6*x2^4 + 21547*x0^6*x2^3 + 435*x0^6*x2^2 + 9399*x0^6*x2 - 5009*x0^6 - 4258*x0^5*x1^5 - 11620*x0^5*x1^4*x2 + 1708*x0^5*x1^4 - 15055*x0^5*x1^3*x2^2 - 12645*x0^5*x1^3*x2 + 231*x0^5*x1^3 + 8778*x0^5*x1^2*x2^3 - 3068*x0^5*x1^2*x2^2 - 15690*x0^5*x1^2*x2 + 14055*x0^5*x1^2 - 24854*x0^5*x1*x2^4 + 30832*x0^5*x1*x2^3 - 23391*x0^5*x1*x2^2 + 17050*x0^5*x1*x2 - 1096*x0^5*x1 - 14110*x0^5*x2^5 + 6600*x0^5*x2^4 - 8936*x0^5*x2^3 - 7033*x0^5*x2^2 - 8071*x0^5*x2 + 21196*x0^5 - 13415*x0^4*x1^6 + 20643*x0^4*x1^5*x2 + 6814*x0^4*x1^5 + 9035*x0^4*x1^4*x2^2 + 20587*x0^4*x1^4*x2 + 11689*x0^4*x1^4 - 13280*x0^4*x1^3*x2^3 - 1159*x0^4*x1^3*x2^2 + 15949*x0^4*x1^3*x2 + 16134*x0^4*x1^3 - 9365*x0^4*x1^2*x2^4 - 9864*x0^4*x1^2*x2^3 - 11264*x0^4*x1^2*x2^2 - 24526*x0^4*x1^2*x2 + 14991*x0^4*x1^2 - 14233*x0^4*x1*x2^5 - 2828*x0^4*x1*x2^4 - 17244*x0^4*x1*x2^3 + 26566*x0^4*x1*x2^2 - 21630*x0^4*x1*x2 + 10696*x0^4*x1 + 16790*x0^4*x2^6 - 27170*x0^4*x2^5 - 19687*x0^4*x2^4 + 30516*x0^4*x2^3 - 15247*x0^4*x2^2 + 4822*x0^4*x2 - 24577*x0^4 + 19441*x0^3*x1^7 - 27787*x0^3*x1^6*x2 - 20949*x0^3*x1^6 - 7927*x0^3*x1^5*x2^2 - 25468*x0^3*x1^5*x2 + 25624*x0^3*x1^5 - 26381*x0^3*x1^4*x2^3 + 27566*x0^3*x1^4*x2^2 + 3355*x0^3*x1^4*x2 + 28664*x0^3*x1^4 - 4638*x0^3*x1^3*x2^4 - 1026*x0^3*x1^3*x2^3 - 8190*x0^3*x1^3*x2^2 - 19895*x0^3*x1^3*x2 - 2187*x0^3*x1^3 - 21085*x0^3*x1^2*x2^5 + 14863*x0^3*x1^2*x2^4 + 23569*x0^3*x1^2*x2^3 - 27100*x0^3*x1^2*x2^2 + 31234*x0^3*x1^2*x2 + 24836*x0^3*x1^2 - 13430*x0^3*x1*x2^6 - 12541*x0^3*x1*x2^5 - 17927*x0^3*x1*x2^4 + 13006*x0^3*x1*x2^3 - 12249*x0^3*x1*x2^2 + 15696*x0^3*x1*x2 - 24365*x0^3*x1 - 31062*x0^3*x2^7 - 28518*x0^3*x2^6 + 26428*x0^3*x2^5 + 21232*x0^3*x2^4 + 17969*x0^3*x2^3 - 32631*x0^3*x2^2 + 14508*x0^3*x2 - 11382*x0^3 - 11845*x0^2*x1^8 - 26914*x0^2*x1^7*x2 + 23672*x0^2*x1^7 + 563*x0^2*x1^6*x2^2 + 9165*x0^2*x1^6*x2 + 7970*x0^2*x1^6 + 2204*x0^2*x1^5*x2^3 + 6043*x0^2*x1^5*x2^2 + 8283*x0^2*x1^5*x2 + 23269*x0^2*x1^5 + 20833*x0^2*x1^4*x2^4 - 187*x0^2*x1^4*x2^3 + 729*x0^2*x1^4*x2^2 - 15501*x0^2*x1^4*x2 - 20667*x0^2*x1^4 - 12060*x0^2*x1^3*x2^5 + 28003*x0^2*x1^3*x2^4 + 237*x0^2*x1^3*x2^3 + 32324*x0^2*x1^3*x2^2 - 9184*x0^2*x1^3*x2 + 8809*x0^2*x1^3 + 24534*x0^2*x1^2*x2^6 + 17614*x0^2*x1^2*x2^5 - 7685*x0^2*x1^2*x2^4 - 13137*x0^2*x1^2*x2^3 + 19789*x0^2*x1^2*x2^2 - 3744*x0^2*x1^2*x2 + 6669*x0^2*x1^2 + 4362*x0^2*x1*x2^7 + 31277*x0^2*x1*x2^6 + 6984*x0^2*x1*x2^5 + 16128*x0^2*x1*x2^4 + 22022*x0^2*x1*x2^3 - 31193*x0^2*x1*x2^2 + 28348*x0^2*x1*x2 - 31953*x0^2*x1 + 13573*x0^2*x2^8 - 5127*x0^2*x2^7 + 24509*x0^2*x2^6 - 17586*x0^2*x2^5 + 289*x0^2*x2^4 + 17747*x0^2*x2^3 + 11732*x0^2*x2^2 - 14891*x0^2*x2 + 805*x0^2 + 19790*x0*x1^9 - 14004*x0*x1^8*x2 + 16069*x0*x1^8 - 23773*x0*x1^7*x2^2 + 4169*x0*x1^7*x2 - 25444*x0*x1^7 + 31414*x0*x1^6*x2^3 + 3875*x0*x1^6*x2^2 - 27714*x0*x1^6*x2 + 26581*x0*x1^6 + 13725*x0*x1^5*x2^4 - 7679*x0*x1^5*x2^3 - 14855*x0*x1^5*x2^2 + 24755*x0*x1^5*x2 - 4318*x0*x1^5 - 11885*x0*x1^4*x2^5 - 29405*x0*x1^4*x2^4 + 30174*x0*x1^4*x2^3 - 23009*x0*x1^4*x2^2 - 27409*x0*x1^4*x2 + 26829*x0*x1^4 + 7619*x0*x1^3*x2^6 + 24890*x0*x1^3*x2^5 + 16938*x0*x1^3*x2^4 - 25781*x0*x1^3*x2^3 - 5180*x0*x1^3*x2^2 + 16548*x0*x1^3*x2 + 27094*x0*x1^3 + 10386*x0*x1^2*x2^7 + 14375*x0*x1^2*x2^6 + 21157*x0*x1^2*x2^5 - 22735*x0*x1^2*x2^4 + 9230*x0*x1^2*x2^3 - 8414*x0*x1^2*x2^2 + 28800*x0*x1^2*x2 - 28811*x0*x1^2 - 1953*x0*x1*x2^8 - 24810*x0*x1*x2^7 - 5727*x0*x1*x2^6 + 6039*x0*x1*x2^5 - 9384*x0*x1*x2^4 + 141*x0*x1*x2^3 - 27617*x0*x1*x2^2 - 28642*x0*x1*x2 - 11349*x0*x1 + 9018*x0*x2^9 + 7196*x0*x2^8 - 32672*x0*x2^7 + 31797*x0*x2^6 - 30514*x0*x2^5 + 3674*x0*x2^4 + 12664*x0*x2^3 + 1315*x0*x2^2 - 15400*x0*x2 - 10377*x0 + 17409*x1^10 + 29236*x1^9*x2 + 8943*x1^9 + 30324*x1^8*x2^2 - 16496*x1^8*x2 - 9133*x1^8 - 4565*x1^7*x2^3 + 22670*x1^7*x2^2 + 20607*x1^7*x2 - 15841*x1^7 - 31724*x1^6*x2^4 - 1416*x1^6*x2^3 - 26573*x1^6*x2^2 + 17293*x1^6*x2 - 22674*x1^6 - 8178*x1^5*x2^5 + 7721*x1^5*x2^4 + 18321*x1^5*x2^3 + 15522*x1^5*x2^2 + 8493*x1^5*x2 - 18778*x1^5 - 31487*x1^4*x2^6 + 31976*x1^4*x2^5 + 27691*x1^4*x2^4 - 19576*x1^4*x2^3 + 19584*x1^4*x2^2 + 8873*x1^4*x2 + 11371*x1^4 - 405*x1^3*x2^7 + 32232*x1^3*x2^6 + 4058*x1^3*x2^5 - 13275*x1^3*x2^4 + 14420*x1^3*x2^3 - 6407*x1^3*x2^2 - 19181*x1^3*x2 + 16162*x1^3 - 28447*x1^2*x2^8 + 6824*x1^2*x2^7 - 5219*x1^2*x2^6 + 20882*x1^2*x2^5 - 23636*x1^2*x2^4 + 20807*x1^2*x2^3 + 22716*x1^2*x2^2 + 1837*x1^2*x2 + 30529*x1^2 - 32172*x1*x2^9 + 28631*x1*x2^8 + 20067*x1*x2^7 + 9972*x1*x2^6 - 739*x1*x2^5 - 19730*x1*x2^4 + 10766*x1*x2^3 - 12503*x1*x2^2 - 24447*x1*x2 + 9192*x1 - 9242*x2^10 + 14775*x2^9 + 26527*x2^8 - 4075*x2^7 + 9151*x2^6 + 7083*x2^5 + 29792*x2^4 + 29540*x2^3 + 24868*x2^2 - 31206*x2 - 27194",
 "-14283*x0^10 + 2597*x0^9*x1 + 24750*x0^9*x2 + 7067*x0^9 + 11106*x0^8*x1^2 - 21260*x0^8*x1*x2 + 14164*x0^8*x1 - 25247*x0^8*x2^2 - 27377*x0^8*x2 + 15299*x0^8 - 1440*x0^7*x1^3 + 1104*x0^7*x1^2*x2 + 2024*x0^7*x1^2 + 13267*x0^7*x1*x2^2 - 20398*x0^7*x1*x2 + 2762*x0^7*x1 - 146*x0^7*x2^3 + 12678*x0^7*x2^2 - 15182*x0^7*x2 - 24263*x0^7 - 24027*x0^6*x1^4 + 11117*x0^6*x1^3*x2 + 13947*x0^6*x1^3 + 15790*x0^6*x1^2*x2^2 + 1081*x0^6*x1^2*x2 + 26977*x0^6*x1^2 - 25976*x0^6*x1*x2^3 - 29227*x0^6*x1*x2^2 + 24640*x0^6*x1*x2 - 19409*x0^6*x1 - 28420*x0^6*x2^4 - 26302*x0^6*x2^3 + 7860*x0^6*x2^2 - 21114*x0^6*x2 + 542*x0^6 - 25230*x0^5*x1^5 + 18182*x0^5*x1^4*x2 - 30560*x0^5*x1^4 - 26328*x0^5*x1^3*x2^2 + 22389*x0^5*x1^3*x2 + 29798*x0^5*x1^3 + 929*x0^5*x1^2*x2^3 - 6674*x0^5*x1^2*x2^2 - 1001*x0^5*x1^2*x2 - 6683*x0^5*x1^2 - 32580*x0^5*x1*x2^4 + 24542*x0^5*x1*x2^3 + 21718*x0^5*x1*x2^2 - 9936*x0^5*x1*x2 - 21487*x0^5*x1 - 24160*x0^5*x2^5 - 26914*x0^5*x2^4 + 20738*x0^5*x2^3 + 20564*x0^5*x2^2 + 12778*x0^5*x2 - 19593*x0^5 - 5769*x0^4*x1^6 + 25967*x0^4*x1^5*x2 - 11013*x0^4*x1^5 + 29947*x0^4*x1^4*x2^2 + 10492*x0^4*x1^4*x2 + 28196*x0^4*x1^4 - 13363*x0^4*x1^3*x2^3 + 25552*x0^4*x1^3*x2^2 - 18606*x0^4*x1^3*x2 + 557*x0^4*x1^3 + 7512*x0^4*x1^2*x2^4 - 8014*x0^4*x1^2*x2^3 + 19499*x0^4*x1^2*x2^2 + 14450*x0^4*x1^2*x2 + 14319*x0^4*x1^2 + 16672*x0^4*x1*x2^5 + 226*x0^4*x1*x2^4 - 16110*x0^4*x1*x2^3 + 17135*x0^4*x1*x2^2 + 2094*x0^4*x1*x2 + 1989*x0^4*x1 - 13723*x0^4*x2^6 - 8501*x0^4*x2^5 - 828*x0^4*x2^4 - 2615*x0^4*x2^3 - 9802*x0^4*x2^2 - 18269*x0^4*x2 - 3859*x0^4 - 2224*x0^3*x1^7 - 16563*x0^3*x1^6*x2 + 23702*x0^3*x1^6 + 15645*x0^3*x1^5*x2^2 - 22814*x0^3*x1^5*x2 - 2684*x0^3*x1^5 - 24764*x0^3*x1^4*x2^3 - 20142*x0^3*x1^4*x2^2 + 6715*x0^3*x1^4*x2 - 10625*x0^3*x1^4 + 16448*x0^3*x1^3*x2^4 + 32547*x0^3*x1^3*x2^3 - 6320*x0^3*x1^3*x2^2 + 23628*x0^3*x1^3*x2 - 19012*x0^3*x1^3 - 8949*x0^3*x1^2*x2^5 - 24772*x0^3*x1^2*x2^4 - 20838*x0^3*x1^2*x2^3 - 28493*x0^3*x1^2*x2^2 + 19265*x0^3*x1^2*x2 - 22154*x0^3*x1^2 - 18471*x0^3*x1*x2^6 - 8479*x0^3*x1*x2^5 - 24177*x0^3*x1*x2^4 + 3046*x0^3*x1*x2^3 + 885*x0^3*x1*x2^2 + 231*x0^3*x1*x2 + 24507*x0^3*x1 + 5123*x0^3*x2^7 + 15608*x0^3*x2^6 - 3828*x0^3*x2^5 - 20300*x0^3*x2^4 + 27245*x0^3*x2^3 - 21014*x0^3*x2^2 - 2651*x0^3*x2 - 22399*x0^3 + 17160*x0^2*x1^8 - 24659*x0^2*x1^7*x2 - 14227*x0^2*x1^7 + 22483*x0^2*x1^6*x2^2 - 8413*x0^2*x1^6*x2 + 17066*x0^2*x1^6 - 891*x0^2*x1^5*x2^3 + 8781*x0^2*x1^5*x2^2 + 20554*x0^2*x1^5*x2 - 23661*x0^2*x1^5 - 10389*x0^2*x1^4*x2^4 - 6304*x0^2*x1^4*x2^3 + 30610*x0^2*x1^4*x2^2 + 290*x0^2*x1^4*x2 - 22651*x0^2*x1^4 + 24730*x0^2*x1^3*x2^5 - 28380*x0^2*x1^3*x2^4 - 9180*x0^2*x1^3*x2^3 - 17861*x0^2*x1^3*x2^2 - 13326*x0^2*x1^3*x2 + 9741*x0^2*x1^3 + 18053*x0^2*x1^2*x2^6 - 28363*x0^2*x1^2*x2^5 + 10195*x0^2*x1^2*x2^4 - 1521*x0^2*x1^2*x2^3 - 9472*x0^2*x1^2*x2^2 + 3947*x0^2*x1^2*x2 - 29191*x0^2*x1^2 + 9786*x0^2*x1*x2^7 - 17785*x0^2*x1*x2^6 + 10416*x0^2*x1*x2^5 - 14440*x0^2*x1*x2^4 + 12661*x0^2*x1*x2^3 - 24841*x0^2*x1*x2^2 + 8454*x0^2*x1*x2 + 23446*x0^2*x1 - 13548*x0^2*x2^8 + 4110*x0^2*x2^7 + 20706*x0^2*x2^6 - 24596*x0^2*x2^5 + 27329*x0^2*x2^4 - 3123*x0^2*x2^3 + 24272*x0^2*x2^2 - 1753*x0^2*x2 - 25022*x0^2 - 9385*x0*x1^9 + 5300*x0*x1^8*x2 - 2398*x0*x1^8 + 30814*x0*x1^7*x2^2 - 22444*x0*x1^7*x2 - 6860*x0*x1^7 + 11433*x0*x1^6*x2^3 + 21979*x0*x1^6*x2^2 + 10259*x0*x1^6*x2 - 30360*x0*x1^6 - 13932*x0*x1^5*x2^4 + 19310*x0*x1^5*x2^3 + 16484*x0*x1^5*x2^2 + 28657*x0*x1^5*x2 + 6787*x0*x1^5 + 23467*x0*x1^4*x2^5 + 6171*x0*x1^4*x2^4 + 15386*x0*x1^4*x2^3 + 2163*x0*x1^4*x2^2 - 17328*x0*x1^4*x2 + 17868*x0*x1^4 - 29042*x0*x1^3*x2^6 - 20048*x0*x1^3*x2^5 + 585*x0*x1^3*x2^4 + 1048*x0*x1^3*x2^3 + 31766*x0*x1^3*x2^2 - 5861*x0*x1^3*x2 + 28758*x0*x1^3 + 10577*x0*x1^2*x2^7 + 28934*x0*x1^2*x2^6 + 24056*x0*x1^2*x2^5 - 14590*x0*x1^2*x2^4 + 17591*x0*x1^2*x2^3 + 25190*x0*x1^2*x2^2 - 592*x0*x1^2*x2 + 27190*x0*x1^2 + 19036*x0*x1*x2^8 - 2095*x0*x1*x2^7 - 32436*x0*x1*x2^6 - 2463*x0*x1*x2^5 + 6735*x0*x1*x2^4 + 24666*x0*x1*x2^3 - 24912*x0*x1*x2^2 + 7785*x0*x1*x2 + 28214*x0*x1 - 20511*x0*x2^9 + 10359*x0*x2^8 + 17595*x0*x2^7 - 14376*x0*x2^6 + 605*x0*x2^5 + 7992*x0*x2^4 + 31051*x0*x2^3 - 29613*x0*x2^2 + 21731*x0*x2 - 29254*x0 + 30992*x1^10 + 15796*x1^9*x2 - 27227*x1^9 + 18349*x1^8*x2^2 - 16469*x1^8*x2 - 13132*x1^8 + 10991*x1^7*x2^3 + 1300*x1^7*x2^2 + 28340*x1^7*x2 - 29975*x1^7 - 24042*x1^6*x2^4 - 1289*x1^6*x2^3 - 27573*x1^6*x2^2 - 6854*x1^6*x2 - 10376*x1^6 - 23020*x1^5*x2^5 + 30961*x1^5*x2^4 - 30405*x1^5*x2^3 + 12610*x1^5*x2^2 - 21493*x1^5*x2 - 8783*x1^5 - 14575*x1^4*x2^6 - 13242*x1^4*x2^5 - 14272*x1^4*x2^4 + 7683*x1^4*x2^3 + 10762*x1^4*x2^2 + 19022*x1^4*x2 + 105*x1^4 + 8693*x1^3*x2^7 + 27311*x1^3*x2^6 + 13345*x1^3*x2^5 - 1042*x1^3*x2^4 - 12930*x1^3*x2^3 - 25375*x1^3*x2^2 - 5494*x1^3*x2 - 8132*x1^3 - 16347*x1^2*x2^8 + 10842*x1^2*x2^7 - 14805*x1^2*x2^6 - 29012*x1^2*x2^5 + 5557*x1^2*x2^4 + 20662*x1^2*x2^3 - 14187*x1^2*x2^2 + 14763*x1^2*x2 + 15032*x1^2 - 29357*x1*x2^9 + 5708*x1*x2^8 + 32720*x1*x2^7 + 11455*x1*x2^6 + 13074*x1*x2^5 + 1871*x1*x2^4 + 15816*x1*x2^3 + 22464*x1*x2^2 - 7607*x1*x2 + 31640*x1 + 10663*x2^10 - 21591*x2^9 - 31540*x2^8 - 24729*x2^7 + 19048*x2^6 + 1038*x2^5 + 30524*x2^4 - 31193*x2^3 + 25485*x2^2 - 31020*x2 + 28777"        };
        const char *vars[] = {"x1", "x2"};
        
        compute_dixon_resultant_string(polys, 3, vars, 2, 65537);
    }
    
    test_case_1();  // No parameters
    test_case_2();
    return 0;
}
