// dixon_test.c

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include "dixon_flint.h"
#include "dixon_interface_flint.h"
#include "fq_mvpoly.h"
#include "fq_unified_interface.h"
#include "unified_mpoly_resultant.h"
#include "dixon_with_ideal_reduction.h"
#include "resultant_with_ideal_reduction.h"
#include "dixon_complexity.h"
#include "polynomial_system_solver.h"

#define UNWRAP(...) __VA_ARGS__
#define COUNT_ARGS(...) (sizeof((const char*[]){__VA_ARGS__}) / sizeof(const char*))

#define DIXON(polys, vars) \
    dixon((const char*[]){UNWRAP polys}, COUNT_ARGS polys, \
          (const char*[]){UNWRAP vars}, COUNT_ARGS vars, ctx)

#define DIXON_WITH_IDEAL(polys, vars) \
    dixon_with_ideal((const char*[]){UNWRAP polys}, COUNT_ARGS polys, \
          (const char*[]){UNWRAP vars}, COUNT_ARGS vars, ideal, ctx)

#define GET_ARG1(a, b) a
#define GET_ARG2(a, b) b
#define APPLY(macro, args) macro args
#define RESULTANT(polys, var) \
    bivariate_resultant(APPLY(GET_ARG1, (UNWRAP polys)), \
                       APPLY(GET_ARG2, (UNWRAP polys)), \
                       UNWRAP(UNWRAP var), ctx)

#define DIXON_COMPLEXITY(polys, vars) \
    dixon_complexity_auto((const char*[]){UNWRAP polys}, COUNT_ARGS polys, \
          (const char*[]){UNWRAP vars}, COUNT_ARGS vars, ctx)
#define MAX_COMPLEXITY(...) \
    extract_max_complexity((const char*[]){__VA_ARGS__}, COUNT_ARGS(__VA_ARGS__))


// ============= Math Utility Functions =============

// Calculate binomial coefficient C(n, k)
slong binomial_coefficient(slong n, slong k) {
    if (k > n || k < 0) return 0;
    if (k == 0 || k == n) return 1;
    
    // Use symmetry to optimize calculation
    if (k > n - k) k = n -  k;
    
    slong result = 1;
    for (slong i = 0; i < k; i++) {
        result = result * (n - i) / (i + 1);
    }
    return result;
}

// Calculate total number of possible monomials for given variables, parameters and degree
slong count_possible_monomials(slong nvars, slong npars, slong max_degree) {
    slong total_indeterminates = nvars + npars;
    return binomial_coefficient(total_indeterminates + max_degree, max_degree);
}

// ============= Polynomial Generation =============

typedef struct {
    slong *exponents;  // Combined exponents for all indeterminates
    slong total_degree;
} monomial_t;

// Generate all possible monomials with degree <= max_degree
void enumerate_all_monomials(monomial_t **monomials, slong *count, 
                            slong total_indeterminates, slong max_degree) {
    // Calculate maximum possible monomials
    slong max_possible = 1;
    for (slong d = 0; d <= max_degree; d++) {
        max_possible += binomial_coefficient(total_indeterminates + d - 1, d);
    }
    
    *monomials = (monomial_t*) malloc(max_possible * sizeof(monomial_t));
    *count = 0;
    
    // Recursive function to generate monomials
    void generate_recursive(slong *current_exp, slong pos, slong remaining_degree) {
        if (pos == total_indeterminates) {
            if (remaining_degree == 0) {
                // Valid monomial found
                (*monomials)[*count].exponents = (slong*) malloc(total_indeterminates * sizeof(slong));
                memcpy((*monomials)[*count].exponents, current_exp, total_indeterminates * sizeof(slong));
                
                // Calculate total degree
                slong total_deg = 0;
                for (slong i = 0; i < total_indeterminates; i++) {
                    total_deg += current_exp[i];
                }
                (*monomials)[*count].total_degree = total_deg;
                (*count)++;
            }
            return;
        }
        
        // Try all possible degrees for current position
        for (slong deg = 0; deg <= remaining_degree; deg++) {
            current_exp[pos] = deg;
            generate_recursive(current_exp, pos + 1, remaining_degree - deg);
        }
    }
    
    // Generate monomials for each total degree
    slong *temp_exp = (slong*) calloc(total_indeterminates, sizeof(slong));
    for (slong degree = 0; degree <= max_degree; degree++) {
        generate_recursive(temp_exp, 0, degree);
    }
    free(temp_exp);
}

// Improved polynomial generation with better density control
void generate_random_polynomial(fq_mvpoly_t *poly, slong nvars, slong npars,
                               slong max_degree, double density_ratio,
                               const fq_nmod_ctx_t ctx, flint_rand_t state) {
    
    fq_mvpoly_init(poly, nvars, npars, ctx);
    
    slong total_indeterminates = nvars + npars;
    
    // Generate all possible monomials
    monomial_t *all_monomials;
    slong total_monomials;
    enumerate_all_monomials(&all_monomials, &total_monomials, total_indeterminates, max_degree);
    
    printf("    Total possible monomials: %ld\n", total_monomials);
    
    // Calculate target number of terms
    slong target_terms = (slong)(density_ratio * total_monomials);
    if (target_terms < 1) target_terms = 1;
    if (target_terms > total_monomials) target_terms = total_monomials;
    
    printf("    Target terms: %ld (%.1f%% density)\n", target_terms, (double)target_terms / total_monomials * 100);
    
    // Shuffle the monomial array to randomize selection
    for (slong i = total_monomials - 1; i > 0; i--) {
        slong j = n_randint(state, i + 1);
        // Swap monomials[i] and monomials[j]
        monomial_t temp = all_monomials[i];
        all_monomials[i] = all_monomials[j];
        all_monomials[j] = temp;
    }
    
    // Select first target_terms monomials and add them to polynomial
    for (slong i = 0; i < target_terms; i++) {
        monomial_t *mon = &all_monomials[i];
        
        // Split exponents into variable and parameter parts
        slong *var_exp = NULL;
        slong *par_exp = NULL;
        
        // Extract variable exponents (first nvars positions)
        if (nvars > 0) {
            int has_var = 0;
            for (slong j = 0; j < nvars; j++) {
                if (mon->exponents[j] > 0) {
                    has_var = 1;
                    break;
                }
            }
            if (has_var || (nvars > 0 && npars == 0)) {  // Include even if all zeros for pure variable case
                var_exp = (slong*) malloc(nvars * sizeof(slong));
                memcpy(var_exp, mon->exponents, nvars * sizeof(slong));
            }
        }
        
        // Extract parameter exponents (last npars positions)
        if (npars > 0) {
            int has_par = 0;
            for (slong j = nvars; j < total_indeterminates; j++) {
                if (mon->exponents[j] > 0) {
                    has_par = 1;
                    break;
                }
            }
            if (has_par || (npars > 0 && nvars == 0)) {  // Include even if all zeros for pure parameter case
                par_exp = (slong*) malloc(npars * sizeof(slong));
                memcpy(par_exp, mon->exponents + nvars, npars * sizeof(slong));
            }
        }
        
        // Generate random non-zero coefficient
        fq_nmod_t coeff;
        fq_nmod_init(coeff, ctx);
        do {
            fq_nmod_randtest(coeff, state, ctx);
        } while (fq_nmod_is_zero(coeff, ctx));
        
        // Add the term to polynomial
        fq_mvpoly_add_term(poly, var_exp, par_exp, coeff);
        
        // Cleanup
        fq_nmod_clear(coeff, ctx);
        if (var_exp) free(var_exp);
        if (par_exp) free(par_exp);
    }
    
    // Cleanup monomial list
    for (slong i = 0; i < total_monomials; i++) {
        free(all_monomials[i].exponents);
    }
    free(all_monomials);
    
    printf("    Generated polynomial with %ld terms (target: %ld, achieved density: %.1f%%)\n", 
           poly->nterms, target_terms, (double)poly->nterms / total_monomials * 100);
}

// Generate polynomial system with specified degrees and density
void generate_polynomial_system(fq_mvpoly_t **polys, slong nvars, slong npolys, 
                               slong npars, const slong *degrees,
                               double density_ratio,
                               const fq_nmod_ctx_t ctx, flint_rand_t state) {
    
    if (degrees == NULL) {
        printf("Error: degrees array cannot be NULL\n");
        return;
    }
    
    *polys = (fq_mvpoly_t*) malloc(npolys * sizeof(fq_mvpoly_t));
    if (!*polys) {
        printf("[ERROR] Failed to allocate polynomial array\n");
        return;
    }
    
    for (slong i = 0; i < npolys; i++) {
        
        slong max_degree = degrees[i];
        if (max_degree <= 0) {
            max_degree = 2;  // Default value
        }
        
        // First few polynomials include parameters
        slong poly_npars = npars;
        
        generate_random_polynomial(&(*polys)[i], nvars, poly_npars, 
                                  max_degree, density_ratio, ctx, state);
        
    }
    
}

// ============= Test Functions =============

void test_dixon_system(const char *test_name, slong nvars, slong npars,
                      ulong p, slong field_degree, const slong *degrees,
                      slong npolys, double density_ratio, flint_rand_t state) {
    printf("\n=== %s ===\n", test_name);
    printf("Field: GF(%lu^%ld), Variables: %ld, Parameters: %ld, Density: %.1f%%\n", 
           p, field_degree, nvars, npars, density_ratio * 100);
    
    // Print degree information
    printf("Polynomial degrees: [");
    for (slong i = 0; i < npolys; i++) {
        printf("%ld", degrees[i]);
        if (i < npolys - 1) printf(", ");
    }
    printf("]\n");
    
    // Calculate theoretical complexity
    slong total_theoretical_terms = 0;
    for (slong i = 0; i < npolys; i++) {
        slong possible = count_possible_monomials(nvars, npars, degrees[i]);
        printf("Polynomial %ld: max possible terms = %ld\n", i, possible);
        total_theoretical_terms += possible;
    }
    printf("System theoretical total terms: %ld\n", total_theoretical_terms);
    
    // Initialize field
    fq_nmod_ctx_t ctx;
    fmpz_t p_fmpz;
    fmpz_init(p_fmpz);
    fmpz_set_ui(p_fmpz, p);
    fq_nmod_ctx_init(ctx, p_fmpz, field_degree, "t");
    fmpz_clear(p_fmpz);
    
    // Generate polynomial system
    fq_mvpoly_t *polys;
    generate_polynomial_system(&polys, nvars, npolys, npars, 
                              degrees, density_ratio, ctx, state);
    
    // Print system information
    printf("\nPolynomial system:\n");
    slong actual_total_terms = 0;
    for (slong i = 0; i < npolys; i++) {
        printf("  p%ld (degree %ld): %ld terms", i, degrees[i], polys[i].nterms);
        if (polys[i].npars > 0) {
            printf(" (with %ld parameters)", polys[i].npars);
        }
        printf("\n");
        actual_total_terms += polys[i].nterms;
    }
    printf("System actual total terms: %ld (density: %.1f%%)\n", 
           actual_total_terms, (double)actual_total_terms / total_theoretical_terms * 100);
    
    // Compute Dixon resultant
    fq_mvpoly_t result;
    
    clock_t start = clock();
    fq_dixon_resultant(&result, polys, nvars, npars);
    clock_t end = clock();
    
    
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Computation time: %.3f seconds\n", elapsed);
    printf("Resultant: %ld terms", result.nterms);
    if (result.npars > 0) {
        printf(" with %ld parameters", result.npars);
    }
    printf("\n");
    
    if (result.nterms > 0 && result.nterms <= 10) {
        printf("  ");
        fq_mvpoly_print(&result, "R");
        printf("\n");
    }
    
    // Cleanup
    fq_mvpoly_clear(&result);
    for (slong i = 0; i < npolys; i++) {
        fq_mvpoly_clear(&polys[i]);
    }
    free(polys);
    fq_nmod_ctx_clear(ctx);
}

void test_xhash() {
    printf("=== Prime Field Polynomial System Dixon Resultant Implementation ===\n\n");    
    const char *ideal = "x7^3 = -19561*x4^3 + 2061*x4^2*x5 + 25073*x4^2*x6 + 19787*x4^2 + 779*x4*x5^2 + 31516*x4*x5*x6 - 17049*x4*x5 + 9065*x4*x6^2 - 14413*x4*x6 - 9964*x4 + 24021*x5^3 - 31858*x5^2*x6 - 20667*x5^2 - 19224*x5*x6^2 - 15430*x5*x6 + 19731*x5 + 31617*x6^3 + 6653*x6^2 - 13781*x6 - 15782,\
    x8^3 = 31617*x4^3 + 9065*x4^2*x5 - 19224*x4^2*x6 + 6653*x4^2 + 25073*x4*x5^2 + 31516*x4*x5*x6 - 14413*x4*x5 - 31858*x4*x6^2 - 15430*x4*x6 - 13781*x4 - 19561*x5^3 + 2061*x5^2*x6 + 19787*x5^2 + 779*x5*x6^2 - 17049*x5*x6 - 9964*x5 + 24021*x6^3 - 20667*x6^2 + 19731*x6 - 15782,\
    x9^3 = 24021*x4^3 - 31858*x4^2*x5 + 779*x4^2*x6 - 20667*x4^2 - 19224*x4*x5^2 + 31516*x4*x5*x6 - 15430*x4*x5 + 2061*x4*x6^2 - 17049*x4*x6 + 19731*x4 + 31617*x5^3 + 9065*x5^2*x6 + 6653*x5^2 + 25073*x5*x6^2 - 14413*x5*x6 - 13781*x5 - 19561*x6^3 + 19787*x6^2 - 9964*x6 - 15782";
    
    const char* f1 = "6653*x0^2 - 13781*x0 - x2^2 - 15782";
    const char* f2 = "20667*x0^2 + 19731*x0 - x3^2 - 15782";
    const char* f3 = "3*x1^3 + 16*x1^2*x2 - 24563*x1^2 - 15*x1*x2^2 + 3*x1*x2*x3 + 8202*x1*x2 + 16*x1*x3^2 - 8170*x1*x3 - 6753*x1 - 16*x2^3 + 8*x2^2*x3 - 24637*x2^2 - 9*x2*x3^2 + 27861*x2*x3 - 19866*x2 - 13*x3^3 + 26199*x3^2 - 26963*x3 - x4 - 32551";
    const char* f4 = "-4*x1^3 - 4*x1^2*x2 - 15*x1^2*x3 - 11483*x1^2 + 15*x1*x2^2 - 12*x1*x2*x3 + 16410*x1*x2 - x1*x3^2 + 19631*x1*x3 + 28842*x1 + 11*x2^3 + 6*x2^2*x3 - 3245*x2^2 - 4*x2*x3^2 - 26213*x2*x3 - 5476*x2 + 12*x3^3 + 11475*x3^2 + 12887*x3 - x5 - 8584";
    const char* f5 = "16*x1^3 - 11*x1^2*x2 + 15*x1^2*x3 + 6597*x1^2 + 4*x1*x2^2 - 15*x1*x2*x3 - 21332*x1*x2 + 14*x1*x3^2 + 18051*x1*x3 - 22387*x1 + 2*x2^3 - 5*x2^2*x3 + 8180*x2^2 - 8*x2*x3^2 - 27885*x2*x3 - 15205*x2 + 15*x3^3 - 21257*x3^2 - 3092*x3 - x6 - 25478";
    const char* f6 = "-19561*x4^3 + 2061*x4^2*x5 + 25073*x4^2*x6 + 19787*x4^2 + 779*x4*x5^2 + 31516*x4*x5*x6 - 17049*x4*x5 + 9065*x4*x6^2 - 14413*x4*x6 - 9964*x4 + 24021*x5^3 - 31858*x5^2*x6 - 20667*x5^2 - 19224*x5*x6^2 - 15430*x5*x6 + 19731*x5 + 31617*x6^3 + 6653*x6^2 - 13781*x6 - x7^3 - 15782";
    const char* f7 = "31617*x4^3 + 9065*x4^2*x5 - 19224*x4^2*x6 + 6653*x4^2 + 25073*x4*x5^2 + 31516*x4*x5*x6 - 14413*x4*x5 - 31858*x4*x6^2 - 15430*x4*x6 - 13781*x4 - 19561*x5^3 + 2061*x5^2*x6 + 19787*x5^2 + 779*x5*x6^2 - 17049*x5*x6 - 9964*x5 + 24021*x6^3 - 20667*x6^2 + 19731*x6 - x8^3 - 15782";
    const char* f8 = "24021*x4^3 - 31858*x4^2*x5 + 779*x4^2*x6 - 20667*x4^2 - 19224*x4*x5^2 + 31516*x4*x5*x6 - 15430*x4*x5 + 2061*x4*x6^2 - 17049*x4*x6 + 19731*x4 + 31617*x5^3 + 9065*x5^2*x6 + 6653*x5^2 + 25073*x5*x6^2 - 14413*x5*x6 - 13781*x5 - 19561*x6^3 + 19787*x6^2 - 9964*x6 - x9^3 - 15782";
    const char* f9 = "3*x7^3 + 16*x7^2*x8 - 24563*x7^2 - 15*x7*x8^2 + 3*x7*x8*x9 + 8202*x7*x8 + 16*x7*x9^2 - 8170*x7*x9 - 6753*x7 - 16*x8^3 + 8*x8^2*x9 - 24637*x8^2 - 9*x8*x9^2 + 27861*x8*x9 - 19866*x8 - 13*x9^3 + 26199*x9^2 - 26963*x9 - 32551";

    fq_nmod_ctx_t ctx;
    mp_limb_t prime = 65537;
    fmpz_t p;
    fmpz_init_set_ui(p, prime);
    fq_nmod_ctx_init(ctx, p, 1, "t");    
    
    printf("Prime field: GF(p) where p = 2^63 - 25 = %llu\n", p);
    printf("\n");
    
    clock_t total_start = clock();
    
    char* r1 = RESULTANT((f1, f2), ("x0"));
    char* r2 = DIXON((r1, f3, f4, f5), ("x1", "x2", "x3"));
    char* r3 = DIXON_WITH_IDEAL((f9, f8), ("x9"));
    char* r4 = DIXON_WITH_IDEAL((r3, f7), ("x8"));
    char* r5 = DIXON_WITH_IDEAL((r4, f6), ("x7"));
    char* r6 = DIXON_COMPLEXITY((r2, r5), ("x4"));

    printf("%s",r6);
    clock_t total_end = clock();
    double total_time = (double)(total_end - total_start) / CLOCKS_PER_SEC;
    
    printf("Total computation time: %.3f seconds\n", total_time);
    
    free(r1);
    free(r2);
    free(r3);
    free(r4);
    free(r5);
    free(r6);
    
    fq_nmod_ctx_clear(ctx);
    printf("\n=== Computation Complete ===\n");
} 

// Main test function
int test_dixon_resultant() {
    // Initialize random state
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, time(NULL) + getpid(), time(NULL) * getpid());
    
    printf("Dixon Resultant Test Suite - Comprehensive Testing\n");
    printf("=================================================\n");

    // ============= 1. BASIC FUNCTIONALITY TESTS =============
    printf("\n1. BASIC FUNCTIONALITY TESTS\n");
    printf("-----------------------------\n");
    
    // Simple quadratic tests
    slong deg_basic1[] = {2, 2, 2};
    test_dixon_system("Basic quadratic system GF(7)", 2, 1, 7, 1, deg_basic1, 3, 0.5, state);
    
    slong deg_basic2[] = {1, 1, 1};
    test_dixon_system("Linear system GF(11)", 2, 1, 11, 1, deg_basic2, 3, 1.0, state);
    
    slong deg_basic3[] = {3, 3, 3};
    test_dixon_system("Cubic system GF(13)", 2, 1, 13, 1, deg_basic3, 3, 0.3, state);

    // ============= 2. FIELD CHARACTERISTIC TESTS =============
    printf("\n2. FIELD CHARACTERISTIC TESTS\n");
    printf("------------------------------\n");
    
    // Small prime fields
    slong deg_field1[] = {2, 2, 2};
    test_dixon_system("Small prime GF(2)", 2, 1, 2, 1, deg_field1, 3, 0.8, state);
    test_dixon_system("Small prime GF(3)", 2, 1, 3, 1, deg_field1, 3, 0.8, state);
    test_dixon_system("Small prime GF(5)", 2, 1, 5, 1, deg_field1, 3, 0.8, state);
    
    // Extension fields of different degrees
    slong deg_ext[] = {2, 2, 2};
    test_dixon_system("Extension field GF(2^4)", 2, 1, 2, 4, deg_ext, 3, 0.7, state);
    test_dixon_system("Extension field GF(3^3)", 2, 1, 3, 3, deg_ext, 3, 0.7, state);
    test_dixon_system("Extension field GF(5^2)", 2, 1, 5, 2, deg_ext, 3, 0.7, state);
    
    // Large prime field
    slong deg_large[] = {2, 2, 2};
    test_dixon_system("Large prime GF(65537)", 2, 1, 65537, 1, deg_large, 3, 0.6, state);

    // ============= 3. POLYNOMIAL STRUCTURE TESTS =============
    printf("\n3. POLYNOMIAL STRUCTURE TESTS\n");
    printf("------------------------------\n");
    
    // Dense vs sparse polynomials
    slong deg_density[] = {3, 3, 3};
    test_dixon_system("Dense polynomials (90%)", 2, 1, 17, 1, deg_density, 3, 0.9, state);
    test_dixon_system("Medium density (50%)", 2, 1, 17, 1, deg_density, 3, 0.5, state);
    test_dixon_system("Sparse polynomials (10%)", 2, 1, 17, 1, deg_density, 3, 0.1, state);
    
    // Mixed degree systems
    slong deg_mixed1[] = {1, 2, 3};
    test_dixon_system("Mixed degrees [1,2,3]", 2, 1, 19, 1, deg_mixed1, 3, 0.7, state);
    
    slong deg_mixed2[] = {1, 3, 5};
    test_dixon_system("Mixed degrees [1,3,5]", 2, 1, 23, 1, deg_mixed2, 3, 0.6, state);
    
    slong deg_mixed3[] = {2, 4, 6};
    test_dixon_system("Mixed degrees [2,4,6]", 2, 1, 29, 1, deg_mixed3, 3, 0.4, state);

    // ============= 4. VARIABLE AND PARAMETER TESTS =============
    printf("\n4. VARIABLE AND PARAMETER TESTS\n");
    printf("--------------------------------\n");
    
    // Different numbers of variables
    slong deg_var1[] = {2, 2};
    test_dixon_system("Single variable system", 1, 2, 31, 1, deg_var1, 2, 0.8, state);
    
    slong deg_var3[] = {2, 2, 2, 2};
    test_dixon_system("Three variable system", 3, 1, 37, 1, deg_var3, 4, 0.5, state);
    
    // Different numbers of parameters
    slong deg_par0[] = {2, 2, 2};
    test_dixon_system("No parameters", 2, 0, 41, 1, deg_par0, 3, 0.7, state);
    
    slong deg_par2[] = {2, 2, 2};
    test_dixon_system("Two parameters", 2, 2, 43, 1, deg_par2, 3, 0.6, state);
    
    slong deg_par3[] = {2, 2, 2};
    test_dixon_system("Three parameters", 2, 3, 47, 1, deg_par3, 3, 0.5, state);

    // ============= 5. EDGE CASE TESTS =============
    printf("\n5. EDGE CASE TESTS\n");
    printf("------------------\n");
    
    // Minimal degree cases
    slong deg_min[] = {1, 1, 1};
    test_dixon_system("All linear (minimal case)", 2, 1, 53, 1, deg_min, 3, 1.0, state);
    
    // High degree but sparse
    slong deg_high_sparse[] = {8, 8, 8};
    test_dixon_system("High degree sparse", 2, 1, 59, 1, deg_high_sparse, 3, 0.05, state);
    
    // Unbalanced degrees
    slong deg_unbal[] = {1, 1, 6};
    test_dixon_system("Unbalanced degrees [1,1,6]", 2, 1, 61, 1, deg_unbal, 3, 0.7, state);

    // ============= 6. PERFORMANCE STRESS TESTS =============
    printf("\n6. PERFORMANCE STRESS TESTS\n");
    printf("----------------------------\n");
    
    // Medium complexity
    slong deg_med[] = {4, 4, 4};
    test_dixon_system("Medium complexity [4,4,4]", 2, 2, 67, 1, deg_med, 3, 0.3, state);
    
    // Higher complexity (be careful with computation time)
    slong deg_stress[] = {3, 4, 5};
    test_dixon_system("Stress test [3,4,5]", 2, 2, 71, 1, deg_stress, 3, 0.2, state);
    
    // Large field stress test
    slong deg_large_field[] = {3, 3, 3};
    test_dixon_system("Large field stress GF(2^12)", 2, 1, 2, 12, deg_large_field, 3, 0.4, state);

    // ============= 7. SPECIAL STRUCTURE TESTS =============
    printf("\n7. SPECIAL STRUCTURE TESTS\n");
    printf("---------------------------\n");
    
    // Binary field tests (important for cryptographic applications)
    slong deg_bin[] = {3, 3, 3};
    test_dixon_system("Binary field GF(2^8)", 2, 1, 2, 8, deg_bin, 3, 0.6, state);
    test_dixon_system("Binary field GF(2^16)", 2, 1, 2, 16, deg_bin, 3, 0.4, state);
    
    // Characteristic 3 tests
    slong deg_char3[] = {2, 3, 3};
    test_dixon_system("Characteristic 3 field", 2, 1, 3, 5, deg_char3, 3, 0.5, state);

    // ============= 8. REGRESSION TESTS =============
    printf("\n8. REGRESSION TESTS\n");
    printf("-------------------\n");
    
    // Previously problematic cases (if any were discovered)
    slong deg_reg1[] = {2, 2, 3};
    test_dixon_system("Regression case 1", 2, 1, 73, 1, deg_reg1, 3, 0.8, state);
    
    slong deg_reg2[] = {1, 4, 4};
    test_dixon_system("Regression case 2", 2, 2, 79, 1, deg_reg2, 3, 0.6, state);

    // ============= 9. BOUNDARY CONDITION TESTS =============
    printf("\n9. BOUNDARY CONDITION TESTS\n");
    printf("----------------------------\n");
    
    // Minimal polynomial system
    slong deg_boundary1[] = {1, 1};
    test_dixon_system("Minimal bivariate", 1, 1, 83, 1, deg_boundary1, 2, 1.0, state);
    
    // Maximum reasonable complexity for automated testing
    slong deg_boundary2[] = {2, 3, 4};
    test_dixon_system("Upper boundary test", 2, 1, 89, 1, deg_boundary2, 3, 0.3, state);

    // ============= 10. MULTI-FIELD COMPARISON =============
    printf("\n10. MULTI-FIELD COMPARISON\n");
    printf("--------------------------\n");
    
    // Same polynomial structure across different fields
    slong deg_comp[] = {2, 2, 3};
    printf("Comparing same structure across different fields:\n");
    test_dixon_system("  Field GF(97)", 2, 1, 97, 1, deg_comp, 3, 0.6, state);
    test_dixon_system("  Field GF(101)", 2, 1, 101, 1, deg_comp, 3, 0.6, state);
    test_dixon_system("  Field GF(2^7)", 2, 1, 2, 7, deg_comp, 3, 0.6, state);

    // ============= 11. SCALABILITY TESTS =============
    printf("\n11. SCALABILITY TESTS\n");
    printf("---------------------\n");
    
    // Test how the algorithm scales with increasing complexity
    printf("Scalability analysis:\n");
    
    slong deg_scale1[] = {2, 2, 2};
    test_dixon_system("  Scale test: degree 2", 2, 1, 103, 1, deg_scale1, 3, 0.5, state);
    
    slong deg_scale2[] = {3, 3, 3};
    test_dixon_system("  Scale test: degree 3", 2, 1, 103, 1, deg_scale2, 3, 0.5, state);
    
    slong deg_scale3[] = {4, 4, 4};
    test_dixon_system("  Scale test: degree 4", 2, 1, 103, 1, deg_scale3, 3, 0.3, state);

    // ============= 12. CORNER CASES =============
    printf("\n12. CORNER CASES\n");
    printf("----------------\n");
    
    // Very sparse high-degree polynomials
    slong deg_corner1[] = {10, 10, 10};
    test_dixon_system("Very sparse high degree", 2, 1, 107, 1, deg_corner1, 3, 0.01, state);
    
    // Many parameters, few variables
    slong deg_corner2[] = {2, 2, 2};
    test_dixon_system("Many parameters", 2, 5, 109, 1, deg_corner2, 3, 0.4, state);

    // Cleanup
    flint_rand_clear(state);
    flint_cleanup();
    
    return 0;
}

void test_dixon(){

    //test_dixon_resultant();
    //test_xhash();
    //test_gf28_conversion();
    //test_iterative_elimination();
    //test_polynomial_solver();
    run_unified_comparison();
}
