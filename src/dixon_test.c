// dixon_test.c

#include "dixon_test.h"

// ============= DEGREE CHECKING STRUCTURES =============

typedef struct {
    char test_name[256];
    slong nvars;
    slong npars;
    slong npolys;
    slong *input_degrees;
    slong resultant_degree;
    slong expected_degree_product;
    int degree_satisfies_bound;
} degree_check_result_t;

// Global storage for degree check results
static degree_check_result_t *degree_results = NULL;
static slong num_degree_results = 0;
static slong degree_results_capacity = 0;

// ============= DEGREE CHECKING FUNCTIONS =============

// Get maximum total degree of a polynomial
slong fq_mvpoly_max_total_degree(const fq_mvpoly_t *poly) {
    slong max_deg = 0;
    
    for (slong i = 0; i < poly->nterms; i++) {
        slong total_deg = 0;
        
        // Sum variable exponents
        if (poly->terms[i].var_exp) {
            for (slong j = 0; j < poly->nvars; j++) {
                total_deg += poly->terms[i].var_exp[j];
            }
        }
        
        // Sum parameter exponents
        if (poly->terms[i].par_exp) {
            for (slong j = 0; j < poly->npars; j++) {
                total_deg += poly->terms[i].par_exp[j];
            }
        }
        
        if (total_deg > max_deg) {
            max_deg = total_deg;
        }
    }
    
    return max_deg;
}

// Add degree check result to global storage
void add_degree_check_result(const char *test_name, slong nvars, slong npars, slong npolys,
                             const slong *input_degrees, slong resultant_degree) {
    // Ensure capacity
    if (num_degree_results >= degree_results_capacity) {
        degree_results_capacity = degree_results_capacity == 0 ? 10 : degree_results_capacity * 2;
        degree_results = (degree_check_result_t*) realloc(degree_results, 
                                                          degree_results_capacity * sizeof(degree_check_result_t));
    }
    
    degree_check_result_t *result = &degree_results[num_degree_results];
    
    // Copy test name
    strncpy(result->test_name, test_name, 255);
    result->test_name[255] = '\0';
    
    result->nvars = nvars;
    result->npars = npars;
    result->npolys = npolys;
    
    // Copy input degrees
    result->input_degrees = (slong*) malloc(npolys * sizeof(slong));
    memcpy(result->input_degrees, input_degrees, npolys * sizeof(slong));
    
    result->resultant_degree = resultant_degree;
    
    // Calculate expected degree product
    result->expected_degree_product = 1;
    for (slong i = 0; i < npolys; i++) {
        result->expected_degree_product *= input_degrees[i];
    }
    
    // Check if degree bound is satisfied
    result->degree_satisfies_bound = (resultant_degree <= result->expected_degree_product);
    
    num_degree_results++;
}

// Print comprehensive summary of all degree checks
void print_degree_check_summary(void) {
    printf("\n");
    printf("╔════════════════════════════════════════════════════════════════════════════╗\n");
    printf("║                    DIXON RESULTANT DEGREE ANALYSIS SUMMARY                 ║\n");
    printf("╚════════════════════════════════════════════════════════════════════════════╝\n");
    printf("\n");
    
    if (num_degree_results == 0) {
        printf("No degree check results collected.\n");
        return;
    }
    
    slong satisfied_count = 0;
    slong violated_count = 0;
    
    for (slong i = 0; i < num_degree_results; i++) {
        degree_check_result_t *res = &degree_results[i];
        
        printf("─────────────────────────────────────────────────────────────────────────────\n");
        printf("Test %ld: %s\n", i + 1, res->test_name);
        printf("  System: %ld vars, %ld params, %ld polynomials\n", 
               res->nvars, res->npars, res->npolys);
        
        printf("  Input degrees: [");
        for (slong j = 0; j < res->npolys; j++) {
            printf("%ld", res->input_degrees[j]);
            if (j < res->npolys - 1) printf(", ");
        }
        printf("]\n");
        
        printf("  Expected bound (∏degrees): %ld\n", res->expected_degree_product);
        printf("  Resultant degree: %ld\n", res->resultant_degree);
        
        double ratio = res->expected_degree_product > 0 ? 
                      (double)res->resultant_degree / res->expected_degree_product : 0.0;
        printf("  Ratio (actual/bound): %.3f\n", ratio);
        
        if (res->degree_satisfies_bound) {
            printf("  Status: ✓ PASSED (deg ≤ bound)\n");
            satisfied_count++;
        } else {
            printf("  Status: ✗ VIOLATED (deg > bound)\n");
            violated_count++;
        }
    }
    
    printf("═════════════════════════════════════════════════════════════════════════════\n");
    printf("OVERALL SUMMARY:\n");
    printf("  Total tests: %ld\n", num_degree_results);
    printf("  Degree bound satisfied: %ld (%.1f%%)\n", 
           satisfied_count, 100.0 * satisfied_count / num_degree_results);
    printf("  Degree bound violated: %ld (%.1f%%)\n", 
           violated_count, 100.0 * violated_count / num_degree_results);
    
    // Calculate statistics
    slong min_resultant_deg = LONG_MAX;
    slong max_resultant_deg = 0;
    double avg_ratio = 0.0;
    
    for (slong i = 0; i < num_degree_results; i++) {
        if (degree_results[i].resultant_degree < min_resultant_deg) {
            min_resultant_deg = degree_results[i].resultant_degree;
        }
        if (degree_results[i].resultant_degree > max_resultant_deg) {
            max_resultant_deg = degree_results[i].resultant_degree;
        }
        if (degree_results[i].expected_degree_product > 0) {
            avg_ratio += (double)degree_results[i].resultant_degree / 
                        degree_results[i].expected_degree_product;
        }
    }
    avg_ratio /= num_degree_results;
    
    printf("\n  Statistics:\n");
    printf("    Min resultant degree: %ld\n", min_resultant_deg);
    printf("    Max resultant degree: %ld\n", max_resultant_deg);
    printf("    Average ratio (deg/bound): %.3f\n", avg_ratio);
    printf("═════════════════════════════════════════════════════════════════════════════\n");
}

// Cleanup degree check results
void cleanup_degree_check_results(void) {
    if (degree_results) {
        for (slong i = 0; i < num_degree_results; i++) {
            if (degree_results[i].input_degrees) {
                free(degree_results[i].input_degrees);
            }
        }
        free(degree_results);
        degree_results = NULL;
    }
    num_degree_results = 0;
    degree_results_capacity = 0;
}


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

// ========== DEGREE CHECKING CODE ==========
    
    // Calculate resultant degree
    slong resultant_degree = fq_mvpoly_max_total_degree(&result);
    printf("Resultant total degree: %ld\n", resultant_degree);
    
    // Calculate expected degree bound (product of input degrees)
    slong expected_product = 1;
    for (slong i = 0; i < npolys; i++) {
        expected_product *= degrees[i];
    }
    printf("Expected degree bound (product): %ld\n", expected_product);
    
    // Check if bound is satisfied
    if (resultant_degree <= expected_product) {
        printf("✓ Degree check PASSED: deg(resultant) ≤ ∏deg(inputs)\n");
    } else {
        printf("✗ Degree check FAILED: deg(resultant) > ∏deg(inputs)\n");
    }
    
    // Store result for summary
    add_degree_check_result(test_name, nvars, npars, npolys, degrees, resultant_degree);
    
    // ========== END DEGREE CHECKING CODE ==========
    
    if (result.nterms > 0 && result.nterms <= 10) {
        printf("  ");
        fq_mvpoly_print(&result, "R");
        printf("\n");
    }
    
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
    printf("\n1. BASIC FUNCTIONALITY TESTS (20 repetitions each)\n");
    printf("--------------------------------------------------\n");
    
    for (int rep = 0; rep < 20; rep++) {
        printf("\n--- Repetition %d/20 ---\n", rep + 1);
        
        slong deg_basic1[] = {3, 4, 5};
        test_dixon_system("Mixed degree (3,4,5) system", 2, 1, 65537, 1, deg_basic1, 3, 0.9, state);
        
        slong deg_basic2[] = {5, 6, 6};
        test_dixon_system("Quintic-sextic (5,6,6) system", 2, 1, 65537, 1, deg_basic2, 3, 1.0, state);
        
        slong deg_basic3[] = {5, 5, 5};
        test_dixon_system("Uniform quintic (5,5,5) system", 2, 1, 65537, 1, deg_basic3, 3, 0.6, state);
        
        slong deg_corner1[] = {2, 2, 2, 2, 2, 2};
        test_dixon_system("Six quadratic equations (2,2,2,2,2,2)", 5, 1, 65537, 1, deg_corner1, 6, 0.7, state);
        
        slong deg_corner2[] = {2, 2, 3, 3};
        test_dixon_system("Mixed quadratic-cubic with parameters (2,2,3,3)", 3, 2, 65537, 1, deg_corner2, 4, 0.8, state);
    }
    


    // Cleanup
    flint_rand_clear(state);
    flint_cleanup();
    
    return 0;
}

void test_dixon(){

    test_dixon_resultant();
    print_degree_check_summary();
    cleanup_degree_check_results();
    
    //test_xhash();
    //test_gf28_conversion();
    //test_iterative_elimination();
    //test_polynomial_solver();
    //run_unified_comparison();
}
