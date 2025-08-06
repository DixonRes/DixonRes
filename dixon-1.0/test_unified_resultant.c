/*
 * test_unified_resultant.c - Clean version with timeout support
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <signal.h>
#include <setjmp.h>
#include <unistd.h>  /* For alarm() */
#include <flint/flint.h>
#include <flint/nmod_mpoly.h>
#include <flint/fq_nmod_mpoly.h>
#include <flint/fmpz.h>
#include <flint/mpoly.h>
#include "unified_mpoly_resultant.h"

/* Global timeout handling */
static jmp_buf timeout_env;
static volatile int timeout_occurred = 0;

void timeout_handler(int sig) {
    timeout_occurred = 1;
    longjmp(timeout_env, 1);
}

/* Set timeout in seconds */
void set_timeout(int seconds) {
    timeout_occurred = 0;
    signal(SIGALRM, timeout_handler);
    alarm(seconds);
}

void clear_timeout() {
    alarm(0);
    signal(SIGALRM, SIG_DFL);
}

/* Test configuration for dense polynomials - renamed to avoid conflict */
typedef struct {
    const char *name;
    slong nvars;
    slong length1;
    slong length2;
    slong exp_bound;
    int var_index;
} resultant_test_case_t;

/* Generate random fq_nmod_mpoly polynomial */
void fq_nmod_mpoly_randtest_dense(fq_nmod_mpoly_t poly, flint_rand_t state,
                                 slong length, slong exp_bound,
                                 const fq_nmod_mpoly_ctx_t ctx) {
    fq_nmod_mpoly_zero(poly, ctx);
    
    slong nvars = ctx->minfo->nvars;
    ulong *exp = (ulong *)calloc(nvars, sizeof(ulong));
    fq_nmod_t coeff;
    fq_nmod_init(coeff, ctx->fqctx);
    
    for (slong i = 0; i < length; i++) {
        /* Random exponents - keep them relatively small for dense polynomials */
        for (slong j = 0; j < nvars; j++) {
            exp[j] = n_randint(state, exp_bound);
        }
        
        /* Random non-zero coefficient */
        fq_nmod_randtest_not_zero(coeff, state, ctx->fqctx);
        fq_nmod_mpoly_set_coeff_fq_nmod_ui(poly, coeff, exp, ctx);
    }
    
    free(exp);
    fq_nmod_clear(coeff, ctx->fqctx);
}

/* Test resultant with timeout */
int test_resultant_with_timeout(const char *field_desc, const fq_nmod_ctx_t fq_ctx,
                               resultant_test_case_t *test_case, flint_rand_t state) {
    printf("Testing %s: %s (nvars=%ld, len=%ldx%ld, exp_bound=%ld, elim_var=%d)... ",
           field_desc, test_case->name, test_case->nvars, 
           test_case->length1, test_case->length2, 
           test_case->exp_bound, test_case->var_index);
    fflush(stdout);
    
    /* Initialize contexts */
    fq_nmod_mpoly_ctx_t fq_mpoly_ctx;
    fq_nmod_mpoly_ctx_init(fq_mpoly_ctx, test_case->nvars, ORD_LEX, fq_ctx);
    
    field_ctx_t field_ctx;
    field_ctx_init(&field_ctx, fq_ctx);
    
    unified_mpoly_ctx_t unified_ctx = unified_mpoly_ctx_init(test_case->nvars, 
                                                            ORD_LEX, &field_ctx);
    
    /* Generate random polynomials */
    fq_nmod_mpoly_t A_fq, B_fq, R_fq_expected;
    fq_nmod_mpoly_init(A_fq, fq_mpoly_ctx);
    fq_nmod_mpoly_init(B_fq, fq_mpoly_ctx);
    fq_nmod_mpoly_init(R_fq_expected, fq_mpoly_ctx);
    
    fq_nmod_mpoly_randtest_dense(A_fq, state, test_case->length1, 
                                test_case->exp_bound, fq_mpoly_ctx);
    fq_nmod_mpoly_randtest_dense(B_fq, state, test_case->length2, 
                                test_case->exp_bound, fq_mpoly_ctx);
    
    /* Convert to unified format */
    unified_mpoly_t A_unified = unified_mpoly_init(unified_ctx);
    unified_mpoly_t B_unified = unified_mpoly_init(unified_ctx);
    unified_mpoly_t R_unified = unified_mpoly_init(unified_ctx);
    
    /* Copy data based on field type */
    if (field_ctx.field_id == FIELD_ID_NMOD) {
        /* For prime fields, convert from fq_nmod to nmod */
        nmod_mpoly_t A_nmod, B_nmod;
        nmod_mpoly_init(A_nmod, &unified_ctx->ctx.nmod_ctx);
        nmod_mpoly_init(B_nmod, &unified_ctx->ctx.nmod_ctx);
        
        /* Convert fq_nmod_mpoly to nmod_mpoly */
        for (slong i = 0; i < fq_nmod_mpoly_length(A_fq, fq_mpoly_ctx); i++) {
            fq_nmod_t coeff;
            ulong *exp = (ulong *)calloc(test_case->nvars, sizeof(ulong));
            fq_nmod_init(coeff, fq_ctx);
            
            fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, A_fq, i, fq_mpoly_ctx);
            fq_nmod_mpoly_get_term_exp_ui(exp, A_fq, i, fq_mpoly_ctx);
            
            ulong c = nmod_poly_get_coeff_ui(coeff, 0);
            nmod_mpoly_set_coeff_ui_ui(A_nmod, c, exp, &unified_ctx->ctx.nmod_ctx);
            
            fq_nmod_clear(coeff, fq_ctx);
            free(exp);
        }
        
        for (slong i = 0; i < fq_nmod_mpoly_length(B_fq, fq_mpoly_ctx); i++) {
            fq_nmod_t coeff;
            ulong *exp = (ulong *)calloc(test_case->nvars, sizeof(ulong));
            fq_nmod_init(coeff, fq_ctx);
            
            fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, B_fq, i, fq_mpoly_ctx);
            fq_nmod_mpoly_get_term_exp_ui(exp, B_fq, i, fq_mpoly_ctx);
            
            ulong c = nmod_poly_get_coeff_ui(coeff, 0);
            nmod_mpoly_set_coeff_ui_ui(B_nmod, c, exp, &unified_ctx->ctx.nmod_ctx);
            
            fq_nmod_clear(coeff, fq_ctx);
            free(exp);
        }
        
        nmod_mpoly_set(GET_NMOD_POLY(A_unified), A_nmod, &unified_ctx->ctx.nmod_ctx);
        nmod_mpoly_set(GET_NMOD_POLY(B_unified), B_nmod, &unified_ctx->ctx.nmod_ctx);
        
        nmod_mpoly_clear(A_nmod, &unified_ctx->ctx.nmod_ctx);
        nmod_mpoly_clear(B_nmod, &unified_ctx->ctx.nmod_ctx);
    } else {
        /* For extension fields, copy directly */
        fq_nmod_mpoly_set(GET_FQ_POLY(A_unified), A_fq, &unified_ctx->ctx.fq_ctx);
        fq_nmod_mpoly_set(GET_FQ_POLY(B_unified), B_fq, &unified_ctx->ctx.fq_ctx);
    }
    
    /* Set up timeout and compute resultants */
    int success_flint = 0, success_unified = 0;
    double time_flint = 0, time_unified = 0;
    
    /* Try FLINT resultant with timeout */
    if (setjmp(timeout_env) == 0) {
        set_timeout(20);
        clock_t start = clock();
        success_flint = fq_nmod_mpoly_resultant(R_fq_expected, A_fq, B_fq, 
                                               test_case->var_index, fq_mpoly_ctx);
        clock_t end = clock();
        time_flint = ((double)(end - start)) / CLOCKS_PER_SEC;
        clear_timeout();
    } else {
        clear_timeout();
        printf("TIMEOUT (FLINT)\n");
        goto cleanup;
    }
    
    if (!success_flint) {
        printf("FAILED (FLINT)\n");
        goto cleanup;
    }
    
    /* Try unified resultant with timeout */
    if (setjmp(timeout_env) == 0) {
        set_timeout(20);
        clock_t start = clock();
        success_unified = unified_mpoly_resultant(R_unified, A_unified, B_unified,
                                                test_case->var_index, unified_ctx);
        clock_t end = clock();
        time_unified = ((double)(end - start)) / CLOCKS_PER_SEC;
        clear_timeout();
    } else {
        clear_timeout();
        printf("TIMEOUT (Unified)\n");
        goto cleanup;
    }
    
    if (!success_unified) {
        printf("FAILED (Unified)\n");
        goto cleanup;
    }
    
    /* Compare results */
    int equal = 0;
    
    if (field_ctx.field_id == FIELD_ID_NMOD) {
        /* Convert nmod result to fq_nmod for comparison */
        fq_nmod_mpoly_t R_unified_fq;
        fq_nmod_mpoly_init(R_unified_fq, fq_mpoly_ctx);
        
        nmod_mpoly_struct *R_nmod = GET_NMOD_POLY(R_unified);
        for (slong i = 0; i < nmod_mpoly_length(R_nmod, &unified_ctx->ctx.nmod_ctx); i++) {
            ulong c = nmod_mpoly_get_term_coeff_ui(R_nmod, i, &unified_ctx->ctx.nmod_ctx);
            ulong *exp = (ulong *)calloc(test_case->nvars, sizeof(ulong));
            nmod_mpoly_get_term_exp_ui(exp, R_nmod, i, &unified_ctx->ctx.nmod_ctx);
            
            fq_nmod_t coeff;
            fq_nmod_init(coeff, fq_ctx);
            fq_nmod_set_ui(coeff, c, fq_ctx);
            fq_nmod_mpoly_set_coeff_fq_nmod_ui(R_unified_fq, coeff, exp, fq_mpoly_ctx);
            fq_nmod_clear(coeff, fq_ctx);
            free(exp);
        }
        
        equal = fq_nmod_mpoly_equal(R_fq_expected, R_unified_fq, fq_mpoly_ctx);
        fq_nmod_mpoly_clear(R_unified_fq, fq_mpoly_ctx);
    } else {
        equal = fq_nmod_mpoly_equal(R_fq_expected, GET_FQ_POLY(R_unified), 
                                   fq_mpoly_ctx);
    }
    
    if (equal) {
        printf("PASS (FLINT: %.3fs, Unified: %.3fs, Ratio: %.2fx)\n", 
               time_flint, time_unified, time_unified / time_flint);
    } else {
        printf("FAIL (results differ)\n");
    }
    
cleanup:
    /* Clean up */
    fq_nmod_mpoly_clear(A_fq, fq_mpoly_ctx);
    fq_nmod_mpoly_clear(B_fq, fq_mpoly_ctx);
    fq_nmod_mpoly_clear(R_fq_expected, fq_mpoly_ctx);
    fq_nmod_mpoly_ctx_clear(fq_mpoly_ctx);
    
    unified_mpoly_clear(A_unified);
    unified_mpoly_clear(B_unified);
    unified_mpoly_clear(R_unified);
    unified_mpoly_ctx_clear(unified_ctx);
    
    return equal;
}

/* Performance comparison table */
void performance_comparison_table(const char *field_desc, const fq_nmod_ctx_t fq_ctx,
                                 flint_rand_t state) {
    printf("\n%s Performance Comparison:\n", field_desc);
    printf("%-20s %-10s %-10s %-10s\n", "Test Case", "FLINT(s)", "Unified(s)", "Ratio");
    printf("%-20s %-10s %-10s %-10s\n", "---------", "--------", "----------", "-----");
    
    /* Dense polynomial test cases (≤4 variables) */
    resultant_test_case_t perf_cases[] = {
        {"2-var small", 2, 20, 20, 5, 0},
        {"2-var medium", 2, 50, 50, 8, 0},
        {"2-var large", 2, 100, 100, 10, 1},
        {"3-var small", 3, 15, 15, 4, 1},
        {"3-var medium", 3, 30, 30, 6, 0},
        {"3-var large", 3, 50, 50, 8, 2},
        {"4-var small", 4, 10, 10, 3, 2},
        {"4-var medium", 4, 20, 20, 4, 1},
    };
    
    /* Initialize contexts */
    fq_nmod_mpoly_ctx_t fq_mpoly_ctx;
    field_ctx_t field_ctx;
    field_ctx_init(&field_ctx, fq_ctx);
    
    for (int i = 0; i < sizeof(perf_cases)/sizeof(perf_cases[0]); i++) {
        resultant_test_case_t *tc = &perf_cases[i];
        
        fq_nmod_mpoly_ctx_init(fq_mpoly_ctx, tc->nvars, ORD_LEX, fq_ctx);
        unified_mpoly_ctx_t unified_ctx = unified_mpoly_ctx_init(tc->nvars, ORD_LEX, &field_ctx);
        
        /* Generate test polynomials */
        fq_nmod_mpoly_t A, B, R_flint;
        fq_nmod_mpoly_init(A, fq_mpoly_ctx);
        fq_nmod_mpoly_init(B, fq_mpoly_ctx);
        fq_nmod_mpoly_init(R_flint, fq_mpoly_ctx);
        
        fq_nmod_mpoly_randtest_dense(A, state, tc->length1, tc->exp_bound, fq_mpoly_ctx);
        fq_nmod_mpoly_randtest_dense(B, state, tc->length2, tc->exp_bound, fq_mpoly_ctx);
        
        /* Time FLINT */
        double time_flint = 0;
        if (setjmp(timeout_env) == 0) {
            set_timeout(20);
            clock_t start = clock();
            int success = fq_nmod_mpoly_resultant(R_flint, A, B, tc->var_index, fq_mpoly_ctx);
            clock_t end = clock();
            clear_timeout();
            if (success) {
                time_flint = ((double)(end - start)) / CLOCKS_PER_SEC;
            }
        } else {
            clear_timeout();
            time_flint = -1; // timeout
        }
        
        /* Prepare unified polynomials */
        unified_mpoly_t A_u = unified_mpoly_init(unified_ctx);
        unified_mpoly_t B_u = unified_mpoly_init(unified_ctx);
        unified_mpoly_t R_u = unified_mpoly_init(unified_ctx);
        
        if (field_ctx.field_id == FIELD_ID_NMOD) {
            /* Convert for prime fields */
            nmod_mpoly_t A_nmod, B_nmod;
            nmod_mpoly_init(A_nmod, &unified_ctx->ctx.nmod_ctx);
            nmod_mpoly_init(B_nmod, &unified_ctx->ctx.nmod_ctx);
            
            for (slong j = 0; j < fq_nmod_mpoly_length(A, fq_mpoly_ctx); j++) {
                fq_nmod_t coeff;
                ulong *exp = (ulong *)calloc(tc->nvars, sizeof(ulong));
                fq_nmod_init(coeff, fq_ctx);
                
                fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, A, j, fq_mpoly_ctx);
                fq_nmod_mpoly_get_term_exp_ui(exp, A, j, fq_mpoly_ctx);
                
                ulong c = nmod_poly_get_coeff_ui(coeff, 0);
                nmod_mpoly_set_coeff_ui_ui(A_nmod, c, exp, &unified_ctx->ctx.nmod_ctx);
                
                fq_nmod_clear(coeff, fq_ctx);
                free(exp);
            }
            
            for (slong j = 0; j < fq_nmod_mpoly_length(B, fq_mpoly_ctx); j++) {
                fq_nmod_t coeff;
                ulong *exp = (ulong *)calloc(tc->nvars, sizeof(ulong));
                fq_nmod_init(coeff, fq_ctx);
                
                fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, B, j, fq_mpoly_ctx);
                fq_nmod_mpoly_get_term_exp_ui(exp, B, j, fq_mpoly_ctx);
                
                ulong c = nmod_poly_get_coeff_ui(coeff, 0);
                nmod_mpoly_set_coeff_ui_ui(B_nmod, c, exp, &unified_ctx->ctx.nmod_ctx);
                
                fq_nmod_clear(coeff, fq_ctx);
                free(exp);
            }
            
            nmod_mpoly_set(GET_NMOD_POLY(A_u), A_nmod, &unified_ctx->ctx.nmod_ctx);
            nmod_mpoly_set(GET_NMOD_POLY(B_u), B_nmod, &unified_ctx->ctx.nmod_ctx);
            
            nmod_mpoly_clear(A_nmod, &unified_ctx->ctx.nmod_ctx);
            nmod_mpoly_clear(B_nmod, &unified_ctx->ctx.nmod_ctx);
        } else {
            fq_nmod_mpoly_set(GET_FQ_POLY(A_u), A, &unified_ctx->ctx.fq_ctx);
            fq_nmod_mpoly_set(GET_FQ_POLY(B_u), B, &unified_ctx->ctx.fq_ctx);
        }
        
        /* Time unified */
        double time_unified = 0;
        if (setjmp(timeout_env) == 0) {
            set_timeout(20);
            clock_t start = clock();
            int success = unified_mpoly_resultant(R_u, A_u, B_u, tc->var_index, unified_ctx);
            clock_t end = clock();
            clear_timeout();
            if (success) {
                time_unified = ((double)(end - start)) / CLOCKS_PER_SEC;
            }
        } else {
            clear_timeout();
            time_unified = -1; // timeout
        }
        
        /* Print results */
        printf("%-20s ", tc->name);
        if (time_flint < 0) {
            printf("%-10s ", "TIMEOUT");
        } else {
            printf("%-10.3f ", time_flint);
        }
        
        if (time_unified < 0) {
            printf("%-10s ", "TIMEOUT");
        } else {
            printf("%-10.3f ", time_unified);
        }
        
        if (time_flint > 0 && time_unified > 0) {
            printf("%-10.2f\n", time_unified / time_flint);
        } else {
            printf("%-10s\n", "N/A");
        }
        
        /* Clean up */
        fq_nmod_mpoly_clear(A, fq_mpoly_ctx);
        fq_nmod_mpoly_clear(B, fq_mpoly_ctx);
        fq_nmod_mpoly_clear(R_flint, fq_mpoly_ctx);
        
        unified_mpoly_clear(A_u);
        unified_mpoly_clear(B_u);
        unified_mpoly_clear(R_u);
        
        fq_nmod_mpoly_ctx_clear(fq_mpoly_ctx);
        unified_mpoly_ctx_clear(unified_ctx);
    }
}

int main() {
    printf("=== Unified Multivariate Polynomial Resultant Test ===\n");
    printf("(Dense polynomials, ≤4 variables, 20s timeout)\n\n");
    
    /* Initialize random state */
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, time(NULL), time(NULL) + 123456);
    
    /* Test cases for dense polynomials */
    resultant_test_case_t test_cases[] = {
        {"2-var tiny", 2, 5, 5, 3, 0},
        {"2-var small", 2, 10, 10, 5, 1},
        {"2-var medium", 2, 30, 30, 8, 0},
        {"3-var tiny", 3, 5, 5, 3, 1},
        {"3-var small", 3, 15, 15, 5, 0},
        {"3-var medium", 3, 25, 25, 6, 2},
        {"4-var tiny", 4, 5, 5, 2, 2},
        {"4-var small", 4, 10, 10, 3, 1},
    };
    int num_tests = sizeof(test_cases) / sizeof(test_cases[0]);
    
    int total_passed = 0;
    int total_tests = 0;
    
    /* Test 1: Prime field Z/p */
    printf("========== PRIME FIELD Z/p TESTS ==========\n");
    {
        fmpz_t p;
        fmpz_init_set_ui(p, 1073741827);  /* Large prime */
        fq_nmod_ctx_t fq_ctx;
        fq_nmod_ctx_init(fq_ctx, p, 1, "x");
        
        for (int i = 0; i < num_tests; i++) {
            if (test_resultant_with_timeout("Z/p", fq_ctx, &test_cases[i], state)) {
                total_passed++;
            }
            total_tests++;
        }
        
        performance_comparison_table("Z/p", fq_ctx, state);
        
        fq_nmod_ctx_clear(fq_ctx);
        fmpz_clear(p);
    }
    
    /* Test 2: GF(2^8) */
    printf("\n========== GF(2^8) TESTS ==========\n");
    {
        fmpz_t p;
        fmpz_init_set_ui(p, 2);
        fq_nmod_ctx_t fq_ctx;
        fq_nmod_ctx_init(fq_ctx, p, 8, "a");
        
        unified_mpoly_enable_all_optimizations(FIELD_ID_GF28);
        
        for (int i = 0; i < num_tests; i++) {
            if (test_resultant_with_timeout("GF(2^8)", fq_ctx, &test_cases[i], state)) {
                total_passed++;
            }
            total_tests++;
        }
        
        performance_comparison_table("GF(2^8)", fq_ctx, state);
        
        fq_nmod_ctx_clear(fq_ctx);
        fmpz_clear(p);
    }
    
    /* Test 3: GF(2^128) - smaller test cases only */
    printf("\n========== GF(2^128) TESTS ==========\n");
    {
        fmpz_t p;
        fmpz_init_set_ui(p, 2);
        
        /* Initialize with custom modulus x^128 + x^7 + x^2 + x + 1 */
        nmod_poly_t modulus;
        nmod_poly_init(modulus, 2);
        nmod_poly_set_coeff_ui(modulus, 128, 1);
        nmod_poly_set_coeff_ui(modulus, 7, 1);
        nmod_poly_set_coeff_ui(modulus, 2, 1);
        nmod_poly_set_coeff_ui(modulus, 1, 1);
        nmod_poly_set_coeff_ui(modulus, 0, 1);
        
        fq_nmod_ctx_t fq_ctx;
        fq_nmod_ctx_init_modulus(fq_ctx, modulus, "a");
        
        unified_mpoly_enable_all_optimizations(FIELD_ID_GF2128);
        
        /* Only test small cases for GF(2^128) */
        resultant_test_case_t gf2128_cases[] = {
            {"2-var tiny", 2, 5, 5, 3, 0},
            {"2-var small", 2, 10, 10, 4, 1},
            {"3-var tiny", 3, 5, 5, 3, 0},
        };
        
        for (int i = 0; i < 3; i++) {
            if (test_resultant_with_timeout("GF(2^128)", fq_ctx, &gf2128_cases[i], state)) {
                total_passed++;
            }
            total_tests++;
        }
        
        performance_comparison_table("GF(2^128)", fq_ctx, state);
        
        fq_nmod_ctx_clear(fq_ctx);
        nmod_poly_clear(modulus);
        fmpz_clear(p);
    }
    
    /* Test 4: General extension field GF(11^3) */
    printf("\n========== GF(11^3) TESTS ==========\n");
    {
        fmpz_t p;
        fmpz_init_set_ui(p, 11);
        fq_nmod_ctx_t fq_ctx;
        fq_nmod_ctx_init(fq_ctx, p, 3, "a");
        
        for (int i = 0; i < num_tests; i++) {
            if (test_resultant_with_timeout("GF(11^3)", fq_ctx, &test_cases[i], state)) {
                total_passed++;
            }
            total_tests++;
        }
        
        performance_comparison_table("GF(11^3)", fq_ctx, state);
        
        fq_nmod_ctx_clear(fq_ctx);
        fmpz_clear(p);
    }
    
    /* Summary */
    printf("\n========== TEST SUMMARY ==========\n");
    printf("Total tests: %d\n", total_tests);
    printf("Passed: %d\n", total_passed);
    printf("Failed: %d\n", total_tests - total_passed);
    
    if (total_passed == total_tests) {
        printf("\n✓ All tests passed!\n");
    } else {
        printf("\n✗ Some tests failed.\n");
    }
    
    /* Clean up */
    flint_rand_clear(state);
    
    return (total_passed == total_tests) ? 0 : 1;
}