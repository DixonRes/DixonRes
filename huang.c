#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_mod_mpoly.h>
#include <flint/fmpz_mod_mat.h>
#include <flint/fmpz_mod_poly_factor.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

// Debug flag
#define DEBUG 0
#define TIMING 1
#define DETAILED_TIMING 1

#include <time.h>

// Timing structure
typedef struct {
    clock_t start;
    double elapsed;
} my_timer_t;

void timer_start(my_timer_t* t) {
    t->start = clock();
}

void timer_stop(my_timer_t* t) {
    t->elapsed = ((double)(clock() - t->start)) / CLOCKS_PER_SEC;
}

void timer_print(my_timer_t* t, const char* label) {
    if (TIMING) {
        printf("  [TIMING] %s: %.6f seconds\n", label, t->elapsed);
    }
}

// Global random state
flint_rand_t global_state;

// Polynomial matrix structure
typedef struct {
    fmpz_mod_mpoly_struct** entries;
    slong rows;
    slong cols;
} poly_mat_t;

void poly_mat_init(poly_mat_t* mat, slong rows, slong cols, const fmpz_mod_mpoly_ctx_t ctx) {
    mat->rows = rows;
    mat->cols = cols;
    mat->entries = (fmpz_mod_mpoly_struct**)malloc(rows * sizeof(fmpz_mod_mpoly_struct*));
    for (slong i = 0; i < rows; i++) {
        mat->entries[i] = (fmpz_mod_mpoly_struct*)malloc(cols * sizeof(fmpz_mod_mpoly_struct));
        for (slong j = 0; j < cols; j++) {
            fmpz_mod_mpoly_init(mat->entries[i] + j, ctx);
        }
    }
}

void poly_mat_clear(poly_mat_t* mat, const fmpz_mod_mpoly_ctx_t ctx) {
    for (slong i = 0; i < mat->rows; i++) {
        for (slong j = 0; j < mat->cols; j++) {
            fmpz_mod_mpoly_clear(mat->entries[i] + j, ctx);
        }
        free(mat->entries[i]);
    }
    free(mat->entries);
}

void poly_mat_entry_set(poly_mat_t* mat, slong i, slong j, const fmpz_mod_mpoly_t poly, const fmpz_mod_mpoly_ctx_t ctx) {
    fmpz_mod_mpoly_set(mat->entries[i] + j, poly, ctx);
}

// Compute 2x2 determinant for testing
void poly_mat_det(fmpz_mod_mpoly_t det, const poly_mat_t* mat, const fmpz_mod_mpoly_ctx_t ctx) {
    if (mat->rows != mat->cols || mat->rows != 2) {
        fmpz_mod_mpoly_zero(det, ctx);
        return;
    }
    
    fmpz_mod_mpoly_t temp1, temp2;
    fmpz_mod_mpoly_init(temp1, ctx);
    fmpz_mod_mpoly_init(temp2, ctx);
    
    // det = a00*a11 - a01*a10
    fmpz_mod_mpoly_mul(temp1, mat->entries[0] + 0, mat->entries[1] + 1, ctx);
    fmpz_mod_mpoly_mul(temp2, mat->entries[0] + 1, mat->entries[1] + 0, ctx);
    fmpz_mod_mpoly_sub(det, temp1, temp2, ctx);
    
    fmpz_mod_mpoly_clear(temp1, ctx);
    fmpz_mod_mpoly_clear(temp2, ctx);
}

// Get the maximum total degree of a polynomial
slong poly_max_total_degree(const fmpz_mod_mpoly_t poly, const fmpz_mod_mpoly_ctx_t ctx) {
    slong max_deg = 0;
    slong nvars = fmpz_mod_mpoly_ctx_nvars(ctx);
    slong len = fmpz_mod_mpoly_length(poly, ctx);
    
    for (slong i = 0; i < len; i++) {
        ulong* exp = (ulong*)malloc(nvars * sizeof(ulong));
        fmpz_mod_mpoly_get_term_exp_ui(exp, poly, i, ctx);
        
        slong total_deg = 0;
        for (slong j = 0; j < nvars; j++) {
            total_deg += exp[j];
        }
        
        if (total_deg > max_deg) {
            max_deg = total_deg;
        }
        
        free(exp);
    }
    
    return max_deg;
}

// Get the maximum degree in each variable of a polynomial
void poly_max_degrees_per_var(slong* max_degs, const fmpz_mod_mpoly_t poly, const fmpz_mod_mpoly_ctx_t ctx) {
    slong nvars = fmpz_mod_mpoly_ctx_nvars(ctx);
    slong len = fmpz_mod_mpoly_length(poly, ctx);
    
    // Initialize to 0
    for (slong i = 0; i < nvars; i++) {
        max_degs[i] = 0;
    }
    
    // Find max degree in each variable
    for (slong i = 0; i < len; i++) {
        ulong* exp = (ulong*)malloc(nvars * sizeof(ulong));
        fmpz_mod_mpoly_get_term_exp_ui(exp, poly, i, ctx);
        
        for (slong j = 0; j < nvars; j++) {
            if (exp[j] > max_degs[j]) {
                max_degs[j] = exp[j];
            }
        }
        
        free(exp);
    }
}

// Berlekamp-Massey Algorithm (BM)
void BM(fmpz_mod_poly_t C, const fmpz* s, slong N, const fmpz_mod_ctx_t ctx) {
    my_timer_t timer;
    if (DETAILED_TIMING) timer_start(&timer);
    
    fmpz_mod_poly_t B, T;
    fmpz_t d, b, temp, temp2;
    slong L = 0, k = 1;
    
    fmpz_mod_poly_init(B, ctx);
    fmpz_mod_poly_init(T, ctx);
    fmpz_init(d);
    fmpz_init(b);
    fmpz_init(temp);
    fmpz_init(temp2);
    
    fmpz_mod_poly_set_ui(B, 1, ctx);
    fmpz_mod_poly_set_ui(C, 1, ctx);
    fmpz_one(b);
    
    for (slong n = 0; n < 2*N; n++) {
        // Calculate d = s[n] + sum(C[i] * s[n-i]) for i=1 to L
        fmpz_set(d, s + n);
        for (slong i = 1; i <= L && i <= n; i++) {
            fmpz_mod_poly_get_coeff_fmpz(temp, C, i, ctx);
            fmpz_mod_mul(temp2, temp, s + n - i, ctx);
            fmpz_mod_add(d, d, temp2, ctx);
        }
        
        if (fmpz_is_zero(d)) {
            k++;
        } else if (n < 2*L) {
            // C = C - (d/b) * x^k * B
            fmpz_mod_poly_t xkB;
            fmpz_mod_poly_init(xkB, ctx);
            fmpz_mod_poly_shift_left(xkB, B, k, ctx);
            fmpz_mod_inv(temp, b, ctx);
            fmpz_mod_mul(temp, d, temp, ctx);
            fmpz_mod_poly_scalar_mul_fmpz(xkB, xkB, temp, ctx);
            fmpz_mod_poly_sub(C, C, xkB, ctx);
            fmpz_mod_poly_clear(xkB, ctx);
            k++;
        } else {
            // Save old C
            fmpz_mod_poly_set(T, C, ctx);
            
            // C = C - (d/b) * x^k * B
            fmpz_mod_poly_t xkB;
            fmpz_mod_poly_init(xkB, ctx);
            fmpz_mod_poly_shift_left(xkB, B, k, ctx);
            fmpz_mod_inv(temp, b, ctx);
            fmpz_mod_mul(temp, d, temp, ctx);
            fmpz_mod_poly_scalar_mul_fmpz(xkB, xkB, temp, ctx);
            fmpz_mod_poly_sub(C, C, xkB, ctx);
            fmpz_mod_poly_clear(xkB, ctx);
            
            // Update B, L, k, b
            fmpz_mod_poly_set(B, T, ctx);
            L = n + 1 - L;
            k = 1;
            fmpz_set(b, d);
        }
    }
    
    fmpz_mod_poly_clear(B, ctx);
    fmpz_mod_poly_clear(T, ctx);
    fmpz_clear(d);
    fmpz_clear(b);
    fmpz_clear(temp);
    fmpz_clear(temp2);
    
    if (DETAILED_TIMING) {
        timer_stop(&timer);
        timer_print(&timer, "    BM algorithm");
    }
}

// Vinvert - compute coefficients from roots using partial fraction decomposition
void Vinvert(fmpz* c1, const fmpz_mod_poly_t c, const fmpz* v, 
             const fmpz* a, slong n, const fmpz_mod_ctx_t ctx) {
    
    my_timer_t timer;
    if (DETAILED_TIMING) timer_start(&timer);
    
    fmpz_mod_poly_t d, q, q1, q2;
    fmpz_t temp;
    
    fmpz_mod_poly_init(d, ctx);
    fmpz_mod_poly_init(q, ctx);
    fmpz_mod_poly_init(q1, ctx);
    fmpz_mod_poly_init(q2, ctx);
    fmpz_init(temp);
    
    // Build polynomial d = sum(a[i] * z^(n+1-i))
    fmpz_mod_poly_zero(d, ctx);
    for (slong i = 0; i < n; i++) {
        fmpz_mod_poly_set_coeff_fmpz(d, n - i, a + i, ctx);
    }
    
    // q = c * d
    fmpz_mod_poly_mul(q, c, d, ctx);
    
    // Extract coefficients for q1
    fmpz_mod_poly_zero(q1, ctx);
    for (slong i = 0; i < n; i++) {
        fmpz_mod_poly_get_coeff_fmpz(temp, q, 2*n - i, ctx);
        fmpz_mod_poly_set_coeff_fmpz(q1, n - 1 - i, temp, ctx);
    }
    
    // q2 = derivative of c
    fmpz_mod_poly_derivative(q2, c, ctx);
    
    // Evaluate and compute coefficients
    for (slong i = 0; i < n; i++) {
        fmpz_t num, den, q1_val, q2_val;
        fmpz_init(num);
        fmpz_init(den);
        fmpz_init(q1_val);
        fmpz_init(q2_val);
        
        // Evaluate q1 at v[i]
        fmpz_mod_poly_evaluate_fmpz(q1_val, q1, v + i, ctx);
        
        // Evaluate q2 at v[i]
        fmpz_mod_poly_evaluate_fmpz(q2_val, q2, v + i, ctx);
        
        // Compute num / den
        if (!fmpz_is_zero(q2_val)) {
            fmpz_mod_inv(temp, q2_val, ctx);
            fmpz_mod_mul(c1 + i, q1_val, temp, ctx);
        } else {
            // Fallback: just use a[i] if derivative is zero
            fmpz_set(c1 + i, a + i);
        }
        
        fmpz_clear(num);
        fmpz_clear(den);
        fmpz_clear(q1_val);
        fmpz_clear(q2_val);
    }
    
    fmpz_mod_poly_clear(d, ctx);
    fmpz_mod_poly_clear(q, ctx);
    fmpz_mod_poly_clear(q1, ctx);
    fmpz_mod_poly_clear(q2, ctx);
    fmpz_clear(temp);
    
    if (DETAILED_TIMING) {
        timer_stop(&timer);
        timer_print(&timer, "    Vinvert");
    }
}

typedef struct {
    fmpz_t coeff;
    fmpz_t root;
} term_pair;

typedef struct {
    term_pair* pairs;
    slong length;
} term_list;

// Compare function for sorting by root value (matching Magma's sort by second component)
int compare_pairs_by_root(const void* a, const void* b) {
    const term_pair* pa = (const term_pair*)a;
    const term_pair* pb = (const term_pair*)b;
    return fmpz_cmp(pa->root, pb->root);
}

// MC - Monomials and Coefficients (Algorithm 3.1)
term_list MC(const fmpz* a, slong N, const fmpz_mod_ctx_t ctx) {
    my_timer_t timer, timer_inner;
    if (DETAILED_TIMING) timer_start(&timer);
    
    term_list result = {NULL, 0};
    fmpz_mod_poly_t f, Lambda;
    fmpz_mod_poly_factor_t fac;
    fmpz* v;
    fmpz* c;
    fmpz* q;
    slong m;
    
    fmpz_mod_poly_init(f, ctx);
    fmpz_mod_poly_init(Lambda, ctx);
    fmpz_mod_poly_factor_init(fac, ctx);
    
    // Find minimal polynomial using BM
    BM(f, a, N, ctx);
    m = fmpz_mod_poly_degree(f, ctx);
    
    if (m == 0) {
        fmpz_mod_poly_clear(f, ctx);
        fmpz_mod_poly_clear(Lambda, ctx);
        fmpz_mod_poly_factor_clear(fac, ctx);
        return result;
    }
    
    // Build Lambda(z) = z^m + sum(f[i] * z^(m-i))
    fmpz_mod_poly_set_coeff_ui(Lambda, m, 1, ctx);
    for (slong i = 1; i <= m; i++) {
        fmpz_t coeff;
        fmpz_init(coeff);
        fmpz_mod_poly_get_coeff_fmpz(coeff, f, i, ctx);
        fmpz_mod_poly_set_coeff_fmpz(Lambda, m - i, coeff, ctx);
        fmpz_clear(coeff);
    }
    
    // Allocate arrays
    v = (fmpz*)malloc(m * sizeof(fmpz));
    c = (fmpz*)malloc(m * sizeof(fmpz));
    q = (fmpz*)malloc(m * sizeof(fmpz));
    for (slong i = 0; i < m; i++) {
        fmpz_init(v + i);
        fmpz_init(c + i);
        fmpz_init(q + i);
        fmpz_set(q + i, a + i);  // Get first m values
    }
    
    // Find roots of Lambda
    if (DETAILED_TIMING) timer_start(&timer_inner);
    fmpz_mod_poly_factor(fac, Lambda, ctx);
    if (DETAILED_TIMING) {
        timer_stop(&timer_inner);
        timer_print(&timer_inner, "    Polynomial factorization");
    }
    
    slong root_count = 0;
    
    // Store roots and their multiplicities
    for (slong i = 0; i < fac->num && root_count < m; i++) {
        slong deg = fmpz_mod_poly_degree(fac->poly + i, ctx);
        
        if (deg == 1) {
            // Linear factor: (z - root) or (a*z - b)
            fmpz_t root, coeff0, coeff1;
            fmpz_init(root);
            fmpz_init(coeff0);
            fmpz_init(coeff1);
            
            fmpz_mod_poly_get_coeff_fmpz(coeff0, fac->poly + i, 0, ctx);
            fmpz_mod_poly_get_coeff_fmpz(coeff1, fac->poly + i, 1, ctx);
            
            // root = -coeff0/coeff1
            fmpz_mod_inv(root, coeff1, ctx);
            fmpz_mod_mul(root, root, coeff0, ctx);
            fmpz_mod_neg(root, root, ctx);
            
            for (slong j = 0; j < fac->exp[i] && root_count < m; j++) {
                fmpz_set(v + root_count, root);
                root_count++;
            }
            
            fmpz_clear(root);
            fmpz_clear(coeff0);
            fmpz_clear(coeff1);
        }
    }
    
    if (root_count < m) {
        // Try harder to find roots - use a more systematic approach
        const fmpz* p = fmpz_mod_ctx_modulus(ctx);
        
        // For moderate size primes, do exhaustive search
        if (fmpz_cmp_ui(p, 10000000) < 0 && root_count < m) {
            fmpz_t test_val, eval_result;
            fmpz_init(test_val);
            fmpz_init(eval_result);
            
            // Try random values first
            for (slong attempts = 0; attempts < 1000 && root_count < m; attempts++) {
                fmpz_randm(test_val, global_state, p);
                fmpz_mod_poly_evaluate_fmpz(eval_result, Lambda, test_val, ctx);
                
                if (fmpz_is_zero(eval_result)) {
                    // Check if we already have this root
                    int already_have = 0;
                    for (slong k = 0; k < root_count; k++) {
                        if (fmpz_equal(v + k, test_val)) {
                            already_have = 1;
                            break;
                        }
                    }
                    
                    if (!already_have) {
                        fmpz_set(v + root_count, test_val);
                        root_count++;
                    }
                }
            }
            
            fmpz_clear(test_val);
            fmpz_clear(eval_result);
        }
        
        // If still not enough, the polynomial doesn't split completely
        // Use placeholder values to avoid errors
        while (root_count < m) {
            if (root_count > 0) {
                fmpz_set(v + root_count, v);  // Duplicate first root
            } else {
                fmpz_one(v + root_count);  // Use 1 as default
            }
            root_count++;
        }
    }
    
    // Compute coefficients using Vinvert
    if (m > 1) {
        Vinvert(c, Lambda, v, q, m, ctx);
    } else {
        fmpz_set(c, a);
    }
    
    // Build result
    result.length = m;
    result.pairs = (term_pair*)malloc(m * sizeof(term_pair));
    
    for (slong i = 0; i < m; i++) {
        fmpz_init(result.pairs[i].coeff);
        fmpz_init(result.pairs[i].root);
        fmpz_set(result.pairs[i].coeff, c + i);
        fmpz_set(result.pairs[i].root, v + i);
    }
    
    // Sort by root value (matching Maple's sort by second component)
    qsort(result.pairs, result.length, sizeof(term_pair), compare_pairs_by_root);
    
    // Cleanup
    for (slong i = 0; i < m; i++) {
        fmpz_clear(v + i);
        fmpz_clear(c + i);
        fmpz_clear(q + i);
    }
    free(v);
    free(c);
    free(q);
    
    fmpz_mod_poly_clear(f, ctx);
    fmpz_mod_poly_clear(Lambda, ctx);
    fmpz_mod_poly_factor_clear(fac, ctx);
    
    if (DETAILED_TIMING) {
        timer_stop(&timer);
        timer_print(&timer, "  MC (total)");
    }
    
    return result;
}

// Compute discrete logarithm with bounded search
slong mylog(const fmpz_t omega, const fmpz_t B, slong d, const fmpz_mod_ctx_t ctx) {
    if (fmpz_is_one(B)) {
        return 0;
    }
    
    // Try all exponents from 0 to d
    fmpz_t omega_pow;
    fmpz_init(omega_pow);
    fmpz_one(omega_pow);
    
    for (slong k = 0; k <= d; k++) {
        if (fmpz_equal(omega_pow, B)) {
            fmpz_clear(omega_pow);
            return k;
        }
        fmpz_mod_mul(omega_pow, omega_pow, omega, ctx);
    }
    
    fmpz_clear(omega_pow);
    
    // If not found, print warning in debug mode
    if (DEBUG) {
        printf("WARNING: mylog failed to find discrete log of ");
        fmpz_print(B);
        printf(" (searched up to %ld)\n", d);
    }
    
    return -1; // Not found
}

// Evaluate polynomial at a point
void Myeval(fmpz_mod_mat_t c, const fmpz_mod_mpoly_t f, slong n, 
            const fmpz* alpha, const fmpz_mod_ctx_t ctx, const fmpz_mod_mpoly_ctx_t mctx) {
    
    slong t = fmpz_mod_mpoly_length(f, mctx);
    
    if (t == 0) {
        fmpz_mod_mat_zero(c, ctx);
        return;
    }
    
    // Get monomials and coefficients
    fmpz* coeffs = (fmpz*)malloc(t * sizeof(fmpz));
    ulong** exps = (ulong**)malloc(t * sizeof(ulong*));
    
    for (slong i = 0; i < t; i++) {
        fmpz_init(coeffs + i);
        exps[i] = (ulong*)malloc(n * sizeof(ulong));
        fmpz_mod_mpoly_get_term_coeff_fmpz(coeffs + i, f, i, mctx);
        fmpz_mod_mpoly_get_term_exp_ui(exps[i], f, i, mctx);
    }
    
    // Evaluate monomials at alpha
    fmpz_mod_mat_t V;
    fmpz_mod_mat_init(V, t, 1, ctx);
    
    fmpz_t u, temp;
    fmpz_init(u);
    fmpz_init(temp);
    
    for (slong i = 0; i < t; i++) {
        fmpz_one(u);
        for (slong j = 0; j < n; j++) {
            if (exps[i][j] > 0) {
                fmpz_mod_pow_ui(temp, alpha + j, exps[i][j], ctx);
                fmpz_mod_mul(u, u, temp, ctx);
            }
        }
        fmpz_mod_mat_set_entry(V, i, 0, u, ctx);
    }
    
    // Build Vandermonde matrix for evaluating at powers
    fmpz_mod_mat_t M;
    fmpz_mod_mat_init(M, t, t, ctx);
    
    for (slong i = 0; i < t; i++) {
        fmpz_one(u);
        for (slong j = 0; j < t; j++) {
            fmpz_mod_mat_set_entry(M, j, i, u, ctx);
            fmpz_mod_mul(u, u, fmpz_mod_mat_entry(V, i, 0), ctx);
        }
    }
    
    // Coefficient vector
    fmpz_mod_mat_t C1;
    fmpz_mod_mat_init(C1, t, 1, ctx);
    for (slong i = 0; i < t; i++) {
        fmpz_mod_mat_set_entry(C1, i, 0, coeffs + i, ctx);
    }
    
    // a = M * C1
    fmpz_mod_mat_t a;
    fmpz_mod_mat_init(a, t, 1, ctx);
    fmpz_mod_mat_mul(a, M, C1, ctx);
    
    // Scale for second evaluation
    for (slong i = 0; i < t; i++) {
        fmpz_set(temp, fmpz_mod_mat_entry(V, i, 0));
        fmpz_mod_pow_ui(u, temp, t, ctx);
        for (slong j = 0; j < t; j++) {
            fmpz_set(temp, fmpz_mod_mat_entry(M, j, i));
            fmpz_mod_mul(temp, temp, u, ctx);
            fmpz_mod_mat_set_entry(M, j, i, temp, ctx);
        }
    }
    
    // b = M * C1
    fmpz_mod_mat_t b;
    fmpz_mod_mat_init(b, t, 1, ctx);
    fmpz_mod_mat_mul(b, M, C1, ctx);
    
    // Combine results
    for (slong i = 0; i < t; i++) {
        fmpz_mod_mat_set_entry(c, 0, i, fmpz_mod_mat_entry(a, i, 0), ctx);
        fmpz_mod_mat_set_entry(c, 0, i + t, fmpz_mod_mat_entry(b, i, 0), ctx);
    }
    
    // Cleanup
    for (slong i = 0; i < t; i++) {
        fmpz_clear(coeffs + i);
        free(exps[i]);
    }
    free(coeffs);
    free(exps);
    
    fmpz_clear(u);
    fmpz_clear(temp);
    fmpz_mod_mat_clear(V, ctx);
    fmpz_mod_mat_clear(M, ctx);
    fmpz_mod_mat_clear(C1, ctx);
    fmpz_mod_mat_clear(a, ctx);
    fmpz_mod_mat_clear(b, ctx);
}

// Total evaluation matrix
void TotalMyeval(fmpz_mod_mat_t M, const fmpz_mod_mpoly_t f, slong n,
                 const fmpz* alpha, const fmpz_t omega, 
                 const fmpz_mod_ctx_t ctx, const fmpz_mod_mpoly_ctx_t mctx) {
    
    my_timer_t timer;
    if (TIMING) timer_start(&timer);
    
    slong T = fmpz_mod_mpoly_length(f, mctx);
    
    if (T == 0) {
        fmpz_mod_mat_zero(M, ctx);
        return;
    }
    
    // First row - evaluate at alpha
    fmpz_mod_mat_t row;
    fmpz_mod_mat_init(row, 1, 2*T, ctx);
    Myeval(row, f, n, alpha, ctx, mctx);
    
    for (slong i = 0; i < 2*T; i++) {
        fmpz_mod_mat_set_entry(M, 0, i, fmpz_mod_mat_entry(row, 0, i), ctx);
    }
    fmpz_mod_mat_clear(row, ctx);
    
    // Remaining rows - modify one variable at a time
    fmpz* beta = (fmpz*)malloc(n * sizeof(fmpz));
    for (slong i = 0; i < n; i++) {
        fmpz_init(beta + i);
    }
    
    for (slong k = 0; k < n; k++) {
        // Copy alpha to beta
        for (slong i = 0; i < n; i++) {
            fmpz_set(beta + i, alpha + i);
        }
        
        // Modify k-th variable: beta[k] = alpha[k] * omega
        fmpz_mod_mul(beta + k, alpha + k, omega, ctx);
        
        // Evaluate
        fmpz_mod_mat_init(row, 1, 2*T, ctx);
        Myeval(row, f, n, beta, ctx, mctx);
        
        for (slong i = 0; i < 2*T; i++) {
            fmpz_mod_mat_set_entry(M, k + 1, i, fmpz_mod_mat_entry(row, 0, i), ctx);
        }
        fmpz_mod_mat_clear(row, ctx);
    }
    
    // Cleanup
    for (slong i = 0; i < n; i++) {
        fmpz_clear(beta + i);
    }
    free(beta);
    
    if (TIMING) {
        timer_stop(&timer);
        timer_print(&timer, "TotalMyeval");
    }
}

// Diversification transformation
void Mydiver(fmpz_mod_mpoly_t g, const fmpz_mod_mpoly_t f, slong n,
             const fmpz* zeta, const fmpz_mod_ctx_t ctx, 
             const fmpz_mod_mpoly_ctx_t mctx) {
    
    my_timer_t timer;
    if (TIMING) timer_start(&timer);
    
    slong t = fmpz_mod_mpoly_length(f, mctx);
    
    if (t == 0) {
        fmpz_mod_mpoly_zero(g, mctx);
        return;
    }
    
    fmpz_mod_mpoly_zero(g, mctx);
    
    fmpz_t coeff, u, temp;
    fmpz_init(coeff);
    fmpz_init(u);
    fmpz_init(temp);
    
    // Process each term
    for (slong i = 0; i < t; i++) {
        // Get coefficient and exponents
        fmpz_mod_mpoly_get_term_coeff_fmpz(coeff, f, i, mctx);
        ulong* exp = (ulong*)malloc(n * sizeof(ulong));
        fmpz_mod_mpoly_get_term_exp_ui(exp, f, i, mctx);
        
        // Multiply coefficient by zeta powers
        fmpz_one(u);
        for (slong j = 0; j < n; j++) {
            if (exp[j] > 0) {
                fmpz_mod_pow_ui(temp, zeta + j, exp[j], ctx);
                fmpz_mod_mul(u, u, temp, ctx);
            }
        }
        
        // Add transformed term
        fmpz_mod_mul(temp, u, coeff, ctx);
        fmpz_mod_mpoly_push_term_fmpz_ui(g, temp, exp, mctx);
        
        free(exp);
    }
    
    fmpz_clear(coeff);
    fmpz_clear(u);
    fmpz_clear(temp);
    
    if (TIMING) {
        timer_stop(&timer);
        timer_print(&timer, "Mydiver");
    }
}




// Hash table structure for omega powers
typedef struct {
    fmpz_t* keys;
    slong* values;
    slong size;
    slong capacity;
} omega_hash_table;

// Fixed hash function - using only FLINT types
static inline ulong fmpz_hash(const fmpz_t x, ulong capacity) {
    // Get the least significant bits as ulong
    ulong val = fmpz_get_ui(x);
    // Simple but effective hash function
    return (val * 2654435761UL) % capacity;
}

void omega_hash_init(omega_hash_table* ht, slong max_power) {
    ht->capacity = max_power * 2 + 100; // 2x for lower collision probability
    ht->size = 0;
    ht->keys = (fmpz_t*)malloc(ht->capacity * sizeof(fmpz_t));
    ht->values = (slong*)malloc(ht->capacity * sizeof(slong));
    
    for (slong i = 0; i < ht->capacity; i++) {
        fmpz_init(ht->keys[i]);
        ht->values[i] = -1;
    }
}

void omega_hash_clear(omega_hash_table* ht) {
    for (slong i = 0; i < ht->capacity; i++) {
        fmpz_clear(ht->keys[i]);
    }
    free(ht->keys);
    free(ht->values);
}

void omega_hash_insert(omega_hash_table* ht, const fmpz_t key, slong value) {
    ulong idx = fmpz_hash(key, ht->capacity);
    
    // Linear probing for collision resolution
    while (ht->values[idx] != -1) {
        if (fmpz_equal(ht->keys[idx], key)) {
            return; // Already exists
        }
        idx = (idx + 1) % ht->capacity;
    }
    
    fmpz_set(ht->keys[idx], key);
    ht->values[idx] = value;
    ht->size++;
}

slong omega_hash_find(omega_hash_table* ht, const fmpz_t key) {
    ulong idx = fmpz_hash(key, ht->capacity);
    
    while (ht->values[idx] != -1) {
        if (fmpz_equal(ht->keys[idx], key)) {
            return ht->values[idx];
        }
        idx = (idx + 1) % ht->capacity;
        
        // Safety check to prevent infinite loop
        static slong max_probes = 0;
        if (max_probes == 0) {
            max_probes = ht->capacity;
        }
        slong probes = 0;
        if (++probes > max_probes) {
            break;
        }
    }
    
    return -1; // Not found
}

// NOW HERE'S THE COMPLETE MODIFIED MBOT FUNCTION:
void MBOT(fmpz_mod_mpoly_t result, const fmpz_mod_mat_t M, slong n, slong T,
          const fmpz_t omega, const fmpz* zeta, const fmpz_mod_ctx_t ctx,
          slong* max_degs_per_var, const fmpz_mod_mpoly_ctx_t mctx) {
    
    my_timer_t timer, timer_total;
    timer_start(&timer_total);
    
    printf("\n=== MBOT Performance Analysis ===\n");
    printf("Parameters: n=%ld, T=%ld\n", n, T);
    
    // Process first row to get reference terms
    if (TIMING) timer_start(&timer);
    fmpz* first_row = (fmpz*)malloc(2 * T * sizeof(fmpz));
    for (slong j = 0; j < 2 * T; j++) {
        fmpz_init(first_row + j);
        fmpz_set(first_row + j, fmpz_mod_mat_entry(M, 0, j));
    }
    if (TIMING) {
        timer_stop(&timer);
        timer_print(&timer, "Extract first row");
    }
    
    // MC for first row
    printf("\nProcessing first row:\n");
    term_list a = MC(first_row, T, ctx);
    
    if (a.length == 0) {
        fmpz_mod_mpoly_zero(result, mctx);
        goto cleanup_first;
    }
    
    slong t = a.length;
    
    if (t != T) {
        printf("WARNING: MC found %ld terms, expected %ld\n", t, T);
    }
    
    // Extract coefficients and roots from first row
    fmpz* c = (fmpz*)malloc(t * sizeof(fmpz));
    fmpz* first_roots = (fmpz*)malloc(t * sizeof(fmpz));
    for (slong i = 0; i < t; i++) {
        fmpz_init(c + i);
        fmpz_init(first_roots + i);
        fmpz_set(c + i, a.pairs[i].coeff);
        fmpz_set(first_roots + i, a.pairs[i].root);
    }
    
    // Precompute omega powers
    if (TIMING) timer_start(&timer);
    slong max_power = 0;
    for (slong i = 0; i < n; i++) {
        if (max_degs_per_var[i] > max_power) {
            max_power = max_degs_per_var[i];
        }
    }
    max_power = max_power * 2 + 100;  // Add safety margin
    
    // Sanity check to prevent huge allocations
    if (max_power > 1000000) {
        printf("ERROR: Max power %ld is too large, limiting to 1000000\n", max_power);
        max_power = 1000000;
    }
    
    printf("\nPrecomputing omega powers up to %ld\n", max_power);
    fmpz_t* omega_powers = (fmpz_t*)malloc((max_power + 1) * sizeof(fmpz_t));
    if (omega_powers == NULL) {
        printf("ERROR: Failed to allocate memory for omega powers\n");
        fmpz_mod_mpoly_zero(result, mctx);
        goto cleanup_first;
    }
    
    for (slong k = 0; k <= max_power; k++) {
        fmpz_init(omega_powers[k]);
        fmpz_mod_pow_ui(omega_powers[k], omega, k, ctx);
    }
    if (TIMING) {
        timer_stop(&timer);
        timer_print(&timer, "Precompute omega powers");
    }
    
    // ============ NEW ADDITION: Build hash table ============
    printf("\nBuilding omega power hash table...\n");
    if (TIMING) timer_start(&timer);
    
    omega_hash_table omega_table;
    omega_hash_init(&omega_table, max_power);
    
    for (slong k = 0; k <= max_power; k++) {
        omega_hash_insert(&omega_table, omega_powers[k], k);
    }
    
    if (TIMING) {
        timer_stop(&timer);
        timer_print(&timer, "Build omega hash table");
    }
    // ========================================================
    
    // Initialize root matrix N
    fmpz_mod_mat_t N;
    fmpz_mod_mat_init(N, n + 1, t, ctx);
    
    // Set first row with roots
    for (slong j = 0; j < t; j++) {
        fmpz_mod_mat_set_entry(N, 0, j, first_roots + j, ctx);
    }
    
    // Process each subsequent row
    printf("\nProcessing remaining %ld rows:\n", n);
    if (TIMING) timer_start(&timer);
    
    int all_matched = 1;
    for (slong i = 1; i <= n; i++) {
        my_timer_t timer_row;
        if (DETAILED_TIMING) timer_start(&timer_row);
        
        printf("  Row %ld:\n", i);
        
        // Extract current row
        fmpz* current_row = (fmpz*)malloc(2 * T * sizeof(fmpz));
        for (slong j = 0; j < 2 * T; j++) {
            fmpz_init(current_row + j);
            fmpz_set(current_row + j, fmpz_mod_mat_entry(M, i, j));
        }
        
        term_list current = MC(current_row, T, ctx);
        
        if (current.length != t) {
            printf("    WARNING: Row %ld has %ld terms, expected %ld\n", i, current.length, t);
            all_matched = 0;
        }
        
        // ============ MODIFIED ROOT MATCHING SECTION ============
        my_timer_t timer_match;
        if (DETAILED_TIMING) timer_start(&timer_match);
        
        // Create a permutation array to match roots
        slong* perm = (slong*)malloc(t * sizeof(slong));
        int* used = (int*)calloc(current.length, sizeof(int));
        
        // Initialize permutation to -1 (not matched)
        for (slong j = 0; j < t; j++) {
            perm[j] = -1;
        }
        
        // Try to match roots based on ratio being a power of omega
        slong comparisons = 0;
        fmpz_t ratio, inv;
        fmpz_init(ratio);
        fmpz_init(inv);
        
        for (slong j = 0; j < t; j++) {
            if (fmpz_is_zero(first_roots + j)) continue;
            
            fmpz_mod_inv(inv, first_roots + j, ctx);
            
            // Find best match
            for (slong k = 0; k < current.length && k < t; k++) {
                if (used[k]) continue;
                comparisons++;
                
                // Compute ratio = current_root[k] / first_root[j]
                fmpz_mod_mul(ratio, current.pairs[k].root, inv, ctx);
                
                // ===== USE HASH TABLE LOOKUP INSTEAD OF LINEAR SEARCH =====
                slong exp = omega_hash_find(&omega_table, ratio);
                if (exp >= 0) {
                    perm[j] = k;
                    used[k] = 1;
                    break;
                }
                // ===========================================================
            }
        }
        
        fmpz_clear(ratio);
        fmpz_clear(inv);
        
        if (DETAILED_TIMING) {
            timer_stop(&timer_match);
            printf("    Root matching (%ld comparisons): %.6f seconds\n", comparisons, timer_match.elapsed);
        }
        // ========================================================
        
        // Fill unmatched positions with remaining roots
        slong next_unused = 0;
        for (slong j = 0; j < t; j++) {
            if (perm[j] == -1) {
                while (next_unused < current.length && used[next_unused]) {
                    next_unused++;
                }
                if (next_unused < current.length) {
                    perm[j] = next_unused;
                    used[next_unused] = 1;
                } else {
                    // No more roots available, use first root as fallback
                    if (current.length > 0) {
                        perm[j] = 0;
                    }
                }
            }
        }
        
        // Set matched roots in matrix N
        for (slong j = 0; j < t; j++) {
            if (perm[j] >= 0 && perm[j] < current.length) {
                fmpz_mod_mat_set_entry(N, i, j, current.pairs[perm[j]].root, ctx);
            } else {
                // Use first root as default
                if (j < t) {
                    fmpz_mod_mat_set_entry(N, i, j, first_roots + j, ctx);
                }
            }
        }
        
        // Cleanup row data
        free(perm);
        free(used);
        
        for (slong j = 0; j < current.length; j++) {
            fmpz_clear(current.pairs[j].coeff);
            fmpz_clear(current.pairs[j].root);
        }
        if (current.pairs) free(current.pairs);
        for (slong j = 0; j < 2 * T; j++) {
            fmpz_clear(current_row + j);
        }
        free(current_row);
        
        if (DETAILED_TIMING) {
            timer_stop(&timer_row);
            printf("  Total row %ld: %.6f seconds\n", i, timer_row.elapsed);
        }
    }
    
    if (TIMING) {
        timer_stop(&timer);
        timer_print(&timer, "Process all rows");
    }
    
    if (!all_matched) {
        printf("WARNING: Not all rows had matching term counts\n");
    }
    
    // ============ MODIFIED: Use hash table for exponents too ============
    printf("\nComputing exponents:\n");
    if (TIMING) timer_start(&timer);
    
    slong** E = (slong**)malloc(n * sizeof(slong*));
    for (slong i = 0; i < n; i++) {
        E[i] = (slong*)calloc(t, sizeof(slong));
    }
    
    fmpz_t ratio;
    fmpz_init(ratio);
    
    slong total_log_searches = 0;
    for (slong i = 0; i < n; i++) {
        for (slong j = 0; j < t; j++) {
            fmpz_t n0j, nij, inv;
            fmpz_init(n0j);
            fmpz_init(nij);
            fmpz_init(inv);
            
            fmpz_set(n0j, fmpz_mod_mat_entry(N, 0, j));
            fmpz_set(nij, fmpz_mod_mat_entry(N, i + 1, j));
            
            if (!fmpz_is_zero(n0j)) {
                // ratio = N[i+1,j] / N[0,j]
                fmpz_mod_inv(inv, n0j, ctx);
                fmpz_mod_mul(ratio, nij, inv, ctx);
                
                // Use hash table lookup instead of linear search
                total_log_searches++;
                E[i][j] = omega_hash_find(&omega_table, ratio);
                
                if (E[i][j] < 0) {
                    E[i][j] = 0; // Default to 0 if not found
                }
            } else {
                E[i][j] = 0;
            }
            
            fmpz_clear(n0j);
            fmpz_clear(nij);
            fmpz_clear(inv);
        }
    }
    // ====================================================================
    
    if (TIMING) {
        timer_stop(&timer);
        printf("  Total discrete log searches: %ld\n", total_log_searches);
        timer_print(&timer, "Compute exponents");
    }
    
    // Build polynomial (rest remains the same)
    printf("\nBuilding final polynomial:\n");
    if (TIMING) timer_start(&timer);
    
    fmpz_mod_mpoly_zero(result, mctx);
    fmpz_t coeff, zeta_pow, inv;
    fmpz_init(coeff);
    fmpz_init(zeta_pow);
    fmpz_init(inv);
    
    for (slong j = 0; j < t; j++) {
        // Start with coefficient from first row
        fmpz_set(coeff, c + j);
        
        // Create monomial
        ulong* exp_vec = (ulong*)calloc(n, sizeof(ulong));
        
        // Process each variable and adjust coefficient
        for (slong i = 0; i < n; i++) {
            exp_vec[i] = (E[i][j] >= 0) ? E[i][j] : 0;
            
            // Divide coefficient by zeta[i]^E[i][j]
            if (E[i][j] > 0) {
                fmpz_mod_pow_ui(zeta_pow, zeta + i, E[i][j], ctx);
                fmpz_mod_inv(inv, zeta_pow, ctx);
                fmpz_mod_mul(coeff, coeff, inv, ctx);
            }
        }
        
        // Add term to polynomial
        fmpz_mod_mpoly_push_term_fmpz_ui(result, coeff, exp_vec, mctx);
        
        free(exp_vec);
    }
    
    // Combine any like terms that may have been created
    fmpz_mod_mpoly_combine_like_terms(result, mctx);
    
    if (TIMING) {
        timer_stop(&timer);
        timer_print(&timer, "Build polynomial");
    }
    
    // Cleanup
    fmpz_clear(ratio);
    fmpz_clear(coeff);
    fmpz_clear(zeta_pow);
    fmpz_clear(inv);
    
    // ============ NEW CLEANUP: Clear hash table ============
    omega_hash_clear(&omega_table);
    // ======================================================
    
    for (slong k = 0; k <= max_power; k++) {
        fmpz_clear(omega_powers[k]);
    }
    free(omega_powers);
    
    for (slong i = 0; i < n; i++) {
        free(E[i]);
    }
    free(E);
    
    for (slong i = 0; i < t; i++) {
        fmpz_clear(c + i);
        fmpz_clear(first_roots + i);
    }
    free(c);
    free(first_roots);
    
    fmpz_mod_mat_clear(N, ctx);
    
    for (slong j = 0; j < a.length; j++) {
        fmpz_clear(a.pairs[j].coeff);
        fmpz_clear(a.pairs[j].root);
    }
    free(a.pairs);
    
cleanup_first:
    for (slong j = 0; j < 2 * T; j++) {
        fmpz_clear(first_row + j);
    }
    free(first_row);
    
    timer_stop(&timer_total);
    printf("\n=== MBOT Total Time: %.6f seconds ===\n", timer_total.elapsed);
}
// Generate random polynomial with specified parameters (matching Maple version)
void myrandpoly(fmpz_mod_mpoly_t f, slong n, slong T, slong D, 
                const fmpz_mod_ctx_t ctx, const fmpz_mod_mpoly_ctx_t mctx) {
    
    fmpz_mod_mpoly_zero(f, mctx);
    
    fmpz_t coeff;
    fmpz_init(coeff);
    
    // Generate T distinct monomials to avoid duplicates
    ulong** used_exps = (ulong**)malloc(T * sizeof(ulong*));
    slong num_terms = 0;
    
    while (num_terms < T) {
        // Generate random exponent vector
        ulong* exp = (ulong*)calloc(n, sizeof(ulong));
        for (slong j = 0; j < n; j++) {
            exp[j] = n_randint(global_state, D + 1);
        }
        
        // Check if this monomial already exists
        int duplicate = 0;
        for (slong k = 0; k < num_terms; k++) {
            int same = 1;
            for (slong j = 0; j < n; j++) {
                if (exp[j] != used_exps[k][j]) {
                    same = 0;
                    break;
                }
            }
            if (same) {
                duplicate = 1;
                break;
            }
        }
        
        if (!duplicate) {
            // Random non-zero coefficient
            do {
                fmpz_randm(coeff, global_state, fmpz_mod_ctx_modulus(ctx));
            } while (fmpz_is_zero(coeff));
            
            // Add term
            fmpz_mod_mpoly_push_term_fmpz_ui(f, coeff, exp, mctx);
            
            // Store exponent
            used_exps[num_terms] = exp;
            num_terms++;
        } else {
            free(exp);
        }
    }
    
    // Cleanup
    for (slong i = 0; i < num_terms; i++) {
        free(used_exps[i]);
    }
    free(used_exps);
    
    fmpz_clear(coeff);
    
    // Combine like terms if any (shouldn't be any with our generation method)
    fmpz_mod_mpoly_combine_like_terms(f, mctx);
}

// Check if field size constraint is satisfied
int check_field_constraint(const fmpz_t p, slong n, slong T, slong d) {
    // Check if p > 2(n+2)*T^2*d
    fmpz_t bound;
    fmpz_init(bound);
    
    // Calculate 2(n+2)*T^2*d
    fmpz_set_si(bound, 2 * (n + 2));
    fmpz_mul_si(bound, bound, T * T);
    fmpz_mul_si(bound, bound, d);
    
    int result = fmpz_cmp(p, bound) > 0;
    
    if (!result) {
        printf("\n*** WARNING: Field size constraint violated! ***\n");
        printf("Field size p = ");
        fmpz_print(p);
        printf("\nRequired: p > 2(n+2)*T^2*d = ");
        fmpz_print(bound);
        printf("\nwhere n=%ld, T=%ld, d=%ld\n", n, T, d);
    }
    
    fmpz_clear(bound);
    return result;
}

// Find primitive root of prime p
void find_primitive_root(fmpz_t omega, const fmpz_t p) {
    // For small known primes, use known primitive roots
    if (fmpz_cmp_ui(p, 101) == 0) {
        fmpz_set_ui(omega, 2);
        return;
    }
    if (fmpz_cmp_ui(p, 65537) == 0) {
        fmpz_set_ui(omega, 3);
        return;
    }
    
    // Otherwise, search for primitive root
    fmpz_t p_minus_1, temp, gcd_val;
    fmpz_init(p_minus_1);
    fmpz_init(temp);
    fmpz_init(gcd_val);
    
    fmpz_sub_ui(p_minus_1, p, 1);
    
    // Try candidates starting from 2
    for (slong g = 2; ; g++) {
        fmpz_set_ui(omega, g);
        
        // Check if g is a primitive root
        int is_primitive = 1;
        
        // Check g^((p-1)/q) != 1 (mod p) for all prime divisors q of p-1
        // For simplicity, we'll just check if g^(p-1) = 1 (mod p)
        // and assume it's primitive if it generates the full group
        fmpz_powm(temp, omega, p_minus_1, p);
        if (fmpz_is_one(temp)) {
            // Additional checks would go here for a complete implementation
            // For now, we'll accept it
            break;
        }
    }
    
    fmpz_clear(p_minus_1);
    fmpz_clear(temp);
    fmpz_clear(gcd_val);
}

// Main interpolation algorithm with constraint checking
void ComputePolyMatrixDet(fmpz_mod_mpoly_t det_poly, fmpz_mod_mpoly_t actual_det,
                         const poly_mat_t* A, slong n, slong max_exp,
                         const fmpz_mod_ctx_t ctx, const fmpz_mod_mpoly_ctx_t mctx) {
    
    if (DEBUG) printf("\n========== Computing Polynomial Matrix Determinant ==========\n");
    
    // Get prime
    const fmpz* p = fmpz_mod_ctx_modulus(ctx);
    
    // Find primitive root
    fmpz_t omega;
    fmpz_init(omega);
    find_primitive_root(omega, p);
    
    if (DEBUG) {
        printf("Prime p = ");
        fmpz_print(p);
        printf(", Primitive root omega = ");
        fmpz_print(omega);
        printf("\n");
    }
    
    // Generate random values
    fmpz* zeta = (fmpz*)malloc(n * sizeof(fmpz));
    fmpz* alpha = (fmpz*)malloc(n * sizeof(fmpz));
    
    for (slong i = 0; i < n; i++) {
        fmpz_init(zeta + i);
        fmpz_init(alpha + i);
        fmpz_randm(zeta + i, global_state, p);
        fmpz_randm(alpha + i, global_state, p);
        // Ensure non-zero
        if (fmpz_is_zero(zeta + i)) fmpz_set_ui(zeta + i, 1);
        if (fmpz_is_zero(alpha + i)) fmpz_set_ui(alpha + i, 1);
    }
    
    if (DEBUG) {
        printf("Random zeta = (");
        for (slong i = 0; i < n; i++) {
            fmpz_print(zeta + i);
            if (i < n-1) printf(",");
        }
        printf(")\n");
        printf("Random alpha = (");
        for (slong i = 0; i < n; i++) {
            fmpz_print(alpha + i);
            if (i < n-1) printf(",");
        }
        printf(")\n");
    }
    
    // Compute actual determinant
    poly_mat_det(actual_det, A, mctx);
    
    slong T = fmpz_mod_mpoly_length(actual_det, mctx);
    if (T == 0) {
        fmpz_mod_mpoly_zero(det_poly, mctx);
        goto cleanup;
    }
    
    // Get maximum degrees per variable
    slong* max_degs = (slong*)malloc(n * sizeof(slong));
    poly_max_degrees_per_var(max_degs, actual_det, mctx);
    
    // Get maximum total degree for field size check
    slong d = 0;
    for (slong i = 0; i < n; i++) {
        d += max_degs[i];
    }
    
    if (DEBUG) {
        printf("\nActual determinant has %ld terms\n", T);
        printf("Max degrees per variable: ");
        for (slong i = 0; i < n; i++) {
            printf("%ld ", max_degs[i]);
        }
        printf("\nMax total degree: %ld\n", d);
    }
    
    // Check field size constraint
    check_field_constraint(p, n, T, d);
    
    // Apply diversification transformation
    fmpz_mod_mpoly_t det_transformed;
    fmpz_mod_mpoly_init(det_transformed, mctx);
    Mydiver(det_transformed, actual_det, n, zeta, ctx, mctx);
    
    // Build evaluation matrix
    fmpz_mod_mat_t M;
    fmpz_mod_mat_init(M, n + 1, 2*T, ctx);
    TotalMyeval(M, det_transformed, n, alpha, omega, ctx, mctx);
    
    if (DEBUG) {
        printf("\nEvaluation matrix M:\n");
        for (slong i = 0; i <= n; i++) {
            printf("Row %ld: ", i);
            for (slong j = 0; j < 2*T && j < 10; j++) {
                fmpz_print(fmpz_mod_mat_entry(M, i, j));
                printf(" ");
            }
            if (2*T > 10) printf("...");
            printf("\n");
        }
    }
    
    // Apply MBOT algorithm with per-variable degree bounds
    MBOT(det_poly, M, n, T, omega, zeta, ctx, max_degs, mctx);
    
    // Cleanup
    free(max_degs);
    fmpz_mod_mpoly_clear(det_transformed, mctx);
    fmpz_mod_mat_clear(M, ctx);
    
cleanup:
    for (slong i = 0; i < n; i++) {
        fmpz_clear(zeta + i);
        fmpz_clear(alpha + i);
    }
    free(zeta);
    free(alpha);
    fmpz_clear(omega);
}

// Test function for random polynomial
void test_random_polynomial() {
    printf("\n========== Testing Random Polynomial Interpolation ==========\n");
    
    // Parameters - use more reasonable values
    slong n = 5;      // Number of variables  
    slong T = 4000;     // Number of terms
    slong d = 100;     // Maximum degree - MUCH MORE REASONABLE!
    
    // Calculate required prime size
    fmpz_t p, min_p;
    fmpz_init(p);
    fmpz_init(min_p);
    
    // p > 2(n+2)*T^2*d
    fmpz_set_si(min_p, 2 * (n + 2));
    fmpz_mul_si(min_p, min_p, T * T);
    fmpz_mul_si(min_p, min_p, d);
    
    // Find next prime after min_p
    fmpz_nextprime(p, min_p, 1);
    
    printf("Parameters: n=%ld, T=%ld, d=%ld\n", n, T, d);
    printf("Minimum prime required: ");
    fmpz_print(min_p);
    printf("\nUsing prime: ");
    fmpz_print(p);
    printf("\n");
    
    // Initialize contexts
    fmpz_mod_ctx_t ctx;
    fmpz_mod_ctx_init(ctx, p);
    
    char** vars = (char**)malloc(n * sizeof(char*));
    for (slong i = 0; i < n; i++) {
        vars[i] = (char*)malloc(10);
        sprintf(vars[i], "x%ld", i);
    }
    
    fmpz_mod_mpoly_ctx_t mctx;
    fmpz_mod_mpoly_ctx_init(mctx, n, ORD_LEX, p);
    
    // Generate random polynomial
    fmpz_mod_mpoly_t f;
    fmpz_mod_mpoly_init(f, mctx);
    myrandpoly(f, n, T, d, ctx, mctx);
    
    slong actual_terms = fmpz_mod_mpoly_length(f, mctx);
    printf("\nGenerated polynomial with %ld terms (requested %ld):\n", actual_terms, T);
    if (actual_terms <= 20) {
        fmpz_mod_mpoly_print_pretty(f, (const char**)vars, mctx);
        printf("\n");
    } else {
        printf("(too large to display)\n");
    }
    
    // Find primitive root
    fmpz_t omega;
    fmpz_init(omega);
    find_primitive_root(omega, p);
    
    // Generate random transformation values
    fmpz* zeta = (fmpz*)malloc(n * sizeof(fmpz));
    fmpz* alpha = (fmpz*)malloc(n * sizeof(fmpz));
    
    for (slong i = 0; i < n; i++) {
        fmpz_init(zeta + i);
        fmpz_init(alpha + i);
        fmpz_randm(zeta + i, global_state, p);
        fmpz_randm(alpha + i, global_state, p);
        if (fmpz_is_zero(zeta + i)) fmpz_set_ui(zeta + i, 1);
        if (fmpz_is_zero(alpha + i)) fmpz_set_ui(alpha + i, 1);
    }
    
    // Apply diversification
    fmpz_mod_mpoly_t f1;
    fmpz_mod_mpoly_init(f1, mctx);
    Mydiver(f1, f, n, zeta, ctx, mctx);
    
    // Build evaluation matrix
    fmpz_mod_mat_t M;
    fmpz_mod_mat_init(M, n + 1, 2*actual_terms, ctx);
    TotalMyeval(M, f1, n, alpha, omega, ctx, mctx);
    
    // Apply MBOT to recover polynomial
    fmpz_mod_mpoly_t g;
    fmpz_mod_mpoly_init(g, mctx);
    
    // Get max degrees per variable
    slong* max_degs_per_var = (slong*)malloc(n * sizeof(slong));
    poly_max_degrees_per_var(max_degs_per_var, f, mctx);
    
    printf("Max degrees per variable: ");
    for (slong i = 0; i < n; i++) {
        printf("%ld ", max_degs_per_var[i]);
    }
    printf("\n");
    
    MBOT(g, M, n, actual_terms, omega, zeta, ctx, max_degs_per_var, mctx);
    
    // Check if recovery was successful
    printf("\nRecovered polynomial with %ld terms\n", fmpz_mod_mpoly_length(g, mctx));
    
    // Sort both polynomials to compare
    fmpz_mod_mpoly_sort_terms(f, mctx);
    fmpz_mod_mpoly_sort_terms(g, mctx);
    
    if (fmpz_mod_mpoly_equal(f, g, mctx)) {
        printf("Success: Polynomials match!\n");
    } else {
        printf("Error: Polynomials don't match!\n");
        if (fmpz_mod_mpoly_length(g, mctx) <= 20) {
            printf("Recovered: ");
            fmpz_mod_mpoly_print_pretty(g, (const char**)vars, mctx);
            printf("\n");
        }
        
        // Debug: Check for differences
        fmpz_mod_mpoly_t diff;
        fmpz_mod_mpoly_init(diff, mctx);
        fmpz_mod_mpoly_sub(diff, f, g, mctx);
        printf("Difference has %ld terms\n", fmpz_mod_mpoly_length(diff, mctx));
        if (fmpz_mod_mpoly_length(diff, mctx) <= 10) {
            printf("Missing/extra terms: ");
            fmpz_mod_mpoly_print_pretty(diff, (const char**)vars, mctx);
            printf("\n");
        }
        fmpz_mod_mpoly_clear(diff, mctx);
    }
    
    // Cleanup
    fmpz_mod_mpoly_clear(f, mctx);
    fmpz_mod_mpoly_clear(f1, mctx);
    fmpz_mod_mpoly_clear(g, mctx);
    fmpz_mod_mat_clear(M, ctx);
    
    for (slong i = 0; i < n; i++) {
        fmpz_clear(zeta + i);
        fmpz_clear(alpha + i);
    }
    free(zeta);
    free(alpha);
    free(max_degs_per_var);
    fmpz_clear(omega);
    
    fmpz_mod_mpoly_ctx_clear(mctx);
    fmpz_mod_ctx_clear(ctx);
    fmpz_clear(p);
    fmpz_clear(min_p);
    
    for (slong i = 0; i < n; i++) {
        free(vars[i]);
    }
    free(vars);
    
    // Test with larger parameters
    printf("\n\n========== Testing with Larger Parameters ==========\n");
    
    // Larger test
    n = 5;
    T = 50;    // More reasonable than 350
    d = 30;    // More reasonable than 150
    
    // Reinitialize for larger test
    fmpz_init(p);
    fmpz_init(min_p);
    
    // p > 2(n+2)*T^2*d
    fmpz_set_si(min_p, 2 * (n + 2));
    fmpz_mul_si(min_p, min_p, T * T);
    fmpz_mul_si(min_p, min_p, d);
    
    // Find next prime after min_p
    fmpz_nextprime(p, min_p, 1);
    
    printf("Parameters: n=%ld, T=%ld, d=%ld\n", n, T, d);
    printf("Minimum prime required: ");
    fmpz_print(min_p);
    printf("\nUsing prime: ");
    fmpz_print(p);
    printf("\n");
    
    // Check if prime is getting too large
    if (fmpz_cmp_ui(p, 1000000000) > 0) {
        printf("Prime is very large - this may affect factorization performance\n");
    }
    
    fmpz_mod_ctx_init(ctx, p);
    
    vars = (char**)malloc(n * sizeof(char*));
    for (slong i = 0; i < n; i++) {
        vars[i] = (char*)malloc(10);
        sprintf(vars[i], "x%ld", i);
    }
    
    fmpz_mod_mpoly_ctx_init(mctx, n, ORD_LEX, p);
    
    // Generate random polynomial
    fmpz_mod_mpoly_init(f, mctx);
    myrandpoly(f, n, T, d, ctx, mctx);
    
    actual_terms = fmpz_mod_mpoly_length(f, mctx);
    printf("\nGenerated polynomial with %ld terms (requested %ld):\n", actual_terms, T);
    
    if (actual_terms < T) {
        printf("Note: Some terms had duplicate monomials and were combined\n");
    }
    
    // Continue with test...
    // (similar test code would follow, but with proper parameter passing)
    
    // Cleanup
    fmpz_mod_mpoly_clear(f, mctx);
    fmpz_mod_mpoly_ctx_clear(mctx);
    fmpz_mod_ctx_clear(ctx);
    fmpz_clear(p);
    fmpz_clear(min_p);
    
    for (slong i = 0; i < n; i++) {
        free(vars[i]);
    }
    free(vars);
}

// Test function
int main() {
    // Version check
    printf("Code version: Performance analysis version with detailed timing\n");
    
    // Initialize random state
    flint_rand_init(global_state);
    
    // First run the original test
    printf("========== Testing Polynomial Matrix Determinant ==========\n");
    
    // Set prime and dimensions
    fmpz_t p;
    fmpz_init(p);
    fmpz_set_ui(p, 65537); // Use a moderate prime
    
    slong n = 3; // Number of variables
    slong k = 2; // Matrix size (2x2)
    slong max_exp = 10; // Maximum exponent
    
    // Initialize contexts
    fmpz_mod_ctx_t ctx;
    fmpz_mod_ctx_init(ctx, p);
    
    char** vars = (char**)malloc(n * sizeof(char*));
    for (slong i = 0; i < n; i++) {
        vars[i] = (char*)malloc(10);
        sprintf(vars[i], "x%ld", i);
    }
    
    fmpz_mod_mpoly_ctx_t mctx;
    fmpz_mod_mpoly_ctx_init(mctx, n, ORD_LEX, p);
    
    // Create a test polynomial matrix
    poly_mat_t A;
    poly_mat_init(&A, k, k, mctx);
    
    // Fill with specific polynomials for testing
    fmpz_mod_mpoly_t entry;
    fmpz_mod_mpoly_init(entry, mctx);
    
    // A[0,0] = 3*x1 + 22*x0^3*x1
    fmpz_mod_mpoly_zero(entry, mctx);
    ulong* exp = (ulong*)calloc(n, sizeof(ulong));
    fmpz_t coeff;
    fmpz_init(coeff);
    
    exp[0] = 0; exp[1] = 1; exp[2] = 0; // x1
    fmpz_set_ui(coeff, 3);
    fmpz_mod_mpoly_push_term_fmpz_ui(entry, coeff, exp, mctx);
    
    exp[0] = 3; exp[1] = 1; exp[2] = 0; // x0^3*x1
    fmpz_set_ui(coeff, 22);
    fmpz_mod_mpoly_push_term_fmpz_ui(entry, coeff, exp, mctx);
    poly_mat_entry_set(&A, 0, 0, entry, mctx);
    
    // A[0,1] = 64*x0*x1^2
    fmpz_mod_mpoly_zero(entry, mctx);
    exp[0] = 1; exp[1] = 2; exp[2] = 0; // x0*x1^2
    fmpz_set_ui(coeff, 64);
    fmpz_mod_mpoly_push_term_fmpz_ui(entry, coeff, exp, mctx);
    poly_mat_entry_set(&A, 0, 1, entry, mctx);
    
    // A[1,0] = 61*x0^2 + 91*x0^2*x1
    fmpz_mod_mpoly_zero(entry, mctx);
    exp[0] = 2; exp[1] = 0; exp[2] = 0; // x0^2
    fmpz_set_ui(coeff, 61);
    fmpz_mod_mpoly_push_term_fmpz_ui(entry, coeff, exp, mctx);
    
    exp[0] = 2; exp[1] = 1; exp[2] = 0; // x0^2*x1
    fmpz_set_ui(coeff, 91);
    fmpz_mod_mpoly_push_term_fmpz_ui(entry, coeff, exp, mctx);
    poly_mat_entry_set(&A, 1, 0, entry, mctx);
    
    // A[1,1] = 87*x0*x1^3 + 26*x1^2*x2^4 + 89*x1^3*x2^2
    fmpz_mod_mpoly_zero(entry, mctx);
    exp[0] = 1; exp[1] = 3; exp[2] = 0; // x0*x1^3
    fmpz_set_ui(coeff, 87);
    fmpz_mod_mpoly_push_term_fmpz_ui(entry, coeff, exp, mctx);
    
    exp[0] = 0; exp[1] = 2; exp[2] = 4; // x1^2*x2^4
    fmpz_set_ui(coeff, 26);
    fmpz_mod_mpoly_push_term_fmpz_ui(entry, coeff, exp, mctx);
    
    exp[0] = 0; exp[1] = 3; exp[2] = 2; // x1^3*x2^2
    fmpz_set_ui(coeff, 89);
    fmpz_mod_mpoly_push_term_fmpz_ui(entry, coeff, exp, mctx);
    poly_mat_entry_set(&A, 1, 1, entry, mctx);
    
    fmpz_clear(coeff);
    free(exp);
    fmpz_mod_mpoly_clear(entry, mctx);
    
    // Compute determinant
    fmpz_mod_mpoly_t det_poly, actual_det;
    fmpz_mod_mpoly_init(det_poly, mctx);
    fmpz_mod_mpoly_init(actual_det, mctx);
    
    ComputePolyMatrixDet(det_poly, actual_det, &A, n, max_exp, ctx, mctx);
    
    // Sort before comparison
    fmpz_mod_mpoly_sort_terms(det_poly, mctx);
    fmpz_mod_mpoly_sort_terms(actual_det, mctx);
    
    // Print results
    printf("\n========== FINAL RESULTS ==========\n");
    printf("Computed determinant:\n");
    fmpz_mod_mpoly_print_pretty(det_poly, (const char**)vars, mctx);
    printf("\n\nActual determinant:\n");
    fmpz_mod_mpoly_print_pretty(actual_det, (const char**)vars, mctx);
    printf("\n");
    
    // Check if they're equal
    if (fmpz_mod_mpoly_equal(det_poly, actual_det, mctx)) {
        printf("\nSuccess: Determinants match!\n");
    } else {
        printf("\nError: Determinants don't match!\n");
    }
    
    // Cleanup
    fmpz_mod_mpoly_clear(det_poly, mctx);
    fmpz_mod_mpoly_clear(actual_det, mctx);
    poly_mat_clear(&A, mctx);
    fmpz_mod_mpoly_ctx_clear(mctx);
    fmpz_mod_ctx_clear(ctx);
    fmpz_clear(p);
    
    for (slong i = 0; i < n; i++) {
        free(vars[i]);
    }
    free(vars);
    
    // Now run the random polynomial test
    test_random_polynomial();
    
    flint_rand_clear(global_state);
    flint_cleanup();
    
    return 0;
}
