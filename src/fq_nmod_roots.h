#ifndef FQ_NMOD_ROOTS_H
#define FQ_NMOD_ROOTS_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <flint/flint.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_poly.h>
#include <flint/fq_nmod_poly_factor.h>
#include <flint/fmpz.h>
#include <gmp.h>
#include <flint/ulong_extras.h>

#define PRIME 9223372036854775783ULL
#define SMALL_PRIME 1073741827ULL

// ========== nmod_poly 版本的根存储结构 ==========
typedef struct {
    mp_limb_t *roots;
    slong *mult;
    slong num;
    slong alloc;
} nmod_roots_struct;
typedef nmod_roots_struct nmod_roots_t[1];

void nmod_roots_init(nmod_roots_t roots) {
    roots->alloc = 4;
    roots->roots = (mp_limb_t*)malloc(roots->alloc * sizeof(mp_limb_t));
    roots->mult = (slong*)malloc(roots->alloc * sizeof(slong));
    roots->num = 0;
}

void nmod_roots_clear(nmod_roots_t roots) {
    free(roots->roots);
    free(roots->mult);
    roots->num = 0;
    roots->alloc = 0;
}

void nmod_roots_fit_length(nmod_roots_t roots, slong len) {
    if (len > roots->alloc) {
        roots->alloc = FLINT_MAX(len, 2 * roots->alloc);
        roots->roots = (mp_limb_t*)realloc(roots->roots, roots->alloc * sizeof(mp_limb_t));
        roots->mult = (slong*)realloc(roots->mult, roots->alloc * sizeof(slong));
    }
}

void nmod_roots_add(nmod_roots_t roots, mp_limb_t root, slong mult) {
    nmod_roots_fit_length(roots, roots->num + 1);
    roots->roots[roots->num] = root;
    roots->mult[roots->num] = mult;
    roots->num++;
}

// ========== fq_nmod_poly 版本的根存储结构 ==========
typedef struct {
    fq_nmod_struct *roots;
    slong *mult;
    slong num;
    slong alloc;
    mp_limb_t p;
    slong d;
} fq_nmod_roots_struct;
typedef fq_nmod_roots_struct fq_nmod_roots_t[1];

void fq_nmod_roots_init(fq_nmod_roots_t roots, const fq_nmod_ctx_t ctx) {
    roots->alloc = 4;
    roots->roots = (fq_nmod_struct*)malloc(roots->alloc * sizeof(fq_nmod_struct));
    roots->mult = (slong*)malloc(roots->alloc * sizeof(slong));
    roots->num = 0;
    roots->p = fq_nmod_ctx_prime(ctx);
    roots->d = fq_nmod_ctx_degree(ctx);
    
    // 初始化根元素
    for (int i = 0; i < roots->alloc; i++) {
        fq_nmod_init(roots->roots + i, ctx);
    }
}

void fq_nmod_roots_clear(fq_nmod_roots_t roots, const fq_nmod_ctx_t ctx) {
    for (int i = 0; i < roots->alloc; i++) {
        fq_nmod_clear(roots->roots + i, ctx);
    }
    free(roots->roots);
    free(roots->mult);
    roots->num = 0;
    roots->alloc = 0;
}

void fq_nmod_roots_fit_length(fq_nmod_roots_t roots, slong len, const fq_nmod_ctx_t ctx) {
    if (len > roots->alloc) {
        slong old_alloc = roots->alloc;
        roots->alloc = FLINT_MAX(len, 2 * roots->alloc);
        roots->roots = (fq_nmod_struct*)realloc(roots->roots, roots->alloc * sizeof(fq_nmod_struct));
        roots->mult = (slong*)realloc(roots->mult, roots->alloc * sizeof(slong));
        
        // 初始化新的元素
        for (slong i = old_alloc; i < roots->alloc; i++) {
            fq_nmod_init(roots->roots + i, ctx);
        }
    }
}

void fq_nmod_roots_add(fq_nmod_roots_t roots, const fq_nmod_t root, slong mult, const fq_nmod_ctx_t ctx) {
    fq_nmod_roots_fit_length(roots, roots->num + 1, ctx);
    fq_nmod_set(roots->roots + roots->num, root, ctx);
    roots->mult[roots->num] = mult;
    roots->num++;
}

// ========== 通用工具函数 ==========
double get_time_roots() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1000000.0;
}

// ========== nmod_poly 版本实现 ==========

void nmod_simple_power_x_mod(nmod_poly_t result, mp_limb_t exp, const nmod_poly_t modulus) {
    if (exp == 0) {
        nmod_poly_one(result);
        return;
    }
    
    nmod_poly_t base;
    nmod_poly_init(base, modulus->mod.n);
    
    nmod_poly_set_coeff_ui(base, 1, 1);  // base = x
    nmod_poly_one(result);  // result = 1
    
    mp_limb_t e = exp;
    while (e > 0) {
        if (e & 1) {
            nmod_poly_mulmod(result, result, base, modulus);
        }
        e >>= 1;
        if (e > 0) {
            nmod_poly_mulmod(base, base, base, modulus);
        }
    }
    
    nmod_poly_clear(base);
}

void nmod_extract_linear_factors(nmod_roots_t roots, const nmod_poly_t poly, flint_rand_t state) {
    slong deg = nmod_poly_degree(poly);
    
    if (deg <= 0) return;
    
    if (deg == 1) {
        mp_limb_t a = nmod_poly_get_coeff_ui(poly, 1);
        mp_limb_t b = nmod_poly_get_coeff_ui(poly, 0);
        
        if (a != 0) {
            mp_limb_t inv_a = n_invmod(a, poly->mod.n);
            mp_limb_t root = n_mulmod2_preinv(poly->mod.n - b, inv_a, poly->mod.n, poly->mod.ninv);
            nmod_roots_add(roots, root, 1);
        }
        return;
    }
    
    // 随机化分离
    nmod_poly_t g, h, r;
    nmod_poly_init(g, poly->mod.n);
    nmod_poly_init(h, poly->mod.n);
    nmod_poly_init(r, poly->mod.n);
    
    for (int attempt = 0; attempt < 10; attempt++) {
        nmod_poly_randtest(g, state, deg);
        
        mp_limb_t exp = (poly->mod.n - 1) / 2;
        nmod_poly_powmod_ui_binexp(h, g, exp, poly);
        
        mp_limb_t constant = nmod_poly_get_coeff_ui(h, 0);
        nmod_poly_set_coeff_ui(h, 0, n_submod(constant, 1, poly->mod.n));
        
        nmod_poly_gcd(r, h, poly);
        slong r_deg = nmod_poly_degree(r);
        
        if (r_deg > 0 && r_deg < deg) {
            nmod_extract_linear_factors(roots, r, state);
            
            nmod_poly_t quotient;
            nmod_poly_init(quotient, poly->mod.n);
            nmod_poly_div(quotient, poly, r);
            nmod_extract_linear_factors(roots, quotient, state);
            nmod_poly_clear(quotient);
            break;
        }
    }
    
    nmod_poly_clear(g);
    nmod_poly_clear(h);
    nmod_poly_clear(r);
}

slong our_nmod_poly_roots(nmod_roots_t roots, const nmod_poly_t poly, int with_multiplicity) {
    slong deg = nmod_poly_degree(poly);
    
    if (deg <= 0) return 0;
    
    if (deg == 1) {
        mp_limb_t a = nmod_poly_get_coeff_ui(poly, 1);
        mp_limb_t b = nmod_poly_get_coeff_ui(poly, 0);
        
        if (a != 0) {
            mp_limb_t inv_a = n_invmod(a, poly->mod.n);
            mp_limb_t root = n_mulmod2_preinv(poly->mod.n - b, inv_a, poly->mod.n, poly->mod.ninv);
            nmod_roots_add(roots, root, 1);
            return 1;
        }
        return 0;
    }
    
    // CZ算法核心
    nmod_poly_t x_to_p, x, frobenius, root_poly;
    nmod_poly_init(x_to_p, poly->mod.n);
    nmod_poly_init(x, poly->mod.n);
    nmod_poly_init(frobenius, poly->mod.n);
    nmod_poly_init(root_poly, poly->mod.n);
    
    nmod_simple_power_x_mod(x_to_p, poly->mod.n, poly);
    nmod_poly_set_coeff_ui(x, 1, 1);
    nmod_poly_sub(frobenius, x_to_p, x);
    nmod_poly_gcd(root_poly, poly, frobenius);
    
    if (nmod_poly_degree(root_poly) > 0) {
        flint_rand_t state;
        flint_rand_init(state);
        nmod_extract_linear_factors(roots, root_poly, state);
        flint_rand_clear(state);
    }
    
    nmod_poly_clear(x_to_p);
    nmod_poly_clear(x);
    nmod_poly_clear(frobenius);
    nmod_poly_clear(root_poly);
    
    return roots->num;
}

// ========== fq_nmod_poly 版本实现 ==========

void fq_nmod_simple_power_x_mod(fq_nmod_poly_t result, const fmpz_t exp, const fq_nmod_poly_t modulus, const fq_nmod_ctx_t ctx) {
    if (fmpz_is_zero(exp)) {
        fq_nmod_poly_one(result, ctx);
        return;
    }
    
    fq_nmod_poly_t base;
    fq_nmod_poly_init(base, ctx);
    
    fq_nmod_poly_gen(base, ctx);  // base = x
    fq_nmod_poly_one(result, ctx);  // result = 1
    
    fmpz_t e;
    fmpz_init_set(e, exp);
    
    while (!fmpz_is_zero(e)) {
        if (fmpz_is_odd(e)) {
            fq_nmod_poly_mulmod(result, result, base, modulus, ctx);
        }
        fmpz_fdiv_q_2exp(e, e, 1);
        if (!fmpz_is_zero(e)) {
            fq_nmod_poly_mulmod(base, base, base, modulus, ctx);
        }
    }
    
    fmpz_clear(e);
    fq_nmod_poly_clear(base, ctx);
}

void fq_nmod_extract_linear_factors(fq_nmod_roots_t roots, const fq_nmod_poly_t poly, flint_rand_t state, const fq_nmod_ctx_t ctx) {
    slong deg = fq_nmod_poly_degree(poly, ctx);
    
    if (deg <= 0) return;
    
    if (deg == 1) {
        // 对于线性多项式 ax + b，根是 -b/a
        fq_nmod_t a, b, root;
        fq_nmod_init(a, ctx);
        fq_nmod_init(b, ctx);
        fq_nmod_init(root, ctx);
        
        fq_nmod_poly_get_coeff(a, poly, 1, ctx);
        fq_nmod_poly_get_coeff(b, poly, 0, ctx);
        
        if (!fq_nmod_is_zero(a, ctx)) {
            fq_nmod_inv(root, a, ctx);  // root = 1/a
            fq_nmod_neg(b, b, ctx);     // b = -b
            fq_nmod_mul(root, root, b, ctx);  // root = -b/a
            fq_nmod_roots_add(roots, root, 1, ctx);
        }
        
        fq_nmod_clear(a, ctx);
        fq_nmod_clear(b, ctx);
        fq_nmod_clear(root, ctx);
        return;
    }
    
    // 随机化分离
    fq_nmod_poly_t g, h, r, one_poly;
    fq_nmod_poly_init(g, ctx);
    fq_nmod_poly_init(h, ctx);
    fq_nmod_poly_init(r, ctx);
    fq_nmod_poly_init(one_poly, ctx);
    
    fmpz_t exp;
    fmpz_init(exp);
    fq_nmod_ctx_order(exp, ctx);
    fmpz_sub_ui(exp, exp, 1);
    fmpz_fdiv_q_2exp(exp, exp, 1);  // exp = (q-1)/2
    
    // 设置常数1多项式
    fq_nmod_t one;
    fq_nmod_init(one, ctx);
    fq_nmod_one(one, ctx);
    fq_nmod_poly_set_fq_nmod(one_poly, one, ctx);
    
    for (int attempt = 0; attempt < 10; attempt++) {
        fq_nmod_poly_randtest(g, state, deg, ctx);
        
        fq_nmod_poly_powmod_fmpz_binexp(h, g, exp, poly, ctx);
        
        // h = h - 1
        fq_nmod_poly_sub(h, h, one_poly, ctx);
        
        fq_nmod_poly_gcd(r, h, poly, ctx);
        slong r_deg = fq_nmod_poly_degree(r, ctx);
        
        if (r_deg > 0 && r_deg < deg) {
            fq_nmod_extract_linear_factors(roots, r, state, ctx);
            
            fq_nmod_poly_t quotient;
            fq_nmod_poly_init(quotient, ctx);
            fq_nmod_poly_div(quotient, poly, r, ctx);
            fq_nmod_extract_linear_factors(roots, quotient, state, ctx);
            fq_nmod_poly_clear(quotient, ctx);
            break;
        }
    }
    
    fmpz_clear(exp);
    fq_nmod_clear(one, ctx);
    fq_nmod_poly_clear(g, ctx);
    fq_nmod_poly_clear(h, ctx);
    fq_nmod_poly_clear(r, ctx);
    fq_nmod_poly_clear(one_poly, ctx);
}

slong our_fq_nmod_poly_roots(fq_nmod_roots_t roots, const fq_nmod_poly_t poly, int with_multiplicity, const fq_nmod_ctx_t ctx) {
    slong deg = fq_nmod_poly_degree(poly, ctx);
    
    if (deg <= 0) return 0;
    
    if (deg == 1) {
        fq_nmod_t a, b, root;
        fq_nmod_init(a, ctx);
        fq_nmod_init(b, ctx);
        fq_nmod_init(root, ctx);
        
        fq_nmod_poly_get_coeff(a, poly, 1, ctx);
        fq_nmod_poly_get_coeff(b, poly, 0, ctx);
        
        if (!fq_nmod_is_zero(a, ctx)) {
            fq_nmod_inv(root, a, ctx);
            fq_nmod_neg(b, b, ctx);
            fq_nmod_mul(root, root, b, ctx);
            fq_nmod_roots_add(roots, root, 1, ctx);
            
            fq_nmod_clear(a, ctx);
            fq_nmod_clear(b, ctx);
            fq_nmod_clear(root, ctx);
            return 1;
        }
        
        fq_nmod_clear(a, ctx);
        fq_nmod_clear(b, ctx);
        fq_nmod_clear(root, ctx);
        return 0;
    }
    
    // CZ算法核心
    fq_nmod_poly_t x_to_q, x, frobenius, root_poly;
    fq_nmod_poly_init(x_to_q, ctx);
    fq_nmod_poly_init(x, ctx);
    fq_nmod_poly_init(frobenius, ctx);
    fq_nmod_poly_init(root_poly, ctx);
    
    fmpz_t q;
    fmpz_init(q);
    fq_nmod_ctx_order(q, ctx);  // q = p^n
    
    fq_nmod_simple_power_x_mod(x_to_q, q, poly, ctx);
    fq_nmod_poly_gen(x, ctx);  // x = X
    fq_nmod_poly_sub(frobenius, x_to_q, x, ctx);  // frobenius = x^q - x
    fq_nmod_poly_gcd(root_poly, poly, frobenius, ctx);
    
    if (fq_nmod_poly_degree(root_poly, ctx) > 0) {
        flint_rand_t state;
        flint_rand_init(state);
        fq_nmod_extract_linear_factors(roots, root_poly, state, ctx);
        flint_rand_clear(state);
    }
    
    fmpz_clear(q);
    fq_nmod_poly_clear(x_to_q, ctx);
    fq_nmod_poly_clear(x, ctx);
    fq_nmod_poly_clear(frobenius, ctx);
    fq_nmod_poly_clear(root_poly, ctx);
    
    return roots->num;
}

// ========== 基准测试函数 ==========

void generate_nmod_poly(nmod_poly_t poly, flint_rand_t state, slong degree, mp_limb_t p) {
    nmod_poly_zero(poly);
    mp_limb_t lead_coeff = n_randint(state, p - 1) + 1;
    nmod_poly_set_coeff_ui(poly, degree, lead_coeff);
    
    for (slong i = 0; i < degree; i++) {
        mp_limb_t coeff = n_randint(state, p);
        nmod_poly_set_coeff_ui(poly, i, coeff);
    }
}

void generate_fq_nmod_poly(fq_nmod_poly_t poly, flint_rand_t state, slong degree, const fq_nmod_ctx_t ctx) {
    fq_nmod_poly_randtest_monic(poly, state, degree + 1, ctx);
}

double benchmark_nmod_roots(slong degree, int num_tests) {
    printf("\n=== nmod_poly CZ根查找测试 ===\n");
    
    nmod_poly_t poly;
    nmod_roots_t roots;
    flint_rand_t state;
    
    nmod_poly_init(poly, PRIME);
    nmod_roots_init(roots);
    flint_rand_init(state);
    
    double total_time = 0.0;
    slong total_roots = 0;
    
    for (int i = 0; i < num_tests; i++) {
        generate_nmod_poly(poly, state, degree, PRIME);
        
        roots->num = 0;
        double start = get_time_roots();
        slong num_roots = our_nmod_poly_roots(roots, poly, 1);
        double end = get_time_roots();
        
        printf("nmod测试 %d: %.6f 秒，找到 %ld 个根\n", i + 1, end - start, num_roots);
        total_time += (end - start);
        total_roots += num_roots;
    }
    
    double avg = total_time / num_tests;
    printf("nmod平均: %.6f 秒，平均根数: %.1f\n", avg, (double)total_roots / num_tests);
    
    nmod_poly_clear(poly);
    nmod_roots_clear(roots);
    flint_rand_clear(state);
    
    return avg;
}

double benchmark_fq_nmod_roots(slong degree, slong extension, int num_tests) {
    printf("\n=== fq_nmod_poly CZ根查找测试 (F_{%llu^%ld}) ===\n", SMALL_PRIME, extension);
    
    fq_nmod_ctx_t ctx;
    fq_nmod_poly_t poly;
    fq_nmod_roots_t roots;
    flint_rand_t state;
    
    // 初始化有限域上下文
    fq_nmod_ctx_init_ui(ctx, SMALL_PRIME, extension, "a");
    fq_nmod_poly_init(poly, ctx);
    fq_nmod_roots_init(roots, ctx);
    flint_rand_init(state);
    
    double total_time = 0.0;
    slong total_roots = 0;
    
    for (int i = 0; i < num_tests; i++) {
        generate_fq_nmod_poly(poly, state, degree, ctx);
        
        roots->num = 0;
        double start = get_time_roots();
        slong num_roots = our_fq_nmod_poly_roots(roots, poly, 1, ctx);
        double end = get_time_roots();
        
        printf("fq_nmod测试 %d: %.6f 秒，找到 %ld 个根\n", i + 1, end - start, num_roots);
        total_time += (end - start);
        total_roots += num_roots;
    }
    
    double avg = total_time / num_tests;
    printf("fq_nmod平均: %.6f 秒，平均根数: %.1f\n", avg, (double)total_roots / num_tests);
    
    fq_nmod_poly_clear(poly, ctx);
    fq_nmod_roots_clear(roots, ctx);
    fq_nmod_ctx_clear(ctx);
    flint_rand_clear(state);
    
    return avg;
}

double benchmark_flint_fq_nmod_factor(slong degree, slong extension, int num_tests) {
    printf("\n=== FLINT fq_nmod_poly因式分解测试 (F_{%llu^%ld}) ===\n", SMALL_PRIME, extension);
    
    fq_nmod_ctx_t ctx;
    fq_nmod_poly_t poly;
    fq_nmod_poly_factor_t factors;
    flint_rand_t state;
    
    fq_nmod_ctx_init_ui(ctx, SMALL_PRIME, extension, "a");
    fq_nmod_poly_init(poly, ctx);
    fq_nmod_poly_factor_init(factors, ctx);
    flint_rand_init(state);
    
    double total_time = 0.0;
    slong total_roots = 0;
    
    for (int i = 0; i < num_tests; i++) {
        generate_fq_nmod_poly(poly, state, degree, ctx);
        
        fq_nmod_poly_factor_clear(factors, ctx);
        fq_nmod_poly_factor_init(factors, ctx);
        
        fq_nmod_t lead_coeff;
        fq_nmod_init(lead_coeff, ctx);
        
        double start = get_time_roots();
        fq_nmod_poly_factor(factors, lead_coeff, poly, ctx);  // 正确的参数顺序
        double end = get_time_roots();
        
        slong linear_factors = 0;
        for (slong j = 0; j < factors->num; j++) {
            if (fq_nmod_poly_degree(factors->poly + j, ctx) == 1) {
                linear_factors += factors->exp[j];
            }
        }
        
        printf("FLINT fq_nmod测试 %d: %.6f 秒，找到 %ld 个根\n", i + 1, end - start, linear_factors);
        total_time += (end - start);
        total_roots += linear_factors;
        
        fq_nmod_clear(lead_coeff, ctx);
    }
    
    double avg = total_time / num_tests;
    printf("FLINT fq_nmod平均: %.6f 秒，平均根数: %.1f\n", avg, (double)total_roots / num_tests);
    
    fq_nmod_poly_clear(poly, ctx);
    fq_nmod_poly_factor_clear(factors, ctx);
    fq_nmod_ctx_clear(ctx);
    flint_rand_clear(state);
    
    return avg;
}

void test_fq_nmod_correctness() {
    printf("\n=== fq_nmod_poly 正确性验证测试 ===\n");
    
    fq_nmod_ctx_t ctx;
    fq_nmod_poly_t poly;
    fq_nmod_roots_t our_roots;
    flint_rand_t state;
    
    // 使用 F_5^2
    fq_nmod_ctx_init_ui(ctx, 5, 2, "a");
    fq_nmod_poly_init(poly, ctx);
    fq_nmod_roots_init(our_roots, ctx);
    flint_rand_init(state);
    
    printf("测试域: F_{5^2} = F_25\n");
    
    // 构造一个简单的多项式：x^2 - 1 = (x-1)(x+1)
    fq_nmod_poly_zero(poly, ctx);
    fq_nmod_t one, neg_one;
    fq_nmod_init(one, ctx);
    fq_nmod_init(neg_one, ctx);
    fq_nmod_one(one, ctx);
    fq_nmod_set_ui(neg_one, 4, ctx);  // -1 在 F_5 中是 4
    
    fq_nmod_poly_set_coeff(poly, 2, one, ctx);      // x^2
    fq_nmod_poly_set_coeff(poly, 0, neg_one, ctx);  // -1 (在F_5中，-1 = 4)
    
    printf("测试多项式: x^2 - 1\n");
    printf("期望根: 1, 4 (在F_5中)\n");
    
    our_roots->num = 0;
    slong num_roots = our_fq_nmod_poly_roots(our_roots, poly, 1, ctx);
    
    printf("我们找到的根: ");
    for (slong i = 0; i < our_roots->num; i++) {
        fq_nmod_print_pretty(our_roots->roots + i, ctx);
        printf(" ");
    }
    printf("(共 %ld 个)\n", num_roots);
    
    // 验证根的正确性
    for (slong i = 0; i < our_roots->num; i++) {
        fq_nmod_t value;
        fq_nmod_init(value, ctx);
        fq_nmod_poly_evaluate_fq_nmod(value, poly, our_roots->roots + i, ctx);
        printf("验证 f(");
        fq_nmod_print_pretty(our_roots->roots + i, ctx);
        printf(") = ");
        fq_nmod_print_pretty(value, ctx);
        printf("\n");
        fq_nmod_clear(value, ctx);
    }
    
    fq_nmod_clear(one, ctx);
    fq_nmod_clear(neg_one, ctx);
    fq_nmod_poly_clear(poly, ctx);
    fq_nmod_roots_clear(our_roots, ctx);
    fq_nmod_ctx_clear(ctx);
    flint_rand_clear(state);
}

void run_unified_comparison() {
    printf("=== 修复版统一的CZ根查找算法测试 ===\n");
    printf("比较 nmod_poly 和 fq_nmod_poly 版本\n");
    printf("====================================\n");
    
    test_fq_nmod_correctness();
    
    slong degrees[] = {100, 500, 1000};
    int num_degrees = sizeof(degrees) / sizeof(degrees[0]);
    int num_tests = 3;
    
    printf("\n=== nmod_poly 测试 (F_%llu) ===\n", PRIME);
    for (int i = 0; i < num_degrees; i++) {
        printf("\n--- 次数 %ld ---\n", degrees[i]);
        benchmark_nmod_roots(degrees[i], num_tests);
    }
    
    printf("\n=== fq_nmod_poly 测试 (F_{%llu^2}) ===\n", SMALL_PRIME);
    for (int i = 0; i < num_degrees; i++) {
        printf("\n--- 次数 %ld ---\n", degrees[i]);
        double our_time = benchmark_fq_nmod_roots(degrees[i], 2, num_tests);
        double flint_time = benchmark_flint_fq_nmod_factor(degrees[i], 2, num_tests);
        
        double ratio = (flint_time > 0) ? our_time / flint_time : 0;
        printf("比率 (我们的/FLINT): %.2fx\n", ratio);
    }
}
/*
int test_roots() {
    printf("=== 修复版统一的CZ根查找算法 ===\n");
    printf("支持 nmod_poly 和 fq_nmod_poly\n\n");
    
    run_unified_comparison();
    
    return 0;
}
*/
#endif