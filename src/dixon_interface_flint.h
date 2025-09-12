/*
 * 完整修复的Dixon结果式字符串接口
 */

#ifndef DIXON_INTERFACE_FLINT_H
#define DIXON_INTERFACE_FLINT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <flint/flint.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_poly.h>
#include <flint/fq_nmod_mat.h>
#include <flint/fmpz.h>
#include <flint/fq_nmod_poly_factor.h>  // 确保包含因式分解头文件
#include "unified_mpoly_resultant.h"
#include "fq_nmod_roots.h"

// 调试开关
#define DEBUG_PARSER 0


#if DEBUG_PARSER
#define DEBUG_PRINT(fmt, ...) printf("[PARSER] " fmt, ##__VA_ARGS__)
#else
#define DEBUG_PRINT(fmt, ...)
#endif

void fq_dixon_resultant(fq_mvpoly_t *result, fq_mvpoly_t *polys, 
                       slong nvars, slong npars);


// ============= Token类型定义 =============

typedef enum {
    TOK_NUMBER,      // 数字
    TOK_VARIABLE,    // 变量 (x, y, z等)
    TOK_PARAMETER,   // 参数 (a, b, c等)
    TOK_GENERATOR,   // 扩域生成元 (t)
    TOK_PLUS,        // +
    TOK_MINUS,       // -
    TOK_MULT,        // *
    TOK_POWER,       // ^
    TOK_LPAREN,      // (
    TOK_RPAREN,      // )
    TOK_EOF          // 结束
} token_type_t;

typedef struct {
    token_type_t type;
    char *str;
    fq_nmod_t value;
    slong int_value;
    const fq_nmod_ctx_struct *ctx;
} token_t;

typedef struct {
    const char *input;
    size_t pos;
    size_t len;
    token_t current;
    
    char **var_names;
    slong nvars;
    char **par_names;
    slong npars;
    slong max_pars;
    
    const fq_nmod_ctx_struct *ctx;
    char *generator_name;
} parser_state_t;

// ============= 修复的字符串解析器 =============

static int at_end(parser_state_t *state) {
    return state->pos >= state->len;
}

static char peek(parser_state_t *state) {
    if (at_end(state)) return '\0';
    return state->input[state->pos];
}

static char advance(parser_state_t *state) {
    if (at_end(state)) return '\0';
    return state->input[state->pos++];
}

static void skip_whitespace(parser_state_t *state) {
    while (!at_end(state) && isspace(peek(state))) {
        advance(state);
    }
}

static void parse_number(parser_state_t *state) {
    size_t start = state->pos;
    
    while (!at_end(state) && isdigit(peek(state))) {
        advance(state);
    }
    
    size_t len = state->pos - start;
    state->current.str = (char*) malloc(len + 1);
    strncpy(state->current.str, state->input + start, len);
    state->current.str[len] = '\0';
    
    state->current.int_value = atol(state->current.str);
    fq_nmod_set_ui(state->current.value, state->current.int_value, state->ctx);
    state->current.type = TOK_NUMBER;
    
    DEBUG_PRINT("Number: %s (%ld)\n", state->current.str, state->current.int_value);
}

static void parse_identifier(parser_state_t *state) {
    size_t start = state->pos;
    
    while (!at_end(state) && (isalnum(peek(state)) || peek(state) == '_')) {
        advance(state);
    }
    
    size_t len = state->pos - start;
    state->current.str = (char*) malloc(len + 1);
    strncpy(state->current.str, state->input + start, len);
    state->current.str[len] = '\0';
    
    if (state->generator_name && strcmp(state->current.str, state->generator_name) == 0) {
        state->current.type = TOK_GENERATOR;
        fq_nmod_gen(state->current.value, state->ctx);
        DEBUG_PRINT("Generator: %s\n", state->current.str);
    } else {
        state->current.type = TOK_VARIABLE;
        DEBUG_PRINT("Identifier: %s\n", state->current.str);
    }
}

static void next_token(parser_state_t *state) {
    skip_whitespace(state);
    
    if (state->current.str) {
        free(state->current.str);
        state->current.str = NULL;
    }
    
    if (at_end(state)) {
        state->current.type = TOK_EOF;
        return;
    }
    
    char c = peek(state);
    
    switch (c) {
        case '+': advance(state); state->current.type = TOK_PLUS; return;
        case '-': advance(state); state->current.type = TOK_MINUS; return;
        case '*': advance(state); state->current.type = TOK_MULT; return;
        case '^': advance(state); state->current.type = TOK_POWER; return;
        case '(': advance(state); state->current.type = TOK_LPAREN; return;
        case ')': advance(state); state->current.type = TOK_RPAREN; return;
    }
    
    if (isdigit(c)) {
        parse_number(state);
    } else if (isalpha(c) || c == '_') {
        parse_identifier(state);
    } else {
        advance(state);
        next_token(state);
    }
}

// 前向声明
static void parse_expression(parser_state_t *state, fq_mvpoly_t *poly);
static void parse_term(parser_state_t *state, fq_mvpoly_t *poly);
static void parse_factor(parser_state_t *state, fq_mvpoly_t *poly);
static void parse_primary(parser_state_t *state, fq_mvpoly_t *poly);

static slong find_or_add_parameter(parser_state_t *state, const char *name) {
    // 检查是否是变量
    for (slong i = 0; i < state->nvars; i++) {
        if (strcmp(state->var_names[i], name) == 0) {
            return -1;
        }
    }
    
    // 检查现有参数
    for (slong i = 0; i < state->npars; i++) {
        if (strcmp(state->par_names[i], name) == 0) {
            return i;
        }
    }
    
    // 添加新参数
    if (state->npars >= state->max_pars) {
        state->max_pars *= 2;
        state->par_names = (char**) realloc(state->par_names, 
                                           state->max_pars * sizeof(char*));
    }
    
    state->par_names[state->npars] = strdup(name);
    DEBUG_PRINT("New parameter: %s (index %ld)\n", name, state->npars);
    return state->npars++;
}

static slong get_variable_index(parser_state_t *state, const char *name) {
    for (slong i = 0; i < state->nvars; i++) {
        if (strcmp(state->var_names[i], name) == 0) {
            return i;
        }
    }
    return -1;
}

static void parse_primary(parser_state_t *state, fq_mvpoly_t *poly) {
    if (state->current.type == TOK_NUMBER) {
        fq_mvpoly_add_term(poly, NULL, NULL, state->current.value);
        next_token(state);
        
    } else if (state->current.type == TOK_GENERATOR) {
        fq_mvpoly_add_term(poly, NULL, NULL, state->current.value);
        next_token(state);
        
    } else if (state->current.type == TOK_VARIABLE) {
        char *name = strdup(state->current.str);
        next_token(state);
        
        slong var_idx = get_variable_index(state, name);
        if (var_idx >= 0) {
            slong *var_exp = (slong*) calloc(state->nvars, sizeof(slong));
            var_exp[var_idx] = 1;
            
            fq_nmod_t one;
            fq_nmod_init(one, state->ctx);
            fq_nmod_one(one, state->ctx);
            fq_mvpoly_add_term(poly, var_exp, NULL, one);
            fq_nmod_clear(one, state->ctx);
            free(var_exp);
        } else {
            slong par_idx = find_or_add_parameter(state, name);
            if (par_idx >= 0) {
                slong *par_exp = (slong*) calloc(state->max_pars, sizeof(slong));
                par_exp[par_idx] = 1;
                
                fq_nmod_t one;
                fq_nmod_init(one, state->ctx);
                fq_nmod_one(one, state->ctx);
                fq_mvpoly_add_term(poly, NULL, par_exp, one);
                fq_nmod_clear(one, state->ctx);
                free(par_exp);
            }
        }
        free(name);
        
    } else if (state->current.type == TOK_LPAREN) {
        next_token(state);
        parse_expression(state, poly);
        if (state->current.type == TOK_RPAREN) {
            next_token(state);
        }
        
    } else if (state->current.type == TOK_MINUS) {
        next_token(state);
        fq_mvpoly_t temp;
        fq_mvpoly_init(&temp, state->nvars, state->max_pars, state->ctx);
        parse_primary(state, &temp);
        
        for (slong i = 0; i < temp.nterms; i++) {
            fq_nmod_t neg_coeff;
            fq_nmod_init(neg_coeff, state->ctx);
            fq_nmod_neg(neg_coeff, temp.terms[i].coeff, state->ctx);
            fq_mvpoly_add_term(poly, temp.terms[i].var_exp, temp.terms[i].par_exp, neg_coeff);
            fq_nmod_clear(neg_coeff, state->ctx);
        }
        fq_mvpoly_clear(&temp);
    }
}

static void parse_factor(parser_state_t *state, fq_mvpoly_t *poly) {
    fq_mvpoly_t base;
    fq_mvpoly_init(&base, state->nvars, state->max_pars, state->ctx);
    parse_primary(state, &base);
    
    if (state->current.type == TOK_POWER) {
        next_token(state);
        if (state->current.type == TOK_NUMBER) {
            slong exp = state->current.int_value;
            next_token(state);
            
            fq_mvpoly_t result;
            fq_mvpoly_pow(&result, &base, exp);
            
            for (slong i = 0; i < result.nterms; i++) {
                fq_mvpoly_add_term(poly, result.terms[i].var_exp, result.terms[i].par_exp, result.terms[i].coeff);
            }
            fq_mvpoly_clear(&result);
        }
    } else {
        for (slong i = 0; i < base.nterms; i++) {
            fq_mvpoly_add_term(poly, base.terms[i].var_exp, base.terms[i].par_exp, base.terms[i].coeff);
        }
    }
    
    fq_mvpoly_clear(&base);
}

static void parse_term(parser_state_t *state, fq_mvpoly_t *poly) {
    fq_mvpoly_t result;
    fq_mvpoly_init(&result, state->nvars, state->max_pars, state->ctx);
    
    parse_factor(state, &result);
    
    while (state->current.type == TOK_MULT) {
        next_token(state);
        
        fq_mvpoly_t factor;
        fq_mvpoly_init(&factor, state->nvars, state->max_pars, state->ctx);
        parse_factor(state, &factor);
        
        fq_mvpoly_t temp;
        fq_mvpoly_mul(&temp, &result, &factor);
        
        fq_mvpoly_clear(&result);
        fq_mvpoly_clear(&factor);
        fq_mvpoly_copy(&result, &temp);
        fq_mvpoly_clear(&temp);
    }
    
    for (slong i = 0; i < result.nterms; i++) {
        fq_mvpoly_add_term_fast(poly, result.terms[i].var_exp, result.terms[i].par_exp, result.terms[i].coeff);
    }
    fq_mvpoly_clear(&result);
}

static void parse_expression(parser_state_t *state, fq_mvpoly_t *poly) {
    int negate = 0;
    if (state->current.type == TOK_MINUS) {
        negate = 1;
        next_token(state);
    } else if (state->current.type == TOK_PLUS) {
        next_token(state);
    }
    
    fq_mvpoly_t first_term;
    fq_mvpoly_init(&first_term, state->nvars, state->max_pars, state->ctx);
    parse_term(state, &first_term);
    
    if (negate) {
        for (slong i = 0; i < first_term.nterms; i++) {
            fq_nmod_t neg_coeff;
            fq_nmod_init(neg_coeff, state->ctx);
            fq_nmod_neg(neg_coeff, first_term.terms[i].coeff, state->ctx);
            fq_mvpoly_add_term_fast(poly, first_term.terms[i].var_exp, first_term.terms[i].par_exp, neg_coeff);
            fq_nmod_clear(neg_coeff, state->ctx);
        }
    } else {
        for (slong i = 0; i < first_term.nterms; i++) {
            fq_mvpoly_add_term_fast(poly, first_term.terms[i].var_exp, first_term.terms[i].par_exp, 
                           first_term.terms[i].coeff);
        }
    }
    fq_mvpoly_clear(&first_term);
    
    while (state->current.type == TOK_PLUS || state->current.type == TOK_MINUS) {
        int subtract = (state->current.type == TOK_MINUS);
        next_token(state);
        
        fq_mvpoly_t term;
        fq_mvpoly_init(&term, state->nvars, state->max_pars, state->ctx);
        parse_term(state, &term);
        
        for (slong i = 0; i < term.nterms; i++) {
            if (subtract) {
                fq_nmod_t neg_coeff;
                fq_nmod_init(neg_coeff, state->ctx);
                fq_nmod_neg(neg_coeff, term.terms[i].coeff, state->ctx);
                fq_mvpoly_add_term_fast(poly, term.terms[i].var_exp, term.terms[i].par_exp, neg_coeff);
                fq_nmod_clear(neg_coeff, state->ctx);
            } else {
                fq_mvpoly_add_term_fast(poly, term.terms[i].var_exp, term.terms[i].par_exp, term.terms[i].coeff);
            }
        }
        fq_mvpoly_clear(&term);
    }
}

// ============= 输出函数 =============

void fq_nmod_print_pretty_enhanced(const fq_nmod_t a, const fq_nmod_ctx_t ctx) {
    if (fq_nmod_is_zero(a, ctx)) {
        printf("0");
        return;
    }
    
    slong degree = fq_nmod_ctx_degree(ctx);
    
    if (degree == 1) {
        nmod_poly_t poly;
        nmod_poly_init(poly, fq_nmod_ctx_prime(ctx));
        fq_nmod_get_nmod_poly(poly, a, ctx);
        
        if (nmod_poly_degree(poly) >= 0) {
            printf("%lu", nmod_poly_get_coeff_ui(poly, 0));
        } else {
            printf("0");
        }
        nmod_poly_clear(poly);
    } else {
        nmod_poly_t poly;
        nmod_poly_init(poly, fq_nmod_ctx_prime(ctx));
        fq_nmod_get_nmod_poly(poly, a, ctx);
        
        slong deg = nmod_poly_degree(poly);
        int first_term = 1;
        
        for (slong i = deg; i >= 0; i--) {
            mp_limb_t coeff = nmod_poly_get_coeff_ui(poly, i);
            if (coeff != 0) {
                if (!first_term) {
                    printf(" + ");
                }
                first_term = 0;
                
                if (i == 0) {
                    printf("%lu", coeff);
                } else if (i == 1) {
                    if (coeff == 1) {
                        printf("t");
                    } else {
                        printf("%lu*t", coeff);
                    }
                } else {
                    if (coeff == 1) {
                        printf("t^%ld", i);
                    } else {
                        printf("%lu*t^%ld", coeff, i);
                    }
                }
            }
        }
        
        if (first_term) {
            printf("0");
        }
        
        nmod_poly_clear(poly);
    }
}

void fq_mvpoly_print_enhanced(const fq_mvpoly_t *p, const char *name) {
    printf("%s", name);
    if (strlen(name) > 0) printf(" = ");
    
    if (p->nterms == 0) {
        printf("0\n");
        return;
    }
    
    char var_names[] = {'x', 'y', 'z', 'w', 'v', 'u'};
    char par_names[] = {'a', 'b', 'c', 'd'};
    
    for (slong i = 0; i < p->nterms; i++) {
        if (i > 0) printf(" + ");
        
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
        
        fq_nmod_t one;
        fq_nmod_init(one, p->ctx);
        fq_nmod_one(one, p->ctx);
        
        int printed_something = 0;  // 跟踪是否已经打印了任何内容
        
        if (fq_nmod_is_one(p->terms[i].coeff, p->ctx)) {
            if (!has_vars) {
                printf("1");
                printed_something = 1;
            }
            // 系数为1且有变量时不打印系数，printed_something保持为0
        } else {
            printf("(");
            fq_nmod_print_pretty_enhanced(p->terms[i].coeff, p->ctx);
            printf(")");
            printed_something = 1;
        }
        
        fq_nmod_clear(one, p->ctx);
        
        // 打印变量
        for (slong j = 0; j < p->nvars; j++) {
            if (p->terms[i].var_exp && p->terms[i].var_exp[j] > 0) {
                if (printed_something) printf("*");  // 如果之前打印过任何内容，先打印乘号
                
                if (j < 6) printf("%c", var_names[j]);
                else printf("x_%ld", j);
                
                if (p->terms[i].var_exp[j] > 1) {
                    printf("^%ld", p->terms[i].var_exp[j]);
                }
                
                printed_something = 1;  // 标记已经打印了内容
            }
        }
        
        // 打印参数
        for (slong j = 0; j < p->npars; j++) {
            if (p->terms[i].par_exp && p->terms[i].par_exp[j] > 0) {
                if (printed_something) printf("*");  // 如果之前打印过任何内容，先打印乘号
                
                if (j < 4) printf("%c", par_names[j]);
                else printf("p_%ld", j);
                
                if (p->terms[i].par_exp[j] > 1) {
                    printf("^%ld", p->terms[i].par_exp[j]);
                }
                
                printed_something = 1;  // 标记已经打印了内容
            }
        }
    }
    printf("\n");
}

void find_and_print_roots_of_univariate_resultant(const fq_mvpoly_t *result, parser_state_t *state) {
    // 检查是否有消除变量
    if (result->nvars != 0) {
        return;  // 还有消除变量，不是最终结果
    }
    
    // 检查实际使用的参数个数
    int *par_used = (int*) calloc(result->npars, sizeof(int));
    slong actual_par_count = 0;
    slong main_par_idx = -1;
    
    // 遍历所有项，统计实际使用的参数
    for (slong t = 0; t < result->nterms; t++) {
        if (result->terms[t].par_exp) {
            for (slong p = 0; p < result->npars; p++) {
                if (result->terms[t].par_exp[p] > 0 && !par_used[p]) {
                    par_used[p] = 1;
                    main_par_idx = p;  // 记住最后一个使用的参数
                    actual_par_count++;
                }
            }
        }
    }
    
    printf("\n=== Polynomial Analysis ===\n");
    printf("Defined parameters: %ld\n", result->npars);
    printf("Actually used parameters: %ld\n", actual_par_count);
    
    if (actual_par_count > 1) {
        printf("Multiple parameters used: ");
        int first = 1;
        for (slong p = 0; p < result->npars; p++) {
            if (par_used[p]) {
                if (!first) printf(", ");
                if (state->par_names && state->par_names[p]) {
                    printf("%s", state->par_names[p]);
                } else {
                    printf("p_%ld", p);
                }
                first = 0;
            }
        }
        printf("\n");
        free(par_used);
        return;
    }
    
    if (actual_par_count == 0) {
        printf("Polynomial is constant.\n");
        free(par_used);
        return;
    }
    
    // actual_par_count == 1，进行求根
    printf("Single parameter detected! Finding roots...\n");
    
    // 将 fq_mvpoly_t 转换为 fq_nmod_poly_t
    fq_nmod_poly_t poly;
    fq_nmod_poly_init(poly, result->ctx);
    
    // 转换：将参数多项式转为单变量多项式
    for (slong i = 0; i < result->nterms; i++) {
        slong degree = 0;
        if (result->terms[i].par_exp && result->terms[i].par_exp[main_par_idx] > 0) {
            degree = result->terms[i].par_exp[main_par_idx];
        }
        fq_nmod_poly_set_coeff(poly, degree, result->terms[i].coeff, result->ctx);
    }
    
    const char *var_name = "unknown";
    if (state->par_names && state->par_names[main_par_idx]) {
        var_name = state->par_names[main_par_idx];
    }
    
    printf("\n=== Finding Roots of Univariate Resultant ===\n");
    printf("Univariate polynomial in %s:\n", var_name);
    printf("  Degree: %ld\n", fq_nmod_poly_degree(poly, result->ctx));
    
    slong degree = fq_nmod_poly_degree(poly, result->ctx);
    if (degree <= 0) {
        printf("  Polynomial is constant or zero, no roots to find.\n");
        fq_nmod_poly_clear(poly, result->ctx);
        free(par_used);
        return;
    }
    
// 使用FLINT的求根算法 - 优化版本，检测域类型
printf("\nFinding roots using appropriate algorithm...\n");

degree = fq_nmod_poly_degree(poly, result->ctx);
if (degree <= 0) {
    printf("  Polynomial is constant or zero, no roots to find.\n");
    fq_nmod_poly_clear(poly, result->ctx);
    free(par_used);
    return;
}

// 检测域类型
slong field_degree = fq_nmod_ctx_degree(result->ctx);
printf("Field degree: %ld\n", field_degree);
slong root_count = 0;
slong total_multiplicity = 0;

if (field_degree == 1) {
    // 素域情况 - 使用更高效的 nmod_poly_roots
    printf("Prime field detected, using nmod_poly_roots...\n");
    
    mp_limb_t prime = fq_nmod_ctx_prime(result->ctx);
    
    // 转换为 nmod_poly_t
    nmod_poly_t nmod_poly;
    nmod_poly_init(nmod_poly, prime);
    
    // 复制系数
    for (slong i = 0; i <= degree; i++) {
        fq_nmod_t coeff;
        fq_nmod_init(coeff, result->ctx);
        fq_nmod_poly_get_coeff(coeff, poly, i, result->ctx);
        
        // 从 fq_nmod_t 提取作为 mp_limb_t
        nmod_poly_t temp_poly;
        nmod_poly_init(temp_poly, prime);
        fq_nmod_get_nmod_poly(temp_poly, coeff, result->ctx);
        
        mp_limb_t coeff_ui = 0;
        if (nmod_poly_degree(temp_poly) >= 0) {
            coeff_ui = nmod_poly_get_coeff_ui(temp_poly, 0);
        }
        
        nmod_poly_set_coeff_ui(nmod_poly, i, coeff_ui);
        
        fq_nmod_clear(coeff, result->ctx);
        nmod_poly_clear(temp_poly);
    }
    
    // 使用 nmod_poly_roots 求根
    //nmod_poly_factor_t nmod_roots;
    //nmod_poly_factor_init(nmod_roots);
    nmod_roots_t nmod_roots;
    nmod_roots_init(nmod_roots);
    slong num_roots = our_nmod_poly_roots(nmod_roots, nmod_poly, 1);  // with_multiplicity = 1
    
    // 输出找到的根
    printf("\nRoots found:\n");
    
    printf("Find %ld roots:\n", num_roots);
    for (slong i = 0; i < nmod_roots->num; i++) {
        printf("  Root %ld: %lu (Multiplicity: %ld)\n", i + 1, 
               nmod_roots->roots[i], nmod_roots->mult[i]);
    }
    
    // 清理 nmod 相关结构
    nmod_poly_clear(nmod_poly);
    //nmod_poly_factor_clear(nmod_roots);
    
} else {
    // 扩域情况 - 使用原来的 fq_nmod_poly_roots
    printf("Extension field detected, using fq_nmod_poly_roots...\n");
    
    //fq_nmod_poly_factor_t roots;
    //fq_nmod_poly_factor_init(roots, result->ctx);
    fq_nmod_roots_t roots;
    fq_nmod_roots_init(roots, result->ctx);
    slong num_roots = our_fq_nmod_poly_roots(roots, poly, 1, result->ctx);
    
    printf("Find %ld roots:\n", num_roots);
    for (slong i = 0; i < roots->num; i++) {
        printf("  root %ld: ", i + 1);
        fq_nmod_print_pretty(roots->roots + i, result->ctx);
        printf(" (Multiplicity: %ld)\n", roots->mult[i]);
    }
    
    // 清理
    //fq_nmod_poly_factor_clear(roots, result->ctx);
}

    // 清理
    fq_nmod_poly_clear(poly, result->ctx);
    //fq_nmod_poly_factor_clear(roots, result->ctx);
    free(par_used);
}
// ============= 数据结构定义 =============

// Dixon结果结构
typedef struct {
    char *poly_string;
    char **remaining_vars;
    slong num_remaining_vars;
    const fq_nmod_ctx_struct *ctx;  // 改为指针类型
    int is_constant;
    fq_nmod_t constant_value;
    char *generator_name;
} dixon_result_t;

// ============= 辅助函数 =============

// 获取域生成元名称
char* get_generator_name(const fq_nmod_ctx_t ctx) {
    return strdup(ctx->var);
}

// Convert fq_nmod_t to string with field generator
char* fq_nmod_to_string_with_gen(const fq_nmod_t elem, const fq_nmod_ctx_t ctx, const char *gen_name) {
    // 初始分配较大的缓冲区
    size_t capacity = 256;
    char *buffer = (char*) malloc(capacity);
    buffer[0] = '\0';
    
    if (fq_nmod_is_zero(elem, ctx)) {
        strcpy(buffer, "0");
        return buffer;
    }
    
    if (fq_nmod_ctx_degree(ctx) == 1) {
        // Prime field element
        nmod_poly_t poly;
        nmod_poly_init(poly, fq_nmod_ctx_prime(ctx));
        fq_nmod_get_nmod_poly(poly, elem, ctx);
        
        if (nmod_poly_degree(poly) >= 0) {
            sprintf(buffer, "%lu", nmod_poly_get_coeff_ui(poly, 0));
        } else {
            strcpy(buffer, "0");
        }
        nmod_poly_clear(poly);
    } else {
        // Extension field element
        nmod_poly_t poly;
        nmod_poly_init(poly, fq_nmod_ctx_prime(ctx));
        fq_nmod_get_nmod_poly(poly, elem, ctx);
        
        slong deg = nmod_poly_degree(poly);
        int first = 1;
        size_t len = 0;
        
        // 动态构建字符串
        strcat(buffer, "(");
        len = 1;
        
        for (slong i = deg; i >= 0; i--) {
            mp_limb_t coeff = nmod_poly_get_coeff_ui(poly, i);
            if (coeff != 0) {
                char temp[128];
                
                if (!first) {
                    sprintf(temp, " + ");
                } else {
                    temp[0] = '\0';
                }
                first = 0;
                
                if (i == 0) {
                    sprintf(temp + strlen(temp), "%lu", coeff);
                } else if (i == 1) {
                    if (coeff == 1) {
                        sprintf(temp + strlen(temp), "%s", gen_name);
                    } else {
                        sprintf(temp + strlen(temp), "%lu*%s", coeff, gen_name);
                    }
                } else {
                    if (coeff == 1) {
                        sprintf(temp + strlen(temp), "%s^%ld", gen_name, i);
                    } else {
                        sprintf(temp + strlen(temp), "%lu*%s^%ld", coeff, gen_name, i);
                    }
                }
                
                // 检查是否需要扩展缓冲区
                size_t temp_len = strlen(temp);
                if (len + temp_len + 2 >= capacity) {
                    capacity = capacity * 2 + temp_len + 100;
                    char *new_buffer = realloc(buffer, capacity);
                    if (!new_buffer) {
                        free(buffer);
                        nmod_poly_clear(poly);
                        return NULL;
                    }
                    buffer = new_buffer;
                }
                
                strcat(buffer, temp);
                len += temp_len;
            }
        }
        
        strcat(buffer, ")");
        
        if (first) {
            strcpy(buffer, "(0)");
        }
        
        nmod_poly_clear(poly);
    }
    
    return buffer;
}

// Convert fq_mvpoly to string
char* fq_mvpoly_to_string_old(const fq_mvpoly_t *poly, char **var_names, const char *gen_name) {
    if (poly->nterms == 0) {
        return strdup("0");
    }
    
    size_t buf_size = poly->nterms * 500 + 4096;
    char *buffer = (char*) malloc(buf_size);
    buffer[0] = '\0';
    
    for (slong i = 0; i < poly->nterms; i++) {
        if (i > 0) strcat(buffer, " + ");
        
        // Add coefficient
        char *coeff_str = fq_nmod_to_string_with_gen(poly->terms[i].coeff, poly->ctx, gen_name);
        strcat(buffer, coeff_str);
        
        // Check if there are variables
        int has_vars = 0;
        if (poly->npars > 0 && poly->terms[i].par_exp) {
            for (slong j = 0; j < poly->npars; j++) {
                if (poly->terms[i].par_exp[j] > 0) {
                    has_vars = 1;
                    break;
                }
            }
        }
        
        // Add variables
        if (has_vars) {
            strcat(buffer, "*");
            int first_var = 1;
            for (slong j = 0; j < poly->npars; j++) {
                if (poly->terms[i].par_exp && poly->terms[i].par_exp[j] > 0) {
                    if (!first_var) strcat(buffer, "*");
                    first_var = 0;
                    
                    strcat(buffer, var_names[j]);
                    if (poly->terms[i].par_exp[j] > 1) {
                        char exp_str[32];
                        sprintf(exp_str, "^%ld", poly->terms[i].par_exp[j]);
                        strcat(buffer, exp_str);
                    }
                }
            }
        }
        
        free(coeff_str);
    }
    
    return buffer;
}

char* fq_mvpoly_to_string_old2(const fq_mvpoly_t *poly, char **var_names, const char *gen_name) {
    if (poly->nterms == 0) {
        return strdup("0");
    }
    
    // 动态分配初始缓冲区
    size_t capacity = poly->nterms * 200 + 4096;
    size_t length = 0;
    char *buffer = (char*) malloc(capacity);
    if (!buffer) return NULL;
    buffer[0] = '\0';
    
    for (slong i = 0; i < poly->nterms; i++) {
        char term_buffer[2048];  // 临时缓冲区for每一项
        term_buffer[0] = '\0';
        
        if (i > 0) strcat(term_buffer, " + ");
        
        // Add coefficient
        char *coeff_str = fq_nmod_to_string_with_gen(poly->terms[i].coeff, poly->ctx, gen_name);
        strcat(term_buffer, coeff_str);
        
        // Check if there are variables or parameters
        int has_vars = 0;
        
        // Check variables (var_exp)
        if (poly->nvars > 0 && poly->terms[i].var_exp) {
            for (slong j = 0; j < poly->nvars; j++) {
                if (poly->terms[i].var_exp[j] > 0) {
                    has_vars = 1;
                    break;
                }
            }
        }
        
        // Check parameters (par_exp) 
        if (!has_vars && poly->npars > 0 && poly->terms[i].par_exp) {
            for (slong j = 0; j < poly->npars; j++) {
                if (poly->terms[i].par_exp[j] > 0) {
                    has_vars = 1;
                    break;
                }
            }
        }
        
        // Add variables first (var_exp)
        if (poly->nvars > 0 && poly->terms[i].var_exp) {
            int printed_something = 0;
            for (slong j = 0; j < poly->nvars; j++) {
                if (poly->terms[i].var_exp[j] > 0) {
                    if (!printed_something && strcmp(coeff_str, "(1)") != 0) {
                        strcat(term_buffer, "*");
                    } else if (printed_something) {
                        strcat(term_buffer, "*");
                    }
                    
                    // Use provided var_names for variables
                    if (var_names && var_names[j]) {
                        strcat(term_buffer, var_names[j]);
                    } else {
                        char temp[32];
                        sprintf(temp, "x_%ld", j);
                        strcat(term_buffer, temp);
                    }
                    
                    if (poly->terms[i].var_exp[j] > 1) {
                        char exp_str[32];
                        sprintf(exp_str, "^%ld", poly->terms[i].var_exp[j]);
                        strcat(term_buffer, exp_str);
                    }
                    printed_something = 1;
                }
            }
            has_vars = printed_something;
        }
        
        // Then add parameters (par_exp)
        if (poly->npars > 0 && poly->terms[i].par_exp) {
            for (slong j = 0; j < poly->npars; j++) {
                if (poly->terms[i].par_exp[j] > 0) {
                    if (has_vars || strcmp(coeff_str, "(1)") != 0) {
                        strcat(term_buffer, "*");
                    }
                    
                    // Use var_names if provided, otherwise use p_j notation
                    if (var_names && var_names[j]) {
                        strcat(term_buffer, var_names[j]);
                    } else {
                        char temp[32];
                        sprintf(temp, "p_%ld", j);
                        strcat(term_buffer, temp);
                    }
                    
                    if (poly->terms[i].par_exp[j] > 1) {
                        char exp_str[32];
                        sprintf(exp_str, "^%ld", poly->terms[i].par_exp[j]);
                        strcat(term_buffer, exp_str);
                    }
                    has_vars = 1;
                }
            }
        }
        
        free(coeff_str);
        
        // 检查缓冲区大小，必要时扩展
        size_t term_len = strlen(term_buffer);
        if (length + term_len + 1 >= capacity) {
            capacity = capacity * 2 + term_len + 1000;
            char *new_buffer = realloc(buffer, capacity);
            if (!new_buffer) {
                free(buffer);
                return NULL;
            }
            buffer = new_buffer;
        }
        
        strcat(buffer, term_buffer);
        length += term_len;
    }
    
    return buffer;
}

// String builder 结构用于高效字符串构建
typedef struct {
    char *buffer;
    size_t length;      // 当前字符串长度
    size_t capacity;    // 总分配容量
} string_builder_t;

// String builder 辅助函数
static void sb_init(string_builder_t *sb, size_t initial_capacity) {
    sb->capacity = initial_capacity < 1024 ? 1024 : initial_capacity;
    sb->buffer = (char*) malloc(sb->capacity);
    if (!sb->buffer) {
        sb->capacity = 0;
        sb->length = 0;
        return;
    }
    sb->buffer[0] = '\0';
    sb->length = 0;
}

static int sb_ensure_capacity(string_builder_t *sb, size_t additional) {
    size_t required = sb->length + additional + 1;
    
    if (required > sb->capacity) {
        size_t new_capacity = sb->capacity * 2;
        while (new_capacity < required) {
            new_capacity *= 2;
        }
        
        char *new_buffer = (char*) realloc(sb->buffer, new_capacity);
        if (!new_buffer) {
            return 0;
        }
        
        sb->buffer = new_buffer;
        sb->capacity = new_capacity;
    }
    return 1;
}

static void sb_append(string_builder_t *sb, const char *str) {
    size_t len = strlen(str);
    if (sb_ensure_capacity(sb, len)) {
        memcpy(sb->buffer + sb->length, str, len + 1);
        sb->length += len;
    }
}

static void sb_append_char(string_builder_t *sb, char c) {
    if (sb_ensure_capacity(sb, 1)) {
        sb->buffer[sb->length++] = c;
        sb->buffer[sb->length] = '\0';
    }
}

static void sb_append_long(string_builder_t *sb, long value) {
    char temp[32];
    sprintf(temp, "%ld", value);
    sb_append(sb, temp);
}

static void sb_append_ulong(string_builder_t *sb, unsigned long value) {
    char temp[32];
    sprintf(temp, "%lu", value);
    sb_append(sb, temp);
}

static char* sb_finalize(string_builder_t *sb) {
    char *result = sb->buffer;
    sb->buffer = NULL;
    sb->capacity = 0;
    sb->length = 0;
    return result;
}

// ============= 优化的系数转字符串函数 =============
static void fq_nmod_to_string_builder(string_builder_t *sb, const fq_nmod_t elem, 
                                      const fq_nmod_ctx_t ctx, const char *gen_name) {
    if (fq_nmod_is_zero(elem, ctx)) {
        sb_append(sb, "0");
        return;
    }
    
    if (fq_nmod_ctx_degree(ctx) == 1) {
        // 素域元素
        nmod_poly_t poly;
        nmod_poly_init(poly, fq_nmod_ctx_prime(ctx));
        fq_nmod_get_nmod_poly(poly, elem, ctx);
        
        if (nmod_poly_degree(poly) >= 0) {
            sb_append_ulong(sb, nmod_poly_get_coeff_ui(poly, 0));
        } else {
            sb_append(sb, "0");
        }
        nmod_poly_clear(poly);
    } else {
        // 扩域元素
        nmod_poly_t poly;
        nmod_poly_init(poly, fq_nmod_ctx_prime(ctx));
        fq_nmod_get_nmod_poly(poly, elem, ctx);
        
        slong deg = nmod_poly_degree(poly);
        
        sb_append_char(sb, '(');
        
        int first = 1;
        for (slong i = deg; i >= 0; i--) {
            mp_limb_t coeff = nmod_poly_get_coeff_ui(poly, i);
            if (coeff != 0) {
                if (!first) {
                    sb_append(sb, " + ");
                }
                first = 0;
                
                if (i == 0) {
                    sb_append_ulong(sb, coeff);
                } else if (i == 1) {
                    if (coeff == 1) {
                        sb_append(sb, gen_name);
                    } else {
                        sb_append_ulong(sb, coeff);
                        sb_append_char(sb, '*');
                        sb_append(sb, gen_name);
                    }
                } else {
                    if (coeff == 1) {
                        sb_append(sb, gen_name);
                        sb_append_char(sb, '^');
                        sb_append_long(sb, i);
                    } else {
                        sb_append_ulong(sb, coeff);
                        sb_append_char(sb, '*');
                        sb_append(sb, gen_name);
                        sb_append_char(sb, '^');
                        sb_append_long(sb, i);
                    }
                }
            }
        }
        
        if (first) {
            sb_append(sb, "0");
        }
        
        sb_append_char(sb, ')');
        
        nmod_poly_clear(poly);
    }
}

// ============= 新的优化版本函数 - 直接替换原来的 fq_mvpoly_to_string =============
char* fq_mvpoly_to_string(const fq_mvpoly_t *poly, char **var_names, const char *gen_name) {
    if (poly->nterms == 0) {
        return strdup("0");
    }
    
    // 估算初始容量：每项大约需要 50 + 10*(变量数+参数数) 字符
    size_t estimated_size = poly->nterms * (50 + 10 * (poly->nvars + poly->npars));
    if (estimated_size < 1024) estimated_size = 1024;
    
    string_builder_t sb;
    sb_init(&sb, estimated_size);
    
    for (slong i = 0; i < poly->nterms; i++) {
        if (i > 0) {
            sb_append(&sb, " + ");
        }
        
        // 添加系数
        fq_nmod_to_string_builder(&sb, poly->terms[i].coeff, poly->ctx, gen_name);
        
        // 检查是否有变量或参数
        int has_vars = 0;
        
        // 检查变量 (var_exp)
        if (poly->nvars > 0 && poly->terms[i].var_exp) {
            for (slong j = 0; j < poly->nvars; j++) {
                if (poly->terms[i].var_exp[j] > 0) {
                    has_vars = 1;
                    break;
                }
            }
        }
        
        // 检查参数 (par_exp)
        if (!has_vars && poly->npars > 0 && poly->terms[i].par_exp) {
            for (slong j = 0; j < poly->npars; j++) {
                if (poly->terms[i].par_exp[j] > 0) {
                    has_vars = 1;
                    break;
                }
            }
        }
        
        // 添加变量 (var_exp)
        if (poly->nvars > 0 && poly->terms[i].var_exp) {
            for (slong j = 0; j < poly->nvars; j++) {
                if (poly->terms[i].var_exp[j] > 0) {
                    sb_append_char(&sb, '*');
                    
                    if (var_names && var_names[j]) {
                        sb_append(&sb, var_names[j]);
                    } else {
                        sb_append(&sb, "x_");
                        sb_append_long(&sb, j);
                    }
                    
                    if (poly->terms[i].var_exp[j] > 1) {
                        sb_append_char(&sb, '^');
                        sb_append_long(&sb, poly->terms[i].var_exp[j]);
                    }
                }
            }
        }
        
        // 添加参数 (par_exp)
        if (poly->npars > 0 && poly->terms[i].par_exp) {
            for (slong j = 0; j < poly->npars; j++) {
                if (poly->terms[i].par_exp[j] > 0) {
                    sb_append_char(&sb, '*');
                    
                    if (var_names && var_names[j]) {
                        sb_append(&sb, var_names[j]);
                    } else {
                        sb_append(&sb, "p_");
                        sb_append_long(&sb, j);
                    }
                    
                    if (poly->terms[i].par_exp[j] > 1) {
                        sb_append_char(&sb, '^');
                        sb_append_long(&sb, poly->terms[i].par_exp[j]);
                    }
                }
            }
        }
    }
    
    return sb_finalize(&sb);
}

// Print remaining variables info
void print_remaining_vars(char **var_names, slong nvars) {
    printf("Remaining variables (%ld): ", nvars);
    if (nvars == 0) {
        printf("none");
    } else {
        for (slong i = 0; i < nvars; i++) {
            if (i > 0) printf(", ");
            printf("%s", var_names[i]);
        }
    }
    printf("\n");
}

// concat polynomials helper function:
char* concat_polynomials(const char* poly1, const char* poly2, 
                              char** buffer, size_t* buffer_size) {
    size_t len1 = strlen(poly1);
    size_t len2 = strlen(poly2);
    size_t needed = len1 + len2 + 2;
    
    if (needed > *buffer_size) {
        *buffer_size = needed;
        *buffer = (char*) realloc(*buffer, *buffer_size);
        if (!*buffer) {
            printf("ERROR: Failed to allocate %zu bytes\n", *buffer_size);
            exit(1);
        }
    }
    
    sprintf(*buffer, "%s,%s", poly1, poly2);
    return *buffer;
}

// ============= Core Dixon Function =============

// Internal computation function
char* compute_dixon_internal(const char **poly_strings, slong npoly_strings,
                           const char **var_names, slong nvars,
                           const fq_nmod_ctx_t ctx,
                           char ***remaining_vars, slong *num_remaining) {
    
    if (npoly_strings != nvars + 1) {
        fprintf(stderr, "Error: Need exactly %ld polynomials for %ld variables\n",
                nvars + 1, nvars);
        *remaining_vars = NULL;
        *num_remaining = 0;
        return strdup("0");
    }
    
    // Get generator name
    char *gen_name = get_generator_name(ctx);
    
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
    state.ctx = ctx;
    state.current.str = NULL;
    fq_nmod_init(state.current.value, ctx);
    state.generator_name = strdup(gen_name);
    
    // First pass: identify parameters
    for (slong i = 0; i < npoly_strings; i++) {
        fq_mvpoly_t temp;
        fq_mvpoly_init(&temp, nvars, state.max_pars, ctx);
        
        state.input = poly_strings[i];
        state.pos = 0;
        state.len = strlen(poly_strings[i]);
        next_token(&state);
        
        parse_expression(&state, &temp);
        fq_mvpoly_clear(&temp);
    }
    
    // Save parameters as remaining variables
    *num_remaining = state.npars;
    if (state.npars > 0) {
        *remaining_vars = (char**) malloc(state.npars * sizeof(char*));
        for (slong i = 0; i < state.npars; i++) {
            (*remaining_vars)[i] = strdup(state.par_names[i]);
        }
    } else {
        *remaining_vars = NULL;
    }
    
    // Second pass: parse polynomials
    fq_mvpoly_t *polys = (fq_mvpoly_t*) malloc(npoly_strings * sizeof(fq_mvpoly_t));
    
    for (slong i = 0; i < npoly_strings; i++) {
        fq_mvpoly_init(&polys[i], nvars, state.npars, ctx);
        
        state.input = poly_strings[i];
        state.pos = 0;
        state.len = strlen(poly_strings[i]);
        if (state.current.str) {
            free(state.current.str);
            state.current.str = NULL;
        }
        next_token(&state);
        
        parse_expression(&state, &polys[i]);
    }
    
    // Compute Dixon resultant
    fq_mvpoly_t dixon_result_poly;
    fq_dixon_resultant(&dixon_result_poly, polys, nvars, state.npars);

    find_and_print_roots_of_univariate_resultant(&dixon_result_poly, &state);
    
    // Convert result to string
    char *result_string;
    if (dixon_result_poly.nterms == 0) {
        result_string = strdup("0");
    } else {
        result_string = fq_mvpoly_to_string(&dixon_result_poly, state.par_names, gen_name);
    }
    
    // Cleanup
    fq_mvpoly_clear(&dixon_result_poly);
    for (slong i = 0; i < npoly_strings; i++) {
        fq_mvpoly_clear(&polys[i]);
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
    
    if (state.generator_name) {
        free(state.generator_name);
    }
    
    fq_nmod_clear(state.current.value, ctx);
    if (state.current.str) {
        free(state.current.str);
    }
    
    free(gen_name);
    
    return result_string;
}

// 辅助函数：去除字符串两端的空格
static char* trim_whitespace(char *str) {
    if (!str) return NULL;
    
    // 去除开头空格
    while (isspace(*str)) str++;
    
    if (*str == '\0') return str;
    
    // 去除结尾空格
    char *end = str + strlen(str) - 1;
    while (end > str && isspace(*end)) end--;
    *(end + 1) = '\0';
    
    return str;
}

// 辅助函数：分割字符串（按逗号）
static char** split_string_old(const char *input, slong *count) {
    if (!input || strlen(input) == 0) {
        *count = 0;
        return NULL;
    }
    
    // 复制输入字符串以便修改
    char *work_str = strdup(input);
    
    // 第一遍：计算有多少个元素
    slong num_elements = 1;
    for (const char *p = input; *p; p++) {
        if (*p == ',') num_elements++;
    }
    
    // 分配结果数组
    char **result = (char**) malloc(num_elements * sizeof(char*));
    
    // 第二遍：分割字符串
    slong idx = 0;
    char *token = strtok(work_str, ",");
    while (token != NULL) {
        token = trim_whitespace(token);
        result[idx++] = strdup(token);
        token = strtok(NULL, ",");
    }
    
    *count = idx;
    free(work_str);
    return result;
}

static char** split_string(const char *input, slong *count) {
    if (!input || strlen(input) == 0) {
        *count = 0;
        return NULL;
    }
    
    size_t input_len = strlen(input);
    //printf("  [DEBUG] Input string length: %zu\n", input_len);
    
    // 复制输入字符串以便修改 - 对于大字符串使用动态分配
    char *work_str = (char*) malloc(input_len + 1);
    if (!work_str) {
        //printf("  [ERROR] Failed to allocate memory for work_str\n");
        *count = 0;
        return NULL;
    }
    memcpy(work_str, input, input_len + 1);
    
    // 第一遍：计算有多少个元素（查找逗号）
    slong num_elements = 1;
    int in_parentheses = 0;
    for (size_t i = 0; i < input_len; i++) {
        if (input[i] == '(') in_parentheses++;
        else if (input[i] == ')') in_parentheses--;
        else if (input[i] == ',' && in_parentheses == 0) {
            num_elements++;
        }
    }
    
    //printf("  [DEBUG] Found %ld polynomials\n", num_elements);
    
    // 分配结果数组
    char **result = (char**) malloc(num_elements * sizeof(char*));
    if (!result) {
        //printf("  [ERROR] Failed to allocate result array\n");
        free(work_str);
        *count = 0;
        return NULL;
    }
    
    // 第二遍：分割字符串（手动处理，不用 strtok）
    slong idx = 0;
    size_t start = 0;
    in_parentheses = 0;
    
    for (size_t i = 0; i <= input_len; i++) {
        if (i < input_len) {
            if (input[i] == '(') in_parentheses++;
            else if (input[i] == ')') in_parentheses--;
        }
        
        if ((i == input_len || (input[i] == ',' && in_parentheses == 0)) && i > start) {
            size_t poly_len = i - start;
            char *poly = (char*) malloc(poly_len + 1);
            memcpy(poly, input + start, poly_len);
            poly[poly_len] = '\0';
            
            // Trim whitespace
            char *trimmed = trim_whitespace(poly);
            result[idx++] = strdup(trimmed);
            if (poly != trimmed) free(poly);
            
            start = i + 1;
        }
    }
    
    *count = idx;
    free(work_str);
    
    //printf("  [DEBUG] Successfully split into %ld polynomials\n", *count);
    return result;
}

// 释放分割后的字符串数组
static void free_split_strings(char **strings, slong count) {
    if (!strings) return;
    for (slong i = 0; i < count; i++) {
        if (strings[i]) free(strings[i]);
    }
    free(strings);
}

// ============= Main Dixon Interface =============

// Unified Dixon function - accepts array of polynomial strings
// Returns result as string, optionally outputs remaining variables
char* dixon(const char **poly_strings, slong num_polys, 
            const char **elim_vars, slong num_elim_vars,
            const fq_nmod_ctx_t ctx) {
    
    printf("\n=== Dixon Computation ===\n");
    clock_t start = clock();
    printf("Eliminating variables: ");
    for (slong i = 0; i < num_elim_vars; i++) {
        if (i > 0) printf(", ");
        printf("%s", elim_vars[i]);
    }
    printf("\n");
    
    //printf("Input polynomials:\n");
    for (slong i = 0; i < num_polys; i++) {
        //printf("  p%ld: %s\n", i, poly_strings[i]);
    }
    
    // Compute Dixon resultant
    char **remaining_vars = NULL;
    slong num_remaining = 0;
    
    char *result = compute_dixon_internal(poly_strings, num_polys, 
                                         elim_vars, num_elim_vars, ctx,
                                         &remaining_vars, &num_remaining);
    
    // printf("Result: %s\n", result); print_remaining_vars(remaining_vars, num_remaining);
    
    // Cleanup remaining vars
    if (remaining_vars) {
        for (slong i = 0; i < num_remaining; i++) {
            free(remaining_vars[i]);
        }
        free(remaining_vars);
    }

    clock_t end = clock();
    printf("Time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("resultant length: %zu characters\n\n", strlen(result));
    
    return result;
}

// Implementation using unified_mpoly_resultant
char* bivariate_resultant(const char *poly1_str, const char *poly2_str,
                                         const char *elim_var, const fq_nmod_ctx_t ctx) {
    clock_t total_start = clock();
    // 获取生成元名称
    char *gen_name = get_generator_name(ctx);
        
    char **remaining_vars = NULL;
    slong num_remaining = 0;
    
    // 初始化解析器状态
    parser_state_t state;
    state.var_names = (char**) malloc(1 * sizeof(char*));
    state.var_names[0] = strdup(elim_var);
    state.nvars = 1;
    state.npars = 0;
    state.max_pars = 16;
    state.par_names = (char**) malloc(state.max_pars * sizeof(char*));
    state.ctx = ctx;
    state.current.str = NULL;
    fq_nmod_init(state.current.value, ctx);
    state.generator_name = strdup(gen_name);
    
    // 第一遍：解析以识别参数
    fq_mvpoly_t temp1, temp2;
    fq_mvpoly_init(&temp1, 1, state.max_pars, ctx);
    fq_mvpoly_init(&temp2, 1, state.max_pars, ctx);

    //printf("saving pars1\n");
    state.input = poly1_str;
    state.pos = 0;
    state.len = strlen(poly1_str);
    next_token(&state);
    parse_expression(&state, &temp1);

    //printf("saving pars2\n");
    state.input = poly2_str;
    state.pos = 0;
    state.len = strlen(poly2_str);
    if (state.current.str) {
        free(state.current.str);
        state.current.str = NULL;
    }
    next_token(&state);
    parse_expression(&state, &temp2);
    
    fq_mvpoly_clear(&temp1);
    fq_mvpoly_clear(&temp2);
    //printf("saving...");
    // 保存参数
    num_remaining = state.npars;
    if (state.npars > 0) {
        remaining_vars = (char**) malloc(state.npars * sizeof(char*));
        for (slong i = 0; i < state.npars; i++) {
            (remaining_vars)[i] = strdup(state.par_names[i]);
        }
    } else {
        remaining_vars = NULL;
    }

    //printf("parsing...\n");
    // 第二遍：正式解析
    fq_mvpoly_t poly1, poly2;
    fq_mvpoly_init(&poly1, 1, state.npars, ctx);
    fq_mvpoly_init(&poly2, 1, state.npars, ctx);
    
    state.input = poly1_str;
    state.pos = 0;
    state.len = strlen(poly1_str);
    if (state.current.str) {
        free(state.current.str);
        state.current.str = NULL;
    }
    next_token(&state);
    parse_expression(&state, &poly1);
    
    state.input = poly2_str;
    state.pos = 0;
    state.len = strlen(poly2_str);
    if (state.current.str) {
        free(state.current.str);
        state.current.str = NULL;
    }
    next_token(&state);
    parse_expression(&state, &poly2);
    
    // 初始化统一字段上下文
    field_ctx_t field_ctx;
    field_ctx_init(&field_ctx, ctx);
    
    // 创建统一多项式上下文
    slong total_vars = 1 + state.npars;  // 消去变量 + 参数
    unified_mpoly_ctx_t unified_ctx = unified_mpoly_ctx_init(total_vars, ORD_LEX, &field_ctx);
    
    // 初始化统一多项式
    unified_mpoly_t A = unified_mpoly_init(unified_ctx);
    unified_mpoly_t B = unified_mpoly_init(unified_ctx);
    unified_mpoly_t R = unified_mpoly_init(unified_ctx);

    //printf("conversion first poly...\n");
    // 转换第一个多项式
    for (slong i = 0; i < poly1.nterms; i++) {
        field_elem_u coeff;
        ulong *exp = (ulong*) calloc(total_vars, sizeof(ulong));
        
        // 转换系数
        fq_nmod_to_field_elem(&coeff, poly1.terms[i].coeff, &field_ctx);
        
        // 设置指数
        if (poly1.terms[i].var_exp) {
            exp[0] = poly1.terms[i].var_exp[0];
        }
        if (poly1.terms[i].par_exp) {
            for (slong j = 0; j < state.npars; j++) {
                exp[1 + j] = poly1.terms[i].par_exp[j];
            }
        }
        
        unified_mpoly_set_coeff_ui(A, &coeff, exp);
        free(exp);
    }

    //printf("conversion second poly...\n");
    // 转换第二个多项式
    for (slong i = 0; i < poly2.nterms; i++) {
        field_elem_u coeff;
        ulong *exp = (ulong*) calloc(total_vars, sizeof(ulong));
        
        // 转换系数
        fq_nmod_to_field_elem(&coeff, poly2.terms[i].coeff, &field_ctx);
        
        // 设置指数
        if (poly2.terms[i].var_exp) {
            exp[0] = poly2.terms[i].var_exp[0];
        }
        if (poly2.terms[i].par_exp) {
            for (slong j = 0; j < state.npars; j++) {
                exp[1 + j] = poly2.terms[i].par_exp[j];
            }
        }
        
        unified_mpoly_set_coeff_ui(B, &coeff, exp);
        free(exp);
    }
    
    printf("Computing resultant using unified interface w.r.t. %s\n", elim_var);
    //printf("  Length of poly1: %ld terms\n", unified_mpoly_length(A));
    //printf("  Length of poly2: %ld terms\n", unified_mpoly_length(B));
    
    // 计算结式
    clock_t start = clock();
    int success = unified_mpoly_resultant(R, A, B, 0, unified_ctx);
    clock_t end = clock();
    
    if (!success) {
        printf("Unified resultant computation failed!\n");
        unified_mpoly_clear(A);
        unified_mpoly_clear(B);
        unified_mpoly_clear(R);
        unified_mpoly_ctx_clear(unified_ctx);
        
        // 清理其他资源
        for (slong i = 0; i < state.nvars; i++) {
            free(state.var_names[i]);
        }
        free(state.var_names);
        for (slong i = 0; i < state.npars; i++) {
            free(state.par_names[i]);
        }
        free(state.par_names);
        if (state.generator_name) free(state.generator_name);
        fq_nmod_clear(state.current.value, ctx);
        if (state.current.str) free(state.current.str);
        free(gen_name);
        fq_mvpoly_clear(&poly1);
        fq_mvpoly_clear(&poly2);
        
        return strdup("0");
    }
    
    printf("Unified resultant computation time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("Resultant has %ld terms\n", unified_mpoly_length(R));
    
    // 将结果转换回fq_mvpoly格式
    fq_mvpoly_t result_mvpoly;
    fq_mvpoly_init(&result_mvpoly, 0, state.npars, ctx);
    
    // 从unified格式转换回来
    slong result_len = unified_mpoly_length(R);
    if (result_len > 0) {
        // 需要遍历结果的所有项
        // 这需要访问unified_mpoly的内部结构
        // 由于unified_mpoly是对nmod_mpoly或fq_nmod_mpoly的封装
        // 我们需要根据field_id来处理
        
        if (field_ctx.field_id == FIELD_ID_NMOD) {
            // 处理nmod情况
            nmod_mpoly_struct *nmod_res = GET_NMOD_POLY(R);
            nmod_mpoly_ctx_struct *nmod_ctx = &(unified_ctx->ctx.nmod_ctx);
            
            for (slong i = 0; i < nmod_mpoly_length(nmod_res, nmod_ctx); i++) {
                ulong coeff_ui = nmod_mpoly_get_term_coeff_ui(nmod_res, i, nmod_ctx);
                ulong *exp = (ulong*) calloc(total_vars, sizeof(ulong));
                nmod_mpoly_get_term_exp_ui(exp, nmod_res, i, nmod_ctx);
                
                // 转换系数
                fq_nmod_t coeff;
                fq_nmod_init(coeff, ctx);
                fq_nmod_set_ui(coeff, coeff_ui, ctx);
                
                // 提取参数指数
                slong *par_exp = NULL;
                if (state.npars > 0) {
                    par_exp = (slong*) calloc(state.npars, sizeof(slong));
                    for (slong j = 0; j < state.npars; j++) {
                        par_exp[j] = exp[1 + j];
                    }
                }
                
                fq_mvpoly_add_term_fast(&result_mvpoly, NULL, par_exp, coeff);
                
                fq_nmod_clear(coeff, ctx);
                free(exp);
                if (par_exp) free(par_exp);
            }
        } else if (field_ctx.field_id == FIELD_ID_FQ_ZECH) {
        // 处理fq_zech情况
        fq_zech_mpoly_struct *zech_res = GET_ZECH_POLY(R);
        fq_zech_mpoly_ctx_struct *zech_ctx = &(unified_ctx->ctx.zech_ctx);
        
        for (slong i = 0; i < fq_zech_mpoly_length(zech_res, zech_ctx); i++) {
            fq_zech_t zech_coeff;
            fq_zech_init(zech_coeff, zech_ctx->fqctx);
            fq_zech_mpoly_get_term_coeff_fq_zech(zech_coeff, zech_res, i, zech_ctx);
            
            // 转换为fq_nmod
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ctx);
            fq_zech_get_fq_nmod(coeff, zech_coeff, zech_ctx->fqctx);
            
            ulong *exp = (ulong*) calloc(total_vars, sizeof(ulong));
            fq_zech_mpoly_get_term_exp_ui(exp, zech_res, i, zech_ctx);
            
            // 提取参数指数
            slong *par_exp = NULL;
            if (state.npars > 0) {
                par_exp = (slong*) calloc(state.npars, sizeof(slong));
                for (slong j = 0; j < state.npars; j++) {
                    par_exp[j] = exp[1 + j];
                }
            }
            
            fq_mvpoly_add_term_fast(&result_mvpoly, NULL, par_exp, coeff);
            
            fq_nmod_clear(coeff, ctx);
            fq_zech_clear(zech_coeff, zech_ctx->fqctx);
            free(exp);
            if (par_exp) free(par_exp);
        }
    } else {
            // 处理fq_nmod情况
            fq_nmod_mpoly_struct *fq_res = GET_FQ_POLY(R);
            fq_nmod_mpoly_ctx_struct *fq_ctx = &(unified_ctx->ctx.fq_ctx);
            printf("fq_mvpoly_add_term_fast\n");
            for (slong i = 0; i < fq_nmod_mpoly_length(fq_res, fq_ctx); i++) {
                fq_nmod_t coeff;
                fq_nmod_init(coeff, ctx);
                fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, fq_res, i, fq_ctx);
                
                ulong *exp = (ulong*) calloc(total_vars, sizeof(ulong));
                fq_nmod_mpoly_get_term_exp_ui(exp, fq_res, i, fq_ctx);
                
                // 提取参数指数
                slong *par_exp = NULL;
                if (state.npars > 0) {
                    par_exp = (slong*) calloc(state.npars, sizeof(slong));
                    for (slong j = 0; j < state.npars; j++) {
                        par_exp[j] = exp[1 + j];
                    }
                }
                
                fq_mvpoly_add_term_fast(&result_mvpoly, NULL, par_exp, coeff);
                
                fq_nmod_clear(coeff, ctx);
                free(exp);
                if (par_exp) free(par_exp);
            }
        }
    }
    fq_mvpoly_make_monic(&result_mvpoly);
    // 如果是单变量多项式，尝试找根
    find_and_print_roots_of_univariate_resultant(&result_mvpoly, &state);
    printf("fq_mvpoly_to_string\n");
    // 转换结果为字符串
    char *result_string;
    if (result_mvpoly.nterms == 0) {
        result_string = strdup("0");
    } else {
        result_string = fq_mvpoly_to_string(&result_mvpoly, state.par_names, gen_name);
    }
    printf("clean up\n");

    // 输出剩余变量信息
    if (num_remaining > 0) {
        printf("Remaining variables (%ld): ", num_remaining);
        for (slong i = 0; i < num_remaining; i++) {
            if (i > 0) printf(", ");
            printf("%s", remaining_vars[i]);
            free(remaining_vars[i]);
        }
        printf("\n");
        free(remaining_vars);
    }
    // 清理
    unified_mpoly_clear(A);
    unified_mpoly_clear(B);
    unified_mpoly_clear(R);
    unified_mpoly_ctx_clear(unified_ctx);
    fq_mvpoly_clear(&poly1);
    fq_mvpoly_clear(&poly2);
    fq_mvpoly_clear(&result_mvpoly);
    
    for (slong i = 0; i < state.nvars; i++) {
        free(state.var_names[i]);
    }
    free(state.var_names);
    for (slong i = 0; i < state.npars; i++) {
        free(state.par_names[i]);
    }
    free(state.par_names);
    if (state.generator_name) free(state.generator_name);
    fq_nmod_clear(state.current.value, ctx);
    if (state.current.str) free(state.current.str);
    free(gen_name);

    clock_t total_end = clock();
    printf("Time: %.3f seconds\n", (double)(total_end - total_start) / CLOCKS_PER_SEC);
    printf("resultant length: %zu characters\n\n", strlen(result_string));
    
    return result_string;
}

char* dixon_str(const char *poly_string,    // 逗号分隔的多项式
                const char *vars_string,     // 逗号分隔的变量
                const fq_nmod_ctx_t ctx) {
    
    printf("\n=== Dixon/Resultant Computation (String Interface) ===\n");
    
    // 分割输入字符串
    slong num_polys, num_vars;
    char **poly_array = split_string(poly_string, &num_polys);
    char **vars_array = split_string(vars_string, &num_vars);
    
    char *result = NULL;
    
    // 检查是否为双变量情况
    if (num_polys == 2 && num_vars == 1) {
        printf("Using unified bivariate resultant for 2 polynomials...\n");

        
        // 调用统一接口的双变量结式计算
        result = bivariate_resultant(poly_array[0], poly_array[1], 
                                                    vars_array[0], ctx);        
    } else {
        // 使用原始Dixon方法
        printf("Using Dixon resultant for %ld polynomials...\n", num_polys);
        
        // 转换为const char**
        const char **poly_strings = (const char**) malloc(num_polys * sizeof(char*));
        const char **elim_vars = (const char**) malloc(num_vars * sizeof(char*));
        
        for (slong i = 0; i < num_polys; i++) {
            poly_strings[i] = poly_array[i];
        }
        for (slong i = 0; i < num_vars; i++) {
            elim_vars[i] = vars_array[i];
        }
        
        // 调用原始dixon函数
        result = dixon(poly_strings, num_polys, elim_vars, num_vars, ctx);
        
        // 清理
        free(poly_strings);
        free(elim_vars);
    }
    
    // 清理分割的字符串
    free_split_strings(poly_array, num_polys);
    free_split_strings(vars_array, num_vars);
    
    return result;
}


#endif // COMPLETE_FIXED_DIXON_INTERFACE_H