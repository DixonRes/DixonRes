/*
 * 完整修复的Dixon结果式字符串接口和插值算法
 * 
 * 修复内容：
 * 1. 字符串解析器的Token识别和运算符优先级
 * 2. 插值算法的度数界计算和多元插值
 * 3. 有限域上的插值点生成
 * 4. Dixon算法中dual variables的正确处理
 */

#ifndef COMPLETE_FIXED_DIXON_INTERFACE_H
#define COMPLETE_FIXED_DIXON_INTERFACE_H

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
// 调试开关
#define DEBUG_PARSER 0


#if DEBUG_PARSER
#define DEBUG_PRINT(fmt, ...) printf("[PARSER] " fmt, ##__VA_ARGS__)
#else
#define DEBUG_PRINT(fmt, ...)
#endif



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
        fq_mvpoly_add_term(poly, result.terms[i].var_exp, result.terms[i].par_exp, result.terms[i].coeff);
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
            fq_mvpoly_add_term(poly, first_term.terms[i].var_exp, first_term.terms[i].par_exp, neg_coeff);
            fq_nmod_clear(neg_coeff, state->ctx);
        }
    } else {
        for (slong i = 0; i < first_term.nterms; i++) {
            fq_mvpoly_add_term(poly, first_term.terms[i].var_exp, first_term.terms[i].par_exp, 
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
                fq_mvpoly_add_term(poly, term.terms[i].var_exp, term.terms[i].par_exp, neg_coeff);
                fq_nmod_clear(neg_coeff, state->ctx);
            } else {
                fq_mvpoly_add_term(poly, term.terms[i].var_exp, term.terms[i].par_exp, term.terms[i].coeff);
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
    // 检查是否为单变量多项式（只有一个参数）
    if (result->npars != 1 || result->nvars != 0) {
        return;  // 不是单变量多项式，直接返回
    }
    
    printf("\n=== Finding Roots of Univariate Resultant ===\n");
    
    // 将 fq_mvpoly_t 转换为 fq_nmod_poly_t
    fq_nmod_poly_t poly;
    fq_nmod_poly_init(poly, result->ctx);
    
    // 转换：将参数多项式转为单变量多项式
    for (slong i = 0; i < result->nterms; i++) {
        slong degree = 0;
        if (result->terms[i].par_exp) {
            degree = result->terms[i].par_exp[0];
        }
        fq_nmod_poly_set_coeff(poly, degree, result->terms[i].coeff, result->ctx);
    }
    
    printf("Univariate polynomial in %s:\n", state->par_names[0]);
    printf("  Degree: %ld\n", fq_nmod_poly_degree(poly, result->ctx));
    printf("  ");
    fq_nmod_poly_print_pretty(poly, state->par_names[0], result->ctx);
    printf("\n");
    
    slong degree = fq_nmod_poly_degree(poly, result->ctx);
    if (degree <= 0) {
        printf("  Polynomial is constant or zero, no roots to find.\n");
        fq_nmod_poly_clear(poly, result->ctx);
        return;
    }
    
    // 使用FLINT的专门求根函数
    printf("\nFinding roots using Rabin's algorithm...\n");
    
    fq_nmod_poly_factor_t roots;
    fq_nmod_poly_factor_init(roots, result->ctx);
    
    // 使用 fq_nmod_poly_roots 直接找根
    // 第三个参数为1表示计算重数
    fq_nmod_poly_roots(roots, poly, 1, result->ctx);
    
    // 输出找到的根
    printf("\nRoots found:\n");
    slong root_count = 0;
    slong total_multiplicity = 0;
    
    for (slong i = 0; i < roots->num; i++) {
        // 每个根都以线性因子 (x - root) 的形式存储
        if (fq_nmod_poly_degree(roots->poly + i, result->ctx) == 1) {
            fq_nmod_t a, b, root;
            fq_nmod_init(a, result->ctx);
            fq_nmod_init(b, result->ctx);
            fq_nmod_init(root, result->ctx);
            
            fq_nmod_poly_get_coeff(a, roots->poly + i, 1, result->ctx);
            fq_nmod_poly_get_coeff(b, roots->poly + i, 0, result->ctx);
            
            if (!fq_nmod_is_zero(a, result->ctx)) {
                // 根 = -b/a
                fq_nmod_neg(root, b, result->ctx);
                fq_nmod_div(root, root, a, result->ctx);
                
                printf("  %s = ", state->par_names[0]);
                fq_nmod_print_pretty_enhanced(root, result->ctx);
                
                if (roots->exp[i] > 1) {
                    printf(" (multiplicity %ld)", roots->exp[i]);
                }
                printf("\n");
                
                root_count++;
                total_multiplicity += roots->exp[i];
            }
            
            fq_nmod_clear(a, result->ctx);
            fq_nmod_clear(b, result->ctx);
            fq_nmod_clear(root, result->ctx);
        } else if (fq_nmod_poly_degree(roots->poly + i, result->ctx) == 0) {
            // 常数因子，跳过
            continue;
        } else {
            // 这不应该发生，因为 roots 函数应该只返回线性因子
            printf("  Warning: non-linear factor of degree %ld found\n", 
                   fq_nmod_poly_degree(roots->poly + i, result->ctx));
        }
    }
    
    if (root_count == 0) {
        printf("  No roots found in the field F_%ld^%ld\n", 
               fq_nmod_ctx_prime(result->ctx), 
               fq_nmod_ctx_degree(result->ctx));
        printf("  The polynomial might be irreducible over this field.\n");
    } else {
        printf("\nTotal roots found: %ld", root_count);
        if (total_multiplicity != root_count) {
            printf(" (total multiplicity: %ld)", total_multiplicity);
        }
        printf("\n");
        
        if (total_multiplicity < degree) {
            printf("Note: Only %ld out of %ld degree accounted for.\n", 
                   total_multiplicity, degree);
            printf("The polynomial has irreducible factors of degree > 1.\n");
        }
    }
    
    // 清理
    fq_nmod_poly_clear(poly, result->ctx);
    fq_nmod_poly_factor_clear(roots, result->ctx);
}

void find_and_print_roots_of_univariate_resultant_factor(const fq_mvpoly_t *result, parser_state_t *state) {
    // 检查是否为单变量多项式（只有一个参数）
    if (result->npars != 1 || result->nvars != 0) {
        return;  // 不是单变量多项式，直接返回
    }
    
    printf("\n=== Finding Roots of Univariate Resultant ===\n");
    
    // 将 fq_mvpoly_t 转换为 fq_nmod_poly_t
    fq_nmod_poly_t poly;
    fq_nmod_poly_init(poly, result->ctx);
    
    // 转换：将参数多项式转为单变量多项式
    for (slong i = 0; i < result->nterms; i++) {
        slong degree = 0;
        if (result->terms[i].par_exp) {
            degree = result->terms[i].par_exp[0];
        }
        fq_nmod_poly_set_coeff(poly, degree, result->terms[i].coeff, result->ctx);
    }
    
    printf("Univariate polynomial in %s:\n", state->par_names[0]);
    printf("  ");
    fq_nmod_poly_print_pretty(poly, state->par_names[0], result->ctx);
    printf("\n");
    
    slong degree = fq_nmod_poly_degree(poly, result->ctx);
    if (degree <= 0) {
        printf("  Polynomial is constant or zero, no roots to find.\n");
        fq_nmod_poly_clear(poly, result->ctx);
        return;
    }
    
    // 使用FLINT内置的因式分解函数
    printf("\nFactoring polynomial...\n");
    
    fq_nmod_poly_factor_t fac;
    fq_nmod_poly_factor_init(fac, result->ctx);

    // 使用FLINT的Cantor-Zassenhaus算法进行因式分解
    printf("\nfactor_kaltofen_shoup: ");
    fq_nmod_poly_factor_kaltofen_shoup(fac, poly, result->ctx); //squarefree_part fq_nmod_poly_factor_cantor_zassenhaus fq_nmod_poly_factor_kaltofen_shoup

    // 输出因式分解结果
    printf("\nFactorization: ");
    
    // 输出各个因式（不处理常数因子，因为结构体可能没有这个成员）
    for (slong i = 0; i < fac->num; i++) {
        if (i > 0) printf(" * ");
        
        printf("(");
        fq_nmod_poly_print_pretty(fac->poly + i, state->par_names[0], result->ctx);
        printf(")");
        
        if (fac->exp[i] > 1) {
            printf("^%ld", fac->exp[i]);
        }
    }
    printf("\n");
    
    // 找出所有根（度数为1的因式）
    printf("\nRoots:\n");
    slong root_count = 0;
    
    for (slong i = 0; i < fac->num; i++) {
        if (fq_nmod_poly_degree(fac->poly + i, result->ctx) == 1) {
            // 线性因式 ax + b = 0，根为 -b/a
            fq_nmod_t a, b, root;
            fq_nmod_init(a, result->ctx);
            fq_nmod_init(b, result->ctx);
            fq_nmod_init(root, result->ctx);
            
            fq_nmod_poly_get_coeff(a, fac->poly + i, 1, result->ctx);
            fq_nmod_poly_get_coeff(b, fac->poly + i, 0, result->ctx);
            
            if (!fq_nmod_is_zero(a, result->ctx)) {
                fq_nmod_neg(root, b, result->ctx);
                fq_nmod_div(root, root, a, result->ctx);
                
                printf("  %s = ", state->par_names[0]);
                fq_nmod_print_pretty_enhanced(root, result->ctx);
                
                // 检查重数
                slong total_multiplicity = fac->exp[i];

                if (total_multiplicity > 1) {
                    printf(" (multiplicity %ld)", total_multiplicity);
                }
                printf("\n");
                root_count++;
            }
            
            fq_nmod_clear(a, result->ctx);
            fq_nmod_clear(b, result->ctx);
            fq_nmod_clear(root, result->ctx);
        }
    }
    
    if (root_count == 0) {
        printf("  No linear factors found (no roots in this field)\n");
        
        // 显示非线性因子的信息
        printf("\nNon-linear irreducible factors:\n");
        for (slong i = 0; i < fac->num; i++) {
            if (fq_nmod_poly_degree(fac->poly + i, result->ctx) > 1) {
                printf("  Factor of degree %ld: ", fq_nmod_poly_degree(fac->poly + i, result->ctx));
                fq_nmod_poly_print_pretty(fac->poly + i, state->par_names[0], result->ctx);
                if (fac->exp[i] > 1) {
                    printf(" (multiplicity %ld)", fac->exp[i]);
                }
                printf("\n");
            }
        }
    } else {
        printf("\nTotal roots found: %ld\n", root_count);
    }
    
    // 清理
    fq_nmod_poly_clear(poly, result->ctx);
    fq_nmod_poly_factor_clear(fac, result->ctx);
}

void compute_fq_dixon_resultant_string_enhanced(const char **poly_strings, slong npoly_strings,
                                               const char **var_names, slong nvars,
                                               const fq_nmod_ctx_t ctx) {
    printf("\n=== Enhanced Dixon Resultant String Interface over F_{p^d} ===\n");
    printf("Field: F_%ld^%ld\n", fq_nmod_ctx_prime(ctx), fq_nmod_ctx_degree(ctx));
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
    
    // 初始化解析器状态
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
    
    if (fq_nmod_ctx_degree(ctx) > 1) {
        state.generator_name = strdup("t");
    } else {
        state.generator_name = NULL;
    }
    
    // 第一遍：识别参数
    printf("\nFirst pass: identifying parameters...\n");
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
    
    printf("Detected parameters (%ld): ", state.npars);
    if (state.npars == 0) {
        printf("none");
    } else {
        for (slong i = 0; i < state.npars; i++) {
            if (i > 0) printf(", ");
            printf("%s", state.par_names[i]);
        }
    }
    printf("\n");
    
    // 第二遍：解析多项式
    fq_mvpoly_t *polys = (fq_mvpoly_t*) malloc(npoly_strings * sizeof(fq_mvpoly_t));
    printf("\nParsing polynomials:\n");
    
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
        
        printf("  p%ld = %s => ", i, poly_strings[i]);
        fq_mvpoly_print_enhanced(&polys[i], "");
    }
    
    // 计算Dixon结果式
    printf("\nComputing Dixon resultant...\n");
    fq_mvpoly_t result;
    
    clock_t start = clock();
    fq_dixon_resultant(&result, polys, nvars, state.npars);
    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    printf("  Dixon Resultant computation time: %.6f seconds\n", elapsed);
    
    printf("\n=== Final Result ===\n");
    printf("Dixon Resultant: ");
    //fq_mvpoly_print_enhanced(&result, "");
    printf("\n");
    find_and_print_roots_of_univariate_resultant(&result, &state);

    // 清理
    fq_mvpoly_clear(&result);
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
}

// ============= 测试函数 =============

// ============= 测试函数 =============


int test_fq_string_interface_enhanced(void) {
    printf("=== Enhanced String Interface Tests ===\n");
    
    // 先运行基础插值测试
    //test_interpolation_basic();
    
    fmpz_t p;
    fmpz_init_set_ui(p, 7); //10000000000037 65539
    
    fq_nmod_ctx_t ctx;
    fq_nmod_ctx_init(ctx, p, 4, "t");
    
    //test_gf2n(ctx);
    
    // 测试1：简单线性多项式（期望非平凡结果）
    {
        printf("\n=== Test 1: Simple Linear Polynomials ===\n");
        const char *polys[] = {"x + 1", "x + 2"};  // 简化：不使用扩域元素
        const char *vars[] = {"x"};
        
        //compute_fq_dixon_resultant_string_enhanced(polys, 2, vars, 1, ctx);
    }
    
    // 测试2：扩域元素
    {
        printf("\n=== Test 2: Extension Field Elements ===\n");
        const char *polys[] = {"x^2 + t", "x + t + 1"};
        const char *vars[] = {"x"};
        
       //compute_fq_dixon_resultant_string_enhanced(polys, 2, vars, 1, ctx);
    }
    
    // 测试3：带参数
    {
        printf("\n=== Test 3: With Parameters ===\n");
        const char *polys[] = {"x + a", "x + b"};
        const char *vars[] = {"x"};
        
        //compute_fq_dixon_resultant_string_enhanced(polys, 2, vars, 1, ctx);
    }
    
    // 测试4：二次多项式
    {
        printf("\n=== Test 4: Quadratic Polynomials ===\n");
        const char *polys[] = {"x^2 + 1", "x^2 + 2"};
        const char *vars[] = {"x"};
        
        //compute_fq_dixon_resultant_string_enhanced(polys, 2, vars, 1, ctx);
    }
    // 测试5：扩域多参数
    {
        printf("\n=== Test 5: Multivariate Case ===\n");
        const char *polys[] = {"x^5 + b*t + 1", "x*y - y^4 + 1 + t", "2*x^3 + y + y^2 "};
        const char *vars[] = {"x", "y"};
        
        //compute_fq_dixon_resultant_string_enhanced(polys, 3, vars, 2, ctx);
    }

    {
        printf("\n\n=== Example 4: Complex Expressions ===\n");
        const char *polys[] = { "-585*x0^4 + 14923*x0^3*x1 - 3928*x0^3*x2 - 8464*x0^3*x3 + 5329*x0^3 + 27412*x0^2*x1^2 + 1094*x0^2*x1*x2 - 2005*x0^2*x1*x3 + 28944*x0^2*x1 + 22840*x0^2*x2^2 - 10535*x0^2*x2*x3 + 3457*x0^2*x2 + 32670*x0^2*x3^2 + 11702*x0^2*x3 - 17263*x0^2 + 27814*x0*x1^3 - 19297*x0*x1^2*x2 - 29347*x0*x1^2*x3 - 10223*x0*x1^2 - 14947*x0*x1*x2^2 - 1048*x0*x1*x2*x3 - 11306*x0*x1*x2 + 14780*x0*x1*x3^2 - 20216*x0*x1*x3 + 27719*x0*x1 + 21721*x0*x2^3 - 1464*x0*x2^2*x3 - 24140*x0*x2^2 + 21754*x0*x2*x3^2 - 20227*x0*x2*x3 - 31448*x0*x2 + 10143*x0*x3^3 + 30903*x0*x3^2 - 18862*x0*x3 - 5172*x0 - 1185*x1^4 - 22600*x1^3*x2 + 9472*x1^3*x3 + 21169*x1^3 - 20168*x1^2*x2^2 - 5782*x1^2*x2*x3 - 31797*x1^2*x2 - 17539*x1^2*x3^2 - 29206*x1^2*x3 - 7310*x1^2 + 9715*x1*x2^3 + 10112*x1*x2^2*x3 + 13469*x1*x2^2 - 6933*x1*x2*x3^2 - 14440*x1*x2*x3 + 8089*x1*x2 - 6386*x1*x3^3 + 15349*x1*x3^2 - 3285*x1*x3 + 31344*x1 + 17380*x2^4 - 31115*x2^3*x3 - 28546*x2^3 - 7940*x2^2*x3^2 + 27435*x2^2*x3 - 22563*x2^2 - 19834*x2*x3^3 + 5916*x2*x3^2 - 29294*x2*x3 - 30168*x2 - 27241*x3^4 - 26171*x3^3 - 19939*x3^2 - 7237*x3 + 24868",
 "-28771*x0^4 - 20161*x0^3*x1 - 31702*x0^3*x2 + 25024*x0^3*x3 - 17797*x0^3 + 27327*x0^2*x1^2 - 4494*x0^2*x1*x2 + 3178*x0^2*x1*x3 + 4381*x0^2*x1 + 10497*x0^2*x2^2 - 20930*x0^2*x2*x3 + 12893*x0^2*x2 - 2534*x0^2*x3^2 - 27089*x0^2*x3 - 16990*x0^2 - 12798*x0*x1^3 + 9172*x0*x1^2*x2 + 25601*x0*x1^2*x3 - 21230*x0*x1^2 - 15073*x0*x1*x2^2 + 23790*x0*x1*x2*x3 + 22627*x0*x1*x2 + 28469*x0*x1*x3^2 - 23855*x0*x1*x3 + 28368*x0*x1 - 18116*x0*x2^3 + 6054*x0*x2^2*x3 - 31329*x0*x2^2 - 16016*x0*x2*x3^2 + 3757*x0*x2*x3 - 11228*x0*x2 + 11945*x0*x3^3 + 2315*x0*x3^2 - 12289*x0*x3 + 27541*x0 - 27221*x1^4 - 32213*x1^3*x2 - 9703*x1^3*x3 - 14193*x1^3 - 32686*x1^2*x2^2 + 11760*x1^2*x2*x3 - 31417*x1^2*x2 - 10266*x1^2*x3^2 - 7288*x1^2*x3 - 32441*x1^2 + 10689*x1*x2^3 - 5585*x1*x2^2*x3 + 28776*x1*x2^2 + 19462*x1*x2*x3^2 + 28483*x1*x2*x3 - 4065*x1*x2 - 1271*x1*x3^3 + 10434*x1*x3^2 + 22460*x1*x3 + 10467*x1 + 24206*x2^4 - 15606*x2^3*x3 + 25632*x2^3 + 25054*x2^2*x3^2 - 29792*x2^2*x3 - 63*x2^2 + 26866*x2*x3^3 - 4898*x2*x3^2 - 7018*x2*x3 + 397*x2 - 19783*x3^4 - 13037*x3^3 - 19552*x3^2 - 10090*x3 - 18178",
 "19485*x0^4 - 12119*x0^3*x1 - 20508*x0^3*x2 + 15895*x0^3*x3 + 10627*x0^3 + 26409*x0^2*x1^2 + 27632*x0^2*x1*x2 - 29208*x0^2*x1*x3 + 31438*x0^2*x1 + 357*x0^2*x2^2 - 30108*x0^2*x2*x3 + 7794*x0^2*x2 - 2160*x0^2*x3^2 - 13209*x0^2*x3 + 7692*x0^2 - 11962*x0*x1^3 - 6701*x0*x1^2*x2 - 30933*x0*x1^2*x3 - 16765*x0*x1^2 + 21230*x0*x1*x2^2 + 6622*x0*x1*x2*x3 + 14832*x0*x1*x2 - 22700*x0*x1*x3^2 + 2813*x0*x1*x3 - 9763*x0*x1 - 6631*x0*x2^3 + 9297*x0*x2^2*x3 - 12132*x0*x2^2 - 21306*x0*x2*x3^2 - 32139*x0*x2*x3 - 1285*x0*x2 + 30342*x0*x3^3 - 30902*x0*x3^2 - 32615*x0*x3 + 6274*x0 + 3727*x1^4 - 7901*x1^3*x2 - 30549*x1^3*x3 + 26464*x1^3 + 14564*x1^2*x2^2 - 22875*x1^2*x2*x3 + 3195*x1^2*x2 + 23231*x1^2*x3^2 - 27313*x1^2*x3 + 10907*x1^2 + 29154*x1*x2^3 - 20699*x1*x2^2*x3 - 14759*x1*x2^2 + 14422*x1*x2*x3^2 + 22147*x1*x2*x3 - 16909*x1*x2 + 19270*x1*x3^3 - 16004*x1*x3^2 - 7608*x1*x3 + 6670*x1 + 28122*x2^4 + 21341*x2^3*x3 - 24516*x2^3 + 22543*x2^2*x3^2 + 12598*x2^2*x3 - 32358*x2^2 + 18980*x2*x3^3 + 11624*x2*x3^2 - 2740*x2*x3 + 4929*x2 + 18145*x3^4 - 13753*x3^3 - 18253*x3^2 - 11884*x3 + 4498",
 "4795*x0^4 - 17971*x0^3*x1 + 5765*x0^3*x2 - 1562*x0^3*x3 - 15820*x0^3 - 31691*x0^2*x1^2 - 16304*x0^2*x1*x2 - 5072*x0^2*x1*x3 + 4490*x0^2*x1 - 4364*x0^2*x2^2 - 27706*x0^2*x2*x3 - 11034*x0^2*x2 + 25489*x0^2*x3^2 + 26272*x0^2*x3 - 31536*x0^2 + 25525*x0*x1^3 - 7619*x0*x1^2*x2 - 8305*x0*x1^2*x3 + 15103*x0*x1^2 - 16695*x0*x1*x2^2 + 9293*x0*x1*x2*x3 + 6242*x0*x1*x2 + 29229*x0*x1*x3^2 - 27009*x0*x1*x3 + 19627*x0*x1 + 18591*x0*x2^3 + 32765*x0*x2^2*x3 - 11819*x0*x2^2 - 26584*x0*x2*x3^2 - 9711*x0*x2*x3 + 30840*x0*x2 + 5366*x0*x3^3 - 24400*x0*x3^2 + 11899*x0*x3 - 4895*x0 - 3125*x1^4 - 842*x1^3*x2 - 28348*x1^3*x3 - 10207*x1^3 - 14544*x1^2*x2^2 - 6393*x1^2*x2*x3 + 14210*x1^2*x2 + 16809*x1^2*x3^2 + 8830*x1^2*x3 + 8114*x1^2 + 13278*x1*x2^3 + 774*x1*x2^2*x3 + 22721*x1*x2^2 + 30690*x1*x2*x3^2 + 2942*x1*x2*x3 - 3914*x1*x2 - 6040*x1*x3^3 - 31907*x1*x3^2 - 2565*x1*x3 - 27219*x1 - 32605*x2^4 - 24512*x2^3*x3 - 412*x2^3 + 5143*x2^2*x3^2 + 12702*x2^2*x3 - 7353*x2^2 - 17992*x2*x3^3 - 15995*x2*x3^2 - 9729*x2*x3 - 30394*x2 - 29315*x3^4 + 12932*x3^3 + 6271*x3^2 + 4796*x3 - 26221 + 1"};
        const char *vars[] = {"x1", "x2", "x3"};
        
        //compute_fq_dixon_resultant_string_enhanced(polys, 4, vars, 3, ctx);
    }
        printf("\n\n=== Example 5: Complex Expressions ===\n");
        const char *polys[] = {           "6961056222999*x1^2 + 492241114309*x1*x2 + 1274654167662*x1*x3 + 8248733682592*x1*x4 + 7088180658272*x1*x5 + 6773861618471*x1*x6 + 7339937066764*x1*x7 + 7311080247713*x1 + 8784848395267*x2^2 + 7471681288516*x2*x3 + 2475968771181*x2*x4 +         9327810548116*x2*x5 + 593958634874*x2*x6 + 7129174689130*x2*x7 + 6049004074731*x2 + 4625604968565*x3^2 + 6253185820602*x3*x4 + 4550991172366*x3*x5 + 5378080585168*x3*x6 + 5482271208701*x3*x7 + 6948043327275*x3 + 5985345136653*x4^2 + 1521244182100*x4*x5        + 5427442402500*x4*x6 + 6436319561526*x4*x7 + 1538672672701*x4 + 9446438301283*x5^2 + 2070362544260*x5*x6 + 5597710908209*x5*x7 + 45067938677*x5 + 6235923167381*x6^2 + 9888660942683*x6*x7 + 6429624798156*x6 + 9439301034271*x7^2 + 51404897377*x7 +         4328023707296",
    "2348031235192*x1^2 + 8665717997373*x1*x2 + 5246007641560*x1*x3 + 9906985446320*x1*x4 + 1039739533794*x1*x5 + 3853637975857*x1*x6 + 3650900250017*x1*x7 + 5617699165163*x1 + 2486685992099*x2^2 + 6516539225264*x2*x3 + 3796214611184*x2*x4 +         6976017446629*x2*x5 + 7027850677083*x2*x6 + 8610707260921*x2*x7 + 9470137480304*x2 + 4765282525408*x3^2 + 9951682888238*x3*x4 + 5731542754562*x3*x5 + 2106154516382*x3*x6 + 1365708569491*x3*x7 + 3802616301549*x3 + 8390225575205*x4^2 +         5842857649957*x4*x5 + 8435359496033*x4*x6 + 1765210986854*x4*x7 + 2058094945956*x4 + 1446485795263*x5^2 + 8233330984410*x5*x6 + 6819040096454*x5*x7 + 548979904995*x5 + 6154573408206*x6^2 + 372431578074*x6*x7 + 95541340093*x6 + 9154025983328*x7^2 +         9236826549547*x7 + 7573449439705",
    "1376242065561*x1^2 + 8281510400554*x1*x2 + 8275213650884*x1*x3 + 2991100950551*x1*x4 + 9844687839744*x1*x5 + 5065687809850*x1*x6 + 3401904569535*x1*x7 + 6668096791669*x1 + 8095313386300*x2^2 + 3951774123601*x2*x3 + 4609356106739*x2*x4 +         7185063247060*x2*x5 + 5976818254726*x2*x6 + 6039118347590*x2*x7 + 6004800535866*x2 + 4930214006646*x3^2 + 2191696971195*x3*x4 + 8254465276952*x3*x5 + 7035814644556*x3*x6 + 4522137159631*x3*x7 + 3731215527237*x3 + 1752510625401*x4^2 +         9609370892494*x4*x5 + 1708026545088*x4*x6 + 3561067996038*x4*x7 + 6292460346011*x4 + 8036365373799*x5^2 + 4816943952435*x5*x6 + 2793898191778*x5*x7 + 8316491495713*x5 + 2175832242583*x6^2 + 4145446018746*x6*x7 + 2145438328560*x6 + 6445662523382*x7^2 +         8877850745983*x7 + 7788550250222",
    "3559720601762*x1^2 + 4080977166341*x1*x2 + 1755793730130*x1*x3 + 4475742260260*x1*x4 + 2545535073536*x1*x5 + 4671895765873*x1*x6 + 453927826882*x1*x7 + 5759562037718*x1 + 7466580129504*x2^2 + 9828672885574*x2*x3 + 7266471059649*x2*x4 +         7382385170429*x2*x5 + 9797161901557*x2*x6 + 1475716872709*x2*x7 + 2582683189886*x2 + 964584327707*x3^2 + 6274061111252*x3*x4 + 5477416323441*x3*x5 + 2769316928042*x3*x6 + 5790043274201*x3*x7 + 961126433202*x3 + 6318241497741*x4^2 + 1140014781892*x4*x5         + 9082144622035*x4*x6 + 4963555971326*x4*x7 + 7010224962159*x4 + 523561394929*x5^2 + 7301801285728*x5*x6 + 3571091639441*x5*x7 + 1906730880464*x5 + 9491899909283*x6^2 + 6553740071787*x6*x7 + 1480990517066*x6 + 9386126648881*x7^2 + 1537502847455*x7 +         3944014457456",
    "3713624134692*x1^2 + 2541938465999*x1*x2 + 4395223983599*x1*x3 + 5416153866146*x1*x4 + 7775613115521*x1*x5 + 3152397964382*x1*x6 + 3119776432181*x1*x7 + 2513215122872*x1 + 4090819413318*x2^2 + 960531507917*x2*x3 + 49395406568*x2*x4 + 6569386367431*x2*x5         + 7682851588529*x2*x6 + 912607515854*x2*x7 + 1057763180753*x2 + 9781545003119*x3^2 + 4729709595794*x3*x4 + 5505234634701*x3*x5 + 3524492443516*x3*x6 + 5383050753103*x3*x7 + 4908712572370*x3 + 8159329278496*x4^2 + 4086211279042*x4*x5 +         1654452506304*x4*x6 + 8908811858439*x4*x7 + 8473746374670*x4 + 7248881767004*x5^2 + 3997291390320*x5*x6 + 7859230109045*x5*x7 + 8530727006017*x5 + 8773615892647*x6^2 + 9786603168054*x6*x7 + 4692227527267*x6 + 3841624285980*x7^2 + 5312928950427*x7 +         8492077339101",
    "2413452089182*x1^2 + 8402051508691*x1*x2 + 9685333107887*x1*x3 + 1155556599672*x1*x4 + 5720001126783*x1*x5 + 1677225499737*x1*x6 + 2071280909477*x1*x7 + 7182221207528*x1 + 2943712649644*x2^2 + 1527726386466*x2*x3 + 5736882321523*x2*x4 +         2188609936636*x2*x5 + 6355204981854*x2*x6 + 7936037963359*x2*x7 + 1540732408408*x2 + 9861793599163*x3^2 + 5709760553348*x3*x4 + 5159829063856*x3*x5 + 234366520065*x3*x6 + 5374614605580*x3*x7 + 7715243783381*x3 + 5887697928670*x4^2 + 5144499328599*x4*x5        + 2180767997056*x4*x6 + 2297919378031*x4*x7 + 1398776345086*x4 + 1300702289754*x5^2 + 7528811943040*x5*x6 + 5957499952889*x5*x7 + 7037143709573*x5 + 8263496801454*x6^2 + 935602898177*x6*x7 + 7176823835641*x6 + 7452987933870*x7^2 + 427517994874*x7 +         9714937211255 + 1"};
        const char *vars[] = {"x1", "x2", "x3", "x4", "x5"};
        
        //compute_fq_dixon_resultant_string_enhanced(polys, 6, vars, 5, ctx);

        {
        printf("\n\n=== Example 6: Complex Expressions ===\n");
        const char *polys[] = {           "    2825464172702*x1^2 + 6870255586359*x1*x2 + 9242344247242*x1*x3 + 9738893996219*x1*x4 + 3250249561545*x1*x5 + 9525611131482*x1*x6 + 6624654987995*x1*x7 + 7280452979381*x1*x8 + 3469220815334*x1 + 5327157877962*x2^2 + 8286524947129*x2*x3 +         1058141132087*x2*x4 + 9210451493676*x2*x5 + 6971821502432*x2*x6 + 9234503490900*x2*x7 + 4899732520675*x2*x8 + 9518921950936*x2 + 190109980441*x3^2 + 654155853484*x3*x4 + 9664270732993*x3*x5 + 2289821953266*x3*x6 + 3616402558089*x3*x7 +         1408852894422*x3*x8 + 7174482283625*x3 + 8707119809487*x4^2 + 2898772030489*x4*x5 + 4200281154609*x4*x6 + 741595753864*x4*x7 + 700823542070*x4*x8 + 1752434375863*x4 + 7847308255215*x5^2 + 6419860964222*x5*x6 + 5651689978877*x5*x7 + 2339141013263*x5*x8         + 5940417360423*x5 + 4499649400703*x6^2 + 2633732084679*x6*x7 + 8200329123430*x6*x8 + 5159144564069*x6 + 1397610490984*x7^2 + 7798200877297*x7*x8 + 8955846586946*x7 + 288760584336*x8^2 + 6080933409197*x8 + 6735266520966",
    "5835533935855*x1^2 + 7190870748051*x1*x2 + 3981511722653*x1*x3 + 4711141595502*x1*x4 + 5283646768807*x1*x5 + 5090484753004*x1*x6 + 9727329407658*x1*x7 + 7676897236099*x1*x8 + 1095521545446*x1 + 9030131398056*x2^2 + 1228556216776*x2*x3 +         1082060907476*x2*x4 + 7565393825677*x2*x5 + 8181448304987*x2*x6 + 2468900654559*x2*x7 + 1790464282990*x2*x8 + 5620690702578*x2 + 9306376346681*x3^2 + 590253977759*x3*x4 + 9436261711115*x3*x5 + 3217198110443*x3*x6 + 1913365893105*x3*x7 +         5788534695471*x3*x8 + 3706178435907*x3 + 8898630171164*x4^2 + 5321571094527*x4*x5 + 814710992764*x4*x6 + 60232570337*x4*x7 + 1020747180347*x4*x8 + 4392435142779*x4 + 507792742442*x5^2 + 3426361236375*x5*x6 + 8972989051544*x5*x7 + 5734720702450*x5*x8 +         1623847996559*x5 + 9284787330376*x6^2 + 1682063505165*x6*x7 + 8988655169395*x6*x8 + 1667312996778*x6 + 4887944787224*x7^2 + 555933645052*x7*x8 + 8060295445452*x7 + 7885565295815*x8^2 + 485181628702*x8 + 9473906894406",
    "5274986957202*x1^2 + 8465732171893*x1*x2 + 9625664501464*x1*x3 + 7401359807190*x1*x4 + 3284532756844*x1*x5 + 1013724577334*x1*x6 + 5418859038721*x1*x7 + 6557246859364*x1*x8 + 6956705192074*x1 + 5218532968633*x2^2 + 6475525409808*x2*x3 +         4681689498109*x2*x4 + 5943282216151*x2*x5 + 494842691368*x2*x6 + 8027497191060*x2*x7 + 2037978156767*x2*x8 + 1055149664660*x2 + 262719020723*x3^2 + 8735348447139*x3*x4 + 9648773525863*x3*x5 + 4242848209530*x3*x6 + 6614379091009*x3*x7 +         3082331593259*x3*x8 + 3446430951914*x3 + 8992914372154*x4^2 + 5962851415420*x4*x5 + 8897703676560*x4*x6 + 7724518923982*x4*x7 + 4711893432105*x4*x8 + 8336439258390*x4 + 4588836740201*x5^2 + 2747653864188*x5*x6 + 1195316076375*x5*x7 +         2795884884979*x5*x8 + 5299365602726*x5 + 5007417727497*x6^2 + 4057239920175*x6*x7 + 2013495471890*x6*x8 + 1806683157665*x6 + 2132114442837*x7^2 + 321460169046*x7*x8 + 3222563566662*x7 + 851474122609*x8^2 + 1301379919483*x8 + 6349397142229",
    "230431862101*x1^2 + 6024429517696*x1*x2 + 5986304206279*x1*x3 + 2442082472494*x1*x4 + 943366133531*x1*x5 + 7614997581981*x1*x6 + 51547929154*x1*x7 + 2978734415598*x1*x8 + 9732980789798*x1 + 2990007770583*x2^2 + 5743133201747*x2*x3 + 8482078203917*x2*x4 +        7783928957666*x2*x5 + 8509262577094*x2*x6 + 5853968560763*x2*x7 + 5113354887819*x2*x8 + 6283841863666*x2 + 5093801541240*x3^2 + 7256194201363*x3*x4 + 1339541041795*x3*x5 + 9079678021369*x3*x6 + 2234753446622*x3*x7 + 6427035317217*x3*x8 +         4925138296173*x3 + 9424131144102*x4^2 + 3616761343772*x4*x5 + 1619347419218*x4*x6 + 9388417322529*x4*x7 + 3108128708938*x4*x8 + 5814791994811*x4 + 5007923669650*x5^2 + 808194948332*x5*x6 + 4603202389751*x5*x7 + 5868020501415*x5*x8 + 370765682338*x5 +         6277044020592*x6^2 + 8094710804601*x6*x7 + 4861168049527*x6*x8 + 2029772637281*x6 + 2483857686360*x7^2 + 8710010398430*x7*x8 + 430423593764*x7 + 8855247769739*x8^2 + 2380493001782*x8 + 6733776138518",
    "2905203928094*x1^2 + 855070481295*x1*x2 + 9650290408720*x1*x3 + 4394482635089*x1*x4 + 5386738931686*x1*x5 + 3204165176498*x1*x6 + 6365650300252*x1*x7 + 1981020407152*x1*x8 + 8751293744519*x1 + 6406221448748*x2^2 + 8269828921854*x2*x3 +         4304189635474*x2*x4 + 2750242559288*x2*x5 + 3218880851847*x2*x6 + 4590551559197*x2*x7 + 8374092191398*x2*x8 + 5609446611792*x2 + 2347958877904*x3^2 + 1228592569818*x3*x4 + 9384221591007*x3*x5 + 2646222497077*x3*x6 + 2837034523107*x3*x7 +         2512475318883*x3*x8 + 722297287982*x3 + 5449514040295*x4^2 + 1136948405779*x4*x5 + 2221274862805*x4*x6 + 7356948792359*x4*x7 + 9536939737291*x4*x8 + 7696532562298*x4 + 8683078581698*x5^2 + 4261374172522*x5*x6 + 2405306252891*x5*x7 + 5804405898597*x5*x8        + 6813298704716*x5 + 197504542053*x6^2 + 1262827401637*x6*x7 + 1709212322176*x6*x8 + 2564727310759*x6 + 2460021384481*x7^2 + 9694085434494*x7*x8 + 6525850640077*x7 + 50495568255*x8^2 + 283017692156*x8 + 679321993640",
    "909579613031*x1^2 + 6293354409992*x1*x2 + 3761384787102*x1*x3 + 5096232797092*x1*x4 + 9594334269382*x1*x5 + 3920843730762*x1*x6 + 8450677899081*x1*x7 + 3287207720085*x1*x8 + 4446104370981*x1 + 9375105788242*x2^2 + 9583521355092*x2*x3 +         7104187635989*x2*x4 + 488346971573*x2*x5 + 4925660362084*x2*x6 + 7246166492675*x2*x7 + 3531935574352*x2*x8 + 6356461426130*x2 + 3536151818477*x3^2 + 1006680877124*x3*x4 + 3781574971988*x3*x5 + 3443197113789*x3*x6 + 8132304929146*x3*x7 +         6627940471504*x3*x8 + 7191509915126*x3 + 5247780048922*x4^2 + 6054447524333*x4*x5 + 6736569712002*x4*x6 + 3622998025051*x4*x7 + 4451123969205*x4*x8 + 5292115148294*x4 + 7919143060202*x5^2 + 2313458921694*x5*x6 + 1076004527854*x5*x7 +         7549935104066*x5*x8 + 6486140207015*x5 + 6941034171343*x6^2 + 9744626435321*x6*x7 + 2079736930839*x6*x8 + 5490315141398*x6 + 9209173885792*x7^2 + 6351265018083*x7*x8 + 3911183763980*x7 + 8685068527617*x8^2 + 8671184186697*x8 + 1406939367478",
    "5578710292472*x1^2 + 4318153353493*x1*x2 + 3434989559513*x1*x3 + 9294906577800*x1*x4 + 7199577272552*x1*x5 + 5365969511254*x1*x6 + 1529369398128*x1*x7 + 6011808063726*x1*x8 + 6125819752160*x1 + 1820631496016*x2^2 + 3659943799209*x2*x3 +         7684096389882*x2*x4 + 5674883037443*x2*x5 + 7874128659215*x2*x6 + 8516623700049*x2*x7 + 8413559627727*x2*x8 + 2111272444325*x2 + 2515694213787*x3^2 + 2548117540516*x3*x4 + 6352554835873*x3*x5 + 7476314734071*x3*x6 + 5994644182512*x3*x7 +         1889828472090*x3*x8 + 4153142629736*x3 + 7833549614966*x4^2 + 7598954585883*x4*x5 + 9903304726748*x4*x6 + 3717752476402*x4*x7 + 5805665406385*x4*x8 + 3105543762902*x4 + 4489115636473*x5^2 + 6642338825019*x5*x6 + 8566851088982*x5*x7 +         6855463649546*x5*x8 + 9810918269302*x5 + 5268640782626*x6^2 + 3700621395197*x6*x7 + 1630791521772*x6*x8 + 2559051398973*x6 + 4494855164328*x7^2 + 5543234649871*x7*x8 + 7111511250184*x7 + 2447882925116*x8^2 + 8340309672103*x8 + 4270832115855",
    "8461226579276*x1^2 + 2672707933238*x1*x2 + 9447719049176*x1*x3 + 2813809456541*x1*x4 + 9677784751625*x1*x5 + 9548268195115*x1*x6 + 8837060277717*x1*x7 + 4567870605927*x1*x8 + 6866902004667*x1 + 4346731802030*x2^2 + 9372589566311*x2*x3 +         5401523791564*x2*x4 + 6786641321102*x2*x5 + 2875008932728*x2*x6 + 7385669963682*x2*x7 + 7433266642424*x2*x8 + 4107139949295*x2 + 6408695516940*x3^2 + 9973704190457*x3*x4 + 3937584973585*x3*x5 + 889172193736*x3*x6 + 7462683108502*x3*x7 +         6970483744939*x3*x8 + 6262731667755*x3 + 5951802966481*x4^2 + 4087836083783*x4*x5 + 383254098112*x4*x6 + 3541621082488*x4*x7 + 2244118418080*x4*x8 + 2494998792929*x4 + 2601486130518*x5^2 + 6754644554077*x5*x6 + 5316575544986*x5*x7 + 7836392979460*x5*x8        + 3291563847335*x5 + 1200712009965*x6^2 + 3391378889449*x6*x7 + 8331897914304*x6*x8 + 6305681643035*x6 + 1729468153857*x7^2 + 8269722316039*x7*x8 + 1322252188605*x7 + 646792563643*x8^2 + 1524008296638*x8 + 2694776922533 + 1"};
        const char *vars[] = {"x1", "x2", "x3", "x4", "x5", "x6", "x7"};
        
        //compute_fq_dixon_resultant_string_enhanced(polys, 8, vars, 7, ctx);
    }
         // Example 6: With more complex expressions
    {
        printf("\n\n=== Example 6: Complex Expressions ===\n");
        const char *polys[] = {  "-4*x0^10 - 2*x0^9*x1 - 3*x0^9*x2 + 5*x0^9 + 5*x0^8*x1^2 - 2*x0^8*x1*x2 + 4*x0^8*x1 - 3*x0^8*x2^2 + 5*x0^8*x2 + 3*x0^7*x1^3 - 3*x0^7*x1^2*x2 - 2*x0^7*x1^2 + 5*x0^7*x1*x2^2 + 2*x0^7*x1*x2 + 3*x0^7*x1 - 2*x0^7*x2^3 - 4*x0^7*x2^2 + 4*x0^7*x2 - 5*x0^7 - 3*x0^6*x1^4 - 2*x0^6*x1^3*x2 - 5*x0^6*x1^3 + 5*x0^6*x1^2*x2^2 + 3*x0^6*x1*x2^3 - x0^6*x1*x2^2 - 3*x0^6*x1*x2 - 2*x0^6*x1 + 4*x0^6*x2^4 + x0^6*x2^3 - 4*x0^6*x2^2 + 3*x0^6*x2 + 3*x0^6 + 4*x0^5*x1^5 + 3*x0^5*x1^4*x2 + 5*x0^5*x1^4 - 2*x0^5*x1^3*x2^2 - 2*x0^5*x1^3 - 2*x0^5*x1^2*x2^3 + 2*x0^5*x1^2*x2^2 + 4*x0^5*x1^2 - 3*x0^5*x1*x2^4 - 4*x0^5*x1*x2^3 + 5*x0^5*x1*x2^2 + 4*x0^5*x1*x2 - 3*x0^5*x1 + 2*x0^5*x2^5 + 3*x0^5*x2^4 - 2*x0^5*x2^3 - 5*x0^5*x2^2 + 3*x0^5*x2 - 5*x0^5 - 5*x0^4*x1^5*x2 + 3*x0^4*x1^4*x2^2 + 3*x0^4*x1^4*x2 - 3*x0^4*x1^3*x2^3 - x0^4*x1^3*x2^2 - 5*x0^4*x1^3*x2 + 5*x0^4*x1^3 - 4*x0^4*x1^2*x2^4 - x0^4*x1^2*x2^3 - 4*x0^4*x1^2*x2^2 + 2*x0^4*x1^2*x2 - 4*x0^4*x1^2 - 3*x0^4*x1*x2^5 - x0^4*x1*x2^4 + 3*x0^4*x1*x2^3 + 5*x0^4*x1*x2^2 - 5*x0^4*x1*x2 - 4*x0^4*x1 + 2*x0^4*x2^6 + 2*x0^4*x2^5 - 2*x0^4*x2^4 + x0^4*x2^3 - x0^4*x2^2 - 5*x0^4*x2 - 2*x0^4 - 4*x0^3*x1^7 + 2*x0^3*x1^6*x2 + 5*x0^3*x1^6 + 3*x0^3*x1^5*x2^2 + 5*x0^3*x1^5*x2 - 5*x0^3*x1^5 + 2*x0^3*x1^4*x2^3 - 2*x0^3*x1^4*x2^2 - 2*x0^3*x1^4*x2 + 5*x0^3*x1^4 + x0^3*x1^3*x2^4 + x0^3*x1^3*x2^3 - 3*x0^3*x1^3*x2^2 + 2*x0^3*x1^3*x2 - 2*x0^3*x1^3 + 3*x0^3*x1^2*x2^5 + 2*x0^3*x1^2*x2^4 - x0^3*x1^2*x2^3 + 4*x0^3*x1^2*x2^2 + 3*x0^3*x1^2*x2 - 5*x0^3*x1^2 + 2*x0^3*x1*x2^6 - 2*x0^3*x1*x2^5 + 3*x0^3*x1*x2^4 - 5*x0^3*x1*x2^2 - 2*x0^3*x1*x2 - 4*x0^3*x1 + x0^3*x2^7 - 5*x0^3*x2^6 - x0^3*x2^5 + 4*x0^3*x2^4 + 5*x0^3*x2^3 - x0^3*x2^2 - 3*x0^3*x2 - 4*x0^3 + 5*x0^2*x1^8 - 3*x0^2*x1^7*x2 + 2*x0^2*x1^7 - x0^2*x1^6*x2 - 4*x0^2*x1^6 - x0^2*x1^5*x2^3 - 5*x0^2*x1^5*x2^2 + 3*x0^2*x1^5*x2 + 4*x0^2*x1^5 + 5*x0^2*x1^4*x2^4 + x0^2*x1^4*x2^3 - 5*x0^2*x1^4*x2 + 4*x0^2*x1^4 + 4*x0^2*x1^3*x2^5 + 4*x0^2*x1^3*x2^4 - 4*x0^2*x1^3*x2^3 + 2*x0^2*x1^3*x2^2 + 3*x0^2*x1^3 + 4*x0^2*x1^2*x2^6 - x0^2*x1^2*x2^5 + 3*x0^2*x1^2*x2^4 - 5*x0^2*x1^2*x2^2 + 2*x0^2*x1^2*x2 - 3*x0^2*x1^2 - 4*x0^2*x1*x2^7 + x0^2*x1*x2^6 + 2*x0^2*x1*x2^5 + 5*x0^2*x1*x2^4 + x0^2*x1*x2^3 - 5*x0^2*x1*x2^2 - x0^2*x1 - 5*x0^2*x2^8 - 5*x0^2*x2^7 + 4*x0^2*x2^6 - 5*x0^2*x2^5 - 4*x0^2*x2^4 + 5*x0^2*x2^3 + 3*x0^2*x2^2 - 5*x0^2*x2 + 3*x0^2 - 5*x0*x1^9 + 4*x0*x1^8*x2 - 4*x0*x1^7*x2^2 - 3*x0*x1^7*x2 + 5*x0*x1^7 - x0*x1^6*x2^3 + x0*x1^6*x2^2 - x0*x1^6*x2 - 5*x0*x1^5*x2^4 + 4*x0*x1^5*x2^3 - 4*x0*x1^5*x2 + 3*x0*x1^4*x2^5 - x0*x1^4*x2^4 + 5*x0*x1^4*x2^3 + 2*x0*x1^4*x2^2 + x0*x1^4*x2 - 4*x0*x1^4 + 2*x0*x1^3*x2^6 - 2*x0*x1^3*x2^5 - 2*x0*x1^3*x2^4 - 3*x0*x1^3*x2 - 2*x0*x1^3 - 4*x0*x1^2*x2^7 + 4*x0*x1^2*x2^6 + 4*x0*x1^2*x2^5 - 3*x0*x1^2*x2^4 + 5*x0*x1^2*x2^3 + 5*x0*x1^2*x2^2 + 2*x0*x1^2*x2 - x0*x1^2 - x0*x1*x2^7 - 5*x0*x1*x2^5 - 4*x0*x1*x2^4 + 2*x0*x1*x2^3 + 5*x0*x1*x2^2 + 4*x0*x1*x2 + x0*x1 - 3*x0*x2^9 + 5*x0*x2^8 - 3*x0*x2^6 - 2*x0*x2^5 - x0*x2^4 - 4*x0*x2^3 + 2*x0*x2^2 + 5*x0*x2 - 2*x0 + x1^10 - 4*x1^9*x2 - 3*x1^9 - 4*x1^8*x2^2 - x1^8*x2 - 2*x1^8 + 5*x1^7*x2^3 - 3*x1^7*x2^2 - 5*x1^7*x2 + 5*x1^7 - 2*x1^6*x2^4 + 5*x1^6*x2^3 - 5*x1^6*x2^2 + x1^6*x2 - 4*x1^6 - 2*x1^5*x2^5 - x1^5*x2^4 - 4*x1^5*x2^3 + 2*x1^5*x2^2 + 2*x1^5*x2 - 4*x1^5 + x1^4*x2^6 + x1^4*x2^5 + 4*x1^4*x2^4 - 5*x1^4*x2^3 - x1^4*x2^2 + 4*x1^4*x2 + 5*x1^4 + 3*x1^3*x2^7 - 5*x1^3*x2^6 + 2*x1^3*x2^5 + 5*x1^3*x2^4 + x1^3*x2^2 - x1^3*x2 + 2*x1^3 + 2*x1^2*x2^8 - 5*x1^2*x2^6 - 3*x1^2*x2^5 - 5*x1^2*x2^3 + x1^2*x2^2 + x1^2 - 4*x1*x2^9 + 5*x1*x2^8 - 5*x1*x2^7 - 3*x1*x2^6 + 3*x1*x2^4 + 5*x1*x2^3 - x1*x2^2 - 2*x1*x2 + 2*x1 - 3*x2^10 - 3*x2^9 - 3*x2^8 + x2^7 + 4*x2^6 - 2*x2^5 - 4*x2^4 + 2*x2^3 + 4*x2^2 + 4*x2 + 2",
 "-4*x0^10 + x0^9*x2 - x0^9 + 4*x0^8*x1^2 + 3*x0^8*x1*x2 + 5*x0^8*x1 - 4*x0^8*x2^2 + x0^8*x2 + 5*x0^8 + 2*x0^7*x1^3 - 3*x0^7*x1^2*x2 - 4*x0^7*x1^2 - 3*x0^7*x1*x2^2 + 4*x0^7*x1*x2 + 4*x0^7*x1 + 2*x0^7*x2^3 - 4*x0^7*x2^2 + 5*x0^7*x2 + 5*x0^7 - 3*x0^6*x1^4 + 5*x0^6*x1^3*x2 + 2*x0^6*x1^3 - 3*x0^6*x1^2*x2^2 + x0^6*x1^2*x2 - 2*x0^6*x1*x2^3 - 3*x0^6*x1*x2^2 - 5*x0^6*x1*x2 + 3*x0^6*x1 - x0^6*x2^4 + x0^6*x2^3 + 3*x0^6*x2^2 + 3*x0^6*x2 + 3*x0^6 + x0^5*x1^5 + x0^5*x1^4*x2 + 5*x0^5*x1^4 + x0^5*x1^3*x2^2 - 2*x0^5*x1^3*x2 + 3*x0^5*x1^3 + 3*x0^5*x1^2*x2^3 + 3*x0^5*x1^2*x2^2 + x0^5*x1^2*x2 + 2*x0^5*x1^2 - x0^5*x1*x2^4 - 3*x0^5*x1*x2^2 + 4*x0^5*x1*x2 - 5*x0^5*x1 + 4*x0^5*x2^5 + 4*x0^5*x2^4 - 4*x0^5*x2^3 - 2*x0^5*x2^2 - 4*x0^5*x2 - 4*x0^5 + 4*x0^4*x1^6 - 3*x0^4*x1^5*x2 + 2*x0^4*x1^5 + 5*x0^4*x1^4*x2^2 - 4*x0^4*x1^4*x2 + 3*x0^4*x1^4 + 2*x0^4*x1^3*x2^3 - 5*x0^4*x1^3*x2^2 - 5*x0^4*x1^3*x2 + 4*x0^4*x1^3 - 2*x0^4*x1^2*x2^4 + 5*x0^4*x1^2*x2^3 - x0^4*x1^2*x2^2 - x0^4*x1^2*x2 - 2*x0^4*x1^2 - 3*x0^4*x1*x2^5 - 2*x0^4*x1*x2^4 - 4*x0^4*x1*x2^3 - 4*x0^4*x1*x2^2 - x0^4*x1*x2 + 5*x0^4*x1 - 3*x0^4*x2^6 - x0^4*x2^5 + 4*x0^4*x2^4 + 2*x0^4*x2^3 - 2*x0^4*x2^2 + 3*x0^4*x2 + 3*x0^4 + 3*x0^3*x1^7 + 5*x0^3*x1^6*x2 + 2*x0^3*x1^6 + x0^3*x1^5*x2^2 - 3*x0^3*x1^5*x2 - 4*x0^3*x1^5 + x0^3*x1^4*x2^3 + 2*x0^3*x1^4*x2 + 4*x0^3*x1^3*x2^4 + 4*x0^3*x1^3*x2^3 - 2*x0^3*x1^3*x2^2 - 4*x0^3*x1^3*x2 - 5*x0^3*x1^3 + 3*x0^3*x1^2*x2^5 - x0^3*x1^2*x2^4 - 3*x0^3*x1^2*x2^3 + 2*x0^3*x1^2*x2^2 + 3*x0^3*x1^2*x2 + 3*x0^3*x1^2 + 2*x0^3*x1*x2^6 - 4*x0^3*x1*x2^5 - 2*x0^3*x1*x2^4 - 4*x0^3*x1*x2^3 - 2*x0^3*x1*x2^2 + 5*x0^3*x1*x2 + x0^3*x2^6 + 4*x0^3*x2^5 + 3*x0^3*x2^4 + 2*x0^3*x2^3 - 2*x0^3*x2^2 + 2*x0^3*x2 - x0^3 - 2*x0^2*x1^8 - 3*x0^2*x1^7*x2 - x0^2*x1^7 - x0^2*x1^6*x2 + 2*x0^2*x1^6 + 3*x0^2*x1^5*x2^3 + 3*x0^2*x1^5*x2^2 - 3*x0^2*x1^5*x2 - 2*x0^2*x1^5 + x0^2*x1^4*x2^4 - 4*x0^2*x1^4*x2^3 + 4*x0^2*x1^4*x2^2 + x0^2*x1^4*x2 + 4*x0^2*x1^4 - 3*x0^2*x1^3*x2^5 + 2*x0^2*x1^3*x2^4 - x0^2*x1^3*x2^3 + 3*x0^2*x1^3*x2^2 - 4*x0^2*x1^3 - x0^2*x1^2*x2^6 + 3*x0^2*x1^2*x2^5 - x0^2*x1^2*x2^4 + 5*x0^2*x1^2*x2^3 + 5*x0^2*x1^2*x2^2 + 5*x0^2*x1^2 - 4*x0^2*x1*x2^7 + 2*x0^2*x1*x2^6 - 3*x0^2*x1*x2^5 - 5*x0^2*x1*x2^3 + x0^2*x1*x2^2 + 3*x0^2*x1*x2 - 2*x0^2*x1 - 2*x0^2*x2^7 - 2*x0^2*x2^6 - 4*x0^2*x2^5 - 3*x0^2*x2^4 + 4*x0^2*x2^3 - 3*x0^2*x2^2 - x0^2*x2 - 5*x0^2 - x0*x1^9 + 3*x0*x1^8*x2 + x0*x1^8 - 5*x0*x1^7*x2^2 - 2*x0*x1^7*x2 + x0*x1^7 + x0*x1^6*x2^3 - 3*x0*x1^6*x2 + 3*x0*x1^6 - 2*x0*x1^5*x2^4 + 5*x0*x1^5*x2^3 - x0*x1^5*x2^2 - 4*x0*x1^5*x2 + 4*x0*x1^5 - 5*x0*x1^4*x2^5 + x0*x1^4*x2^4 - 2*x0*x1^4*x2^3 - 5*x0*x1^4*x2^2 + 2*x0*x1^4*x2 + 5*x0*x1^4 - 3*x0*x1^3*x2^6 + 2*x0*x1^3*x2^5 - 2*x0*x1^3*x2^4 + x0*x1^3*x2^3 + 2*x0*x1^3*x2^2 + 5*x0*x1^3*x2 + 4*x0*x1^3 - 4*x0*x1^2*x2^7 + 3*x0*x1^2*x2^6 + x0*x1^2*x2^5 - x0*x1^2*x2^4 - 4*x0*x1^2*x2^3 - 2*x0*x1^2*x2^2 + 4*x0*x1^2*x2 + 2*x0*x1^2 + 2*x0*x1*x2^8 + x0*x1*x2^7 - 3*x0*x1*x2^6 + x0*x1*x2^4 - x0*x1*x2^3 - 3*x0*x1*x2^2 - x0*x1*x2 + 2*x0*x1 - 2*x0*x2^9 - 5*x0*x2^8 + 4*x0*x2^7 - 3*x0*x2^6 + 3*x0*x2^5 + 2*x0*x2^4 + 4*x0*x2^3 + 3*x0*x2^2 + 5*x0*x2 + 5*x0 + 5*x1^10 - 2*x1^9*x2 + 5*x1^9 + 3*x1^8*x2^2 - 4*x1^8*x2 - 5*x1^8 + 2*x1^7*x2^3 + 3*x1^7*x2^2 + 3*x1^7*x2 - 4*x1^7 - x1^6*x2^3 + 3*x1^6*x2^2 + 3*x1^6*x2 - 2*x1^6 - 5*x1^5*x2^5 + x1^5*x2^4 + 5*x1^5*x2^3 + 2*x1^5*x2^2 + 5*x1^5*x2 + 2*x1^5 - 3*x1^4*x2^6 + 5*x1^4*x2^5 + 4*x1^4*x2^4 + 2*x1^4*x2^3 - 4*x1^4*x2^2 + 2*x1^4*x2 + 3*x1^4 + 2*x1^3*x2^7 - 4*x1^3*x2^6 + 5*x1^3*x2^5 + 5*x1^3*x2^4 - 4*x1^3*x2^3 - 2*x1^3*x2 + x1^3 + 4*x1^2*x2^8 - 2*x1^2*x2^7 - 5*x1^2*x2^6 - 5*x1^2*x2^5 + 4*x1^2*x2^3 - 3*x1^2*x2 + 3*x1^2 - 5*x1*x2^8 + 4*x1*x2^7 - 4*x1*x2^6 + 5*x1*x2^5 - 3*x1*x2^4 - x1*x2^3 - x1*x2^2 - 4*x1*x2 - 2*x2^10 - 3*x2^9 + x2^8 + x2^7 - 5*x2^6 + 3*x2^5 - 4*x2^4 + 5*x2^3 + 2*x2 - 1",
 "-3*x0^10 - 5*x0^9*x2 + x0^9 - 5*x0^8*x1^2 + 5*x0^8*x1*x2 - 2*x0^8*x1 - 3*x0^8*x2^2 + 4*x0^8 + x0^7*x1^3 - 4*x0^7*x1^2*x2 + 2*x0^7*x1^2 + x0^7*x1*x2^2 + 2*x0^7*x1 + 2*x0^7*x2^3 + 2*x0^7*x2^2 + x0^7*x2 + 4*x0^7 + 4*x0^6*x1^4 + x0^6*x1^3*x2 - 4*x0^6*x1^3 - 2*x0^6*x1^2*x2^2 + 2*x0^6*x1^2*x2 + 3*x0^6*x1^2 + x0^6*x1*x2^2 + 3*x0^6*x1*x2 - 3*x0^6*x1 + 3*x0^6*x2^4 + 5*x0^6*x2^3 - 5*x0^6*x2^2 - x0^6*x2 - x0^6 - 3*x0^5*x1^5 + 2*x0^5*x1^4*x2 + x0^5*x1^4 + x0^5*x1^3*x2^2 + 4*x0^5*x1^3*x2 + 2*x0^5*x1^3 - 3*x0^5*x1^2*x2^3 + 5*x0^5*x1^2*x2^2 - x0^5*x1^2*x2 - 3*x0^5*x1^2 - 3*x0^5*x1*x2^4 - 5*x0^5*x1*x2^3 + 3*x0^5*x1*x2^2 + 3*x0^5*x1*x2 + 5*x0^5*x1 + 5*x0^5*x2^5 - x0^5*x2^4 + 3*x0^5*x2^3 - 5*x0^5*x2^2 + 3*x0^5*x2 - 2*x0^5 + 2*x0^4*x1^6 + 4*x0^4*x1^5*x2 + 2*x0^4*x1^5 - 3*x0^4*x1^4*x2^2 + 3*x0^4*x1^4*x2 + x0^4*x1^4 + 5*x0^4*x1^3*x2^3 - 2*x0^4*x1^3*x2^2 + 4*x0^4*x1^3*x2 - x0^4*x1^3 + 3*x0^4*x1^2*x2^4 - 4*x0^4*x1^2*x2^3 + 3*x0^4*x1^2*x2^2 - x0^4*x1^2*x2 - 4*x0^4*x1^2 - x0^4*x1*x2^5 - x0^4*x1*x2^4 + 5*x0^4*x1*x2^3 + 2*x0^4*x1*x2^2 + 3*x0^4*x1*x2 + 5*x0^4*x1 + 4*x0^4*x2^6 + 3*x0^4*x2^5 - 5*x0^4*x2^4 + 4*x0^4*x2^3 - 2*x0^4*x2^2 + 2*x0^4*x2 - x0^4 - 2*x0^3*x1^7 + 3*x0^3*x1^6*x2 - 4*x0^3*x1^6 - 2*x0^3*x1^5*x2^2 - 2*x0^3*x1^5 + 2*x0^3*x1^4*x2^3 + 3*x0^3*x1^4*x2^2 - 3*x0^3*x1^4 - 3*x0^3*x1^3*x2^3 + 2*x0^3*x1^3*x2^2 + 4*x0^3*x1^3*x2 - 5*x0^3*x1^3 - 4*x0^3*x1^2*x2^5 - 2*x0^3*x1^2*x2^4 - 5*x0^3*x1^2*x2^3 - 3*x0^3*x1^2*x2^2 - 5*x0^3*x1^2*x2 + 4*x0^3*x1^2 + 3*x0^3*x1*x2^6 + 3*x0^3*x1*x2^5 + 5*x0^3*x1*x2^3 - 4*x0^3*x1*x2^2 + 2*x0^3*x1*x2 + 4*x0^3*x1 - x0^3*x2^7 + 4*x0^3*x2^6 - 5*x0^3*x2^5 + 3*x0^3*x2^4 + 3*x0^3*x2^3 + 5*x0^3*x2^2 + 5*x0^3*x2 - 4*x0^3 - 3*x0^2*x1^8 + 4*x0^2*x1^7*x2 - 5*x0^2*x1^7 - 2*x0^2*x1^6*x2^2 - 3*x0^2*x1^6*x2 + 4*x0^2*x1^5*x2^3 - 4*x0^2*x1^5*x2^2 - 2*x0^2*x1^5*x2 + 2*x0^2*x1^5 + 2*x0^2*x1^4*x2^4 + x0^2*x1^4*x2^3 + x0^2*x1^4*x2^2 - 3*x0^2*x1^4*x2 - 2*x0^2*x1^4 + 2*x0^2*x1^3*x2^5 + 3*x0^2*x1^3*x2^4 + 2*x0^2*x1^3*x2^3 + 2*x0^2*x1^3*x2 - x0^2*x1^3 + 4*x0^2*x1^2*x2^6 - x0^2*x1^2*x2^5 - 4*x0^2*x1^2*x2^4 + 4*x0^2*x1^2*x2^3 + x0^2*x1^2*x2^2 + 3*x0^2*x1^2*x2 - 4*x0^2*x1^2 + 3*x0^2*x1*x2^7 - 4*x0^2*x1*x2^5 - 4*x0^2*x1*x2^4 + 4*x0^2*x1*x2^2 + 3*x0^2*x1*x2 - 4*x0^2*x1 + 5*x0^2*x2^8 + 3*x0^2*x2^7 - x0^2*x2^6 - 5*x0^2*x2^5 - 2*x0^2*x2^4 - 2*x0^2*x2^3 - 5*x0^2*x2^2 - 2*x0^2*x2 + x0*x1^8 + 4*x0*x1^7*x2^2 + 2*x0*x1^7*x2 + 5*x0*x1^7 - 5*x0*x1^6*x2^3 - x0*x1^6*x2^2 - 3*x0*x1^6 - 5*x0*x1^5*x2^4 + 4*x0*x1^5*x2^3 - 4*x0*x1^5*x2^2 + 3*x0*x1^5*x2 + 4*x0*x1^5 - x0*x1^4*x2^5 - 4*x0*x1^4*x2^4 + 2*x0*x1^4*x2^3 + 4*x0*x1^4*x2^2 + 4*x0*x1^4*x2 + 5*x0*x1^4 - x0*x1^3*x2^6 - 3*x0*x1^3*x2^5 - 5*x0*x1^3*x2^4 + x0*x1^3*x2^3 + 3*x0*x1^3*x2^2 + 2*x0*x1^3*x2 - 5*x0*x1^3 - 5*x0*x1^2*x2^7 + x0*x1^2*x2^5 - 2*x0*x1^2*x2^4 - 3*x0*x1^2*x2^3 + 5*x0*x1^2*x2^2 + 3*x0*x1^2*x2 - 5*x0*x1^2 + 5*x0*x1*x2^8 + 5*x0*x1*x2^7 + x0*x1*x2^6 + 4*x0*x1*x2^5 - 2*x0*x1*x2^4 - 2*x0*x1*x2^3 - 5*x0*x1*x2^2 + 5*x0*x1 - 3*x0*x2^9 + 3*x0*x2^7 - 4*x0*x2^6 + 5*x0*x2^4 - x0*x2^3 - 2*x0*x2^2 - 3*x0*x2 + 3*x0 - 2*x1^10 + 5*x1^9 + x1^8*x2^2 - 3*x1^8*x2 - 5*x1^8 + 3*x1^7*x2^3 - 3*x1^7*x2^2 - 4*x1^7*x2 + 3*x1^7 + 5*x1^6*x2^4 - 5*x1^6*x2^3 + 2*x1^6*x2^2 - 3*x1^6*x2 - 2*x1^5*x2^5 - x1^5*x2^4 - 3*x1^5*x2^3 - 2*x1^5*x2^2 + 4*x1^5*x2 + 4*x1^5 - x1^4*x2^6 + 2*x1^4*x2^5 + x1^4*x2^4 + 5*x1^4*x2^3 + x1^4*x2^2 - 5*x1^4*x2 - 2*x1^4 + x1^3*x2^7 - x1^3*x2^6 - 4*x1^3*x2^5 + x1^3*x2^4 + 4*x1^3*x2^3 - 4*x1^3*x2^2 + 4*x1^3*x2 + 5*x1^3 - 4*x1^2*x2^8 + 4*x1^2*x2^7 - 5*x1^2*x2^6 + x1^2*x2^5 - 5*x1^2*x2^4 - x1^2*x2^3 + 3*x1^2*x2^2 - 4*x1^2*x2 + 3*x1^2 + 2*x1*x2^9 - 3*x1*x2^7 - 5*x1*x2^6 + 2*x1*x2^5 - 3*x1*x2^4 + 5*x1*x2^3 - 4*x1*x2^2 - 3*x1*x2 + 2*x1 - 3*x2^10 - 4*x2^9 + 4*x2^8 - 3*x2^7 + 5*x2^6 + 3*x2^5 + 2*x2^4 - 4*x2^3 - 5*x2^2 - 2*x2 + t"};
        const char *vars[] = {"x1", "x2"};
        compute_fq_dixon_resultant_string_enhanced(polys, 3, vars, 2, ctx);//10000000000037
    }
    fq_nmod_ctx_clear(ctx);
    fmpz_clear(p);
    
    return 0;
}

#endif // COMPLETE_FIXED_DIXON_INTERFACE_H