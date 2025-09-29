// Complete fixed Dixon resultant string interface implementation
#include "dixon_interface_flint.h"

// Fixed string parser implementation

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

void next_token(parser_state_t *state) {
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

static slong find_or_add_parameter(parser_state_t *state, const char *name) {
    // Check if it's a variable
    for (slong i = 0; i < state->nvars; i++) {
        if (strcmp(state->var_names[i], name) == 0) {
            return -1;
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

void parse_primary(parser_state_t *state, fq_mvpoly_t *poly) {
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

void parse_factor(parser_state_t *state, fq_mvpoly_t *poly) {
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

void parse_term(parser_state_t *state, fq_mvpoly_t *poly) {
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

void parse_expression(parser_state_t *state, fq_mvpoly_t *poly) {
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




void find_and_print_roots_of_univariate_resultant(const fq_mvpoly_t *result, parser_state_t *state) {
    // Check if there are elimination variables
    if (result->nvars != 0) {
        return;  // Still has elimination variables, not final result
    }
    
    // Check actual number of parameters used
    int *par_used = (int*) calloc(result->npars, sizeof(int));
    slong actual_par_count = 0;
    slong main_par_idx = -1;
    
    // Iterate through all terms to count actually used parameters
    for (slong t = 0; t < result->nterms; t++) {
        if (result->terms[t].par_exp) {
            for (slong p = 0; p < result->npars; p++) {
                if (result->terms[t].par_exp[p] > 0 && !par_used[p]) {
                    par_used[p] = 1;
                    main_par_idx = p;  // Remember the last used parameter
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
    
    // actual_par_count == 1, proceed with root finding
    printf("Single parameter detected! Finding roots...\n");
    
    // Convert fq_mvpoly_t to fq_nmod_poly_t
    fq_nmod_poly_t poly;
    fq_nmod_poly_init(poly, result->ctx);
    
    // Convert: transform parameter polynomial to univariate polynomial
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
    
// Use FLINT's root finding algorithm - optimized version, detect field type
printf("\nFinding roots using appropriate algorithm...\n");

degree = fq_nmod_poly_degree(poly, result->ctx);
if (degree <= 0) {
    printf("  Polynomial is constant or zero, no roots to find.\n");
    fq_nmod_poly_clear(poly, result->ctx);
    free(par_used);
    return;
}

// Detect field type
slong field_degree = fq_nmod_ctx_degree(result->ctx);
printf("Field degree: %ld\n", field_degree);
slong root_count = 0;
slong total_multiplicity = 0;

if (field_degree == 1) {
    // Prime field case - use more efficient nmod_poly_roots
    printf("Prime field detected, using nmod_poly_roots...\n");
    
    mp_limb_t prime = fq_nmod_ctx_prime(result->ctx);
    
    // Convert to nmod_poly_t
    nmod_poly_t nmod_poly;
    nmod_poly_init(nmod_poly, prime);
    
    // Copy coefficients
    for (slong i = 0; i <= degree; i++) {
        fq_nmod_t coeff;
        fq_nmod_init(coeff, result->ctx);
        fq_nmod_poly_get_coeff(coeff, poly, i, result->ctx);
        
        // Extract from fq_nmod_t as mp_limb_t
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
    
    // Use nmod_poly_roots for root finding
    nmod_roots_t nmod_roots;
    nmod_roots_init(nmod_roots);
    slong num_roots = our_nmod_poly_roots(nmod_roots, nmod_poly, 1);  // with_multiplicity = 1
    
    // Output found roots
    printf("\nRoots found:\n");
    
    printf("Find %ld roots:\n", num_roots);
    for (slong i = 0; i < nmod_roots->num; i++) {
        printf("  Root %ld: %lu (Multiplicity: %ld)\n", i + 1, 
               nmod_roots->roots[i], nmod_roots->mult[i]);
    }
    
    // Clean up nmod related structures
    nmod_poly_clear(nmod_poly);
    
} else {
    // Extension field case - use original fq_nmod_poly_roots
    printf("Extension field detected, using fq_nmod_poly_roots...\n");
    
    fq_nmod_roots_t roots;
    fq_nmod_roots_init(roots, result->ctx);
    slong num_roots = our_fq_nmod_poly_roots(roots, poly, 1, result->ctx);
    
    printf("Find %ld roots:\n", num_roots);
    for (slong i = 0; i < roots->num; i++) {
        printf("  root %ld: ", i + 1);
        fq_nmod_print_pretty(roots->roots + i, result->ctx);
        printf(" (Multiplicity: %ld)\n", roots->mult[i]);
    }
    
    // Clean up
}

    // Clean up
    fq_nmod_poly_clear(poly, result->ctx);
    free(par_used);
}

// Helper functions

// Get field generator name
char* get_generator_name(const fq_nmod_ctx_t ctx) {
    return strdup(ctx->var);
}

// Convert fq_nmod_t to string with field generator
char* fq_nmod_to_string_with_gen(const fq_nmod_t elem, const fq_nmod_ctx_t ctx, const char *gen_name) {
    // Initially allocate a large buffer
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
        
        // Dynamically build string
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
                
                // Check if buffer expansion is needed
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

// String builder helper functions
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

// Optimized coefficient to string function
// Updated string builder version for consistent formatting
void fq_nmod_to_string_builder(string_builder_t *sb, const fq_nmod_t elem, 
                                    const fq_nmod_ctx_t ctx, const char *gen_name) {
    if (fq_nmod_is_zero(elem, ctx)) {
        sb_append(sb, "0");
        return;
    }
    
    if (fq_nmod_ctx_degree(ctx) == 1) {
        // Prime field element - no parentheses needed
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
        // Extension field element
        nmod_poly_t poly;
        nmod_poly_init(poly, fq_nmod_ctx_prime(ctx));
        fq_nmod_get_nmod_poly(poly, elem, ctx);
        
        slong deg = nmod_poly_degree(poly);
        
        // Count non-zero terms
        int term_count = 0;
        for (slong i = deg; i >= 0; i--) {
            if (nmod_poly_get_coeff_ui(poly, i) != 0) {
                term_count++;
            }
        }
        
        // Use parentheses for multi-term expressions
        int use_parens = (term_count > 1);
        
        if (use_parens) sb_append_char(sb, '(');
        
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
                        sb_append(sb, gen_name ? gen_name : "t");
                    } else {
                        sb_append_ulong(sb, coeff);
                        sb_append_char(sb, '*');
                        sb_append(sb, gen_name ? gen_name : "t");
                    }
                } else {
                    if (coeff == 1) {
                        sb_append(sb, gen_name ? gen_name : "t");
                        sb_append_char(sb, '^');
                        sb_append_long(sb, i);
                    } else {
                        sb_append_ulong(sb, coeff);
                        sb_append_char(sb, '*');
                        sb_append(sb, gen_name ? gen_name : "t");
                        sb_append_char(sb, '^');
                        sb_append_long(sb, i);
                    }
                }
            }
        }
        
        if (first) {
            sb_append(sb, "0");
        }
        
        if (use_parens) sb_append_char(sb, ')');
        
        nmod_poly_clear(poly);
    }
}

// New optimized version function - directly replace original fq_mvpoly_to_string
char* fq_mvpoly_to_string(const fq_mvpoly_t *poly, char **par_names, const char *gen_name) {
    if (poly->nterms == 0) {
        return strdup("0");
    }
    
    // Estimate initial capacity
    size_t estimated_size = poly->nterms * (50 + 10 * (poly->nvars + poly->npars));
    if (estimated_size < 1024) estimated_size = 1024;
    
    string_builder_t sb;
    sb_init(&sb, estimated_size);
    
    for (slong i = 0; i < poly->nterms; i++) {
        if (i > 0) {
            sb_append(&sb, " + ");
        }
        
        // Check if we have any variables or parameters for this term
        int has_vars_or_pars = 0;
        
        // Check variables (var_exp)
        if (poly->nvars > 0 && poly->terms[i].var_exp) {
            for (slong j = 0; j < poly->nvars; j++) {
                if (poly->terms[i].var_exp[j] > 0) {
                    has_vars_or_pars = 1;
                    break;
                }
            }
        }
        
        // Check parameters (par_exp)
        if (!has_vars_or_pars && poly->npars > 0 && poly->terms[i].par_exp) {
            for (slong j = 0; j < poly->npars; j++) {
                if (poly->terms[i].par_exp[j] > 0) {
                    has_vars_or_pars = 1;
                    break;
                }
            }
        }
        
        // Handle coefficient printing using fixed function
        fq_nmod_t one;
        fq_nmod_init(one, poly->ctx);
        fq_nmod_one(one, poly->ctx);
        
        if (fq_nmod_is_one(poly->terms[i].coeff, poly->ctx) && has_vars_or_pars) {
            // Coefficient is 1 and we have variables/parameters, don't print coefficient
        } else {
            // Add coefficient with proper parentheses using fixed function
            fq_nmod_to_string_builder(&sb, poly->terms[i].coeff, poly->ctx, gen_name);
            
            if (has_vars_or_pars) {
                sb_append_char(&sb, '*');
            }
        }
        
        fq_nmod_clear(one, poly->ctx);
        
        // Track what has been printed for this term
        int term_has_content = 0;
        
        // Add variables (var_exp) - use default variable names
        if (poly->nvars > 0 && poly->terms[i].var_exp) {
            for (slong j = 0; j < poly->nvars; j++) {
                if (poly->terms[i].var_exp[j] > 0) {
                    // FIXED: Only add * if something was already printed for this term
                    if (term_has_content) {
                        sb_append_char(&sb, '*');
                    }
                    
                    // Use standard variable names x, y, z, etc.
                    if (j < 6) {
                        char var_names[] = {'x', 'y', 'z', 'w', 'v', 'u'};
                        sb_append_char(&sb, var_names[j]);
                    } else {
                        sb_append(&sb, "x_");
                        sb_append_long(&sb, j);
                    }
                    
                    if (poly->terms[i].var_exp[j] > 1) {
                        sb_append_char(&sb, '^');
                        sb_append_long(&sb, poly->terms[i].var_exp[j]);
                    }
                    term_has_content = 1;
                }
            }
        }
        
        // Add parameters using ACTUAL parameter names
        if (poly->npars > 0 && poly->terms[i].par_exp) {
            for (slong j = 0; j < poly->npars; j++) {
                if (poly->terms[i].par_exp[j] > 0) {
                    // FIXED: Only add * if something was already printed for this term
                    if (term_has_content) {
                        sb_append_char(&sb, '*');
                    }
                    
                    // Use actual parameter names if provided
                    if (par_names && par_names[j]) {
                        sb_append(&sb, par_names[j]);
                    } else {
                        // Fallback to default parameter names
                        if (j < 4) {
                            char default_par_names[] = {'a', 'b', 'c', 'd'};
                            sb_append_char(&sb, default_par_names[j]);
                        } else {
                            sb_append(&sb, "p_");
                            sb_append_long(&sb, j);
                        }
                    }
                    
                    if (poly->terms[i].par_exp[j] > 1) {
                        sb_append_char(&sb, '^');
                        sb_append_long(&sb, poly->terms[i].par_exp[j]);
                    }
                    term_has_content = 1;
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

// Concat polynomials helper function:
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

// Core Dixon Function

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
    
    // Initialize parser state with original variable names
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
    
    // Compute Dixon resultant with original names
    fq_mvpoly_t dixon_result_poly;
    fq_dixon_resultant_with_names(&dixon_result_poly, polys, nvars, state.npars,
                                 state.var_names, state.par_names, gen_name);

    // Find roots with proper parameter names
    find_and_print_roots_of_univariate_resultant(&dixon_result_poly, &state);
    
    // Convert result to string with original parameter names
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

// Helper function: remove whitespace from both ends of string
static char* trim_whitespace(char *str) {
    if (!str) return NULL;
    
    // Remove leading whitespace
    while (isspace(*str)) str++;
    
    if (*str == '\0') return str;
    
    // Remove trailing whitespace
    char *end = str + strlen(str) - 1;
    while (end > str && isspace(*end)) end--;
    *(end + 1) = '\0';
    
    return str;
}

// Helper function: split string (by comma)
char** split_string(const char *input, slong *count) {
    if (!input || strlen(input) == 0) {
        *count = 0;
        return NULL;
    }
    
    size_t input_len = strlen(input);
    
    // Copy input string for modification - use dynamic allocation for large strings
    char *work_str = (char*) malloc(input_len + 1);
    if (!work_str) {
        *count = 0;
        return NULL;
    }
    memcpy(work_str, input, input_len + 1);
    
    // First pass: count how many elements (look for commas)
    slong num_elements = 1;
    int in_parentheses = 0;
    for (size_t i = 0; i < input_len; i++) {
        if (input[i] == '(') in_parentheses++;
        else if (input[i] == ')') in_parentheses--;
        else if (input[i] == ',' && in_parentheses == 0) {
            num_elements++;
        }
    }
    
    // Allocate result array
    char **result = (char**) malloc(num_elements * sizeof(char*));
    if (!result) {
        free(work_str);
        *count = 0;
        return NULL;
    }
    
    // Second pass: split string (manual handling, don't use strtok)
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
    
    return result;
}

// Free split string array
void free_split_strings(char **strings, slong count) {
    if (!strings) return;
    for (slong i = 0; i < count; i++) {
        if (strings[i]) free(strings[i]);
    }
    free(strings);
}

// Main Dixon Interface

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
    
    for (slong i = 0; i < num_polys; i++) {
        printf("  p%ld: %s\n", i, poly_strings[i]);
    }
    
    // Compute Dixon resultant
    char **remaining_vars = NULL;
    slong num_remaining = 0;
    
    char *result = compute_dixon_internal(poly_strings, num_polys, 
                                         elim_vars, num_elim_vars, ctx,
                                         &remaining_vars, &num_remaining);
    
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
    // Get generator name
    char *gen_name = get_generator_name(ctx);
        
    char **remaining_vars = NULL;
    slong num_remaining = 0;
    
    // Initialize parser state
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
    
    // First pass: parse to identify parameters
    fq_mvpoly_t temp1, temp2;
    fq_mvpoly_init(&temp1, 1, state.max_pars, ctx);
    fq_mvpoly_init(&temp2, 1, state.max_pars, ctx);

    state.input = poly1_str;
    state.pos = 0;
    state.len = strlen(poly1_str);
    next_token(&state);
    parse_expression(&state, &temp1);

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
    
    // Save parameters
    num_remaining = state.npars;
    if (state.npars > 0) {
        remaining_vars = (char**) malloc(state.npars * sizeof(char*));
        for (slong i = 0; i < state.npars; i++) {
            (remaining_vars)[i] = strdup(state.par_names[i]);
        }
    } else {
        remaining_vars = NULL;
    }

    // Second pass: formal parsing
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
    
    // Initialize unified field context
    field_ctx_t field_ctx;
    field_ctx_init(&field_ctx, ctx);
    
    // Create unified multivariate polynomial context
    slong total_vars = 1 + state.npars;  // elimination variable + parameters
    unified_mpoly_ctx_t unified_ctx = unified_mpoly_ctx_init(total_vars, ORD_LEX, &field_ctx);
    
    // Initialize unified polynomials
    unified_mpoly_t A = unified_mpoly_init(unified_ctx);
    unified_mpoly_t B = unified_mpoly_init(unified_ctx);
    unified_mpoly_t R = unified_mpoly_init(unified_ctx);

    // Convert first polynomial
    for (slong i = 0; i < poly1.nterms; i++) {
        field_elem_u coeff;
        ulong *exp = (ulong*) calloc(total_vars, sizeof(ulong));
        
        // Convert coefficient
        fq_nmod_to_field_elem(&coeff, poly1.terms[i].coeff, &field_ctx);
        
        // Set exponents
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

    // Convert second polynomial
    for (slong i = 0; i < poly2.nterms; i++) {
        field_elem_u coeff;
        ulong *exp = (ulong*) calloc(total_vars, sizeof(ulong));
        
        // Convert coefficient
        fq_nmod_to_field_elem(&coeff, poly2.terms[i].coeff, &field_ctx);
        
        // Set exponents
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
    
    // Compute resultant
    clock_t start = clock();
    int success = unified_mpoly_resultant(R, A, B, 0, unified_ctx);
    clock_t end = clock();
    
    if (!success) {
        printf("Unified resultant computation failed!\n");
        unified_mpoly_clear(A);
        unified_mpoly_clear(B);
        unified_mpoly_clear(R);
        unified_mpoly_ctx_clear(unified_ctx);
        
        // Clean up other resources
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
    
    // Convert result back to fq_mvpoly format
    fq_mvpoly_t result_mvpoly;
    fq_mvpoly_init(&result_mvpoly, 0, state.npars, ctx);
    
    // Convert from unified format back
    slong result_len = unified_mpoly_length(R);
    if (result_len > 0) {
        // Need to iterate through all terms of the result
        // This requires accessing the internal structure of unified_mpoly
        // Since unified_mpoly is a wrapper for nmod_mpoly or fq_nmod_mpoly
        // We need to handle based on field_id
        
        if (field_ctx.field_id == FIELD_ID_NMOD) {
            // Handle nmod case
            nmod_mpoly_struct *nmod_res = GET_NMOD_POLY(R);
            nmod_mpoly_ctx_struct *nmod_ctx = &(unified_ctx->ctx.nmod_ctx);
            
            for (slong i = 0; i < nmod_mpoly_length(nmod_res, nmod_ctx); i++) {
                ulong coeff_ui = nmod_mpoly_get_term_coeff_ui(nmod_res, i, nmod_ctx);
                ulong *exp = (ulong*) calloc(total_vars, sizeof(ulong));
                nmod_mpoly_get_term_exp_ui(exp, nmod_res, i, nmod_ctx);
                
                // Convert coefficient
                fq_nmod_t coeff;
                fq_nmod_init(coeff, ctx);
                fq_nmod_set_ui(coeff, coeff_ui, ctx);
                
                // Extract parameter exponents
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
        // Handle fq_zech case
        fq_zech_mpoly_struct *zech_res = GET_ZECH_POLY(R);
        fq_zech_mpoly_ctx_struct *zech_ctx = &(unified_ctx->ctx.zech_ctx);
        
        for (slong i = 0; i < fq_zech_mpoly_length(zech_res, zech_ctx); i++) {
            fq_zech_t zech_coeff;
            fq_zech_init(zech_coeff, zech_ctx->fqctx);
            fq_zech_mpoly_get_term_coeff_fq_zech(zech_coeff, zech_res, i, zech_ctx);
            
            // Convert to fq_nmod
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ctx);
            fq_zech_get_fq_nmod(coeff, zech_coeff, zech_ctx->fqctx);
            
            ulong *exp = (ulong*) calloc(total_vars, sizeof(ulong));
            fq_zech_mpoly_get_term_exp_ui(exp, zech_res, i, zech_ctx);
            
            // Extract parameter exponents
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
            // Handle fq_nmod case
            fq_nmod_mpoly_struct *fq_res = GET_FQ_POLY(R);
            fq_nmod_mpoly_ctx_struct *fq_ctx = &(unified_ctx->ctx.fq_ctx);
            printf("fq_mvpoly_add_term_fast\n");
            for (slong i = 0; i < fq_nmod_mpoly_length(fq_res, fq_ctx); i++) {
                fq_nmod_t coeff;
                fq_nmod_init(coeff, ctx);
                fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, fq_res, i, fq_ctx);
                
                ulong *exp = (ulong*) calloc(total_vars, sizeof(ulong));
                fq_nmod_mpoly_get_term_exp_ui(exp, fq_res, i, fq_ctx);
                
                // Extract parameter exponents
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
    // If it's a univariate polynomial, try to find roots
    find_and_print_roots_of_univariate_resultant(&result_mvpoly, &state);
    //printf("fq_mvpoly_to_string\n");
    // Convert result to string
    char *result_string;
    if (result_mvpoly.nterms == 0) {
        result_string = strdup("0");
    } else {
        result_string = fq_mvpoly_to_string(&result_mvpoly, state.par_names, gen_name);
    }
    //printf("%s", result_string);
    //printf("clean up\n");

    // Output remaining variable info
    if (num_remaining > 0) {
        printf("Remaining variables (%ld): ", num_remaining);
        for (slong i = 0; i < num_remaining; i++) {
            if (i > 0) printf(", ");
            printf("%s\n", remaining_vars[i]);
            free(remaining_vars[i]);
        }
        printf("\n");
        free(remaining_vars);
    }
    // Clean up
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

char* dixon_str(const char *poly_string,    // comma-separated polynomials
                const char *vars_string,     // comma-separated variables
                const fq_nmod_ctx_t ctx) {
    
    printf("\n=== Dixon/Resultant Computation (String Interface) ===\n");
    
    // Split input strings
    slong num_polys, num_vars;
    char **poly_array = split_string(poly_string, &num_polys);
    char **vars_array = split_string(vars_string, &num_vars);
    
    char *result = NULL;
    
    // Check if it's bivariate case
    if (num_polys == 2 && num_vars == 1) {
        printf("Using unified bivariate resultant for 2 polynomials...\n");

        
        // Call unified interface bivariate resultant computation
        result = bivariate_resultant(poly_array[0], poly_array[1], 
                                                    vars_array[0], ctx);        
    } else {
        // Use original Dixon method
        printf("Using Dixon resultant for %ld polynomials...\n", num_polys);
        
        // Convert to const char**
        const char **poly_strings = (const char**) malloc(num_polys * sizeof(char*));
        const char **elim_vars = (const char**) malloc(num_vars * sizeof(char*));
        
        for (slong i = 0; i < num_polys; i++) {
            poly_strings[i] = poly_array[i];
        }
        for (slong i = 0; i < num_vars; i++) {
            elim_vars[i] = vars_array[i];
        }
        
        // Call original dixon function
        result = dixon(poly_strings, num_polys, elim_vars, num_vars, ctx);
        
        // Clean up
        free(poly_strings);
        free(elim_vars);
    }
    
    // Clean up split strings
    free_split_strings(poly_array, num_polys);
    free_split_strings(vars_array, num_vars);
    
    return result;
}