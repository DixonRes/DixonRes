/* dixon_with_ideal_reduction.h - Dixon resultant with ideal reduction using native FLINT types */

#ifndef DIXON_WITH_IDEAL_REDUCTION_H
#define DIXON_WITH_IDEAL_REDUCTION_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/nmod.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_mpoly.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_mpoly.h>


/* Triangular ideal structure for both nmod and fq_nmod */
typedef struct {
    void **generators;            /* Array of generators (either nmod_mpoly or fq_nmod_mpoly) */
    slong num_gens;              /* Number of generators */
    slong *var_indices;          /* Which variable each generator eliminates */
    char **var_names;            /* Names of variables being reduced */
    slong max_gens;              /* Maximum number of generators */
    int is_prime_field;          /* 1 for prime field, 0 for extension field */
    union {
        nmod_mpoly_ctx_struct *nmod_ctx;  /* Changed to pointer */
        fq_nmod_mpoly_ctx_struct *fq_ctx; /* Changed to pointer */
    } ctx;
    const fq_nmod_ctx_struct *field_ctx;  /* Field context for conversions */
} unified_triangular_ideal_t;

/* Initialize triangular ideal */
void unified_triangular_ideal_init(unified_triangular_ideal_t *ideal, 
                                  slong max_gens, slong nvars, 
                                  const fq_nmod_ctx_t field_ctx) {
    ideal->generators = (void**) flint_malloc(max_gens * sizeof(void*));
    ideal->var_indices = (slong*) flint_malloc(max_gens * sizeof(slong));
    ideal->var_names = (char**) flint_malloc(max_gens * sizeof(char*));
    ideal->num_gens = 0;
    ideal->max_gens = max_gens;
    ideal->field_ctx = field_ctx;
    
    /* Initialize var_names to NULL */
    for (slong i = 0; i < max_gens; i++) {
        ideal->var_names[i] = NULL;
    }
    
    /* Determine if prime field or extension field */
    if (fq_nmod_ctx_degree(field_ctx) == 1) {
        ideal->is_prime_field = 1;
        ideal->ctx.nmod_ctx = (nmod_mpoly_ctx_struct*) flint_malloc(sizeof(nmod_mpoly_ctx_struct));
        nmod_mpoly_ctx_init(ideal->ctx.nmod_ctx, nvars, ORD_LEX, fq_nmod_ctx_prime(field_ctx));
        
        for (slong i = 0; i < max_gens; i++) {
            ideal->generators[i] = flint_malloc(sizeof(nmod_mpoly_struct));
            nmod_mpoly_init((nmod_mpoly_struct*)ideal->generators[i], ideal->ctx.nmod_ctx);
        }
    } else {
        ideal->is_prime_field = 0;
        ideal->ctx.fq_ctx = (fq_nmod_mpoly_ctx_struct*) flint_malloc(sizeof(fq_nmod_mpoly_ctx_struct));
        fq_nmod_mpoly_ctx_init(ideal->ctx.fq_ctx, nvars, ORD_LEX, field_ctx);
        
        for (slong i = 0; i < max_gens; i++) {
            ideal->generators[i] = flint_malloc(sizeof(fq_nmod_mpoly_struct));
            fq_nmod_mpoly_init((fq_nmod_mpoly_struct*)ideal->generators[i], ideal->ctx.fq_ctx);
        }
    }
}

/* Clear triangular ideal */
void unified_triangular_ideal_clear(unified_triangular_ideal_t *ideal) {
    if (ideal->is_prime_field) {
        for (slong i = 0; i < ideal->max_gens; i++) {
            nmod_mpoly_clear((nmod_mpoly_struct*)ideal->generators[i], ideal->ctx.nmod_ctx);
            flint_free(ideal->generators[i]);
        }
        nmod_mpoly_ctx_clear(ideal->ctx.nmod_ctx);
        flint_free(ideal->ctx.nmod_ctx);
    } else {
        for (slong i = 0; i < ideal->max_gens; i++) {
            fq_nmod_mpoly_clear((fq_nmod_mpoly_struct*)ideal->generators[i], ideal->ctx.fq_ctx);
            flint_free(ideal->generators[i]);
        }
        fq_nmod_mpoly_ctx_clear(ideal->ctx.fq_ctx);
        flint_free(ideal->ctx.fq_ctx);
    }
    
    /* Clear variable names */
    for (slong i = 0; i < ideal->max_gens; i++) {
        if (ideal->var_names[i] != NULL) {
            free(ideal->var_names[i]);
        }
    }
    
    flint_free(ideal->generators);
    flint_free(ideal->var_indices);
    flint_free(ideal->var_names);
}

/* Helper function to detect main variable from generator polynomial */
static slong detect_main_variable(const fq_mvpoly_t *gen, char **detected_var_name, 
                                 char **var_names_hint) {
    slong detected_var_idx = -1;
    slong vars_with_degree_3 = 0;
    
    /* Scan all terms to find maximum degree for each variable */
    slong *max_degrees = (slong*) flint_calloc(gen->nvars, sizeof(slong));
    
    for (slong t = 0; t < gen->nterms; t++) {
        if (gen->terms[t].var_exp) {
            for (slong v = 0; v < gen->nvars; v++) {
                if (gen->terms[t].var_exp[v] > max_degrees[v]) {
                    max_degrees[v] = gen->terms[t].var_exp[v];
                }
            }
        }
    }
    
    /* Find variable with degree 3 */
    for (slong v = 0; v < gen->nvars; v++) {
        if (max_degrees[v] == 3) {
            detected_var_idx = v;
            vars_with_degree_3++;
        }
    }
    
    /* Set detected variable name */
    if (detected_var_idx >= 0 && var_names_hint != NULL && var_names_hint[detected_var_idx] != NULL) {
        *detected_var_name = strdup(var_names_hint[detected_var_idx]);
    } else if (detected_var_idx >= 0) {
        /* Generate default name */
        char temp[32];
        sprintf(temp, "x_%ld", detected_var_idx);
        *detected_var_name = strdup(temp);
    }
    
    flint_free(max_degrees);
    
    if (vars_with_degree_3 != 1) {
        printf("Warning: Found %ld variables with degree 3 (expected 1)\n", vars_with_degree_3);
    }
    
    return detected_var_idx;
}

/* Add generator to ideal from fq_mvpoly with automatic variable detection */
void unified_triangular_ideal_add_generator_from_mvpoly_auto(unified_triangular_ideal_t *ideal,
                                                            const fq_mvpoly_t *gen, 
                                                            char **var_names_hint) {
    if (ideal->num_gens >= ideal->max_gens) {
        printf("Error: Maximum number of generators exceeded\n");
        return;
    }
    
    /* Detect main variable */
    char *detected_var_name = NULL;
    slong detected_var_idx = detect_main_variable(gen, &detected_var_name, var_names_hint);
    
    if (detected_var_idx < 0) {
        printf("Error: Could not detect main variable in generator\n");
        return;
    }
    
    printf("DEBUG: Detected main variable '%s' at index %ld for generator %ld\n", 
           detected_var_name, detected_var_idx, ideal->num_gens);
    
    /* Convert and store generator */
    if (ideal->is_prime_field) {
        fq_mvpoly_to_nmod_mpoly((nmod_mpoly_struct*)ideal->generators[ideal->num_gens], gen, 
                               ideal->ctx.nmod_ctx);
    } else {
        fq_mvpoly_to_fq_nmod_mpoly((fq_nmod_mpoly_struct*)ideal->generators[ideal->num_gens], gen,
                                  ideal->ctx.fq_ctx);
    }
    
    ideal->var_indices[ideal->num_gens] = detected_var_idx;
    ideal->var_names[ideal->num_gens] = detected_var_name;
    ideal->num_gens++;
}

/* Find variable index by name in the current context */
static slong find_variable_by_name(const char *var_name, char **current_var_names, slong nvars) {
    for (slong i = 0; i < nvars; i++) {
        if (current_var_names && current_var_names[i] && 
            strcmp(current_var_names[i], var_name) == 0) {
            return i;
        }
    }
    return -1;
}

/* Reduction functions with variable name awareness */
void triangular_ideal_reduce_nmod_mpoly_with_names(nmod_mpoly_t poly, 
                                                   const unified_triangular_ideal_t *ideal,
                                                   char **current_var_names) {
    if (ideal->num_gens == 0) return;
    
    /* Static counter for debug limiting */
    static slong call_count = 0;
    static const slong DEBUG_LIMIT = 5;
    int debug_enabled = 0; //(call_count < DEBUG_LIMIT);
    
    call_count++;
    
    slong nvars = nmod_mpoly_ctx_nvars(ideal->ctx.nmod_ctx);
    
    if (debug_enabled) {
        printf("\n=== DEBUG: triangular_ideal_reduce (call #%ld) ===\n", call_count);
        printf("Polynomial context has %ld variables\n", nvars);
        printf("Current variable names: ");
        if (current_var_names) {
            for (slong i = 0; i < nvars; i++) {
                printf("%s ", current_var_names[i] ? current_var_names[i] : "?");
            }
        } else {
            printf("(none provided)");
        }
        printf("\n");
        
        printf("Ideal has %ld generators:\n", ideal->num_gens);
        for (slong i = 0; i < ideal->num_gens; i++) {
            printf("  Generator %ld: reduces '%s' (stored index %ld)\n", 
                   i, ideal->var_names[i] ? ideal->var_names[i] : "?", 
                   ideal->var_indices[i]);
        }
    }
    
    /* Process each generator */
    for (slong gen_idx = ideal->num_gens - 1; gen_idx >= 0; gen_idx--) {
        /* Find current index of the variable to reduce */
        slong var_idx = -1;
        
        if (ideal->var_names[gen_idx] && current_var_names) {
            /* Use variable name to find current index */
            var_idx = find_variable_by_name(ideal->var_names[gen_idx], 
                                          current_var_names, nvars);
            
            if (debug_enabled) {
                printf("\nGenerator %ld: Looking for variable '%s'... ", 
                       gen_idx, ideal->var_names[gen_idx]);
                if (var_idx >= 0) {
                    printf("found at index %ld\n", var_idx);
                } else {
                    printf("NOT FOUND!\n");
                }
            }
        } else {
            /* Fall back to stored index */
            var_idx = ideal->var_indices[gen_idx];
            
            if (debug_enabled) {
                printf("\nGenerator %ld: Using stored index %ld\n", gen_idx, var_idx);
            }
        }
        
        if (var_idx < 0 || var_idx >= nvars) {
            if (debug_enabled) {
                printf("  Skipping: invalid variable index\n");
            }
            continue;
        }
        
        nmod_mpoly_struct *gen = (nmod_mpoly_struct*)ideal->generators[gen_idx];
        
        /* Check for terms needing reduction */
        if (debug_enabled) {
            slong terms_to_reduce = 0;
            for (slong t = 0; t < nmod_mpoly_length(poly, ideal->ctx.nmod_ctx); t++) {
                ulong *exp = (ulong*) flint_malloc(nvars * sizeof(ulong));
                nmod_mpoly_get_term_exp_ui(exp, poly, t, ideal->ctx.nmod_ctx);
                
                if (exp[var_idx] >= 3) {
                    terms_to_reduce++;
                }
                flint_free(exp);
            }
            printf("  Terms needing reduction: %ld\n", terms_to_reduce);
        }
        
        /* Extract f_i (remove leading term var^3) */
        nmod_mpoly_t f_i;
        nmod_mpoly_init(f_i, ideal->ctx.nmod_ctx);
        
        for (slong t = 0; t < nmod_mpoly_length(gen, ideal->ctx.nmod_ctx); t++) {
            ulong *exp = (ulong*) flint_malloc(nvars * sizeof(ulong));
            nmod_mpoly_get_term_exp_ui(exp, gen, t, ideal->ctx.nmod_ctx);
            
            /* Skip var^3 term */
            int is_var_cubed = (exp[var_idx] == 3);
            for (slong v = 0; v < nvars; v++) {
                if (v != var_idx && exp[v] > 0) {
                    is_var_cubed = 0;
                    break;
                }
            }
            
            if (!is_var_cubed) {
                mp_limb_t coeff = nmod_mpoly_get_term_coeff_ui(gen, t, ideal->ctx.nmod_ctx);
                coeff = nmod_neg(coeff, ideal->ctx.nmod_ctx->mod);
                nmod_mpoly_t temp;
                nmod_mpoly_init(temp, ideal->ctx.nmod_ctx);
                nmod_mpoly_set_coeff_ui_ui(temp, coeff, exp, ideal->ctx.nmod_ctx);
                nmod_mpoly_add(f_i, f_i, temp, ideal->ctx.nmod_ctx);
                nmod_mpoly_clear(temp, ideal->ctx.nmod_ctx);
            }
            flint_free(exp);
        }
        
        /* Reduce all var^k (k≥3) terms */
        nmod_mpoly_t new_poly;
        nmod_mpoly_init(new_poly, ideal->ctx.nmod_ctx);
        
        slong nterms = nmod_mpoly_length(poly, ideal->ctx.nmod_ctx);
        slong reduced_count = 0;
        
        for (slong t = 0; t < nterms; t++) {
            ulong *exp = (ulong*) flint_malloc(nvars * sizeof(ulong));
            nmod_mpoly_get_term_exp_ui(exp, poly, t, ideal->ctx.nmod_ctx);
            mp_limb_t coeff = nmod_mpoly_get_term_coeff_ui(poly, t, ideal->ctx.nmod_ctx);
            
            slong var_deg = exp[var_idx];
            if (var_deg >= 3) {
                reduced_count++;
                
                slong quotient = var_deg / 3;
                slong remainder = var_deg % 3;
                
                exp[var_idx] = remainder;
                
                /* Compute f_i^quotient */
                nmod_mpoly_t f_i_power;
                nmod_mpoly_init(f_i_power, ideal->ctx.nmod_ctx);
                if (quotient == 1) {
                    nmod_mpoly_set(f_i_power, f_i, ideal->ctx.nmod_ctx);
                } else {
                    nmod_mpoly_pow_ui(f_i_power, f_i, quotient, ideal->ctx.nmod_ctx);
                }
                
                /* Create monomial */
                nmod_mpoly_t monom;
                nmod_mpoly_init(monom, ideal->ctx.nmod_ctx);
                nmod_mpoly_set_coeff_ui_ui(monom, coeff, exp, ideal->ctx.nmod_ctx);
                
                /* Multiply and add */
                nmod_mpoly_t product;
                nmod_mpoly_init(product, ideal->ctx.nmod_ctx);
                nmod_mpoly_mul(product, monom, f_i_power, ideal->ctx.nmod_ctx);
                nmod_mpoly_add(new_poly, new_poly, product, ideal->ctx.nmod_ctx);
                
                nmod_mpoly_clear(monom, ideal->ctx.nmod_ctx);
                nmod_mpoly_clear(f_i_power, ideal->ctx.nmod_ctx);
                nmod_mpoly_clear(product, ideal->ctx.nmod_ctx);
            } else {
                /* Degree < 3, copy directly */
                nmod_mpoly_t temp;
                nmod_mpoly_init(temp, ideal->ctx.nmod_ctx);
                nmod_mpoly_set_coeff_ui_ui(temp, coeff, exp, ideal->ctx.nmod_ctx);
                nmod_mpoly_add(new_poly, new_poly, temp, ideal->ctx.nmod_ctx);
                nmod_mpoly_clear(temp, ideal->ctx.nmod_ctx);
            }
            
            flint_free(exp);
        }
        
        if (debug_enabled && reduced_count > 0) {
            printf("  Reduced %ld terms\n", reduced_count);
        }
        
        nmod_mpoly_swap(poly, new_poly, ideal->ctx.nmod_ctx);
        nmod_mpoly_clear(new_poly, ideal->ctx.nmod_ctx);
        nmod_mpoly_clear(f_i, ideal->ctx.nmod_ctx);
    }
    
    if (debug_enabled) {
        printf("\nFinal polynomial has %ld terms\n", 
               nmod_mpoly_length(poly, ideal->ctx.nmod_ctx));
        
        /* Check for remaining high-degree terms */
        printf("Checking for remaining high-degree terms:\n");
        for (slong i = 0; i < ideal->num_gens; i++) {
            if (ideal->var_names[i] && current_var_names) {
                slong var_idx = find_variable_by_name(ideal->var_names[i], 
                                                     current_var_names, nvars);
                if (var_idx >= 0) {
                    slong high_degree_count = 0;
                    slong max_degree = 0;
                    
                    for (slong t = 0; t < nmod_mpoly_length(poly, ideal->ctx.nmod_ctx); t++) {
                        ulong *exp = (ulong*) flint_malloc(nvars * sizeof(ulong));
                        nmod_mpoly_get_term_exp_ui(exp, poly, t, ideal->ctx.nmod_ctx);
                        if (exp[var_idx] >= 3) {
                            high_degree_count++;
                            if (exp[var_idx] > max_degree) {
                                max_degree = exp[var_idx];
                            }
                        }
                        flint_free(exp);
                    }
                    
                    if (high_degree_count > 0) {
                        printf("  WARNING: Variable '%s' still has %ld terms with degree >= 3 (max: %ld)\n", 
                               ideal->var_names[i], high_degree_count, max_degree);
                    } else {
                        printf("  Variable '%s': all degrees < 3 ✓\n", ideal->var_names[i]);
                    }
                }
            }
        }
        
        printf("=== END DEBUG (call #%ld) ===\n\n", call_count);
        
        if (call_count == DEBUG_LIMIT) {
            printf("*** Debug output limit reached. Suppressing further output. ***\n\n");
        }
    }
}

/* Similar function for fq_nmod_mpoly */
void triangular_ideal_reduce_fq_nmod_mpoly_with_names(fq_nmod_mpoly_t poly,
                                                      const unified_triangular_ideal_t *ideal,
                                                      char **current_var_names) {
    if (ideal->num_gens == 0) return;
    
    slong nvars = fq_nmod_mpoly_ctx_nvars(ideal->ctx.fq_ctx);
    
    /* Process each generator from highest to lowest */
    for (slong gen_idx = ideal->num_gens - 1; gen_idx >= 0; gen_idx--) {
        /* Find current index of the variable to reduce */
        slong var_idx = -1;
        
        if (ideal->var_names[gen_idx] && current_var_names) {
            var_idx = find_variable_by_name(ideal->var_names[gen_idx], 
                                          current_var_names, nvars);
        } else {
            var_idx = ideal->var_indices[gen_idx];
        }
        
        if (var_idx < 0 || var_idx >= nvars) continue;
        
        fq_nmod_mpoly_struct *gen = (fq_nmod_mpoly_struct*)ideal->generators[gen_idx];
        
        /* Extract f_i */
        fq_nmod_mpoly_t f_i;
        fq_nmod_mpoly_init(f_i, ideal->ctx.fq_ctx);
        
        slong nterms = fq_nmod_mpoly_length(gen, ideal->ctx.fq_ctx);
        for (slong t = 0; t < nterms; t++) {
            ulong *exp = (ulong*) flint_malloc(nvars * sizeof(ulong));
            fq_nmod_mpoly_get_term_exp_ui(exp, gen, t, ideal->ctx.fq_ctx);
            
            int is_var_cubed = (exp[var_idx] == 3);
            for (slong v = 0; v < nvars; v++) {
                if (v != var_idx && exp[v] > 0) {
                    is_var_cubed = 0;
                    break;
                }
            }
            
            if (!is_var_cubed) {
                fq_nmod_t coeff, neg_coeff;
                fq_nmod_init(coeff, ideal->field_ctx);
                fq_nmod_init(neg_coeff, ideal->field_ctx);
                
                fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, gen, t, ideal->ctx.fq_ctx);
                fq_nmod_neg(neg_coeff, coeff, ideal->field_ctx);
                
                fq_nmod_mpoly_t temp;
                fq_nmod_mpoly_init(temp, ideal->ctx.fq_ctx);
                fq_nmod_mpoly_set_coeff_fq_nmod_ui(temp, neg_coeff, exp, ideal->ctx.fq_ctx);
                fq_nmod_mpoly_add(f_i, f_i, temp, ideal->ctx.fq_ctx);
                fq_nmod_mpoly_clear(temp, ideal->ctx.fq_ctx);
                
                fq_nmod_clear(coeff, ideal->field_ctx);
                fq_nmod_clear(neg_coeff, ideal->field_ctx);
            }
            flint_free(exp);
        }
        
        /* Reduce all high-degree terms */
        fq_nmod_mpoly_t new_poly;
        fq_nmod_mpoly_init(new_poly, ideal->ctx.fq_ctx);
        
        slong poly_nterms = fq_nmod_mpoly_length(poly, ideal->ctx.fq_ctx);
        for (slong t = 0; t < poly_nterms; t++) {
            ulong *exp = (ulong*) flint_malloc(nvars * sizeof(ulong));
            fq_nmod_mpoly_get_term_exp_ui(exp, poly, t, ideal->ctx.fq_ctx);
            
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ideal->field_ctx);
            fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, poly, t, ideal->ctx.fq_ctx);
            
            slong var_deg = exp[var_idx];
            if (var_deg >= 3) {
                slong quotient = var_deg / 3;
                slong remainder = var_deg % 3;
                
                exp[var_idx] = remainder;
                
                fq_nmod_mpoly_t f_i_power;
                fq_nmod_mpoly_init(f_i_power, ideal->ctx.fq_ctx);
                if (quotient == 1) {
                    fq_nmod_mpoly_set(f_i_power, f_i, ideal->ctx.fq_ctx);
                } else {
                    fq_nmod_mpoly_pow_ui(f_i_power, f_i, quotient, ideal->ctx.fq_ctx);
                }
                
                fq_nmod_mpoly_t monom;
                fq_nmod_mpoly_init(monom, ideal->ctx.fq_ctx);
                fq_nmod_mpoly_set_coeff_fq_nmod_ui(monom, coeff, exp, ideal->ctx.fq_ctx);
                
                fq_nmod_mpoly_t product;
                fq_nmod_mpoly_init(product, ideal->ctx.fq_ctx);
                fq_nmod_mpoly_mul(product, monom, f_i_power, ideal->ctx.fq_ctx);
                fq_nmod_mpoly_add(new_poly, new_poly, product, ideal->ctx.fq_ctx);
                
                fq_nmod_mpoly_clear(monom, ideal->ctx.fq_ctx);
                fq_nmod_mpoly_clear(f_i_power, ideal->ctx.fq_ctx);
                fq_nmod_mpoly_clear(product, ideal->ctx.fq_ctx);
            } else {
                fq_nmod_mpoly_t temp;
                fq_nmod_mpoly_init(temp, ideal->ctx.fq_ctx);
                fq_nmod_mpoly_set_coeff_fq_nmod_ui(temp, coeff, exp, ideal->ctx.fq_ctx);
                fq_nmod_mpoly_add(new_poly, new_poly, temp, ideal->ctx.fq_ctx);
                fq_nmod_mpoly_clear(temp, ideal->ctx.fq_ctx);
            }
            
            fq_nmod_clear(coeff, ideal->field_ctx);
            flint_free(exp);
        }
        
        fq_nmod_mpoly_swap(poly, new_poly, ideal->ctx.fq_ctx);
        fq_nmod_mpoly_clear(new_poly, ideal->ctx.fq_ctx);
        fq_nmod_mpoly_clear(f_i, ideal->ctx.fq_ctx);
    }
}

/* Wrapper functions that maintain compatibility */
void triangular_ideal_reduce_nmod_mpoly(nmod_mpoly_t poly, 
                                       const unified_triangular_ideal_t *ideal) {
    triangular_ideal_reduce_nmod_mpoly_with_names(poly, ideal, NULL);
}

void triangular_ideal_reduce_fq_nmod_mpoly(fq_nmod_mpoly_t poly,
                                          const unified_triangular_ideal_t *ideal) {
    triangular_ideal_reduce_fq_nmod_mpoly_with_names(poly, ideal, NULL);
}

/* Compute determinant with reduction for nmod_mpoly matrix - with variable names */
void compute_nmod_det_with_triangular_reduction_with_names(nmod_mpoly_t det,
                                                          nmod_mpoly_struct **matrix,
                                                          slong size,
                                                          const unified_triangular_ideal_t *ideal,
                                                          char **current_var_names) {
    if (size == 0) {
        nmod_mpoly_one(det, ideal->ctx.nmod_ctx);
        return;
    }
    
    if (size == 1) {
        nmod_mpoly_set(det, &matrix[0][0], ideal->ctx.nmod_ctx);
        triangular_ideal_reduce_nmod_mpoly_with_names(det, ideal, current_var_names);
        return;
    }
    
    if (size == 2) {
        /* 2x2 matrix determinant */
        nmod_mpoly_t ad, bc, temp;
        nmod_mpoly_init(ad, ideal->ctx.nmod_ctx);
        nmod_mpoly_init(bc, ideal->ctx.nmod_ctx);
        nmod_mpoly_init(temp, ideal->ctx.nmod_ctx);
        
        nmod_mpoly_mul(ad, &matrix[0][0], &matrix[1][1], ideal->ctx.nmod_ctx);
        nmod_mpoly_mul(bc, &matrix[0][1], &matrix[1][0], ideal->ctx.nmod_ctx);
        nmod_mpoly_sub(det, ad, bc, ideal->ctx.nmod_ctx);
        
        /* Only reduce once at the end */
        triangular_ideal_reduce_nmod_mpoly_with_names(det, ideal, current_var_names);
        
        nmod_mpoly_clear(ad, ideal->ctx.nmod_ctx);
        nmod_mpoly_clear(bc, ideal->ctx.nmod_ctx);
        nmod_mpoly_clear(temp, ideal->ctx.nmod_ctx);
        return;
    }
    
    /* For larger matrices, use cofactor expansion */
    nmod_mpoly_zero(det, ideal->ctx.nmod_ctx);
    
    /* Expand along first row */
    for (slong j = 0; j < size; j++) {
        if (nmod_mpoly_is_zero(&matrix[0][j], ideal->ctx.nmod_ctx)) continue;
        
        /* Build submatrix */
        nmod_mpoly_struct **submatrix = (nmod_mpoly_struct**) 
            flint_malloc((size-1) * sizeof(nmod_mpoly_struct*));
        for (slong i = 0; i < size-1; i++) {
            submatrix[i] = (nmod_mpoly_struct*) 
                flint_malloc((size-1) * sizeof(nmod_mpoly_struct));
            for (slong k = 0; k < size-1; k++) {
                nmod_mpoly_init(&submatrix[i][k], ideal->ctx.nmod_ctx);
            }
        }
        
        /* Copy submatrix */
        for (slong i = 1; i < size; i++) {
            slong col_idx = 0;
            for (slong k = 0; k < size; k++) {
                if (k != j) {
                    nmod_mpoly_set(&submatrix[i-1][col_idx], &matrix[i][k], 
                                  ideal->ctx.nmod_ctx);
                    col_idx++;
                }
            }
        }
        
        /* Compute minor recursively */
        nmod_mpoly_t minor;
        nmod_mpoly_init(minor, ideal->ctx.nmod_ctx);
        compute_nmod_det_with_triangular_reduction_with_names(minor, submatrix, size-1, ideal, current_var_names);
        
        /* Compute contribution */
        nmod_mpoly_t contrib;
        nmod_mpoly_init(contrib, ideal->ctx.nmod_ctx);
        nmod_mpoly_mul(contrib, &matrix[0][j], minor, ideal->ctx.nmod_ctx);
        
        if (j % 2 == 0) {
            nmod_mpoly_add(det, det, contrib, ideal->ctx.nmod_ctx);
        } else {
            nmod_mpoly_sub(det, det, contrib, ideal->ctx.nmod_ctx);
        }
        
        /* Reduce periodically */
        if ((j + 1) % 3 == 0 || j == size - 1) {
            triangular_ideal_reduce_nmod_mpoly_with_names(det, ideal, current_var_names);
        }
        
        /* Cleanup */
        nmod_mpoly_clear(minor, ideal->ctx.nmod_ctx);
        nmod_mpoly_clear(contrib, ideal->ctx.nmod_ctx);
        
        for (slong i = 0; i < size-1; i++) {
            for (slong k = 0; k < size-1; k++) {
                nmod_mpoly_clear(&submatrix[i][k], ideal->ctx.nmod_ctx);
            }
            flint_free(submatrix[i]);
        }
        flint_free(submatrix);
    }
}

/* Similar function for fq_nmod_mpoly - with variable names */
void compute_fq_nmod_det_with_triangular_reduction_with_names(fq_nmod_mpoly_t det,
                                                             fq_nmod_mpoly_struct **matrix,
                                                             slong size,
                                                             const unified_triangular_ideal_t *ideal,
                                                             char **current_var_names) {
    if (size == 0) {
        fq_nmod_mpoly_one(det, ideal->ctx.fq_ctx);
        return;
    }
    
    if (size == 1) {
        fq_nmod_mpoly_set(det, &matrix[0][0], ideal->ctx.fq_ctx);
        triangular_ideal_reduce_fq_nmod_mpoly_with_names(det, ideal, current_var_names);
        return;
    }
    
    if (size == 2) {
        fq_nmod_mpoly_t ad, bc;
        fq_nmod_mpoly_init(ad, ideal->ctx.fq_ctx);
        fq_nmod_mpoly_init(bc, ideal->ctx.fq_ctx);
        
        fq_nmod_mpoly_mul(ad, &matrix[0][0], &matrix[1][1], ideal->ctx.fq_ctx);
        fq_nmod_mpoly_mul(bc, &matrix[0][1], &matrix[1][0], ideal->ctx.fq_ctx);
        fq_nmod_mpoly_sub(det, ad, bc, ideal->ctx.fq_ctx);
        
        triangular_ideal_reduce_fq_nmod_mpoly_with_names(det, ideal, current_var_names);
        
        fq_nmod_mpoly_clear(ad, ideal->ctx.fq_ctx);
        fq_nmod_mpoly_clear(bc, ideal->ctx.fq_ctx);
        return;
    }
    
    fq_nmod_mpoly_zero(det, ideal->ctx.fq_ctx);
    
    for (slong j = 0; j < size; j++) {
        if (fq_nmod_mpoly_is_zero(&matrix[0][j], ideal->ctx.fq_ctx)) continue;
        
        fq_nmod_mpoly_struct **submatrix = (fq_nmod_mpoly_struct**) 
            flint_malloc((size-1) * sizeof(fq_nmod_mpoly_struct*));
        for (slong i = 0; i < size-1; i++) {
            submatrix[i] = (fq_nmod_mpoly_struct*) 
                flint_malloc((size-1) * sizeof(fq_nmod_mpoly_struct));
            for (slong k = 0; k < size-1; k++) {
                fq_nmod_mpoly_init(&submatrix[i][k], ideal->ctx.fq_ctx);
            }
        }
        
        for (slong i = 1; i < size; i++) {
            slong col_idx = 0;
            for (slong k = 0; k < size; k++) {
                if (k != j) {
                    fq_nmod_mpoly_set(&submatrix[i-1][col_idx], &matrix[i][k], 
                                     ideal->ctx.fq_ctx);
                    col_idx++;
                }
            }
        }
        
        fq_nmod_mpoly_t minor;
        fq_nmod_mpoly_init(minor, ideal->ctx.fq_ctx);
        compute_fq_nmod_det_with_triangular_reduction_with_names(minor, submatrix, size-1, ideal, current_var_names);
        
        fq_nmod_mpoly_t contrib;
        fq_nmod_mpoly_init(contrib, ideal->ctx.fq_ctx);
        fq_nmod_mpoly_mul(contrib, &matrix[0][j], minor, ideal->ctx.fq_ctx);
        
        if (j % 2 == 0) {
            fq_nmod_mpoly_add(det, det, contrib, ideal->ctx.fq_ctx);
        } else {
            fq_nmod_mpoly_sub(det, det, contrib, ideal->ctx.fq_ctx);
        }
        
        if ((j + 1) % 3 == 0 || j == size - 1) {
            triangular_ideal_reduce_fq_nmod_mpoly_with_names(det, ideal, current_var_names);
        }
        
        fq_nmod_mpoly_clear(minor, ideal->ctx.fq_ctx);
        fq_nmod_mpoly_clear(contrib, ideal->ctx.fq_ctx);
        
        for (slong i = 0; i < size-1; i++) {
            for (slong k = 0; k < size-1; k++) {
                fq_nmod_mpoly_clear(&submatrix[i][k], ideal->ctx.fq_ctx);
            }
            flint_free(submatrix[i]);
        }
        flint_free(submatrix);
    }
}

/* Old wrapper functions for compatibility */
void compute_nmod_det_with_triangular_reduction(nmod_mpoly_t det,
                                               nmod_mpoly_struct **matrix,
                                               slong size,
                                               const unified_triangular_ideal_t *ideal) {
    compute_nmod_det_with_triangular_reduction_with_names(det, matrix, size, ideal, NULL);
}

void compute_fq_nmod_det_with_triangular_reduction(fq_nmod_mpoly_t det,
                                                  fq_nmod_mpoly_struct **matrix,
                                                  slong size,
                                                  const unified_triangular_ideal_t *ideal) {
    compute_fq_nmod_det_with_triangular_reduction_with_names(det, matrix, size, ideal, NULL);
}

/* Update the old wrappers to maintain compatibility */
void nmod_mpoly_reduce_by_ideal(nmod_mpoly_t poly, 
                               const unified_triangular_ideal_t *ideal,
                               slong max_var) {
    triangular_ideal_reduce_nmod_mpoly(poly, ideal);
}

void fq_nmod_mpoly_reduce_by_ideal(fq_nmod_mpoly_t poly,
                                  const unified_triangular_ideal_t *ideal,
                                  slong max_var) {
    triangular_ideal_reduce_fq_nmod_mpoly(poly, ideal);
}

void compute_nmod_det_with_reduction(nmod_mpoly_t det,
                                   nmod_mpoly_struct **matrix,
                                   slong size,
                                   const unified_triangular_ideal_t *ideal) {
    compute_nmod_det_with_triangular_reduction(det, matrix, size, ideal);
}

void compute_fq_nmod_det_with_reduction(fq_nmod_mpoly_t det,
                                      fq_nmod_mpoly_struct **matrix,
                                      slong size,
                                      const unified_triangular_ideal_t *ideal) {
    compute_fq_nmod_det_with_triangular_reduction(det, matrix, size, ideal);
}

/* Convert matrix and compute determinant with reduction */
void compute_det_with_reduction_from_mvpoly(fq_mvpoly_t *result,
                                           fq_mvpoly_t **matrix,
                                           slong size,
                                           const unified_triangular_ideal_t *ideal,
                                           char **current_var_names) {
    printf("Computing determinant with ideal reduction (matrix size: %ld x %ld)\n", size, size);
    
    if (ideal->is_prime_field) {
        printf("Using nmod_mpoly for prime field computation\n");
        
        /* Convert to nmod_mpoly matrix */
        nmod_mpoly_struct **nmod_matrix = (nmod_mpoly_struct**) flint_malloc(size * sizeof(nmod_mpoly_struct*));
        for (slong i = 0; i < size; i++) {
            nmod_matrix[i] = (nmod_mpoly_struct*) flint_malloc(size * sizeof(nmod_mpoly_struct));
            for (slong j = 0; j < size; j++) {
                nmod_mpoly_init(&nmod_matrix[i][j], ideal->ctx.nmod_ctx);
                fq_mvpoly_to_nmod_mpoly(&nmod_matrix[i][j], &matrix[i][j], ideal->ctx.nmod_ctx);
            }
        }
        
        /* Compute determinant with reduction */
        nmod_mpoly_t det;
        nmod_mpoly_init(det, ideal->ctx.nmod_ctx);
        
        clock_t start = clock();
        /* Use the version with variable names */
        compute_nmod_det_with_triangular_reduction_with_names(det, nmod_matrix, size, ideal, current_var_names);
        clock_t end = clock();
        printf("Determinant computation time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
        
        /* Apply final reduction with variable names */
        if (current_var_names) {
            triangular_ideal_reduce_nmod_mpoly_with_names(det, ideal, current_var_names);
        }
        
        /* Convert result back */
        nmod_mpoly_to_fq_mvpoly(result, det, 0, matrix[0][0].npars, ideal->ctx.nmod_ctx, ideal->field_ctx);
        
        /* Cleanup */
        nmod_mpoly_clear(det, ideal->ctx.nmod_ctx);
        for (slong i = 0; i < size; i++) {
            for (slong j = 0; j < size; j++) {
                nmod_mpoly_clear(&nmod_matrix[i][j], ideal->ctx.nmod_ctx);
            }
            flint_free(nmod_matrix[i]);
        }
        flint_free(nmod_matrix);
        
    } else {
        printf("Using fq_nmod_mpoly for extension field computation\n");
        
        /* Similar for fq_nmod_mpoly */
        fq_nmod_mpoly_struct **fq_matrix = (fq_nmod_mpoly_struct**) flint_malloc(size * sizeof(fq_nmod_mpoly_struct*));
        for (slong i = 0; i < size; i++) {
            fq_matrix[i] = (fq_nmod_mpoly_struct*) flint_malloc(size * sizeof(fq_nmod_mpoly_struct));
            for (slong j = 0; j < size; j++) {
                fq_nmod_mpoly_init(&fq_matrix[i][j], ideal->ctx.fq_ctx);
                fq_mvpoly_to_fq_nmod_mpoly(&fq_matrix[i][j], &matrix[i][j], ideal->ctx.fq_ctx);
            }
        }
        
        fq_nmod_mpoly_t det;
        fq_nmod_mpoly_init(det, ideal->ctx.fq_ctx);
        
        clock_t start = clock();
        /* Use the version with variable names */
        compute_fq_nmod_det_with_triangular_reduction_with_names(det, fq_matrix, size, ideal, current_var_names);
        clock_t end = clock();
        printf("Determinant computation time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
        
        if (current_var_names) {
            triangular_ideal_reduce_fq_nmod_mpoly_with_names(det, ideal, current_var_names);
        }
        
        fq_nmod_mpoly_to_fq_mvpoly(result, det, 0, matrix[0][0].npars, ideal->ctx.fq_ctx, ideal->field_ctx);
        
        fq_nmod_mpoly_clear(det, ideal->ctx.fq_ctx);
        for (slong i = 0; i < size; i++) {
            for (slong j = 0; j < size; j++) {
                fq_nmod_mpoly_clear(&fq_matrix[i][j], ideal->ctx.fq_ctx);
            }
            flint_free(fq_matrix[i]);
        }
        flint_free(fq_matrix);
    }
}

void construct_triangular_ideal_from_strings(unified_triangular_ideal_t *ideal,
                                            const char **ideal_gens,
                                            slong num_gens,
                                            const char **var_names,
                                            slong nvars,
                                            const fq_nmod_ctx_t ctx) {
    printf("Constructing triangular ideal from %ld generators\n", num_gens);
    
    // Initialize the ideal with the correct size
    unified_triangular_ideal_init(ideal, num_gens, nvars, ctx);
    
    // Print the generators being added
    printf("Triangular ideal generators:\n");
    for (slong i = 0; i < num_gens; i++) {
        printf("  g%ld: %s\n", i+1, ideal_gens[i]);
    }
    printf("\n");
    
    // Parse and add each generator
    for (slong i = 0; i < num_gens; i++) {
        // Initialize parser state
        parser_state_t state;
        state.var_names = (char**) malloc(nvars * sizeof(char*));
        for (slong j = 0; j < nvars; j++) {
            state.var_names[j] = strdup(var_names[j]);
        }
        state.nvars = nvars;
        state.npars = 0;
        state.max_pars = 16;
        state.par_names = (char**) malloc(state.max_pars * sizeof(char*));
        state.ctx = ctx;
        state.generator_name = get_generator_name(ctx);
        fq_nmod_init(state.current.value, ctx);
        state.current.str = NULL;
        
        // Parse the generator polynomial
        fq_mvpoly_t gen_poly;
        fq_mvpoly_init(&gen_poly, nvars, 0, ctx);
        
        state.input = ideal_gens[i];
        state.pos = 0;
        state.len = strlen(ideal_gens[i]);
        next_token(&state);
        
        parse_expression(&state, &gen_poly);
        
        // Add to ideal with automatic variable detection
        unified_triangular_ideal_add_generator_from_mvpoly_auto(ideal, &gen_poly, 
                                                               state.var_names);
        
        // Cleanup
        fq_mvpoly_clear(&gen_poly);
        for (slong j = 0; j < nvars; j++) {
            free(state.var_names[j]);
        }
        free(state.var_names);
        for (slong j = 0; j < state.npars; j++) {
            free(state.par_names[j]);
        }
        free(state.par_names);
        free(state.generator_name);
        fq_nmod_clear(state.current.value, ctx);
        if (state.current.str) free(state.current.str);
    }
    
    // Print summary of ideal construction
    printf("Ideal construction complete:\n");
    printf("  Number of generators: %ld\n", ideal->num_gens);
    for (slong i = 0; i < ideal->num_gens; i++) {
        printf("  Generator %ld reduces variable '%s' (index %ld)\n", 
               i, ideal->var_names[i] ? ideal->var_names[i] : "?", 
               ideal->var_indices[i]);
    }
    printf("\n");
}

/* Main Dixon with ideal reduction function */
char* dixon_with_ideal_reduction(const char **poly_strings, slong num_polys,
                                const char **elim_vars, slong num_elim_vars,
                                const fq_nmod_ctx_t ctx,
                                unified_triangular_ideal_t *ideal) {
    printf("\n=== Dixon with Ideal Reduction ===\n");
    printf("Eliminating variables: ");
    for (slong i = 0; i < num_elim_vars; i++) {
        if (i > 0) printf(", ");
        printf("%s", elim_vars[i]);
    }
    printf("\n");
    
    printf("Input polynomials:\n");
    for (slong i = 0; i < num_polys; i++) {
        printf("  p%ld: %s\n", i, poly_strings[i]);
    }
    
    /* Parse polynomials */
    parser_state_t state;
    state.var_names = (char**) malloc(num_elim_vars * sizeof(char*));
    for (slong i = 0; i < num_elim_vars; i++) {
        state.var_names[i] = strdup(elim_vars[i]);
    }
    state.nvars = num_elim_vars;
    state.npars = 0;
    state.max_pars = 16;
    state.par_names = (char**) malloc(state.max_pars * sizeof(char*));
    state.ctx = ctx;
    state.generator_name = get_generator_name(ctx);
    fq_nmod_init(state.current.value, ctx);
    state.current.str = NULL;
    
    /* First pass: identify parameters */
    for (slong i = 0; i < num_polys; i++) {
        fq_mvpoly_t temp;
        fq_mvpoly_init(&temp, num_elim_vars, state.max_pars, ctx);
        
        state.input = poly_strings[i];
        state.pos = 0;
        state.len = strlen(poly_strings[i]);
        next_token(&state);
        
        parse_expression(&state, &temp);
        fq_mvpoly_clear(&temp);
    }
    
    /* DEBUG: Print parameter order */
    printf("\n=== DEBUG: Parameter names after parsing ===\n");
    printf("Number of parameters: %ld\n", state.npars);
    for (slong i = 0; i < state.npars; i++) {
        printf("  Parameter[%ld]: %s\n", i, state.par_names[i]);
    }
    printf("=== END DEBUG ===\n\n");
    
    /* Parse polynomials */
    fq_mvpoly_t *polys = (fq_mvpoly_t*) malloc(num_polys * sizeof(fq_mvpoly_t));
    
    for (slong i = 0; i < num_polys; i++) {
        fq_mvpoly_init(&polys[i], num_elim_vars, state.npars, ctx);
        
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
    
    /* Build cancellation matrix */
    printf("\nStep 1: Build Cancellation Matrix\n");
    fq_mvpoly_t **M_mvpoly;
    build_fq_cancellation_matrix_mvpoly(&M_mvpoly, polys, num_elim_vars, state.npars);
    
    /* Perform row operations */
    printf("\nStep 2: Perform Matrix Row Operations\n");
    fq_mvpoly_t **modified_M_mvpoly;
    perform_fq_matrix_row_operations_mvpoly(&modified_M_mvpoly, &M_mvpoly, num_elim_vars, state.npars);
    
    /* Compute determinant of modified matrix */
    printf("\nStep 3: Compute determinant of modified matrix\n");
    fq_mvpoly_t d_poly;
    compute_fq_cancel_matrix_det(&d_poly, modified_M_mvpoly, num_elim_vars, state.npars, DET_METHOD_RECURSIVE);
    printf("Dixon polynomial has %ld terms\n", d_poly.nterms);
    
    /* Extract coefficient matrix */
    printf("\nStep 4: Extract coefficient matrix\n");
    fq_mvpoly_t **coeff_matrix = NULL;
    slong *row_indices = (slong*) flint_malloc(d_poly.nterms * sizeof(slong));
    slong *col_indices = (slong*) flint_malloc(d_poly.nterms * sizeof(slong));
    slong matrix_size = 0;
    
    extract_fq_coefficient_matrix_from_dixon(&coeff_matrix, row_indices, col_indices,
                                            &matrix_size, &d_poly, num_elim_vars, state.npars);
    
    /* Compute determinant with ideal reduction */
    fq_mvpoly_t result_poly;
    if (matrix_size > 0) {
        printf("\nStep 5: Compute determinant of coefficient matrix with ideal reduction\n");
        printf("Matrix size: %ld x %ld\n", matrix_size, matrix_size);
        
        /* 修复：确保变量名数组大小正确 */
        char **var_names_for_reduction = state.par_names;
        if (ideal) {
            /* 获取理想期望的变量数 */
            slong ideal_nvars = 0;
            if (ideal->is_prime_field && ideal->ctx.nmod_ctx) {
                ideal_nvars = nmod_mpoly_ctx_nvars(ideal->ctx.nmod_ctx);
            } else if (!ideal->is_prime_field && ideal->ctx.fq_ctx) {
                ideal_nvars = fq_nmod_mpoly_ctx_nvars(ideal->ctx.fq_ctx);
            }
            
            /* 如果理想期望更多变量，创建扩展的数组 */
            if (ideal_nvars > state.npars) {
                var_names_for_reduction = (char**) calloc(ideal_nvars, sizeof(char*));
                /* 复制现有的参数名 */
                for (slong i = 0; i < state.npars; i++) {
                    var_names_for_reduction[i] = state.par_names[i];
                }
                /* 其余位置保持为NULL */
            }
        }
        
        /* CRITICAL: Pass var_names_for_reduction to ensure variable names are available during reduction */
        compute_det_with_reduction_from_mvpoly(&result_poly, coeff_matrix, matrix_size, ideal, 
                                             var_names_for_reduction);
        
        /* 如果创建了临时数组，清理它（但不要释放字符串） */
        if (var_names_for_reduction != state.par_names) {
            free(var_names_for_reduction);
        }
        
        printf("Resultant has %ld terms\n", result_poly.nterms);
        
        /* Cleanup coefficient matrix */
        for (slong i = 0; i < matrix_size; i++) {
            for (slong j = 0; j < matrix_size; j++) {
                fq_mvpoly_clear(&coeff_matrix[i][j]);
            }
            flint_free(coeff_matrix[i]);
        }
        flint_free(coeff_matrix);
    } else {
        fq_mvpoly_init(&result_poly, 0, state.npars, ctx);
        printf("Warning: Empty coefficient matrix, resultant is 0\n");
    }
    
    /* Convert result to string */
    char *result = fq_mvpoly_to_string(&result_poly, state.par_names, state.generator_name);
    
    /* Print remaining variables */
    printf("Remaining variables: ");
    if (state.npars == 0) {
        printf("none");
    } else {
        for (slong i = 0; i < state.npars; i++) {
            if (i > 0) printf(", ");
            printf("%s", state.par_names[i]);
        }
    }
    printf("\n");
    
    /* Cleanup */
    fq_mvpoly_clear(&result_poly);
    flint_free(row_indices);
    flint_free(col_indices);
    
    for (slong i = 0; i <= num_elim_vars; i++) {
        for (slong j = 0; j <= num_elim_vars; j++) {
            fq_mvpoly_clear(&M_mvpoly[i][j]);
            fq_mvpoly_clear(&modified_M_mvpoly[i][j]);
        }
        flint_free(M_mvpoly[i]);
        flint_free(modified_M_mvpoly[i]);
    }
    flint_free(M_mvpoly);
    flint_free(modified_M_mvpoly);
    fq_mvpoly_clear(&d_poly);
    
    for (slong i = 0; i < num_polys; i++) {
        fq_mvpoly_clear(&polys[i]);
    }
    free(polys);
    
    for (slong i = 0; i < num_elim_vars; i++) {
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
    
    printf("\n=== Dixon Resultant Computation Complete ===\n");
    
    return result;
}

void construct_triangular_ideal_str(unified_triangular_ideal_t *ideal,
                                           const char *gens_string,   // 逗号分隔的生成元
                                           const char *vars_string,   // 逗号分隔的变量
                                           const fq_nmod_ctx_t ctx) {
    
    slong num_gens, num_vars;
    char **gens_array = split_string(gens_string, &num_gens);
    char **vars_array = split_string(vars_string, &num_vars);
    
    // 转换为const char**
    const char **ideal_gens = (const char**) malloc(num_gens * sizeof(char*));
    const char **var_names = (const char**) malloc(num_vars * sizeof(char*));
    
    for (slong i = 0; i < num_gens; i++) {
        ideal_gens[i] = gens_array[i];
    }
    for (slong i = 0; i < num_vars; i++) {
        var_names[i] = vars_array[i];
    }
    
    // 调用原始函数
    construct_triangular_ideal_from_strings(ideal, ideal_gens, num_gens, 
                                          var_names, num_vars, ctx);
    
    // 清理
    free(ideal_gens);
    free(var_names);
    free_split_strings(gens_array, num_gens);
    free_split_strings(vars_array, num_vars);
}

// 组合接口：Dixon with ideal reduction (字符串版本)
char* dixon_with_ideal_reduction_str(const char *poly_string,
                                     const char *elim_vars_string,
                                     const char *ideal_gens_string,
                                     const char *all_vars_string,
                                     const fq_nmod_ctx_t ctx) {
    
    // 分割所有输入
    slong num_polys, num_elim_vars, num_gens, num_all_vars;
    char **poly_array = split_string(poly_string, &num_polys);
    char **elim_array = split_string(elim_vars_string, &num_elim_vars);
    char **gens_array = split_string(ideal_gens_string, &num_gens);
    char **all_vars_array = split_string(all_vars_string, &num_all_vars);
    
    // 构造理想
    unified_triangular_ideal_t ideal;
    const char **ideal_gens = (const char**) malloc(num_gens * sizeof(char*));
    const char **all_var_names = (const char**) malloc(num_all_vars * sizeof(char*));
    
    for (slong i = 0; i < num_gens; i++) {
        ideal_gens[i] = gens_array[i];
    }
    for (slong i = 0; i < num_all_vars; i++) {
        all_var_names[i] = all_vars_array[i];
    }
    
    construct_triangular_ideal_from_strings(&ideal, ideal_gens, num_gens,
                                          all_var_names, num_all_vars, ctx);
    
    // 准备多项式和变量数组
    const char **poly_strings = (const char**) malloc(num_polys * sizeof(char*));
    const char **elim_vars = (const char**) malloc(num_elim_vars * sizeof(char*));
    
    for (slong i = 0; i < num_polys; i++) {
        poly_strings[i] = poly_array[i];
    }
    for (slong i = 0; i < num_elim_vars; i++) {
        elim_vars[i] = elim_array[i];
    }
    
    // 调用dixon with ideal reduction
    char *result = dixon_with_ideal_reduction(poly_strings, num_polys,
                                            elim_vars, num_elim_vars,
                                            ctx, &ideal);
    
    // 清理
    unified_triangular_ideal_clear(&ideal);
    free(poly_strings);
    free(elim_vars);
    free(ideal_gens);
    free(all_var_names);
    free_split_strings(poly_array, num_polys);
    free_split_strings(elim_array, num_elim_vars);
    free_split_strings(gens_array, num_gens);
    free_split_strings(all_vars_array, num_all_vars);
    
    return result;
}

#endif /* DIXON_WITH_IDEAL_REDUCTION_H */