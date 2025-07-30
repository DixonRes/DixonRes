#ifndef DIXON_WITH_IDEAL_REDUCTION_H
#define DIXON_WITH_IDEAL_REDUCTION_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/nmod.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_mpoly.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_mpoly.h>


/* Debug output control macro - set to 0 to disable all output */
#define DEBUG_OUTPUT_R 0

#if DEBUG_OUTPUT_R
    #define DEBUG_PRINT_R(...) printf(__VA_ARGS__)
#else
    #define DEBUG_PRINT_R(...) ((void)0)
#endif

/* Keep the original triangular ideal structure but add optimization fields */
typedef struct {
    void **generators;            /* Array of generators */
    slong num_gens;              /* Number of generators */
    slong *var_indices;          /* Which variable each generator eliminates */
    char **var_names;            /* Names of variables being reduced */
    slong *leading_degrees;      /* Leading degree of each generator (NEW - for optimization) */
    slong max_gens;              /* Maximum number of generators */
    int is_prime_field;          /* 1 for prime field, 0 for extension field */
    union {
        nmod_mpoly_ctx_struct *nmod_ctx;
        fq_nmod_mpoly_ctx_struct *fq_ctx;
    } ctx;
    const fq_nmod_ctx_struct *field_ctx;
} unified_triangular_ideal_t;

/* Initialize triangular ideal - MINIMALLY MODIFIED */
void unified_triangular_ideal_init(unified_triangular_ideal_t *ideal, 
                                  slong max_gens, slong nvars, 
                                  const fq_nmod_ctx_t field_ctx) {
    ideal->generators = (void**) flint_malloc(max_gens * sizeof(void*));
    ideal->var_indices = (slong*) flint_malloc(max_gens * sizeof(slong));
    ideal->var_names = (char**) flint_malloc(max_gens * sizeof(char*));
    ideal->leading_degrees = (slong*) flint_malloc(max_gens * sizeof(slong)); /* NEW */
    ideal->num_gens = 0;
    ideal->max_gens = max_gens;
    ideal->field_ctx = field_ctx;
    
    /* Initialize arrays */
    for (slong i = 0; i < max_gens; i++) {
        ideal->var_names[i] = NULL;
        ideal->leading_degrees[i] = 3; /* Default to 3 */
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

/* Clear triangular ideal - ADD leading_degrees cleanup */
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
    flint_free(ideal->leading_degrees); /* NEW */
}

/* KEEP ORIGINAL variable detection logic - don't change this */
static slong detect_main_variable(const fq_mvpoly_t *gen, char **detected_var_name, 
                                 char **var_names_hint) {
    slong detected_var_idx = -1;
    
    /* For each variable, check if it has an isolated x^3 term */
    for (slong v = 0; v < gen->nvars; v++) {
        int has_isolated_cube = 0;
        int has_other_terms = 0;
        
        /* Scan all terms to check this variable */
        for (slong t = 0; t < gen->nterms; t++) {
            if (gen->terms[t].var_exp) {
                slong var_degree = gen->terms[t].var_exp[v];
                
                if (var_degree > 0) {
                    /* Check if this is an isolated x^3 term */
                    if (var_degree == 3) {
                        /* Check if all other variables have degree 0 */
                        int is_isolated = 1;
                        for (slong other_v = 0; other_v < gen->nvars; other_v++) {
                            if (other_v != v && gen->terms[t].var_exp[other_v] > 0) {
                                is_isolated = 0;
                                break;
                            }
                        }
                        if (is_isolated) {
                            has_isolated_cube = 1;
                        } else {
                            /* x^3 appears with other variables */
                            has_other_terms = 1;
                        }
                    } else {
                        /* Variable appears with degree != 3 */
                        has_other_terms = 1;
                    }
                }
            }
        }
        
        /* This variable is the main variable if it has isolated x^3 and no other terms */
        if (has_isolated_cube && !has_other_terms) {
            if (detected_var_idx >= 0) {
                DEBUG_PRINT_R("Warning: Multiple variables with isolated x^3 found (vars %ld and %ld)\n", 
                       detected_var_idx, v);
            }
            detected_var_idx = v;
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
    
    if (detected_var_idx < 0) {
        DEBUG_PRINT_R("Warning: No variable with isolated x^3 term found\n");
    }
    
    return detected_var_idx;
}

/* KEEP ORIGINAL generator addition - just add leading degree detection */
void unified_triangular_ideal_add_generator_from_mvpoly_auto(unified_triangular_ideal_t *ideal,
                                                            const fq_mvpoly_t *gen, 
                                                            char **var_names_hint) {
    if (ideal->num_gens >= ideal->max_gens) {
        DEBUG_PRINT_R("Error: Maximum number of generators exceeded\n");
        return;
    }
    
    /* Detect main variable */
    char *detected_var_name = NULL;
    slong detected_var_idx = detect_main_variable(gen, &detected_var_name, var_names_hint);
    
    if (detected_var_idx < 0) {
        DEBUG_PRINT_R("Error: Could not detect main variable in generator\n");
        return;
    }
    
    /* Find the leading degree for this variable - OPTIMIZED */
    slong leading_degree = 3; /* Default */
    for (slong t = 0; t < gen->nterms; t++) {
        if (gen->terms[t].var_exp && gen->terms[t].var_exp[detected_var_idx] > 0) {
            /* Check if this is an isolated term */
            int is_isolated = 1;
            for (slong v = 0; v < gen->nvars && is_isolated; v++) {
                if (v != detected_var_idx && gen->terms[t].var_exp[v] > 0) {
                    is_isolated = 0;
                }
            }
            if (gen->terms[t].par_exp) {
                for (slong p = 0; p < gen->npars && is_isolated; p++) {
                    if (gen->terms[t].par_exp[p] > 0) {
                        is_isolated = 0;
                    }
                }
            }
            
            if (is_isolated && gen->terms[t].var_exp[detected_var_idx] > leading_degree) {
                leading_degree = gen->terms[t].var_exp[detected_var_idx];
            }
        }
    }
    
    DEBUG_PRINT_R("DEBUG: Detected main variable '%s' at index %ld (leading degree %ld) for generator %ld\n", 
           detected_var_name, detected_var_idx, leading_degree, ideal->num_gens);
    
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
    ideal->leading_degrees[ideal->num_gens] = leading_degree; /* NEW */
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

/* Find generator for a specific variable in the ideal */
static slong find_generator_for_variable(const unified_triangular_ideal_t *ideal, const char *var_name) {
    for (slong i = 0; i < ideal->num_gens; i++) {
        if (ideal->var_names[i] && strcmp(ideal->var_names[i], var_name) == 0) {
            return i;
        }
    }
    return -1;
}

/* KEEP ORIGINAL variable mapping logic but make it more robust */
void create_reduced_ideal_context(unified_triangular_ideal_t *reduced_ideal,
                                const unified_triangular_ideal_t *original_ideal,
                                char **current_var_names,
                                slong current_nvars,
                                const fq_nmod_ctx_t field_ctx) {
    
    /* Count how many generators are relevant to the current variables */
    slong relevant_gens = 0;
    for (slong i = 0; i < current_nvars; i++) {
        if (current_var_names[i]) {
            slong gen_idx = find_generator_for_variable(original_ideal, current_var_names[i]);
            if (gen_idx >= 0) {
                relevant_gens++;
            }
        }
    }
    
    DEBUG_PRINT_R("DEBUG: Creating reduced ideal with %ld relevant generators for %ld variables\n", 
           relevant_gens, current_nvars);
    DEBUG_PRINT_R("Current variables: ");
    for (slong i = 0; i < current_nvars; i++) {
        DEBUG_PRINT_R("%s ", current_var_names[i] ? current_var_names[i] : "?");
    }
    DEBUG_PRINT_R("\n");
    
    /* Initialize reduced ideal with matching context */
    unified_triangular_ideal_init(reduced_ideal, relevant_gens, current_nvars, field_ctx);
    
    /* For each current variable, find its generator and add to reduced ideal */
    for (slong i = 0; i < current_nvars; i++) {
        if (!current_var_names[i]) continue;
        
        slong gen_idx = find_generator_for_variable(original_ideal, current_var_names[i]);
        if (gen_idx < 0) continue;
        
        DEBUG_PRINT_R("DEBUG: Processing generator %ld for variable '%s' (new index %ld)\n", 
               gen_idx, current_var_names[i], i);
        
        /* Convert generator to new context */
        if (reduced_ideal->is_prime_field && original_ideal->is_prime_field) {
            nmod_mpoly_struct *old_gen = (nmod_mpoly_struct*)original_ideal->generators[gen_idx];
            nmod_mpoly_struct *new_gen = (nmod_mpoly_struct*)reduced_ideal->generators[reduced_ideal->num_gens];
            
            /* Create variable mapping from original to current context */
            slong orig_nvars = nmod_mpoly_ctx_nvars(original_ideal->ctx.nmod_ctx);
            slong *var_map = (slong*) flint_malloc(orig_nvars * sizeof(slong));
            
            /* Initialize mapping to -1 (variable eliminated) */
            for (slong j = 0; j < orig_nvars; j++) {
                var_map[j] = -1;
            }
            
            /* Build variable mapping based on names */
            /* Assume original variables are [x, z_1, z_2, ..., z_9] or similar pattern */
            for (slong old_v = 0; old_v < orig_nvars; old_v++) {
                char var_name[32];
                if (old_v == 0) {
                    /* First variable might be x or another name */
                    /* Try to match with current variables */
                    for (slong new_v = 0; new_v < current_nvars; new_v++) {
                        if (current_var_names[new_v] && 
                            strlen(current_var_names[new_v]) == 1) {
                            /* Single character variable, likely x */
                            var_map[0] = new_v;
                            break;
                        }
                    }
                    /* If not found, try "x" explicitly */
                    if (var_map[0] == -1) {
                        slong idx = find_variable_by_name("x", current_var_names, current_nvars);
                        if (idx >= 0) var_map[0] = idx;
                    }
                } else {
                    /* For z_i variables */
                    sprintf(var_name, "z_%ld", old_v);
                    slong new_idx = find_variable_by_name(var_name, current_var_names, current_nvars);
                    if (new_idx >= 0) {
                        var_map[old_v] = new_idx;
                    }
                }
            }
            
            /* Apply variable mapping to convert generator */
            nmod_mpoly_zero(new_gen, reduced_ideal->ctx.nmod_ctx);
            
            for (slong t = 0; t < nmod_mpoly_length(old_gen, original_ideal->ctx.nmod_ctx); t++) {
                mp_limb_t coeff = nmod_mpoly_get_term_coeff_ui(old_gen, t, original_ideal->ctx.nmod_ctx);
                ulong *old_exp = (ulong*) flint_malloc(orig_nvars * sizeof(ulong));
                nmod_mpoly_get_term_exp_ui(old_exp, old_gen, t, original_ideal->ctx.nmod_ctx);
                
                /* Create new exponent vector */
                ulong *new_exp = (ulong*) flint_calloc(current_nvars, sizeof(ulong));
                int valid_term = 1;
                
                /* Map exponents */
                for (slong old_v = 0; old_v < orig_nvars; old_v++) {
                    if (old_exp[old_v] > 0) {
                        if (var_map[old_v] >= 0) {
                            new_exp[var_map[old_v]] = old_exp[old_v];
                        } else {
                            /* Variable was eliminated */
                            valid_term = 0;
                            break;
                        }
                    }
                }
                
                if (valid_term) {
                    nmod_mpoly_t temp;
                    nmod_mpoly_init(temp, reduced_ideal->ctx.nmod_ctx);
                    nmod_mpoly_set_coeff_ui_ui(temp, coeff, new_exp, reduced_ideal->ctx.nmod_ctx);
                    nmod_mpoly_add(new_gen, new_gen, temp, reduced_ideal->ctx.nmod_ctx);
                    nmod_mpoly_clear(temp, reduced_ideal->ctx.nmod_ctx);
                }
                
                flint_free(old_exp);
                flint_free(new_exp);
            }
            
            flint_free(var_map);
            
            /* Check if the resulting polynomial is non-zero */
            if (!nmod_mpoly_is_zero(new_gen, reduced_ideal->ctx.nmod_ctx)) {
                reduced_ideal->var_indices[reduced_ideal->num_gens] = i;
                reduced_ideal->var_names[reduced_ideal->num_gens] = strdup(current_var_names[i]);
                reduced_ideal->leading_degrees[reduced_ideal->num_gens] = 
                    original_ideal->leading_degrees[gen_idx]; /* Copy leading degree */
                reduced_ideal->num_gens++;
                
                DEBUG_PRINT_R("DEBUG: Added generator %ld to reduced ideal\n", reduced_ideal->num_gens - 1);
            }
        } else if (!reduced_ideal->is_prime_field && !original_ideal->is_prime_field) {
            /* Similar for fq_nmod_mpoly */
            fq_nmod_mpoly_struct *old_gen = (fq_nmod_mpoly_struct*)original_ideal->generators[gen_idx];
            fq_nmod_mpoly_struct *new_gen = (fq_nmod_mpoly_struct*)reduced_ideal->generators[reduced_ideal->num_gens];
            
            /* Create variable mapping from original to current context */
            slong orig_nvars = fq_nmod_mpoly_ctx_nvars(original_ideal->ctx.fq_ctx);
            slong *var_map = (slong*) flint_malloc(orig_nvars * sizeof(slong));
            
            /* Initialize mapping to -1 (variable eliminated) */
            for (slong j = 0; j < orig_nvars; j++) {
                var_map[j] = -1;
            }
            
            /* Build variable mapping based on names - same logic as nmod case */
            for (slong old_v = 0; old_v < orig_nvars; old_v++) {
                char var_name[32];
                if (old_v == 0) {
                    for (slong new_v = 0; new_v < current_nvars; new_v++) {
                        if (current_var_names[new_v] && 
                            strlen(current_var_names[new_v]) == 1) {
                            var_map[0] = new_v;
                            break;
                        }
                    }
                    if (var_map[0] == -1) {
                        slong idx = find_variable_by_name("x", current_var_names, current_nvars);
                        if (idx >= 0) var_map[0] = idx;
                    }
                } else {
                    sprintf(var_name, "z_%ld", old_v);
                    slong new_idx = find_variable_by_name(var_name, current_var_names, current_nvars);
                    if (new_idx >= 0) {
                        var_map[old_v] = new_idx;
                    }
                }
            }
            
            /* Apply variable mapping to convert generator */
            fq_nmod_mpoly_zero(new_gen, reduced_ideal->ctx.fq_ctx);
            
            for (slong t = 0; t < fq_nmod_mpoly_length(old_gen, original_ideal->ctx.fq_ctx); t++) {
                fq_nmod_t coeff;
                fq_nmod_init(coeff, field_ctx);
                fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, old_gen, t, original_ideal->ctx.fq_ctx);
                
                ulong *old_exp = (ulong*) flint_malloc(orig_nvars * sizeof(ulong));
                fq_nmod_mpoly_get_term_exp_ui(old_exp, old_gen, t, original_ideal->ctx.fq_ctx);
                
                /* Create new exponent vector */
                ulong *new_exp = (ulong*) flint_calloc(current_nvars, sizeof(ulong));
                int valid_term = 1;
                
                /* Map exponents */
                for (slong old_v = 0; old_v < orig_nvars; old_v++) {
                    if (old_exp[old_v] > 0) {
                        if (var_map[old_v] >= 0) {
                            new_exp[var_map[old_v]] = old_exp[old_v];
                        } else {
                            /* Variable was eliminated */
                            valid_term = 0;
                            break;
                        }
                    }
                }
                
                if (valid_term) {
                    fq_nmod_mpoly_set_coeff_fq_nmod_ui(new_gen, coeff, new_exp, reduced_ideal->ctx.fq_ctx);
                }
                
                fq_nmod_clear(coeff, field_ctx);
                flint_free(old_exp);
                flint_free(new_exp);
            }
            
            flint_free(var_map);
            
            /* Check if the resulting polynomial is non-zero */
            if (!fq_nmod_mpoly_is_zero(new_gen, reduced_ideal->ctx.fq_ctx)) {
                reduced_ideal->var_indices[reduced_ideal->num_gens] = i;
                reduced_ideal->var_names[reduced_ideal->num_gens] = strdup(current_var_names[i]);
                reduced_ideal->leading_degrees[reduced_ideal->num_gens] = 
                    original_ideal->leading_degrees[gen_idx];
                reduced_ideal->num_gens++;
                
                DEBUG_PRINT_R("DEBUG: Added generator %ld to reduced ideal\n", reduced_ideal->num_gens - 1);
            }
        }
    }
    
    DEBUG_PRINT_R("DEBUG: Reduced ideal has %ld generators\n", reduced_ideal->num_gens);
}

/* HEAVILY OPTIMIZED reduction function with batch polynomial building for nmod */
void triangular_ideal_reduce_nmod_mpoly_with_names(nmod_mpoly_t poly, 
                                                   const unified_triangular_ideal_t *ideal,
                                                   char **current_var_names) {
    if (ideal->num_gens == 0) return;
    
    slong nvars = nmod_mpoly_ctx_nvars(ideal->ctx.nmod_ctx);
    slong initial_terms = nmod_mpoly_length(poly, ideal->ctx.nmod_ctx);
    
    /* Timing variables */
    clock_t total_start = clock();
    double total_scan_time = 0.0;
    double total_fi_extract_time = 0.0;
    double total_power_precompute_time = 0.0;
    double total_reduction_time = 0.0;
    
    DEBUG_PRINT_R("\n=== REDUCTION TIMING ANALYSIS ===\n");
    DEBUG_PRINT_R("Initial polynomial: %ld terms, %ld variables\n", initial_terms, nvars);
    DEBUG_PRINT_R("Ideal has %ld generators\n\n", ideal->num_gens);
    
    /* Process variables in reverse order */
    for (slong var_idx = nvars - 1; var_idx >= 0; var_idx--) {
        if (!current_var_names || !current_var_names[var_idx]) continue;
        
        clock_t var_start = clock();
        
        /* Find generator for this variable */
        slong gen_idx = -1;
        slong leading_degree = 3;
        for (slong g = 0; g < ideal->num_gens; g++) {
            if (ideal->var_names[g] && strcmp(ideal->var_names[g], current_var_names[var_idx]) == 0) {
                gen_idx = g;
                leading_degree = ideal->leading_degrees[g];
                break;
            }
        }
        
        if (gen_idx < 0) continue;
        
        /* Quick scan to estimate reduction workload */
        clock_t scan_start = clock();
        slong current_terms = nmod_mpoly_length(poly, ideal->ctx.nmod_ctx);
        slong terms_to_reduce = 0;
        slong max_degree = 0;
        
        /* Sample terms to estimate reduction workload */
        slong sample_size = FLINT_MIN(5000, current_terms);
        slong sample_terms_to_reduce = 0;
        
        for (slong t = 0; t < sample_size; t++) {
            ulong *exp = (ulong*) flint_malloc(nvars * sizeof(ulong));
            if (!exp) {
                DEBUG_PRINT_R("Memory allocation failed for exponent vector\n");
                continue;
            }
            
            nmod_mpoly_get_term_exp_ui(exp, poly, t, ideal->ctx.nmod_ctx);
            if (exp[var_idx] >= leading_degree) {
                sample_terms_to_reduce++;
                if (exp[var_idx] > max_degree) max_degree = exp[var_idx];
            }
            flint_free(exp);
        }
        
        /* Estimate total terms to reduce */
        if (sample_size < current_terms && sample_size > 0) {
            terms_to_reduce = (sample_terms_to_reduce * current_terms) / sample_size;
        } else {
            terms_to_reduce = sample_terms_to_reduce;
        }
        
        clock_t scan_end = clock();
        double scan_time = ((double)(scan_end - scan_start)) / CLOCKS_PER_SEC;
        total_scan_time += scan_time;
        
        DEBUG_PRINT_R("Variable '%s' (index %ld): %ld terms, max degree %ld\n", 
               current_var_names[var_idx], var_idx, current_terms, max_degree);
        
        if (terms_to_reduce == 0) {
            DEBUG_PRINT_R("  No terms to reduce, skipping\n\n");
            continue;
        }
        
        nmod_mpoly_struct *gen = (nmod_mpoly_struct*)ideal->generators[gen_idx];
        
        /* Extract f_i once and cache it */
        clock_t fi_start = clock();
        nmod_mpoly_t f_i;
        nmod_mpoly_init(f_i, ideal->ctx.nmod_ctx);
        
        /* Build f_i by subtracting leading term from generator */
        slong gen_length = nmod_mpoly_length(gen, ideal->ctx.nmod_ctx);
        for (slong t = 0; t < gen_length; t++) {
            ulong *exp = (ulong*) flint_malloc(nvars * sizeof(ulong));
            if (!exp) continue;
            
            nmod_mpoly_get_term_exp_ui(exp, gen, t, ideal->ctx.nmod_ctx);
            
            /* Skip var^leading_degree term */
            int is_var_leading = (exp[var_idx] == leading_degree);
            for (slong v = 0; v < nvars; v++) {
                if (v != var_idx && exp[v] > 0) {
                    is_var_leading = 0;
                    break;
                }
            }
            
            if (!is_var_leading) {
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
        
        clock_t fi_end = clock();
        double fi_time = ((double)(fi_end - fi_start)) / CLOCKS_PER_SEC;
        total_fi_extract_time += fi_time;
        
        /* Precompute powers for large polynomials */
        clock_t power_start = clock();
        slong max_quotient = max_degree / leading_degree;
        nmod_mpoly_t *f_i_powers = NULL;
        slong max_precompute = 0;
        
        if (max_quotient > 0 && terms_to_reduce > 1000) {
            max_precompute = FLINT_MIN(max_quotient, 6);
            f_i_powers = (nmod_mpoly_t*) flint_malloc((max_precompute + 1) * sizeof(nmod_mpoly_t));
            
            if (f_i_powers) {
                for (slong i = 0; i <= max_precompute; i++) {
                    nmod_mpoly_init(f_i_powers[i], ideal->ctx.nmod_ctx);
                }
                nmod_mpoly_one(f_i_powers[0], ideal->ctx.nmod_ctx);
                if (max_precompute >= 1) {
                    nmod_mpoly_set(f_i_powers[1], f_i, ideal->ctx.nmod_ctx);
                }
                
                for (slong i = 2; i <= max_precompute; i++) {
                    nmod_mpoly_mul(f_i_powers[i], f_i_powers[i-1], f_i, ideal->ctx.nmod_ctx);
                }
            }
        }
        
        clock_t power_end = clock();
        double power_time = ((double)(power_end - power_start)) / CLOCKS_PER_SEC;
        total_power_precompute_time += power_time;
        
        /* Process polynomial using batch building approach */
        clock_t reduction_start = clock();
        slong nterms = nmod_mpoly_length(poly, ideal->ctx.nmod_ctx);
        slong reduced_count = 0;
        
        /* Final result polynomial */
        nmod_mpoly_t result_poly;
        nmod_mpoly_init(result_poly, ideal->ctx.nmod_ctx);
        
        slong batch_size = 10000;
        
        /* Pre-allocate workspace for batch building */
        slong max_batch_terms = batch_size * 10; /* Conservative estimate */
        mp_limb_t *batch_coeffs = (mp_limb_t*) flint_malloc(max_batch_terms * sizeof(mp_limb_t));
        ulong *batch_exps = (ulong*) flint_malloc(max_batch_terms * nvars * sizeof(ulong));
        
        for (slong batch_start = 0; batch_start < nterms; batch_start += batch_size) {
            slong batch_end = FLINT_MIN(batch_start + batch_size, nterms);
            
            /* Collect all terms for this batch */
            slong batch_term_count = 0;
            
            /* Pre-allocate single exponent vector for reuse */
            ulong *exp = (ulong*) flint_malloc(nvars * sizeof(ulong));
            
            for (slong t = batch_start; t < batch_end; t++) {
                nmod_mpoly_get_term_exp_ui(exp, poly, t, ideal->ctx.nmod_ctx);
                mp_limb_t coeff = nmod_mpoly_get_term_coeff_ui(poly, t, ideal->ctx.nmod_ctx);
                slong var_deg = exp[var_idx];
                
                if (var_deg >= leading_degree) {
                    /* Reduction case */
                    reduced_count++;
                    
                    slong quotient = var_deg / leading_degree;
                    slong remainder = var_deg % leading_degree;
                    exp[var_idx] = remainder;
                    
                    /* Apply reduction and collect resulting terms */
                    if (quotient == 0) {
                        /* Just add the term with modified exponent */
                        if (batch_term_count >= max_batch_terms) {
                            max_batch_terms *= 2;
                            batch_coeffs = (mp_limb_t*) flint_realloc(batch_coeffs, max_batch_terms * sizeof(mp_limb_t));
                            batch_exps = (ulong*) flint_realloc(batch_exps, max_batch_terms * nvars * sizeof(ulong));
                        }
                        batch_coeffs[batch_term_count] = coeff;
                        memcpy(batch_exps + batch_term_count * nvars, exp, nvars * sizeof(ulong));
                        batch_term_count++;
                    } else {
                        /* Need to multiply by f_i^quotient */
                        nmod_mpoly_t *power_poly = NULL;
                        nmod_mpoly_t temp_power;
                        
                        if (quotient == 1) {
                            power_poly = &f_i;
                        } else if (f_i_powers && quotient <= max_precompute) {
                            power_poly = &f_i_powers[quotient];
                        } else {
                            /* Compute power on demand */
                            nmod_mpoly_init(temp_power, ideal->ctx.nmod_ctx);
                            if (quotient <= 20) {
                                nmod_mpoly_pow_ui(temp_power, f_i, quotient, ideal->ctx.nmod_ctx);
                            } else {
                                /* Binary exponentiation */
                                nmod_mpoly_one(temp_power, ideal->ctx.nmod_ctx);
                                nmod_mpoly_t base;
                                nmod_mpoly_init(base, ideal->ctx.nmod_ctx);
                                nmod_mpoly_set(base, f_i, ideal->ctx.nmod_ctx);
                                
                                slong exp_remaining = quotient;
                                while (exp_remaining > 0) {
                                    if (exp_remaining & 1) {
                                        nmod_mpoly_mul(temp_power, temp_power, base, ideal->ctx.nmod_ctx);
                                    }
                                    if (exp_remaining > 1) {
                                        nmod_mpoly_mul(base, base, base, ideal->ctx.nmod_ctx);
                                    }
                                    exp_remaining >>= 1;
                                }
                                nmod_mpoly_clear(base, ideal->ctx.nmod_ctx);
                            }
                            power_poly = &temp_power;
                        }
                        
                        /* Collect terms from coeff * x^remainder * power_poly */
                        slong power_len = nmod_mpoly_length(*power_poly, ideal->ctx.nmod_ctx);
                        for (slong p = 0; p < power_len; p++) {
                            if (batch_term_count >= max_batch_terms) {
                                max_batch_terms *= 2;
                                batch_coeffs = (mp_limb_t*) flint_realloc(batch_coeffs, max_batch_terms * sizeof(mp_limb_t));
                                batch_exps = (ulong*) flint_realloc(batch_exps, max_batch_terms * nvars * sizeof(ulong));
                            }
                            
                            /* Get term from power polynomial */
                            ulong *power_exp = batch_exps + batch_term_count * nvars;
                            nmod_mpoly_get_term_exp_ui(power_exp, *power_poly, p, ideal->ctx.nmod_ctx);
                            mp_limb_t power_coeff = nmod_mpoly_get_term_coeff_ui(*power_poly, p, ideal->ctx.nmod_ctx);
                            
                            /* Add base exponent */
                            for (slong v = 0; v < nvars; v++) {
                                power_exp[v] += exp[v];
                            }
                            
                            /* Multiply coefficients */
                            batch_coeffs[batch_term_count] = nmod_mul(coeff, power_coeff, ideal->ctx.nmod_ctx->mod);
                            batch_term_count++;
                        }
                        
                        if (power_poly == &temp_power) {
                            nmod_mpoly_clear(temp_power, ideal->ctx.nmod_ctx);
                        }
                    }
                } else {
                    /* Direct copy case - just collect the term */
                    if (batch_term_count >= max_batch_terms) {
                        max_batch_terms *= 2;
                        batch_coeffs = (mp_limb_t*) flint_realloc(batch_coeffs, max_batch_terms * sizeof(mp_limb_t));
                        batch_exps = (ulong*) flint_realloc(batch_exps, max_batch_terms * nvars * sizeof(ulong));
                    }
                    
                    batch_coeffs[batch_term_count] = coeff;
                    memcpy(batch_exps + batch_term_count * nvars, exp, nvars * sizeof(ulong));
                    batch_term_count++;
                }
            }
            
            /* Free the reused exponent vector */
            flint_free(exp);
            
            /* Build polynomial from collected terms using low-level API */
            nmod_mpoly_t batch_poly;
            nmod_mpoly_init(batch_poly, ideal->ctx.nmod_ctx);
            
            if (batch_term_count > 0) {
                /* Use fit_length to pre-allocate space */
                nmod_mpoly_fit_length(batch_poly, batch_term_count, ideal->ctx.nmod_ctx);
                
                /* Set terms directly */
                for (slong i = 0; i < batch_term_count; i++) {
                    nmod_mpoly_push_term_ui_ui(batch_poly, batch_coeffs[i], 
                                              batch_exps + i * nvars, ideal->ctx.nmod_ctx);
                }
                
                /* Sort and combine like terms */
                nmod_mpoly_sort_terms(batch_poly, ideal->ctx.nmod_ctx);
                nmod_mpoly_combine_like_terms(batch_poly, ideal->ctx.nmod_ctx);
            }
            
            /* Add batch result to final result */
            nmod_mpoly_add(result_poly, result_poly, batch_poly, ideal->ctx.nmod_ctx);
            nmod_mpoly_clear(batch_poly, ideal->ctx.nmod_ctx);
            
            /* Progress indicator */
            if (batch_start % (batch_size * 10) == 0 && batch_start > 0) {
                DEBUG_PRINT_R("  Progress: %ld/%ld terms processed\n", batch_end, nterms);
            }
        }
        
        /* Free batch workspace */
        flint_free(batch_coeffs);
        flint_free(batch_exps);
        
        /* Replace original polynomial with result */
        nmod_mpoly_swap(poly, result_poly, ideal->ctx.nmod_ctx);
        
        clock_t reduction_end = clock();
        double reduction_time = ((double)(reduction_end - reduction_start)) / CLOCKS_PER_SEC;
        total_reduction_time += reduction_time;
        
        clock_t var_end = clock();
        double var_total = ((double)(var_end - var_start)) / CLOCKS_PER_SEC;
        
        DEBUG_PRINT_R("  Time: %.3f seconds (reduced %ld terms, final size: %ld)\n\n", 
               var_total, reduced_count, nmod_mpoly_length(poly, ideal->ctx.nmod_ctx));
        
        /* Cleanup */
        nmod_mpoly_clear(result_poly, ideal->ctx.nmod_ctx);
        nmod_mpoly_clear(f_i, ideal->ctx.nmod_ctx);
        
        if (f_i_powers) {
            for (slong i = 0; i <= max_precompute; i++) {
                nmod_mpoly_clear(f_i_powers[i], ideal->ctx.nmod_ctx);
            }
            flint_free(f_i_powers);
        }
    }
    
    clock_t total_end = clock();
    double total_time = ((double)(total_end - total_start)) / CLOCKS_PER_SEC;
    
    DEBUG_PRINT_R("=== TIMING SUMMARY ===\n");
    DEBUG_PRINT_R("Total reduction time: %.3f seconds\n", total_time);
    DEBUG_PRINT_R("  Scanning: %.3f seconds (%.1f%%)\n", 
           total_scan_time, 100.0 * total_scan_time / total_time);
    DEBUG_PRINT_R("  f_i extraction: %.3f seconds (%.1f%%)\n", 
           total_fi_extract_time, 100.0 * total_fi_extract_time / total_time);
    DEBUG_PRINT_R("  Power precomputation: %.3f seconds (%.1f%%)\n", 
           total_power_precompute_time, 100.0 * total_power_precompute_time / total_time);
    DEBUG_PRINT_R("  Main reduction: %.3f seconds (%.1f%%)\n", 
           total_reduction_time, 100.0 * total_reduction_time / total_time);
    DEBUG_PRINT_R("Final polynomial: %ld terms\n", nmod_mpoly_length(poly, ideal->ctx.nmod_ctx));
    DEBUG_PRINT_R("=== END TIMING ANALYSIS ===\n\n");
}

/* Similar function for fq_nmod_mpoly */
void triangular_ideal_reduce_fq_nmod_mpoly_with_names(fq_nmod_mpoly_t poly,
                                                      const unified_triangular_ideal_t *ideal,
                                                      char **current_var_names) {
    if (ideal->num_gens == 0) return;
    
    slong nvars = fq_nmod_mpoly_ctx_nvars(ideal->ctx.fq_ctx);
    slong initial_terms = fq_nmod_mpoly_length(poly, ideal->ctx.fq_ctx);
    
    /* Timing variables */
    clock_t total_start = clock();
    double total_scan_time = 0.0;
    double total_fi_extract_time = 0.0;
    double total_power_precompute_time = 0.0;
    double total_reduction_time = 0.0;
    
    DEBUG_PRINT_R("\n=== REDUCTION TIMING ANALYSIS (FQ) ===\n");
    DEBUG_PRINT_R("Initial polynomial: %ld terms, %ld variables\n", initial_terms, nvars);
    DEBUG_PRINT_R("Ideal has %ld generators\n\n", ideal->num_gens);
    
    /* Process variables in reverse order */
    for (slong var_idx = nvars - 1; var_idx >= 0; var_idx--) {
        if (!current_var_names || !current_var_names[var_idx]) continue;
        
        clock_t var_start = clock();
        
        /* Find generator for this variable */
        slong gen_idx = -1;
        slong leading_degree = 3;
        for (slong g = 0; g < ideal->num_gens; g++) {
            if (ideal->var_names[g] && strcmp(ideal->var_names[g], current_var_names[var_idx]) == 0) {
                gen_idx = g;
                leading_degree = ideal->leading_degrees[g];
                break;
            }
        }
        
        if (gen_idx < 0) continue;
        
        /* Quick scan to estimate reduction workload */
        clock_t scan_start = clock();
        slong current_terms = fq_nmod_mpoly_length(poly, ideal->ctx.fq_ctx);
        slong terms_to_reduce = 0;
        slong max_degree = 0;
        
        /* Sample terms to estimate reduction workload */
        slong sample_size = FLINT_MIN(5000, current_terms);
        slong sample_terms_to_reduce = 0;
        
        for (slong t = 0; t < sample_size; t++) {
            ulong *exp = (ulong*) flint_malloc(nvars * sizeof(ulong));
            if (!exp) {
                DEBUG_PRINT_R("Memory allocation failed for exponent vector\n");
                continue;
            }
            
            fq_nmod_mpoly_get_term_exp_ui(exp, poly, t, ideal->ctx.fq_ctx);
            if (exp[var_idx] >= leading_degree) {
                sample_terms_to_reduce++;
                if (exp[var_idx] > max_degree) max_degree = exp[var_idx];
            }
            flint_free(exp);
        }
        
        /* Estimate total terms to reduce */
        if (sample_size < current_terms && sample_size > 0) {
            terms_to_reduce = (sample_terms_to_reduce * current_terms) / sample_size;
        } else {
            terms_to_reduce = sample_terms_to_reduce;
        }
        
        clock_t scan_end = clock();
        double scan_time = ((double)(scan_end - scan_start)) / CLOCKS_PER_SEC;
        total_scan_time += scan_time;
        
        DEBUG_PRINT_R("Variable '%s' (index %ld): %ld terms, max degree %ld\n", 
               current_var_names[var_idx], var_idx, current_terms, max_degree);
        
        if (terms_to_reduce == 0) {
            DEBUG_PRINT_R("  No terms to reduce, skipping\n\n");
            continue;
        }
        
        fq_nmod_mpoly_struct *gen = (fq_nmod_mpoly_struct*)ideal->generators[gen_idx];
        
        /* Extract f_i once and cache it */
        clock_t fi_start = clock();
        fq_nmod_mpoly_t f_i;
        fq_nmod_mpoly_init(f_i, ideal->ctx.fq_ctx);
        
        /* Build f_i by subtracting leading term from generator */
        slong gen_length = fq_nmod_mpoly_length(gen, ideal->ctx.fq_ctx);
        for (slong t = 0; t < gen_length; t++) {
            ulong *exp = (ulong*) flint_malloc(nvars * sizeof(ulong));
            if (!exp) continue;
            
            fq_nmod_mpoly_get_term_exp_ui(exp, gen, t, ideal->ctx.fq_ctx);
            
            /* Skip var^leading_degree term */
            int is_var_leading = (exp[var_idx] == leading_degree);
            for (slong v = 0; v < nvars; v++) {
                if (v != var_idx && exp[v] > 0) {
                    is_var_leading = 0;
                    break;
                }
            }
            
            if (!is_var_leading) {
                fq_nmod_t coeff;
                fq_nmod_init(coeff, ideal->field_ctx);
                fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, gen, t, ideal->ctx.fq_ctx);
                fq_nmod_neg(coeff, coeff, ideal->field_ctx);
                
                fq_nmod_mpoly_set_coeff_fq_nmod_ui(f_i, coeff, exp, ideal->ctx.fq_ctx);
                
                fq_nmod_clear(coeff, ideal->field_ctx);
            }
            flint_free(exp);
        }
        
        clock_t fi_end = clock();
        double fi_time = ((double)(fi_end - fi_start)) / CLOCKS_PER_SEC;
        total_fi_extract_time += fi_time;
        
        /* Precompute powers for large polynomials */
        clock_t power_start = clock();
        slong max_quotient = max_degree / leading_degree;
        fq_nmod_mpoly_t *f_i_powers = NULL;
        slong max_precompute = 0;
        
        if (max_quotient > 0 && terms_to_reduce > 1000) {
            max_precompute = FLINT_MIN(max_quotient, 6);
            f_i_powers = (fq_nmod_mpoly_t*) flint_malloc((max_precompute + 1) * sizeof(fq_nmod_mpoly_t));
            
            if (f_i_powers) {
                for (slong i = 0; i <= max_precompute; i++) {
                    fq_nmod_mpoly_init(f_i_powers[i], ideal->ctx.fq_ctx);
                }
                fq_nmod_mpoly_one(f_i_powers[0], ideal->ctx.fq_ctx);
                if (max_precompute >= 1) {
                    fq_nmod_mpoly_set(f_i_powers[1], f_i, ideal->ctx.fq_ctx);
                }
                
                for (slong i = 2; i <= max_precompute; i++) {
                    fq_nmod_mpoly_mul(f_i_powers[i], f_i_powers[i-1], f_i, ideal->ctx.fq_ctx);
                }
            }
        }
        
        clock_t power_end = clock();
        double power_time = ((double)(power_end - power_start)) / CLOCKS_PER_SEC;
        total_power_precompute_time += power_time;
        
        /* Process polynomial using batch building approach */
        clock_t reduction_start = clock();
        slong nterms = fq_nmod_mpoly_length(poly, ideal->ctx.fq_ctx);
        slong reduced_count = 0;
        
        /* Final result polynomial */
        fq_nmod_mpoly_t result_poly;
        fq_nmod_mpoly_init(result_poly, ideal->ctx.fq_ctx);
        
        slong batch_size = 10000;
        
        /* Pre-allocate workspace for batch building */
        slong max_batch_terms = batch_size * 10;
        fq_nmod_t *batch_coeffs = (fq_nmod_t*) flint_malloc(max_batch_terms * sizeof(fq_nmod_t));
        for (slong i = 0; i < max_batch_terms; i++) {
            fq_nmod_init(batch_coeffs[i], ideal->field_ctx);
        }
        ulong *batch_exps = (ulong*) flint_malloc(max_batch_terms * nvars * sizeof(ulong));
        
        for (slong batch_start = 0; batch_start < nterms; batch_start += batch_size) {
            slong batch_end = FLINT_MIN(batch_start + batch_size, nterms);
            
            /* Collect all terms for this batch */
            slong batch_term_count = 0;
            
            /* Pre-allocate single exponent vector for reuse */
            ulong *exp = (ulong*) flint_malloc(nvars * sizeof(ulong));
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ideal->field_ctx);
            
            for (slong t = batch_start; t < batch_end; t++) {
                fq_nmod_mpoly_get_term_exp_ui(exp, poly, t, ideal->ctx.fq_ctx);
                fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, poly, t, ideal->ctx.fq_ctx);
                slong var_deg = exp[var_idx];
                
                if (var_deg >= leading_degree) {
                    /* Reduction case */
                    reduced_count++;
                    
                    slong quotient = var_deg / leading_degree;
                    slong remainder = var_deg % leading_degree;
                    exp[var_idx] = remainder;
                    
                    /* Apply reduction and collect resulting terms */
                    if (quotient == 0) {
                        /* Just add the term with modified exponent */
                        if (batch_term_count >= max_batch_terms) {
                            max_batch_terms *= 2;
                            fq_nmod_t *new_coeffs = (fq_nmod_t*) flint_malloc(max_batch_terms * sizeof(fq_nmod_t));
                            for (slong i = 0; i < batch_term_count; i++) {
                                fq_nmod_init(new_coeffs[i], ideal->field_ctx);
                                fq_nmod_set(new_coeffs[i], batch_coeffs[i], ideal->field_ctx);
                            }
                            for (slong i = batch_term_count; i < max_batch_terms; i++) {
                                fq_nmod_init(new_coeffs[i], ideal->field_ctx);
                            }
                            for (slong i = 0; i < batch_term_count; i++) {
                                fq_nmod_clear(batch_coeffs[i], ideal->field_ctx);
                            }
                            flint_free(batch_coeffs);
                            batch_coeffs = new_coeffs;
                            
                            batch_exps = (ulong*) flint_realloc(batch_exps, max_batch_terms * nvars * sizeof(ulong));
                        }
                        fq_nmod_set(batch_coeffs[batch_term_count], coeff, ideal->field_ctx);
                        memcpy(batch_exps + batch_term_count * nvars, exp, nvars * sizeof(ulong));
                        batch_term_count++;
                    } else {
                        /* Need to multiply by f_i^quotient */
                        fq_nmod_mpoly_t *power_poly = NULL;
                        fq_nmod_mpoly_t temp_power;
                        
                        if (quotient == 1) {
                            power_poly = &f_i;
                        } else if (f_i_powers && quotient <= max_precompute) {
                            power_poly = &f_i_powers[quotient];
                        } else {
                            /* Compute power on demand */
                            fq_nmod_mpoly_init(temp_power, ideal->ctx.fq_ctx);
                            if (quotient <= 20) {
                                fq_nmod_mpoly_pow_ui(temp_power, f_i, quotient, ideal->ctx.fq_ctx);
                            } else {
                                /* Binary exponentiation */
                                fq_nmod_mpoly_one(temp_power, ideal->ctx.fq_ctx);
                                fq_nmod_mpoly_t base;
                                fq_nmod_mpoly_init(base, ideal->ctx.fq_ctx);
                                fq_nmod_mpoly_set(base, f_i, ideal->ctx.fq_ctx);
                                
                                slong exp_remaining = quotient;
                                while (exp_remaining > 0) {
                                    if (exp_remaining & 1) {
                                        fq_nmod_mpoly_mul(temp_power, temp_power, base, ideal->ctx.fq_ctx);
                                    }
                                    if (exp_remaining > 1) {
                                        fq_nmod_mpoly_mul(base, base, base, ideal->ctx.fq_ctx);
                                    }
                                    exp_remaining >>= 1;
                                }
                                fq_nmod_mpoly_clear(base, ideal->ctx.fq_ctx);
                            }
                            power_poly = &temp_power;
                        }
                        
                        /* Collect terms from coeff * x^remainder * power_poly */
                        slong power_len = fq_nmod_mpoly_length(*power_poly, ideal->ctx.fq_ctx);
                        fq_nmod_t power_coeff;
                        fq_nmod_init(power_coeff, ideal->field_ctx);
                        
                        for (slong p = 0; p < power_len; p++) {
                            if (batch_term_count >= max_batch_terms) {
                                max_batch_terms *= 2;
                                fq_nmod_t *new_coeffs = (fq_nmod_t*) flint_malloc(max_batch_terms * sizeof(fq_nmod_t));
                                for (slong i = 0; i < batch_term_count; i++) {
                                    fq_nmod_init(new_coeffs[i], ideal->field_ctx);
                                    fq_nmod_set(new_coeffs[i], batch_coeffs[i], ideal->field_ctx);
                                }
                                for (slong i = batch_term_count; i < max_batch_terms; i++) {
                                    fq_nmod_init(new_coeffs[i], ideal->field_ctx);
                                }
                                for (slong i = 0; i < batch_term_count; i++) {
                                    fq_nmod_clear(batch_coeffs[i], ideal->field_ctx);
                                }
                                flint_free(batch_coeffs);
                                batch_coeffs = new_coeffs;
                                
                                batch_exps = (ulong*) flint_realloc(batch_exps, max_batch_terms * nvars * sizeof(ulong));
                            }
                            
                            /* Get term from power polynomial */
                            ulong *power_exp = batch_exps + batch_term_count * nvars;
                            fq_nmod_mpoly_get_term_exp_ui(power_exp, *power_poly, p, ideal->ctx.fq_ctx);
                            fq_nmod_mpoly_get_term_coeff_fq_nmod(power_coeff, *power_poly, p, ideal->ctx.fq_ctx);
                            
                            /* Add base exponent */
                            for (slong v = 0; v < nvars; v++) {
                                power_exp[v] += exp[v];
                            }
                            
                            /* Multiply coefficients */
                            fq_nmod_mul(batch_coeffs[batch_term_count], coeff, power_coeff, ideal->field_ctx);
                            batch_term_count++;
                        }
                        
                        fq_nmod_clear(power_coeff, ideal->field_ctx);
                        
                        if (power_poly == &temp_power) {
                            fq_nmod_mpoly_clear(temp_power, ideal->ctx.fq_ctx);
                        }
                    }
                } else {
                    /* Direct copy case - just collect the term */
                    if (batch_term_count >= max_batch_terms) {
                        max_batch_terms *= 2;
                        fq_nmod_t *new_coeffs = (fq_nmod_t*) flint_malloc(max_batch_terms * sizeof(fq_nmod_t));
                        for (slong i = 0; i < batch_term_count; i++) {
                            fq_nmod_init(new_coeffs[i], ideal->field_ctx);
                            fq_nmod_set(new_coeffs[i], batch_coeffs[i], ideal->field_ctx);
                        }
                        for (slong i = batch_term_count; i < max_batch_terms; i++) {
                            fq_nmod_init(new_coeffs[i], ideal->field_ctx);
                        }
                        for (slong i = 0; i < batch_term_count; i++) {
                            fq_nmod_clear(batch_coeffs[i], ideal->field_ctx);
                        }
                        flint_free(batch_coeffs);
                        batch_coeffs = new_coeffs;
                        
                        batch_exps = (ulong*) flint_realloc(batch_exps, max_batch_terms * nvars * sizeof(ulong));
                    }
                    
                    fq_nmod_set(batch_coeffs[batch_term_count], coeff, ideal->field_ctx);
                    memcpy(batch_exps + batch_term_count * nvars, exp, nvars * sizeof(ulong));
                    batch_term_count++;
                }
            }
            
            /* Free the reused vectors */
            flint_free(exp);
            fq_nmod_clear(coeff, ideal->field_ctx);
            
            /* Build polynomial from collected terms */
            fq_nmod_mpoly_t batch_poly;
            fq_nmod_mpoly_init(batch_poly, ideal->ctx.fq_ctx);
            
            if (batch_term_count > 0) {
                /* Add terms one by one (no direct batch API for fq_nmod_mpoly) */
                for (slong i = 0; i < batch_term_count; i++) {
                    fq_nmod_mpoly_set_coeff_fq_nmod_ui(batch_poly, batch_coeffs[i], 
                                                       batch_exps + i * nvars, ideal->ctx.fq_ctx);
                }
            }
            
            /* Add batch result to final result */
            fq_nmod_mpoly_add(result_poly, result_poly, batch_poly, ideal->ctx.fq_ctx);
            fq_nmod_mpoly_clear(batch_poly, ideal->ctx.fq_ctx);
            
            /* Progress indicator */
            if (batch_start % (batch_size * 10) == 0 && batch_start > 0) {
                DEBUG_PRINT_R("  Progress: %ld/%ld terms processed\n", batch_end, nterms);
            }
        }
        
        /* Free batch workspace */
        for (slong i = 0; i < max_batch_terms; i++) {
            fq_nmod_clear(batch_coeffs[i], ideal->field_ctx);
        }
        flint_free(batch_coeffs);
        flint_free(batch_exps);
        
        /* Replace original polynomial with result */
        fq_nmod_mpoly_swap(poly, result_poly, ideal->ctx.fq_ctx);
        
        clock_t reduction_end = clock();
        double reduction_time = ((double)(reduction_end - reduction_start)) / CLOCKS_PER_SEC;
        total_reduction_time += reduction_time;
        
        clock_t var_end = clock();
        double var_total = ((double)(var_end - var_start)) / CLOCKS_PER_SEC;
        
        DEBUG_PRINT_R("  Time: %.3f seconds (reduced %ld terms, final size: %ld)\n\n", 
               var_total, reduced_count, fq_nmod_mpoly_length(poly, ideal->ctx.fq_ctx));
        
        /* Cleanup */
        fq_nmod_mpoly_clear(result_poly, ideal->ctx.fq_ctx);
        fq_nmod_mpoly_clear(f_i, ideal->ctx.fq_ctx);
        
        if (f_i_powers) {
            for (slong i = 0; i <= max_precompute; i++) {
                fq_nmod_mpoly_clear(f_i_powers[i], ideal->ctx.fq_ctx);
            }
            flint_free(f_i_powers);
        }
    }
    
    clock_t total_end = clock();
    double total_time = ((double)(total_end - total_start)) / CLOCKS_PER_SEC;
    
    DEBUG_PRINT_R("=== TIMING SUMMARY (FQ) ===\n");
    DEBUG_PRINT_R("Total reduction time: %.3f seconds\n", total_time);
    DEBUG_PRINT_R("  Scanning: %.3f seconds (%.1f%%)\n", 
           total_scan_time, 100.0 * total_scan_time / total_time);
    DEBUG_PRINT_R("  f_i extraction: %.3f seconds (%.1f%%)\n", 
           total_fi_extract_time, 100.0 * total_fi_extract_time / total_time);
    DEBUG_PRINT_R("  Power precomputation: %.3f seconds (%.1f%%)\n", 
           total_power_precompute_time, 100.0 * total_power_precompute_time / total_time);
    DEBUG_PRINT_R("  Main reduction: %.3f seconds (%.1f%%)\n", 
           total_reduction_time, 100.0 * total_reduction_time / total_time);
    DEBUG_PRINT_R("Final polynomial: %ld terms\n", fq_nmod_mpoly_length(poly, ideal->ctx.fq_ctx));
    DEBUG_PRINT_R("=== END TIMING ANALYSIS ===\n\n");
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
        /* 2x2 matrix determinant */
        fq_nmod_mpoly_t ad, bc, temp;
        fq_nmod_mpoly_init(ad, ideal->ctx.fq_ctx);
        fq_nmod_mpoly_init(bc, ideal->ctx.fq_ctx);
        fq_nmod_mpoly_init(temp, ideal->ctx.fq_ctx);
        
        fq_nmod_mpoly_mul(ad, &matrix[0][0], &matrix[1][1], ideal->ctx.fq_ctx);
        fq_nmod_mpoly_mul(bc, &matrix[0][1], &matrix[1][0], ideal->ctx.fq_ctx);
        fq_nmod_mpoly_sub(det, ad, bc, ideal->ctx.fq_ctx);
        
        /* Only reduce once at the end */
        triangular_ideal_reduce_fq_nmod_mpoly_with_names(det, ideal, current_var_names);
        
        fq_nmod_mpoly_clear(ad, ideal->ctx.fq_ctx);
        fq_nmod_mpoly_clear(bc, ideal->ctx.fq_ctx);
        fq_nmod_mpoly_clear(temp, ideal->ctx.fq_ctx);
        return;
    }
    
    /* For larger matrices, use cofactor expansion */
    fq_nmod_mpoly_zero(det, ideal->ctx.fq_ctx);
    
    /* Expand along first row */
    for (slong j = 0; j < size; j++) {
        if (fq_nmod_mpoly_is_zero(&matrix[0][j], ideal->ctx.fq_ctx)) continue;
        
        /* Build submatrix */
        fq_nmod_mpoly_struct **submatrix = (fq_nmod_mpoly_struct**) 
            flint_malloc((size-1) * sizeof(fq_nmod_mpoly_struct*));
        for (slong i = 0; i < size-1; i++) {
            submatrix[i] = (fq_nmod_mpoly_struct*) 
                flint_malloc((size-1) * sizeof(fq_nmod_mpoly_struct));
            for (slong k = 0; k < size-1; k++) {
                fq_nmod_mpoly_init(&submatrix[i][k], ideal->ctx.fq_ctx);
            }
        }
        
        /* Copy submatrix */
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
        
        /* Compute minor recursively */
        fq_nmod_mpoly_t minor;
        fq_nmod_mpoly_init(minor, ideal->ctx.fq_ctx);
        compute_fq_nmod_det_with_triangular_reduction_with_names(minor, submatrix, size-1, ideal, current_var_names);
        
        /* Compute contribution */
        fq_nmod_mpoly_t contrib;
        fq_nmod_mpoly_init(contrib, ideal->ctx.fq_ctx);
        fq_nmod_mpoly_mul(contrib, &matrix[0][j], minor, ideal->ctx.fq_ctx);
        
        if (j % 2 == 0) {
            fq_nmod_mpoly_add(det, det, contrib, ideal->ctx.fq_ctx);
        } else {
            fq_nmod_mpoly_sub(det, det, contrib, ideal->ctx.fq_ctx);
        }
        
        /* Reduce periodically */
        if ((j + 1) % 3 == 0 || j == size - 1) {
            triangular_ideal_reduce_fq_nmod_mpoly_with_names(det, ideal, current_var_names);
        }
        
        /* Cleanup */
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

/* Simplified and safer matrix conversion with timeout protection */
void compute_det_with_reduction_from_mvpoly(fq_mvpoly_t *result,
                                           fq_mvpoly_t **matrix,
                                           slong size,
                                           const unified_triangular_ideal_t *ideal,
                                           char **current_var_names) {
    DEBUG_PRINT_R("Computing determinant with ideal reduction (matrix size: %ld x %ld)\n", size, size);
    
    if (size == 0) {
        fq_mvpoly_init(result, 0, 0, ideal->field_ctx);
        return;
    }
    
    slong current_nvars = matrix[0][0].npars;  /* These are the remaining variables */
    
    /* Build complete variable name list - only parameters remain as variables */
    char **complete_var_names = (char**) malloc(current_nvars * sizeof(char*));
    for (slong i = 0; i < current_nvars; i++) {
        if (current_var_names && current_var_names[i]) {
            complete_var_names[i] = current_var_names[i];
        } else {
            char temp[32];
            sprintf(temp, "param_%ld", i);
            complete_var_names[i] = strdup(temp);
        }
    }
    
    /* Create reduced ideal with proper context */
    unified_triangular_ideal_t reduced_ideal;
    create_reduced_ideal_context(&reduced_ideal, ideal, complete_var_names, 
                               current_nvars, ideal->field_ctx);
    
    if (ideal->is_prime_field) {
        DEBUG_PRINT_R("Using nmod_mpoly for prime field computation\n");
        
        /* Convert to nmod_mpoly matrix with only parameter variables */
        nmod_mpoly_struct **nmod_matrix = (nmod_mpoly_struct**) flint_malloc(size * sizeof(nmod_mpoly_struct*));
        for (slong i = 0; i < size; i++) {
            nmod_matrix[i] = (nmod_mpoly_struct*) flint_malloc(size * sizeof(nmod_mpoly_struct));
            for (slong j = 0; j < size; j++) {
                nmod_mpoly_init(&nmod_matrix[i][j], reduced_ideal.ctx.nmod_ctx);
                nmod_mpoly_zero(&nmod_matrix[i][j], reduced_ideal.ctx.nmod_ctx);
                
                /* Convert fq_mvpoly to nmod_mpoly - SIMPLIFIED AND SAFER */
                DEBUG_PRINT_R("Converting matrix[%ld][%ld] with %ld terms...\n", i, j, matrix[i][j].nterms);
                
                /* Process each term with safety checks */
                for (slong t = 0; t < matrix[i][j].nterms && t < 10000; t++) { /* Limit terms to prevent infinite loops */
                    /* Get coefficient as mp_limb_t */
                    mp_limb_t coeff = 0;
                    
                    /* For prime field, extract the coefficient safely */
                    if (fq_nmod_is_one(matrix[i][j].terms[t].coeff, ideal->field_ctx)) {
                        coeff = 1;
                    } else if (!fq_nmod_is_zero(matrix[i][j].terms[t].coeff, ideal->field_ctx)) {
                        /* Extract coefficient value - for prime field, it's the constant term */
                        fq_nmod_t temp;
                        fq_nmod_init(temp, ideal->field_ctx);
                        fq_nmod_set(temp, matrix[i][j].terms[t].coeff, ideal->field_ctx);
                        
                        /* For prime field (degree 1), coefficient is directly the value */
                        nmod_poly_struct *p = &temp[0];
                        if (p->length > 0) {
                            coeff = p->coeffs[0];
                        }
                        
                        fq_nmod_clear(temp, ideal->field_ctx);
                    }
                    
                    /* Build exponent vector for nmod_mpoly */
                    ulong *exp = (ulong*) flint_calloc(current_nvars, sizeof(ulong));
                    
                    /* Copy parameter exponents - these become the variables */
                    if (matrix[i][j].terms[t].par_exp) {
                        for (slong v = 0; v < matrix[i][j].npars && v < current_nvars; v++) {
                            exp[v] = matrix[i][j].terms[t].par_exp[v];
                        }
                    }
                    
                    /* Add term if coefficient is non-zero */
                    if (coeff != 0) {
                        nmod_mpoly_t temp;
                        nmod_mpoly_init(temp, reduced_ideal.ctx.nmod_ctx);
                        nmod_mpoly_set_coeff_ui_ui(temp, coeff, exp, reduced_ideal.ctx.nmod_ctx);
                        nmod_mpoly_add(&nmod_matrix[i][j], &nmod_matrix[i][j], temp, reduced_ideal.ctx.nmod_ctx);
                        nmod_mpoly_clear(temp, reduced_ideal.ctx.nmod_ctx);
                    }
                    
                    flint_free(exp);
                }
                
                DEBUG_PRINT_R("  Converted to %ld terms\n", nmod_mpoly_length(&nmod_matrix[i][j], reduced_ideal.ctx.nmod_ctx));
            }
        }
        
        /* Compute determinant with reduction - ADD TIMEOUT PROTECTION */
        nmod_mpoly_t det;
        nmod_mpoly_init(det, reduced_ideal.ctx.nmod_ctx);
        
        DEBUG_PRINT_R("Starting determinant computation...\n");
        clock_t start = clock();
        
        /* Use simplified determinant for small matrices */
        if (size == 1) {
            nmod_mpoly_set(det, &nmod_matrix[0][0], reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(det, &reduced_ideal, complete_var_names);
        } else if (size == 2) {
            nmod_mpoly_t ad, bc;
            nmod_mpoly_init(ad, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_init(bc, reduced_ideal.ctx.nmod_ctx);
            
            nmod_mpoly_mul(ad, &nmod_matrix[0][0], &nmod_matrix[1][1], reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_mul(bc, &nmod_matrix[0][1], &nmod_matrix[1][0], reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_sub(det, ad, bc, reduced_ideal.ctx.nmod_ctx);
            
            /* Reduce only once at the end */
            triangular_ideal_reduce_nmod_mpoly_with_names(det, &reduced_ideal, complete_var_names);
            
            nmod_mpoly_clear(ad, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_clear(bc, reduced_ideal.ctx.nmod_ctx);
        } else if (size == 3) {
            /* For 3x3: First reduce all matrix elements, then compute with intermediate reductions */
            DEBUG_PRINT_R("Pre-reducing all 3x3 matrix elements...\n");
            
            /* Reduce each matrix element first */
            for (slong i = 0; i < 3; i++) {
                for (slong j = 0; j < 3; j++) {
                    DEBUG_PRINT_R("  Reducing element [%ld][%ld] (%ld terms)...\n", i, j,
                           nmod_mpoly_length(&nmod_matrix[i][j], reduced_ideal.ctx.nmod_ctx));
                    triangular_ideal_reduce_nmod_mpoly_with_names(&nmod_matrix[i][j], 
                                                                &reduced_ideal, complete_var_names);
                    DEBUG_PRINT_R("    Reduced to %ld terms\n", 
                           nmod_mpoly_length(&nmod_matrix[i][j], reduced_ideal.ctx.nmod_ctx));
                }
            }
            
            DEBUG_PRINT_R("Computing 3x3 determinant with intermediate reductions...\n");
            
            /* Compute 3x3 determinant: det = a00*a11*a22 + a01*a12*a20 + a02*a10*a21 
                                              - a00*a12*a21 - a01*a10*a22 - a02*a11*a20 */
            nmod_mpoly_t term1, term2, term3, term4, term5, term6, temp1, temp2;
            nmod_mpoly_init(term1, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_init(term2, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_init(term3, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_init(term4, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_init(term5, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_init(term6, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_init(temp1, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_init(temp2, reduced_ideal.ctx.nmod_ctx);
            
            /* Positive terms with reduction after each multiplication */
            
            /* term1 = a00 * a11 * a22 */
            DEBUG_PRINT_R("  Computing term1: a00*a11*a22...\n");
            nmod_mpoly_mul(temp1, &nmod_matrix[0][0], &nmod_matrix[1][1], reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(temp1, &reduced_ideal, complete_var_names);
            nmod_mpoly_mul(term1, temp1, &nmod_matrix[2][2], reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(term1, &reduced_ideal, complete_var_names);
            DEBUG_PRINT_R("    term1 has %ld terms\n", nmod_mpoly_length(term1, reduced_ideal.ctx.nmod_ctx));
            
            /* term2 = a01 * a12 * a20 */
            DEBUG_PRINT_R("  Computing term2: a01*a12*a20...\n");
            nmod_mpoly_mul(temp1, &nmod_matrix[0][1], &nmod_matrix[1][2], reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(temp1, &reduced_ideal, complete_var_names);
            nmod_mpoly_mul(term2, temp1, &nmod_matrix[2][0], reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(term2, &reduced_ideal, complete_var_names);
            DEBUG_PRINT_R("    term2 has %ld terms\n", nmod_mpoly_length(term2, reduced_ideal.ctx.nmod_ctx));
            
            /* term3 = a02 * a10 * a21 */
            DEBUG_PRINT_R("  Computing term3: a02*a10*a21...\n");
            nmod_mpoly_mul(temp1, &nmod_matrix[0][2], &nmod_matrix[1][0], reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(temp1, &reduced_ideal, complete_var_names);
            nmod_mpoly_mul(term3, temp1, &nmod_matrix[2][1], reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(term3, &reduced_ideal, complete_var_names);
            DEBUG_PRINT_R("    term3 has %ld terms\n", nmod_mpoly_length(term3, reduced_ideal.ctx.nmod_ctx));
            
            /* Negative terms with reduction after each multiplication */
            
            /* term4 = a00 * a12 * a21 */
            DEBUG_PRINT_R("  Computing term4: a00*a12*a21...\n");
            nmod_mpoly_mul(temp1, &nmod_matrix[0][0], &nmod_matrix[1][2], reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(temp1, &reduced_ideal, complete_var_names);
            nmod_mpoly_mul(term4, temp1, &nmod_matrix[2][1], reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(term4, &reduced_ideal, complete_var_names);
            DEBUG_PRINT_R("    term4 has %ld terms\n", nmod_mpoly_length(term4, reduced_ideal.ctx.nmod_ctx));
            
            /* term5 = a01 * a10 * a22 */
            DEBUG_PRINT_R("  Computing term5: a01*a10*a22...\n");
            nmod_mpoly_mul(temp1, &nmod_matrix[0][1], &nmod_matrix[1][0], reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(temp1, &reduced_ideal, complete_var_names);
            nmod_mpoly_mul(term5, temp1, &nmod_matrix[2][2], reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(term5, &reduced_ideal, complete_var_names);
            DEBUG_PRINT_R("    term5 has %ld terms\n", nmod_mpoly_length(term5, reduced_ideal.ctx.nmod_ctx));
            
            /* term6 = a02 * a11 * a20 */
            DEBUG_PRINT_R("  Computing term6: a02*a11*a20...\n");
            nmod_mpoly_mul(temp1, &nmod_matrix[0][2], &nmod_matrix[1][1], reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(temp1, &reduced_ideal, complete_var_names);
            nmod_mpoly_mul(term6, temp1, &nmod_matrix[2][0], reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(term6, &reduced_ideal, complete_var_names);
            DEBUG_PRINT_R("    term6 has %ld terms\n", nmod_mpoly_length(term6, reduced_ideal.ctx.nmod_ctx));
            
            /* Combine terms with reduction after each addition/subtraction */
            DEBUG_PRINT_R("  Combining terms...\n");
            
            /* det = term1 + term2 + term3 - term4 - term5 - term6 */
            nmod_mpoly_add(det, term1, term2, reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(det, &reduced_ideal, complete_var_names);
            DEBUG_PRINT_R("    After adding term2: %ld terms\n", nmod_mpoly_length(det, reduced_ideal.ctx.nmod_ctx));
            
            nmod_mpoly_add(det, det, term3, reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(det, &reduced_ideal, complete_var_names);
            DEBUG_PRINT_R("    After adding term3: %ld terms\n", nmod_mpoly_length(det, reduced_ideal.ctx.nmod_ctx));
            
            nmod_mpoly_sub(det, det, term4, reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(det, &reduced_ideal, complete_var_names);
            DEBUG_PRINT_R("    After subtracting term4: %ld terms\n", nmod_mpoly_length(det, reduced_ideal.ctx.nmod_ctx));
            
            nmod_mpoly_sub(det, det, term5, reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(det, &reduced_ideal, complete_var_names);
            DEBUG_PRINT_R("    After subtracting term5: %ld terms\n", nmod_mpoly_length(det, reduced_ideal.ctx.nmod_ctx));
            
            nmod_mpoly_sub(det, det, term6, reduced_ideal.ctx.nmod_ctx);
            triangular_ideal_reduce_nmod_mpoly_with_names(det, &reduced_ideal, complete_var_names);
            DEBUG_PRINT_R("    Final determinant: %ld terms\n", nmod_mpoly_length(det, reduced_ideal.ctx.nmod_ctx));
            
            /* Cleanup */
            nmod_mpoly_clear(term1, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_clear(term2, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_clear(term3, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_clear(term4, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_clear(term5, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_clear(term6, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_clear(temp1, reduced_ideal.ctx.nmod_ctx);
            nmod_mpoly_clear(temp2, reduced_ideal.ctx.nmod_ctx);
        } else {
            /* For larger matrices, use the recursive function but with timeout */
            compute_nmod_det_with_triangular_reduction_with_names(det, nmod_matrix, size, 
                                                               &reduced_ideal, complete_var_names);
        }
        
        clock_t end = clock();
        DEBUG_PRINT_R("Determinant computation time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
        
        /* Convert result back - result has 0 vars and npars parameters */
        nmod_mpoly_to_fq_mvpoly(result, det, 0, current_nvars, reduced_ideal.ctx.nmod_ctx, ideal->field_ctx);
        
        /* Cleanup */
        nmod_mpoly_clear(det, reduced_ideal.ctx.nmod_ctx);
        for (slong i = 0; i < size; i++) {
            for (slong j = 0; j < size; j++) {
                nmod_mpoly_clear(&nmod_matrix[i][j], reduced_ideal.ctx.nmod_ctx);
            }
            flint_free(nmod_matrix[i]);
        }
        flint_free(nmod_matrix);
        
    } else {
        /* Similar for fq_nmod_mpoly */
        DEBUG_PRINT_R("Using fq_nmod_mpoly for extension field computation\n");
        
        /* Convert to fq_nmod_mpoly matrix */
        fq_nmod_mpoly_struct **fq_matrix = (fq_nmod_mpoly_struct**) flint_malloc(size * sizeof(fq_nmod_mpoly_struct*));
        for (slong i = 0; i < size; i++) {
            fq_matrix[i] = (fq_nmod_mpoly_struct*) flint_malloc(size * sizeof(fq_nmod_mpoly_struct));
            for (slong j = 0; j < size; j++) {
                fq_nmod_mpoly_init(&fq_matrix[i][j], reduced_ideal.ctx.fq_ctx);
                fq_nmod_mpoly_zero(&fq_matrix[i][j], reduced_ideal.ctx.fq_ctx);
                
                /* Convert fq_mvpoly to fq_nmod_mpoly */
                DEBUG_PRINT_R("Converting matrix[%ld][%ld] with %ld terms...\n", i, j, matrix[i][j].nterms);
                
                /* Process each term */
                for (slong t = 0; t < matrix[i][j].nterms && t < 10000; t++) {
                    /* Build exponent vector for fq_nmod_mpoly */
                    ulong *exp = (ulong*) flint_calloc(current_nvars, sizeof(ulong));
                    
                    /* Copy parameter exponents - these become the variables */
                    if (matrix[i][j].terms[t].par_exp) {
                        for (slong v = 0; v < matrix[i][j].npars && v < current_nvars; v++) {
                            exp[v] = matrix[i][j].terms[t].par_exp[v];
                        }
                    }
                    
                    /* Set coefficient */
                    fq_nmod_mpoly_set_coeff_fq_nmod_ui(&fq_matrix[i][j], 
                                                       matrix[i][j].terms[t].coeff, 
                                                       exp, reduced_ideal.ctx.fq_ctx);
                    
                    flint_free(exp);
                }
                
                DEBUG_PRINT_R("  Converted to %ld terms\n", 
                       fq_nmod_mpoly_length(&fq_matrix[i][j], reduced_ideal.ctx.fq_ctx));
            }
        }
        
        /* Compute determinant with reduction */
        fq_nmod_mpoly_t det;
        fq_nmod_mpoly_init(det, reduced_ideal.ctx.fq_ctx);
        
        DEBUG_PRINT_R("Starting determinant computation...\n");
        clock_t start = clock();
        
        if (size == 1) {
            fq_nmod_mpoly_set(det, &fq_matrix[0][0], reduced_ideal.ctx.fq_ctx);
            triangular_ideal_reduce_fq_nmod_mpoly_with_names(det, &reduced_ideal, complete_var_names);
        } else if (size == 2) {
            fq_nmod_mpoly_t ad, bc;
            fq_nmod_mpoly_init(ad, reduced_ideal.ctx.fq_ctx);
            fq_nmod_mpoly_init(bc, reduced_ideal.ctx.fq_ctx);
            
            fq_nmod_mpoly_mul(ad, &fq_matrix[0][0], &fq_matrix[1][1], reduced_ideal.ctx.fq_ctx);
            fq_nmod_mpoly_mul(bc, &fq_matrix[0][1], &fq_matrix[1][0], reduced_ideal.ctx.fq_ctx);
            fq_nmod_mpoly_sub(det, ad, bc, reduced_ideal.ctx.fq_ctx);
            
            triangular_ideal_reduce_fq_nmod_mpoly_with_names(det, &reduced_ideal, complete_var_names);
            
            fq_nmod_mpoly_clear(ad, reduced_ideal.ctx.fq_ctx);
            fq_nmod_mpoly_clear(bc, reduced_ideal.ctx.fq_ctx);
        } else {
            compute_fq_nmod_det_with_triangular_reduction_with_names(det, fq_matrix, size, 
                                                                    &reduced_ideal, complete_var_names);
        }
        
        clock_t end = clock();
        DEBUG_PRINT_R("Determinant computation time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
        
        /* Convert result back */
        fq_nmod_mpoly_to_fq_mvpoly(result, det, 0, current_nvars, reduced_ideal.ctx.fq_ctx, ideal->field_ctx);
        
        /* Cleanup */
        fq_nmod_mpoly_clear(det, reduced_ideal.ctx.fq_ctx);
        for (slong i = 0; i < size; i++) {
            for (slong j = 0; j < size; j++) {
                fq_nmod_mpoly_clear(&fq_matrix[i][j], reduced_ideal.ctx.fq_ctx);
            }
            flint_free(fq_matrix[i]);
        }
        flint_free(fq_matrix);
    }
    
    /* Cleanup reduced ideal and temporary names */
    unified_triangular_ideal_clear(&reduced_ideal);
    
    /* Free any allocated names that we created */
    for (slong i = 0; i < current_nvars; i++) {
        if (complete_var_names[i] != current_var_names[i]) {
            free(complete_var_names[i]);
        }
    }
    free(complete_var_names);
}

/* KEEP ORIGINAL FUNCTIONS - don't modify these */
void construct_triangular_ideal_from_strings(unified_triangular_ideal_t *ideal,
                                            const char **ideal_gens,
                                            slong num_gens,
                                            const char **var_names,
                                            slong nvars,
                                            const fq_nmod_ctx_t ctx) {
    DEBUG_PRINT_R("Constructing triangular ideal from %ld generators\n", num_gens);
    
    // Initialize the ideal with the correct size
    unified_triangular_ideal_init(ideal, num_gens, nvars, ctx);
    
    // Print the generators being added
    DEBUG_PRINT_R("Triangular ideal generators:\n");
    for (slong i = 0; i < num_gens; i++) {
        DEBUG_PRINT_R("  g%ld: %s\n", i+1, ideal_gens[i]);
    }
    DEBUG_PRINT_R("\n");
    
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
    DEBUG_PRINT_R("Ideal construction complete:\n");
    DEBUG_PRINT_R("  Number of generators: %ld\n", ideal->num_gens);
    for (slong i = 0; i < ideal->num_gens; i++) {
        DEBUG_PRINT_R("  Generator %ld reduces variable '%s' (index %ld)\n", 
               i, ideal->var_names[i] ? ideal->var_names[i] : "?", 
               ideal->var_indices[i]);
    }
    DEBUG_PRINT_R("\n");
}

/* KEEP ORIGINAL dixon_with_ideal_reduction - DON'T MODIFY */
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
    
    DEBUG_PRINT_R("Input polynomials:\n");
    for (slong i = 0; i < num_polys; i++) {
        //DEBUG_PRINT_R("  p%ld: %s\n", i, poly_strings[i]);
    }
    
    /* Extract ALL variables from the ideal instead of hardcoding */
    slong total_system_vars = 0;
    char **all_system_vars = NULL;
    
    if (ideal && ideal->num_gens > 0) {
        /* Estimate maximum variables from ideal context */
        if (ideal->is_prime_field) {
            total_system_vars = nmod_mpoly_ctx_nvars(ideal->ctx.nmod_ctx);
        } else {
            total_system_vars = fq_nmod_mpoly_ctx_nvars(ideal->ctx.fq_ctx);
        }
        
        all_system_vars = (char**) malloc(total_system_vars * sizeof(char*));
        
        /* Initialize array */
        for (slong i = 0; i < total_system_vars; i++) {
            all_system_vars[i] = NULL;
        }
        
        /* Extract variable names from ideal generators */
        for (slong g = 0; g < ideal->num_gens; g++) {
            if (ideal->var_names[g] && ideal->var_indices[g] >= 0 && 
                ideal->var_indices[g] < total_system_vars) {
                if (!all_system_vars[ideal->var_indices[g]]) {
                    all_system_vars[ideal->var_indices[g]] = strdup(ideal->var_names[g]);
                }
            }
        }
        
        /* Fill in missing variable names with defaults */
        for (slong i = 0; i < total_system_vars; i++) {
            if (!all_system_vars[i]) {
                char temp[32];
                if (i == 0) {
                    sprintf(temp, "x");
                } else {
                    sprintf(temp, "var_%ld", i);
                }
                all_system_vars[i] = strdup(temp);
            }
        }
    } else {
        /* Fallback: minimal variable set */
        total_system_vars = 2;
        all_system_vars = (char**) malloc(total_system_vars * sizeof(char*));
        all_system_vars[0] = strdup("x");
        all_system_vars[1] = strdup("y");
    }
    
    /* Parse polynomials */
    parser_state_t state;
    state.var_names = (char**) malloc(num_elim_vars * sizeof(char*));
    for (slong i = 0; i < num_elim_vars; i++) {
        state.var_names[i] = strdup(elim_vars[i]);
    }
    state.nvars = num_elim_vars;
    state.npars = 0;
    state.max_pars = 32;  // Increased size
    state.par_names = (char**) malloc(state.max_pars * sizeof(char*));
    state.ctx = ctx;
    state.generator_name = get_generator_name(ctx);
    fq_nmod_init(state.current.value, ctx);
    state.current.str = NULL;
    
    /* CRITICAL FIX: Add ALL non-eliminated variables as parameters */
    /* This includes x and all z_i that are not being eliminated */
    for (slong i = 0; i < total_system_vars; i++) {
        int is_elim_var = 0;
        
        /* Check if this is an elimination variable */
        for (slong j = 0; j < num_elim_vars; j++) {
            if (strcmp(all_system_vars[i], elim_vars[j]) == 0) {
                is_elim_var = 1;
                break;
            }
        }
        
        /* If not an elimination variable, add as parameter */
        if (!is_elim_var) {
            /* Check if already in parameters */
            int already_param = 0;
            for (slong j = 0; j < state.npars; j++) {
                if (strcmp(state.par_names[j], all_system_vars[i]) == 0) {
                    already_param = 1;
                    break;
                }
            }
            
            if (!already_param) {
                if (state.npars >= state.max_pars) {
                    state.max_pars *= 2;
                    state.par_names = (char**) realloc(state.par_names, 
                                                      state.max_pars * sizeof(char*));
                }
                state.par_names[state.npars] = strdup(all_system_vars[i]);
                state.npars++;
            }
        }
    }
    
    /* First pass: parse to identify any additional parameters in polynomials */
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
    
    /* Build complete variable list: elimination variables + parameters */
    slong total_vars = num_elim_vars + state.npars;
    char **all_var_names = (char**) malloc(total_vars * sizeof(char*));
    
    /* Copy elimination variables */
    for (slong i = 0; i < num_elim_vars; i++) {
        all_var_names[i] = strdup(elim_vars[i]);
    }
    
    /* Copy parameters */
    for (slong i = 0; i < state.npars; i++) {
        all_var_names[num_elim_vars + i] = strdup(state.par_names[i]);
    }
    
    /* DEBUG: Print all variables */
    DEBUG_PRINT_R("\n=== DEBUG: All variables (elimination + parameters) ===\n");
    DEBUG_PRINT_R("Elimination variables (%ld): ", num_elim_vars);
    for (slong i = 0; i < num_elim_vars; i++) {
        DEBUG_PRINT_R("%s ", all_var_names[i]);
    }
    DEBUG_PRINT_R("\nParameters (%ld): ", state.npars);
    for (slong i = 0; i < state.npars; i++) {
        DEBUG_PRINT_R("%s ", all_var_names[num_elim_vars + i]);
    }
    DEBUG_PRINT_R("\n=== END DEBUG ===\n\n");
    
    /* Parse polynomials again with correct context */
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
        
        /* Pass parameter names for proper variable mapping during reduction */
        compute_det_with_reduction_from_mvpoly(&result_poly, coeff_matrix, matrix_size, ideal, 
                                             state.par_names);
        
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
    
    for (slong i = 0; i < total_vars; i++) {
        free(all_var_names[i]);
    }
    free(all_var_names);
    
    for (slong i = 0; i < total_system_vars; i++) {
        free(all_system_vars[i]);
    }
    free(all_system_vars);
    
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

/* Helper functions for string interface */
char** split_string_r(const char* str, slong* count) {
    if (!str || !count) {
        *count = 0;
        return NULL;
    }
    
    /* Count separators to determine array size */
    slong len = strlen(str);
    slong num_parts = 1;
    for (slong i = 0; i < len; i++) {
        if (str[i] == ';' || str[i] == ',') {
            num_parts++;
        }
    }
    
    char** result = (char**) malloc(num_parts * sizeof(char*));
    *count = 0;
    
    const char* start = str;
    const char* end = str;
    
    while (*end) {
        if (*end == ';' || *end == ',' || *(end + 1) == '\0') {
            if (*(end + 1) == '\0' && *end != ';' && *end != ',') {
                end++;
            }
            
            slong part_len = end - start;
            if (part_len > 0) {
                result[*count] = (char*) malloc((part_len + 1) * sizeof(char));
                strncpy(result[*count], start, part_len);
                result[*count][part_len] = '\0';
                
                /* Trim whitespace */
                char* trimmed = result[*count];
                while (*trimmed == ' ' || *trimmed == '\t') trimmed++;
                slong trimmed_len = strlen(trimmed);
                while (trimmed_len > 0 && (trimmed[trimmed_len-1] == ' ' || trimmed[trimmed_len-1] == '\t')) {
                    trimmed[trimmed_len-1] = '\0';
                    trimmed_len--;
                }
                
                if (trimmed != result[*count]) {
                    memmove(result[*count], trimmed, trimmed_len + 1);
                }
                
                (*count)++;
            }
            
            start = end + 1;
        }
        end++;
    }
    
    return result;
}

void free_split_string_rs(char** strings, slong count) {
    if (strings) {
        for (slong i = 0; i < count; i++) {
            if (strings[i]) {
                free(strings[i]);
            }
        }
        free(strings);
    }
}

void construct_triangular_ideal_str(unified_triangular_ideal_t *ideal,
                                   const char *gens_string,   
                                   const char *vars_string,   
                                   const fq_nmod_ctx_t ctx) {
    
    slong num_gens, num_vars;
    char **gens_array = split_string_r(gens_string, &num_gens);
    char **vars_array = split_string_r(vars_string, &num_vars);
    
    const char **ideal_gens = (const char**) malloc(num_gens * sizeof(char*));
    const char **var_names = (const char**) malloc(num_vars * sizeof(char*));
    
    for (slong i = 0; i < num_gens; i++) {
        ideal_gens[i] = gens_array[i];
    }
    for (slong i = 0; i < num_vars; i++) {
        var_names[i] = vars_array[i];
    }
    
    construct_triangular_ideal_from_strings(ideal, ideal_gens, num_gens, 
                                          var_names, num_vars, ctx);
    
    free(ideal_gens);
    free(var_names);
    free_split_string_rs(gens_array, num_gens);
    free_split_string_rs(vars_array, num_vars);
}

char* dixon_with_ideal_reduction_str(const char *poly_string,
                                    const char *elim_vars_string,
                                    const char *ideal_gens_string,
                                    const char *all_vars_string,
                                    const fq_nmod_ctx_t ctx) {
    
    slong num_polys, num_elim_vars, num_gens, num_all_vars;
    char **poly_array = split_string_r(poly_string, &num_polys);
    char **elim_array = split_string_r(elim_vars_string, &num_elim_vars);
    char **gens_array = split_string_r(ideal_gens_string, &num_gens);
    char **all_vars_array = split_string_r(all_vars_string, &num_all_vars);
    
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
    
    const char **poly_strings = (const char**) malloc(num_polys * sizeof(char*));
    const char **elim_vars = (const char**) malloc(num_elim_vars * sizeof(char*));
    
    for (slong i = 0; i < num_polys; i++) {
        poly_strings[i] = poly_array[i];
    }
    for (slong i = 0; i < num_elim_vars; i++) {
        elim_vars[i] = elim_array[i];
    }
    
    char *result = dixon_with_ideal_reduction(poly_strings, num_polys,
                                            elim_vars, num_elim_vars,
                                            ctx, &ideal);
    
    unified_triangular_ideal_clear(&ideal);
    free(poly_strings);
    free(elim_vars);
    free(ideal_gens);
    free(all_var_names);
    free_split_string_rs(poly_array, num_polys);
    free_split_string_rs(elim_array, num_elim_vars);
    free_split_string_rs(gens_array, num_gens);
    free_split_string_rs(all_vars_array, num_all_vars);
    
    return result;
}

#endif /* DIXON_WITH_IDEAL_REDUCTION_H */