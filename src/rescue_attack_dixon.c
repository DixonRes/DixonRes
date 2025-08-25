// gcc -O3 -march=native -o rescue rescue_attack_dixon.c -lflint -lmpfr -lgmp -lpthread -L/home/suohaohai02/mylinks -lflint -lstdc++ -lpml2 -fopenmp
/* rescue_attack_dixon.c - Build Rescue system and perform elimination */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>
#include <flint/nmod_mpoly.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_poly.h>
#include <flint/fq_nmod_poly_factor.h>

#include "dixon_flint.h"
#include "dixon_with_ideal_reduction.h"
//#include "resultant_with_ideal_reduction.h"
/* Rescue parameters */
#define RESCUE_P      0x64ec6dd0392073ULL
#define RESCUE_D      3
#define RESCUE_N      55
#define RESCUE_T      3
#define RESCUE_R      4

#define RESCUE_ALPHA0 1ULL
#define RESCUE_ALPHA1 28407454060060786ULL
#define RESCUE_ALPHA2 0ULL

#define RESCUE_BETA0  8827379057032300ULL
#define RESCUE_BETA1  11447156669198427ULL
#define RESCUE_BETA2  6228383861343320ULL

/* Round constants */
mp_limb_t round_constants[36] = {
    27498789407612574ULL, 22052048824113086ULL, 2941138345603856ULL, 
    25244655091391380ULL, 27864432703557507ULL, 10655481870000693ULL,
    12029986412025777ULL, 19484968723987120ULL, 13695593575120347ULL,
    11363016928652497ULL, 4848123124202167ULL, 202662769141873ULL,
    10388262623185440ULL, 9627757603546016ULL, 15459910451817652ULL,
    10472020289946983ULL, 18225111667665967ULL, 21110845623706186ULL,
    5962772680402733ULL, 7588371581963569ULL, 24389178909705370ULL,
    7372288933856325ULL, 13064415177158011ULL, 25314800050133692ULL,
    7495881900077152ULL, 27937615950064382ULL, 21671355553166345ULL,
    418504290749794ULL, 1626942147407155ULL, 8523003710549504ULL,
    25482532338301975ULL, 23263909157345227ULL, 19560921910827017ULL,
    18762876693549419ULL, 12826507476828121ULL, 27770298731614453ULL
};

/* MDS matrix */
mp_limb_t MDS[3][3] = {{2, 7, 1}, {1, 2, 7}, {7, 1, 2}};

static slong ctr_cst = 0;

/* Convert polynomial to string - output ALL terms */
char* poly_to_string_complete(const nmod_mpoly_t poly, const char **var_names, 
                             slong nvars, const nmod_mpoly_ctx_t ctx) {
    slong nterms = nmod_mpoly_length(poly, ctx);
    size_t buf_size = nterms * 200 + 1024;  // Allocate enough space
    char *buffer = (char*) malloc(buf_size);
    buffer[0] = '\0';
    
    if (nterms == 0) {
        strcpy(buffer, "0");
        return buffer;
    }
    
    for (slong i = 0; i < nterms; i++) {
        if (i > 0) strcat(buffer, " + ");
        
        mp_limb_t coeff = nmod_mpoly_get_term_coeff_ui(poly, i, ctx);
        ulong *exp = (ulong*) flint_malloc(nvars * sizeof(ulong));
        nmod_mpoly_get_term_exp_ui(exp, poly, i, ctx);
        
        char coeff_str[32];
        sprintf(coeff_str, "%lu", coeff);
        strcat(buffer, coeff_str);
        
        for (slong j = 0; j < nvars; j++) {
            if (exp[j] > 0) {
                strcat(buffer, "*");
                strcat(buffer, var_names[j]);
                if (exp[j] > 1) {
                    char exp_str[32];
                    sprintf(exp_str, "^%lu", exp[j]);
                    strcat(buffer, exp_str);
                }
            }
        }
        
        if (strlen(buffer) > buf_size - 500) {
            buf_size *= 2;
            buffer = (char*) realloc(buffer, buf_size);
        }
        
        flint_free(exp);
    }
    
    return buffer;
}

/* Get degrees of all variables in polynomial */
void get_poly_degrees(const nmod_mpoly_t poly, slong *degrees, slong nvars, 
                     const nmod_mpoly_ctx_t ctx) {
    for (slong i = 0; i < nvars; i++) {
        degrees[i] = 0;
    }
    
    slong nterms = nmod_mpoly_length(poly, ctx);
    for (slong t = 0; t < nterms; t++) {
        ulong *exp = (ulong*) flint_malloc(nvars * sizeof(ulong));
        nmod_mpoly_get_term_exp_ui(exp, poly, t, ctx);
        
        for (slong v = 0; v < nvars; v++) {
            if (exp[v] > degrees[v]) {
                degrees[v] = exp[v];
            }
        }
        
        flint_free(exp);
    }
}

/* Apply MDS matrix transformation */
void apply_mds_nmod(nmod_mpoly_t *state, nmod_mpoly_t *new_state, 
                   const nmod_mpoly_ctx_t ctx) {
    for (slong i = 0; i < RESCUE_T; i++) {
        nmod_mpoly_zero(new_state[i], ctx);
        
        for (slong j = 0; j < RESCUE_T; j++) {
            if (MDS[i][j] != 0) {
                nmod_mpoly_t temp;
                nmod_mpoly_init(temp, ctx);
                nmod_mpoly_scalar_mul_ui(temp, state[j], MDS[i][j], ctx);
                nmod_mpoly_add(new_state[i], new_state[i], temp, ctx);
                nmod_mpoly_clear(temp, ctx);
            }
        }
    }
}

/* Add round constants */
void add_round_constants_nmod(nmod_mpoly_t *state, const nmod_mpoly_ctx_t ctx) {
    for (slong i = 0; i < RESCUE_T; i++) {
        nmod_mpoly_add_ui(state[i], state[i], round_constants[ctr_cst++], ctx);
    }
}

/* Main function */
int main(void) {
    mp_limb_t prime = RESCUE_P;
    
    printf("\n=== Rescue Polynomial System with Elimination ===\n");
    printf("Prime: %lu\n", prime);
    printf("Rounds: %d\n\n", RESCUE_R);
    
    clock_t total_start = clock();
    
    slong total_vars = 1 + 3 * RESCUE_R;  // x + z_1,...,z_9
    
    /* Create variable names */
    const char **var_names = (const char**) malloc(total_vars * sizeof(char*));
    var_names[0] = strdup("z0");
    for (slong i = 1; i < total_vars; i++) {
        char *name = (char*) malloc(32);
        sprintf(name, "z%ld", i);
        var_names[i] = name;
    }
    
    /* Create contexts */
    nmod_mpoly_ctx_t ctx;
    nmod_mpoly_ctx_init(ctx, total_vars, ORD_LEX, prime);
    
    fq_nmod_ctx_t fq_ctx;
    fq_nmod_ctx_init_ui(fq_ctx, prime, 1, "t");
    
    /* Initialize state arrays */
    nmod_mpoly_t state[RESCUE_T];
    nmod_mpoly_t temp_state[RESCUE_T];
    for (slong i = 0; i < RESCUE_T; i++) {
        nmod_mpoly_init(state[i], ctx);
        nmod_mpoly_init(temp_state[i], ctx);
    }
    
    /* Store ideal generators */
    char *ideal_strings[RESCUE_R * RESCUE_T];
    slong ideal_count = 0;
    
    ctr_cst = 2 * RESCUE_T;  // Skip first round constants
    
    /* Initial state */
    ulong *exp = (ulong*) flint_calloc(total_vars, sizeof(ulong));
    
    /* state[0] = ALPHA0*x + BETA0 */
    if (RESCUE_ALPHA0 != 0) {
        exp[0] = 1;
        nmod_mpoly_push_term_ui_ui(state[0], RESCUE_ALPHA0, exp, ctx);
        exp[0] = 0;
    }
    nmod_mpoly_add_ui(state[0], state[0], RESCUE_BETA0, ctx);
    
    /* state[1] = ALPHA1*x + BETA1 */
    if (RESCUE_ALPHA1 != 0) {
        exp[0] = 1;
        nmod_mpoly_push_term_ui_ui(state[1], RESCUE_ALPHA1, exp, ctx);
        exp[0] = 0;
    }
    nmod_mpoly_add_ui(state[1], state[1], RESCUE_BETA1, ctx);
    
    /* state[2] = BETA2 */
    nmod_mpoly_add_ui(state[2], state[2], RESCUE_BETA2, ctx);
    
    flint_free(exp);
    
    /* Process rounds and collect ideal generators */
    for (slong round = 1; round <= RESCUE_R; round++) {
        /* S-box: state^3 */
        for (slong i = 0; i < RESCUE_T; i++) {
            nmod_mpoly_pow_ui(temp_state[i], state[i], RESCUE_D, ctx);
        }
        
        /* First linear layer */
        nmod_mpoly_t after_linear[RESCUE_T];
        for (slong i = 0; i < RESCUE_T; i++) {
            nmod_mpoly_init(after_linear[i], ctx);
        }
        apply_mds_nmod(temp_state, after_linear, ctx);
        
        /* First constants */
        add_round_constants_nmod(after_linear, ctx);
        
        /* Create ideal generators: z_i^3 - expression */
        slong base_z = 1 + 3 * (round - 1);
        for (slong i = 0; i < RESCUE_T; i++) {
            nmod_mpoly_t ideal_gen;
            nmod_mpoly_init(ideal_gen, ctx);
            
            /* z_i^3 term */
            ulong *z_exp = (ulong*) flint_calloc(total_vars, sizeof(ulong));
            z_exp[base_z + i] = 3;
            nmod_mpoly_push_term_ui_ui(ideal_gen, 1, z_exp, ctx);
            flint_free(z_exp);
            
            /* Subtract expression */
            nmod_mpoly_sub(ideal_gen, ideal_gen, after_linear[i], ctx);
            
            ideal_strings[ideal_count] = poly_to_string_complete(ideal_gen, var_names, total_vars, ctx);
            ideal_count++;
            
            nmod_mpoly_clear(ideal_gen, ctx);
        }
        
        /* Replace with z variables for next round */
        for (slong i = 0; i < RESCUE_T; i++) {
            nmod_mpoly_zero(state[i], ctx);
            ulong *z_exp = (ulong*) flint_calloc(total_vars, sizeof(ulong));
            z_exp[base_z + i] = 1;
            nmod_mpoly_push_term_ui_ui(state[i], 1, z_exp, ctx);
            flint_free(z_exp);
        }
        
        /* Second linear layer */
        nmod_mpoly_t after_linear2[RESCUE_T];
        for (slong i = 0; i < RESCUE_T; i++) {
            nmod_mpoly_init(after_linear2[i], ctx);
        }
        apply_mds_nmod(state, after_linear2, ctx);
        
        /* Second constants */
        add_round_constants_nmod(after_linear2, ctx);
        
        /* Update state for next round */
        for (slong i = 0; i < RESCUE_T; i++) {
            nmod_mpoly_set(state[i], after_linear2[i], ctx);
        }
        
        /* Cleanup */
        for (slong i = 0; i < RESCUE_T; i++) {
            nmod_mpoly_clear(after_linear[i], ctx);
            nmod_mpoly_clear(after_linear2[i], ctx);
        }
    }
    
    /* Print the triangular ideal generators */
    printf("[+] Triangular ideal generators:\n");
    for (slong i = 0; i < ideal_count; i++) {
        printf("g_%ld = %s\n", i+1, ideal_strings[i]);
    }
    
    /* Print final polynomial */
    printf("\n[+] Final polynomial (state[2]):\n");
    char *current_poly = poly_to_string_complete(state[2], var_names, total_vars, ctx);
    printf("Y = %s\n", current_poly);
    
    /* Get initial degrees */
    slong *degrees = (slong*) calloc(total_vars, sizeof(slong));
    get_poly_degrees(state[2], degrees, total_vars, ctx);
    printf("\n[+] Initial degrees: ");
    for (slong i = 0; i < total_vars; i++) {
        printf("%ld ", degrees[i]);
    }
    printf("\n");
    free(degrees);
    
    /* Construct triangular ideal */
    printf("\n[+] Constructing triangular ideal...\n");
    unified_triangular_ideal_t ideal;
    construct_triangular_ideal_from_strings(&ideal, 
                                          (const char**)ideal_strings, ideal_count,
                                          var_names, total_vars, fq_ctx);
    
    /* ========== ELIMINATION PHASE ========== */
    printf("\n[+] Starting elimination phase...\n");
    clock_t elim_start = clock();
    
    /* Eliminate variables from z_9 down to z_1 */
    for (slong elim_round = ideal_count; elim_round >= 1; elim_round--) {
        printf("\n=== Round %ld: Eliminating %s ===\n", ideal_count - elim_round + 1, var_names[elim_round]);
        
        /* Get degrees before elimination */
        degrees = (slong*) calloc(total_vars, sizeof(slong));
        
        /* Parse current polynomial to get degrees */
        parser_state_t parse_state;
        parse_state.var_names = (char**) var_names;
        parse_state.nvars = 0;  // We're using parameters
        parse_state.npars = total_vars;
        parse_state.max_pars = total_vars;
        parse_state.par_names = (char**) var_names;
        parse_state.ctx = fq_ctx;
        parse_state.generator_name = get_generator_name(fq_ctx);
        fq_nmod_init(parse_state.current.value, fq_ctx);
        parse_state.current.str = NULL;
        parse_state.input = current_poly;
        parse_state.pos = 0;
        parse_state.len = strlen(current_poly);
        next_token(&parse_state);
        
        fq_mvpoly_t temp_mvpoly;
        fq_mvpoly_init(&temp_mvpoly, 0, total_vars, fq_ctx);
        parse_expression(&parse_state, &temp_mvpoly);
        
        /* Get degrees from mvpoly */
        for (slong t = 0; t < temp_mvpoly.nterms; t++) {
            if (temp_mvpoly.terms[t].par_exp) {
                for (slong v = 0; v < total_vars; v++) {
                    if (temp_mvpoly.terms[t].par_exp[v] > degrees[v]) {
                        degrees[v] = temp_mvpoly.terms[t].par_exp[v];
                    }
                }
            }
        }
        
        printf("[+] Degrees before elimination: ");
        for (slong i = 0; i < total_vars; i++) {
            if (degrees[i] > 0) {
                printf("%s:%ld ", var_names[i], degrees[i]);
            }
        }
        printf("\n");
        printf("[+] Maximum degree in x: %ld\n", degrees[0]);
        printf("[+] Polynomial has %ld terms\n", temp_mvpoly.nterms);
        
        free(degrees);
        fq_mvpoly_clear(&temp_mvpoly);
        fq_nmod_clear(parse_state.current.value, fq_ctx);
        if (parse_state.current.str) free(parse_state.current.str);
        if (parse_state.generator_name) free(parse_state.generator_name);
        
        /* Perform Dixon elimination with ideal reduction */
        clock_t round_start = clock();
        
        /* Prepare polynomials for elimination */
        const char *poly_array[2] = { current_poly, ideal_strings[elim_round - 1] };
        const char *elim_vars[1] = { var_names[elim_round] };
        
        printf("[+] Eliminating %s using Dixon with ideal reduction...\n", var_names[elim_round]);
        
        char *new_poly = dixon_with_ideal_reduction(poly_array, 2, elim_vars, 1, 
                                                   fq_ctx, &ideal); //elimination_with_ideal_reduction dixon_with_ideal_reduction
        
        double round_time = (double)(clock() - round_start) / CLOCKS_PER_SEC;
        printf("[+] Elimination completed in %.3f seconds\n", round_time);
        
        /* Update current polynomial */
        free(current_poly);
        current_poly = new_poly;
        
        /* Check if polynomial became zero */
        if (strcmp(current_poly, "0") == 0) {
            printf("\n[!] Polynomial became zero, stopping elimination\n");
            break;
        }
        
        // Print result preview (first 500 chars)
        printf("[+] Result preview: ");
        if (strlen(current_poly) > 500) {
            printf("%.500s...\n", current_poly);
        } else {
            printf("%s\n", current_poly);
        }
        
    }
    
    double elim_time = (double)(clock() - elim_start) / CLOCKS_PER_SEC;
    printf("\n[+] Elimination complete (%.3f seconds)\n", elim_time);
    
    //printf("\n[+] Final univariate polynomial in x:\n");
    //printf("%s\n", current_poly);
    
    /* ========== ROOT FINDING ========== */
    if (strcmp(current_poly, "0") != 0) {
        printf("\n[+] Finding roots of univariate polynomial...\n");
        
        /* Parse univariate polynomial back to nmod_poly */
        parser_state_t parse_state;
        parse_state.var_names = (char**) malloc(1 * sizeof(char*));
        parse_state.var_names[0] = strdup("z0");
        parse_state.nvars = 1;
        parse_state.npars = 0;
        parse_state.max_pars = 16;
        parse_state.par_names = (char**) malloc(parse_state.max_pars * sizeof(char*));
        parse_state.ctx = fq_ctx;
        parse_state.generator_name = get_generator_name(fq_ctx);
        fq_nmod_init(parse_state.current.value, fq_ctx);
        parse_state.current.str = NULL;
        parse_state.input = current_poly;
        parse_state.pos = 0;
        parse_state.len = strlen(current_poly);
        next_token(&parse_state);
        
        fq_mvpoly_t univar_mvpoly;
        fq_mvpoly_init(&univar_mvpoly, 1, 0, fq_ctx);
        parse_expression(&parse_state, &univar_mvpoly);
        
        /* Convert to nmod_poly */
        nmod_poly_t univar_poly;
        nmod_poly_init(univar_poly, prime);
        
        /* Build polynomial from mvpoly terms */
        for (slong t = 0; t < univar_mvpoly.nterms; t++) {
            mp_limb_t coeff = 0;
            
            /* Extract coefficient */
            if (fq_nmod_is_one(univar_mvpoly.terms[t].coeff, fq_ctx)) {
                coeff = 1;
            } else if (!fq_nmod_is_zero(univar_mvpoly.terms[t].coeff, fq_ctx)) {
                nmod_poly_struct *p = &univar_mvpoly.terms[t].coeff[0];
                if (p->length > 0) {
                    coeff = p->coeffs[0];
                }
            }
            
            /* Get degree */
            slong deg = 0;
            if (univar_mvpoly.terms[t].var_exp) {
                deg = univar_mvpoly.terms[t].var_exp[0];
            }
            
            /* Set coefficient */
            if (coeff != 0) {
                nmod_poly_set_coeff_ui(univar_poly, deg, coeff);
            }
        }
        
        printf("[+] Univariate polynomial degree: %ld\n", nmod_poly_degree(univar_poly));
        
        /* Find roots using FLINT's dedicated root finding function */
        printf("[+] Univariate polynomial degree: %ld\n", nmod_poly_degree(univar_poly));
        
        /* Use FLINT's nmod_poly_roots function */
        nmod_poly_factor_t roots;
        nmod_poly_factor_init(roots);
        
        clock_t roots_start = clock();
        nmod_poly_roots(roots, univar_poly, 1);  /* 1 = with multiplicity */
        double roots_time = (double)(clock() - roots_start) / CLOCKS_PER_SEC;
        
        printf("[+] Root finding completed in %.3f seconds\n", roots_time);
        printf("[+] Number of roots found: %ld\n", roots->num);
        
        /* Print all roots */
        printf("\n[+] Roots in F_%lu:\n", prime);
        
        for (slong i = 0; i < roots->num; i++) {
            /* Each root is stored as a linear polynomial: x - root */
            /* So for polynomial ax + b, the root is -b/a */
            mp_limb_t a = nmod_poly_get_coeff_ui(roots->p + i, 1);
            mp_limb_t b = nmod_poly_get_coeff_ui(roots->p + i, 0);
            
            if (a != 0) {
                /* Compute root = -b/a mod p */
                nmod_t mod = univar_poly->mod;
                mp_limb_t a_inv = nmod_inv(a, mod);
                mp_limb_t neg_b = nmod_neg(b, mod);
                mp_limb_t root = nmod_mul(neg_b, a_inv, mod);
                
                printf("  Root %ld: x = %lu", i + 1, root);
                if (roots->exp[i] > 1) {
                    printf(" (multiplicity %ld)", roots->exp[i]);
                }
                printf("\n");
                
                /* Verify root */
                mp_limb_t eval = nmod_poly_evaluate_nmod(univar_poly, root);
                if (eval == 0) {
                    printf("    Verification: PASSED\n");
                } else {
                    printf("    Verification: FAILED (f(%lu) = %lu)\n", root, eval);
                }
            }
        }
        
        if (roots->num == 0) {
            printf("  No roots found in F_%lu\n", prime);
        } else {
            printf("\n[+] Total roots found: %ld\n", roots->num);
        }
        
        /* Cleanup */
        nmod_poly_factor_clear(roots);
        nmod_poly_clear(univar_poly);
        fq_mvpoly_clear(&univar_mvpoly);
        free(parse_state.var_names[0]);
        free(parse_state.var_names);
        for (slong i = 0; i < parse_state.npars; i++) {
            free(parse_state.par_names[i]);
        }
        free(parse_state.par_names);
        if (parse_state.generator_name) free(parse_state.generator_name);
        fq_nmod_clear(parse_state.current.value, fq_ctx);
        if (parse_state.current.str) free(parse_state.current.str);
    }
    
    double total_time = (double)(clock() - total_start) / CLOCKS_PER_SEC;
    printf("\n[+] Total computation time: %.3f seconds\n", total_time);
    
    /* Cleanup */
    unified_triangular_ideal_clear(&ideal);
    free(current_poly);
    
    for (slong i = 0; i < ideal_count; i++) {
        free(ideal_strings[i]);
    }
    
    for (slong i = 0; i < RESCUE_T; i++) {
        nmod_mpoly_clear(state[i], ctx);
        nmod_mpoly_clear(temp_state[i], ctx);
    }
    
    for (slong i = 0; i < total_vars; i++) {
        free((char*)var_names[i]);
    }
    free(var_names);
    
    nmod_mpoly_ctx_clear(ctx);
    fq_nmod_ctx_clear(fq_ctx);
    
    return 0;
}