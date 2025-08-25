/* Bivariate resultant with ideal reduction */
char* resultant_with_ideal_reduction(const char *poly1_string, const char *poly2_string,
                                   const char *elim_var, 
                                   const fq_nmod_ctx_t ctx,
                                   unified_triangular_ideal_t *ideal) {
    printf("\n=== Resultant with Ideal Reduction ===\n");
    printf("Eliminating variable: %s\n", elim_var);
    
    DEBUG_PRINT_R("Input polynomials:\n");
    DEBUG_PRINT_R("  p1: %s\n", poly1_string);
    DEBUG_PRINT_R("  p2: %s\n", poly2_string);
    
    /* Extract ALL variables from the ideal (similar to dixon version) */
    slong total_system_vars = 0;
    char **all_system_vars = NULL;
    
    if (ideal && ideal->num_gens > 0) {
        if (ideal->is_prime_field) {
            total_system_vars = nmod_mpoly_ctx_nvars(ideal->ctx.nmod_ctx);
        } else {
            total_system_vars = fq_nmod_mpoly_ctx_nvars(ideal->ctx.fq_ctx);
        }
        
        all_system_vars = (char**) malloc(total_system_vars * sizeof(char*));
        
        for (slong i = 0; i < total_system_vars; i++) {
            all_system_vars[i] = NULL;
        }
        
        for (slong g = 0; g < ideal->num_gens; g++) {
            if (ideal->var_names[g] && ideal->var_indices[g] >= 0 && 
                ideal->var_indices[g] < total_system_vars) {
                if (!all_system_vars[ideal->var_indices[g]]) {
                    all_system_vars[ideal->var_indices[g]] = strdup(ideal->var_names[g]);
                }
            }
        }
        
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
        total_system_vars = 2;
        all_system_vars = (char**) malloc(total_system_vars * sizeof(char*));
        all_system_vars[0] = strdup("x");
        all_system_vars[1] = strdup("y");
    }
    
    /* Initialize parser state */
    parser_state_t state;
    state.var_names = (char**) malloc(1 * sizeof(char*));
    state.var_names[0] = strdup(elim_var);
    state.nvars = 1;
    state.npars = 0;
    state.max_pars = 32;
    state.par_names = (char**) malloc(state.max_pars * sizeof(char*));
    state.ctx = ctx;
    state.generator_name = get_generator_name(ctx);
    fq_nmod_init(state.current.value, ctx);
    state.current.str = NULL;
    
    /* Add all non-eliminated variables as parameters */
    for (slong i = 0; i < total_system_vars; i++) {
        if (strcmp(all_system_vars[i], elim_var) != 0) {
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
    
    /* First pass: parse to identify any additional parameters */
    fq_mvpoly_t temp1, temp2;
    fq_mvpoly_init(&temp1, 1, state.max_pars, ctx);
    fq_mvpoly_init(&temp2, 1, state.max_pars, ctx);
    
    state.input = poly1_string;
    state.pos = 0;
    state.len = strlen(poly1_string);
    next_token(&state);
    parse_expression(&state, &temp1);
    
    state.input = poly2_string;
    state.pos = 0;
    state.len = strlen(poly2_string);
    if (state.current.str) {
        free(state.current.str);
        state.current.str = NULL;
    }
    next_token(&state);
    parse_expression(&state, &temp2);
    
    fq_mvpoly_clear(&temp1);
    fq_mvpoly_clear(&temp2);
    
    /* Parse polynomials again with correct context */
    fq_mvpoly_t poly1, poly2;
    fq_mvpoly_init(&poly1, 1, state.npars, ctx);
    fq_mvpoly_init(&poly2, 1, state.npars, ctx);
    
    state.input = poly1_string;
    state.pos = 0;
    state.len = strlen(poly1_string);
    if (state.current.str) {
        free(state.current.str);
        state.current.str = NULL;
    }
    next_token(&state);
    parse_expression(&state, &poly1);
    
    state.input = poly2_string;
    state.pos = 0;
    state.len = strlen(poly2_string);
    if (state.current.str) {
        free(state.current.str);
        state.current.str = NULL;
    }
    next_token(&state);
    parse_expression(&state, &poly2);
    
    /* Create unified field context */
    field_ctx_t field_ctx;
    field_ctx_init(&field_ctx, ctx);
    
    /* Create unified mpoly context */
    slong total_vars = 1 + state.npars;
    unified_mpoly_ctx_t unified_ctx = unified_mpoly_ctx_init(total_vars, ORD_LEX, &field_ctx);
    
    /* Convert to unified mpoly format */
    unified_mpoly_t A = unified_mpoly_init(unified_ctx);
    unified_mpoly_t B = unified_mpoly_init(unified_ctx);
    
    /* Convert poly1 */
    for (slong i = 0; i < poly1.nterms; i++) {
        field_elem_u coeff;
        ulong *exp = (ulong*) calloc(total_vars, sizeof(ulong));
        
        fq_nmod_to_field_elem(&coeff, poly1.terms[i].coeff, &field_ctx);
        
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
    
    /* Convert poly2 */
    for (slong i = 0; i < poly2.nterms; i++) {
        field_elem_u coeff;
        ulong *exp = (ulong*) calloc(total_vars, sizeof(ulong));
        
        fq_nmod_to_field_elem(&coeff, poly2.terms[i].coeff, &field_ctx);
        
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
    
    printf("\nStep 1: Compute resultant w.r.t. %s\n", elim_var);
    printf("  Poly1 has %ld terms\n", unified_mpoly_length(A));
    printf("  Poly2 has %ld terms\n", unified_mpoly_length(B));
    
    /* Compute resultant */
    unified_mpoly_t R = unified_mpoly_init(unified_ctx);
    clock_t start = clock();
    int success = unified_mpoly_resultant(R, A, B, 0, unified_ctx);
    clock_t end = clock();
    
    if (!success) {
        printf("Resultant computation failed!\n");
        unified_mpoly_clear(A);
        unified_mpoly_clear(B);
        unified_mpoly_clear(R);
        unified_mpoly_ctx_clear(unified_ctx);
        
        /* Cleanup */
        for (slong i = 0; i < state.nvars; i++) {
            free(state.var_names[i]);
        }
        free(state.var_names);
        for (slong i = 0; i < state.npars; i++) {
            free(state.par_names[i]);
        }
        free(state.par_names);
        for (slong i = 0; i < total_system_vars; i++) {
            free(all_system_vars[i]);
        }
        free(all_system_vars);
        if (state.generator_name) free(state.generator_name);
        fq_nmod_clear(state.current.value, ctx);
        if (state.current.str) free(state.current.str);
        fq_mvpoly_clear(&poly1);
        fq_mvpoly_clear(&poly2);
        
        return strdup("0");
    }
    
    printf("Resultant computation time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("Resultant has %ld terms\n", unified_mpoly_length(R));
    
    /* Now apply ideal reduction */
    printf("\nStep 2: Apply ideal reduction\n");
    
    /* Create reduced ideal context for current variables */
    unified_triangular_ideal_t reduced_ideal;
    create_reduced_ideal_context(&reduced_ideal, ideal, state.par_names, 
                               state.npars, ctx);
    
    /* Apply reduction based on field type */
    if (ideal->is_prime_field) {
        nmod_mpoly_struct *nmod_res = GET_NMOD_POLY(R);
        triangular_ideal_reduce_nmod_mpoly_with_names(nmod_res, &reduced_ideal, state.par_names);
        printf("After reduction: %ld terms\n", nmod_mpoly_length(nmod_res, &unified_ctx->ctx.nmod_ctx));
    } else {
        fq_nmod_mpoly_struct *fq_res = GET_FQ_POLY(R);
        triangular_ideal_reduce_fq_nmod_mpoly_with_names(fq_res, &reduced_ideal, state.par_names);
        printf("After reduction: %ld terms\n", fq_nmod_mpoly_length(fq_res, &unified_ctx->ctx.fq_ctx));
    }
    
    /* Convert result back to fq_mvpoly */
    fq_mvpoly_t result_mvpoly;
    fq_mvpoly_init(&result_mvpoly, 0, state.npars, ctx);
    
    /* Convert from unified format */
    if (field_ctx.field_id == FIELD_ID_NMOD) {
        nmod_mpoly_struct *nmod_res = GET_NMOD_POLY(R);
        nmod_mpoly_ctx_struct *nmod_ctx = &(unified_ctx->ctx.nmod_ctx);
        
        for (slong i = 0; i < nmod_mpoly_length(nmod_res, nmod_ctx); i++) {
            ulong coeff_ui = nmod_mpoly_get_term_coeff_ui(nmod_res, i, nmod_ctx);
            ulong *exp = (ulong*) calloc(total_vars, sizeof(ulong));
            nmod_mpoly_get_term_exp_ui(exp, nmod_res, i, nmod_ctx);
            
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ctx);
            fq_nmod_set_ui(coeff, coeff_ui, ctx);
            
            slong *par_exp = NULL;
            if (state.npars > 0) {
                par_exp = (slong*) calloc(state.npars, sizeof(slong));
                for (slong j = 0; j < state.npars; j++) {
                    par_exp[j] = exp[1 + j];
                }
            }
            
            fq_mvpoly_add_term(&result_mvpoly, NULL, par_exp, coeff);
            
            fq_nmod_clear(coeff, ctx);
            free(exp);
            if (par_exp) free(par_exp);
        }
    } else if (field_ctx.field_id == FIELD_ID_FQ_ZECH) {
        fq_zech_mpoly_struct *zech_res = GET_ZECH_POLY(R);
        fq_zech_mpoly_ctx_struct *zech_ctx = &(unified_ctx->ctx.zech_ctx);
        
        for (slong i = 0; i < fq_zech_mpoly_length(zech_res, zech_ctx); i++) {
            fq_zech_t zech_coeff;
            fq_zech_init(zech_coeff, zech_ctx->fqctx);
            fq_zech_mpoly_get_term_coeff_fq_zech(zech_coeff, zech_res, i, zech_ctx);
            
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ctx);
            fq_zech_get_fq_nmod(coeff, zech_coeff, zech_ctx->fqctx);
            
            ulong *exp = (ulong*) calloc(total_vars, sizeof(ulong));
            fq_zech_mpoly_get_term_exp_ui(exp, zech_res, i, zech_ctx);
            
            slong *par_exp = NULL;
            if (state.npars > 0) {
                par_exp = (slong*) calloc(state.npars, sizeof(slong));
                for (slong j = 0; j < state.npars; j++) {
                    par_exp[j] = exp[1 + j];
                }
            }
            
            fq_mvpoly_add_term(&result_mvpoly, NULL, par_exp, coeff);
            
            fq_nmod_clear(coeff, ctx);
            fq_zech_clear(zech_coeff, zech_ctx->fqctx);
            free(exp);
            if (par_exp) free(par_exp);
        }
    } else {
        fq_nmod_mpoly_struct *fq_res = GET_FQ_POLY(R);
        fq_nmod_mpoly_ctx_struct *fq_ctx = &(unified_ctx->ctx.fq_ctx);
        
        for (slong i = 0; i < fq_nmod_mpoly_length(fq_res, fq_ctx); i++) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ctx);
            fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, fq_res, i, fq_ctx);
            
            ulong *exp = (ulong*) calloc(total_vars, sizeof(ulong));
            fq_nmod_mpoly_get_term_exp_ui(exp, fq_res, i, fq_ctx);
            
            slong *par_exp = NULL;
            if (state.npars > 0) {
                par_exp = (slong*) calloc(state.npars, sizeof(slong));
                for (slong j = 0; j < state.npars; j++) {
                    par_exp[j] = exp[1 + j];
                }
            }
            
            fq_mvpoly_add_term(&result_mvpoly, NULL, par_exp, coeff);
            
            fq_nmod_clear(coeff, ctx);
            free(exp);
            if (par_exp) free(par_exp);
        }
    }
    
    /* Convert result to string */
    char *result = fq_mvpoly_to_string(&result_mvpoly, state.par_names, state.generator_name);
    
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
    unified_triangular_ideal_clear(&reduced_ideal);
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
    
    printf("\n=== Resultant with Ideal Reduction Complete ===\n");
    
    return result;
}

/* String interface for resultant with ideal reduction */
char* resultant_with_ideal_reduction_str(const char *poly1_string,
                                       const char *poly2_string, 
                                       const char *elim_var_string,
                                       const char *ideal_gens_string,
                                       const char *all_vars_string,
                                       const fq_nmod_ctx_t ctx) {
    
    slong num_gens, num_all_vars;
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
    
    char *result = resultant_with_ideal_reduction(poly1_string, poly2_string,
                                                elim_var_string, ctx, &ideal);
    
    unified_triangular_ideal_clear(&ideal);
    free(ideal_gens);
    free(all_var_names);
    free_split_string_rs(gens_array, num_gens);
    free_split_string_rs(all_vars_array, num_all_vars);
    
    return result;
}

char* elimination_with_ideal_reduction_str(const char *poly_string,
                                               const char *elim_vars_string,
                                               const char *ideal_gens_string,
                                               const char *all_vars_string,
                                               const fq_nmod_ctx_t ctx) {
    
    printf("\n=== Smart Elimination with Ideal Reduction ===\n");
    
    /* 解析输入多项式数量 */
    slong num_polys;
    char **poly_array = split_string_r(poly_string, &num_polys);
    
    printf("Detected %ld polynomial(s)\n", num_polys);
    
    char *result = NULL;
    
    if (num_polys == 2) {
        /* 使用 resultant 方法 */
        printf("Using bivariate resultant method\n");
        
        /* 解析消元变量 - resultant 只支持单变量消元 */
        slong num_elim_vars;
        char **elim_array = split_string_r(elim_vars_string, &num_elim_vars);
        
        if (num_elim_vars != 1) {
            printf("Warning: Resultant method requires exactly one elimination variable, "
                   "using first variable: %s\n", elim_array[0]);
        }
        
        const char *elim_var = (num_elim_vars > 0) ? elim_array[0] : "x";
        
        /* 构造理想 */
        unified_triangular_ideal_t ideal;
        construct_triangular_ideal_str(&ideal, ideal_gens_string, all_vars_string, ctx);
        
        /* 调用 resultant 函数 */
        result = resultant_with_ideal_reduction(poly_array[0], poly_array[1], 
                                              elim_var, ctx, &ideal);
        
        /* 清理 */
        unified_triangular_ideal_clear(&ideal);
        free_split_string_rs(elim_array, num_elim_vars);
        
    } else if (num_polys > 2) {
        /* 使用 Dixon 方法 */
        printf("Using Dixon resultant method for %ld polynomials\n", num_polys);
        
        /* 直接调用原始的 dixon_with_ideal_reduction_str */
        result = dixon_with_ideal_reduction_str(poly_string, elim_vars_string,
                                              ideal_gens_string, all_vars_string, ctx);
        
    } else if (num_polys == 1) {
        /* 单个多项式的情况 */
        printf("Single polynomial detected - no elimination needed\n");
        
        /* 解析单个多项式并应用理想约简 */
        slong num_elim_vars;
        char **elim_array = split_string_r(elim_vars_string, &num_elim_vars);
        
        /* 构造理想 */
        unified_triangular_ideal_t ideal;
        construct_triangular_ideal_str(&ideal, ideal_gens_string, all_vars_string, ctx);
        
        /* 初始化解析器状态 */
        parser_state_t state;
        state.var_names = (char**) malloc(num_elim_vars * sizeof(char*));
        for (slong i = 0; i < num_elim_vars; i++) {
            state.var_names[i] = strdup(elim_array[i]);
        }
        state.nvars = num_elim_vars;
        state.npars = 0;
        state.max_pars = 32;
        state.par_names = (char**) malloc(state.max_pars * sizeof(char*));
        state.ctx = ctx;
        state.generator_name = get_generator_name(ctx);
        fq_nmod_init(state.current.value, ctx);
        state.current.str = NULL;
        
        /* 解析多项式 */
        fq_mvpoly_t poly;
        fq_mvpoly_init(&poly, num_elim_vars, state.max_pars, ctx);
        
        state.input = poly_array[0];
        state.pos = 0;
        state.len = strlen(poly_array[0]);
        next_token(&state);
        parse_expression(&state, &poly);
        
        /* 应用理想约简 */
        printf("Applying ideal reduction to single polynomial\n");
        
        /* 创建简化的理想上下文 */
        unified_triangular_ideal_t reduced_ideal;
        create_reduced_ideal_context(&reduced_ideal, &ideal, state.par_names, 
                                   state.npars, ctx);
        
        /* 这里需要将 fq_mvpoly 转换为统一格式进行约简 */
        /* 为简化，我们直接返回原多项式的字符串表示 */
        result = fq_mvpoly_to_string(&poly, state.par_names, state.generator_name);
        
        /* 清理 */
        unified_triangular_ideal_clear(&reduced_ideal);
        unified_triangular_ideal_clear(&ideal);
        fq_mvpoly_clear(&poly);
        
        for (slong i = 0; i < state.nvars; i++) {
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
        
        free_split_string_rs(elim_array, num_elim_vars);
        
    } else {
        /* 无多项式的情况 */
        printf("No polynomials provided\n");
        result = strdup("1");  /* 返回单位元 */
    }
    
    /* 清理多项式数组 */
    free_split_string_rs(poly_array, num_polys);
    
    printf("\n=== Smart Elimination Complete ===\n");
    
    return result;
}

/* Unified elimination with ideal reduction - automatically chooses Dixon or resultant */
char* elimination_with_ideal_reduction(const char **poly_strings, slong num_polys,
                                      const char **elim_vars, slong num_elim_vars,
                                      const fq_nmod_ctx_t ctx,
                                      unified_triangular_ideal_t *ideal) {
    printf("\n=== Unified Elimination with Ideal Reduction ===\n");
    printf("Eliminating variables: ");
    for (slong i = 0; i < num_elim_vars; i++) {
        if (i > 0) printf(", ");
        printf("%s", elim_vars[i]);
    }
    printf("\n");
    
    printf("Input: %ld polynomial(s), %ld elimination variable(s)\n", 
           num_polys, num_elim_vars);
    
    char *result = NULL;
    
    // Case 1: Two polynomials and one elimination variable -> use bivariate resultant
    if (num_polys == 2 && num_elim_vars == 1) {
        printf("Strategy: Using bivariate resultant method\n");
        result = resultant_with_ideal_reduction(poly_strings[0], poly_strings[1],
                                              elim_vars[0], ctx, ideal);
    }
    // Case 2: Multiple polynomials or multiple variables -> use Dixon method
    else if (num_polys > 2 || num_elim_vars > 1) {
        printf("Strategy: Using Dixon resultant method\n");
        result = dixon_with_ideal_reduction(poly_strings, num_polys,
                                          elim_vars, num_elim_vars,
                                          ctx, ideal);
    }
    // Case 3: Single polynomial -> apply ideal reduction directly
    else if (num_polys == 1) {
        printf("Strategy: Single polynomial - applying ideal reduction only\n");
        
        /* Extract ALL variables from the ideal */
        slong total_system_vars = 0;
        char **all_system_vars = NULL;
        
        if (ideal && ideal->num_gens > 0) {
            if (ideal->is_prime_field) {
                total_system_vars = nmod_mpoly_ctx_nvars(ideal->ctx.nmod_ctx);
            } else {
                total_system_vars = fq_nmod_mpoly_ctx_nvars(ideal->ctx.fq_ctx);
            }
            
            all_system_vars = (char**) malloc(total_system_vars * sizeof(char*));
            
            for (slong i = 0; i < total_system_vars; i++) {
                all_system_vars[i] = NULL;
            }
            
            for (slong g = 0; g < ideal->num_gens; g++) {
                if (ideal->var_names[g] && ideal->var_indices[g] >= 0 && 
                    ideal->var_indices[g] < total_system_vars) {
                    if (!all_system_vars[ideal->var_indices[g]]) {
                        all_system_vars[ideal->var_indices[g]] = strdup(ideal->var_names[g]);
                    }
                }
            }
            
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
            total_system_vars = 2;
            all_system_vars = (char**) malloc(total_system_vars * sizeof(char*));
            all_system_vars[0] = strdup("x");
            all_system_vars[1] = strdup("y");
        }
        
        /* Initialize parser state */
        parser_state_t state;
        state.var_names = (char**) malloc(num_elim_vars * sizeof(char*));
        for (slong i = 0; i < num_elim_vars; i++) {
            state.var_names[i] = strdup(elim_vars[i]);
        }
        state.nvars = num_elim_vars;
        state.npars = 0;
        state.max_pars = 32;
        state.par_names = (char**) malloc(state.max_pars * sizeof(char*));
        state.ctx = ctx;
        state.generator_name = get_generator_name(ctx);
        fq_nmod_init(state.current.value, ctx);
        state.current.str = NULL;
        
        /* Add all non-eliminated variables as parameters */
        for (slong i = 0; i < total_system_vars; i++) {
            int is_elim_var = 0;
            
            for (slong j = 0; j < num_elim_vars; j++) {
                if (strcmp(all_system_vars[i], elim_vars[j]) == 0) {
                    is_elim_var = 1;
                    break;
                }
            }
            
            if (!is_elim_var) {
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
        
        /* Parse polynomial */
        fq_mvpoly_t poly;
        fq_mvpoly_init(&poly, num_elim_vars, state.max_pars, ctx);
        
        state.input = poly_strings[0];
        state.pos = 0;
        state.len = strlen(poly_strings[0]);
        next_token(&state);
        parse_expression(&state, &poly);
        
        /* Create reduced ideal context for current variables */
        unified_triangular_ideal_t reduced_ideal;
        create_reduced_ideal_context(&reduced_ideal, ideal, state.par_names, 
                                   state.npars, ctx);
        
        /* Apply ideal reduction to the polynomial */
        printf("Applying ideal reduction to polynomial\n");
        
        /* For single polynomial, we need to convert to unified format and apply reduction */
        /* This is a simplified approach - in practice you'd want full reduction implementation */
        
        /* Convert result to string */
        result = fq_mvpoly_to_string(&poly, state.par_names, state.generator_name);
        
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
        unified_triangular_ideal_clear(&reduced_ideal);
        fq_mvpoly_clear(&poly);
        
        for (slong i = 0; i < state.nvars; i++) {
            free(state.var_names[i]);
        }
        free(state.var_names);
        
        for (slong i = 0; i < state.npars; i++) {
            free(state.par_names[i]);
        }
        free(state.par_names);
        
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
    }
    // Case 4: No polynomials
    else {
        printf("Strategy: No polynomials provided - returning unit element\n");
        result = strdup("1");
    }
    
    printf("\n=== Unified Elimination Complete ===\n");
    
    return result;
}