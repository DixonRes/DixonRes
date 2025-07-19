/*
 * 使用FLINT fq_nmod_mpoly进行多项式运算的递归行列式计算
 * 
 * Enhanced version with polynomial matrix determinant support
 */

#ifndef FQ_MPOLY_MAT_DET_H
#define FQ_MPOLY_MAT_DET_H

#include <flint/fq_nmod_mpoly.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_mat.h>
#include <flint/fq_nmod_poly.h>
#include "fq_poly_mat_det.h" // Include the polynomial matrix determinant function

// 调试开关
#define DEBUG_FQ_DET 0

#if DEBUG_FQ_DET
#define DET_PRINT(fmt, ...) printf("[FQ_DET] " fmt, ##__VA_ARGS__)
#else
#define DET_PRINT(fmt, ...)
#endif


// ============= 转换函数 =============

// 将fq_mvpoly_t转换为fq_nmod_mpoly_t
void fq_mvpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t mpoly, const fq_mvpoly_t *poly, 
                               fq_nmod_mpoly_ctx_t mpoly_ctx) {
    fq_nmod_mpoly_zero(mpoly, mpoly_ctx);
    
    if (poly->nterms == 0) return;
    
    slong total_vars = poly->nvars + poly->npars;
    
    DET_PRINT("Converting fq_mvpoly (%ld terms) to fq_nmod_mpoly (%ld vars)\n", 
              poly->nterms, total_vars);
    
    for (slong i = 0; i < poly->nterms; i++) {
        // 创建组合指数向量（变量 + 参数）
        ulong *exps = (ulong*) flint_calloc(total_vars, sizeof(ulong));
        
        // 复制变量指数
        if (poly->terms[i].var_exp && poly->nvars > 0) {
            for (slong j = 0; j < poly->nvars; j++) {
                exps[j] = (ulong)poly->terms[i].var_exp[j];
            }
        }
        
        // 复制参数指数
        if (poly->terms[i].par_exp && poly->npars > 0) {
            for (slong j = 0; j < poly->npars; j++) {
                exps[poly->nvars + j] = (ulong)poly->terms[i].par_exp[j];
            }
        }
        
        // 添加项到mpoly
        fq_nmod_mpoly_push_term_fq_nmod_ui(mpoly, poly->terms[i].coeff, exps, mpoly_ctx);
        
        flint_free(exps);
    }
    
    // 整理多项式
    fq_nmod_mpoly_sort_terms(mpoly, mpoly_ctx);
    fq_nmod_mpoly_combine_like_terms(mpoly, mpoly_ctx);
    
    DET_PRINT("Converted to fq_nmod_mpoly with %ld terms\n", 
              fq_nmod_mpoly_length(mpoly, mpoly_ctx));
}

void fq_nmod_mpoly_to_fq_mvpoly(fq_mvpoly_t *poly, const fq_nmod_mpoly_t mpoly,
                                               slong nvars, slong npars, 
                                               fq_nmod_mpoly_ctx_t mpoly_ctx, const fq_nmod_ctx_t ctx) {
    fq_mvpoly_init(poly, nvars, npars, ctx);
    
    slong nterms = fq_nmod_mpoly_length(mpoly, mpoly_ctx);
    if (nterms == 0) return;
    
    DET_PRINT("Safe optimized conversion: %ld terms\n", nterms);
    
    slong total_vars = nvars + npars;
    
    // 预分配主结构
    poly->alloc = nterms;
    poly->terms = (fq_monomial_t*) flint_realloc(poly->terms, poly->alloc * sizeof(fq_monomial_t));
    poly->nterms = nterms;
    
    // 关键修复：每个指数向量仍然单独分配，但优化分配过程
    ulong *exp_buffer = (ulong*) flint_malloc(total_vars * sizeof(ulong));
    
    clock_t start_process = clock();
    
    for (slong i = 0; i < nterms; i++) {
        // 初始化系数
        fq_nmod_init(poly->terms[i].coeff, ctx);
        
        // 一次性获取系数和指数
        fq_nmod_mpoly_get_term_coeff_fq_nmod(poly->terms[i].coeff, mpoly, i, mpoly_ctx);
        fq_nmod_mpoly_get_term_exp_ui(exp_buffer, mpoly, i, mpoly_ctx);
        
        // 为每个项单独分配指数向量（保持兼容性）
        if (nvars > 0) {
            poly->terms[i].var_exp = (slong*) flint_calloc(nvars, sizeof(slong));
            for (slong j = 0; j < nvars; j++) {
                poly->terms[i].var_exp[j] = (slong)exp_buffer[j];
            }
        } else {
            poly->terms[i].var_exp = NULL;
        }
        
        if (npars > 0) {
            poly->terms[i].par_exp = (slong*) flint_calloc(npars, sizeof(slong));
            for (slong j = 0; j < npars; j++) {
                poly->terms[i].par_exp[j] = (slong)exp_buffer[nvars + j];
            }
        } else {
            poly->terms[i].par_exp = NULL;
        }
    }
    
    clock_t end_process = clock();
    double process_time = (double)(end_process - start_process) / CLOCKS_PER_SEC;
    
    flint_free(exp_buffer);
    
    DET_PRINT("Safe optimized conversion time: %.6f seconds\n", process_time);
}

// ============= FLINT-based 多项式矩阵行列式计算 =============

// 函数1：将整个fq_mvpoly_t矩阵转换为fq_nmod_mpoly_t矩阵
void fq_matrix_mvpoly_to_mpoly(fq_nmod_mpoly_t **mpoly_matrix, 
                               fq_mvpoly_t **mvpoly_matrix, 
                               slong size, 
                               fq_nmod_mpoly_ctx_t mpoly_ctx) {
    DET_PRINT("Converting %ld x %ld fq_mvpoly matrix to fq_nmod_mpoly format\n", size, size);
    
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            fq_nmod_mpoly_init(mpoly_matrix[i][j], mpoly_ctx);
            fq_mvpoly_to_fq_nmod_mpoly(mpoly_matrix[i][j], &mvpoly_matrix[i][j], mpoly_ctx);
            
            DET_PRINT("  M[%ld][%ld]: %ld terms -> %ld terms\n", 
                      i, j, mvpoly_matrix[i][j].nterms, 
                      fq_nmod_mpoly_length(mpoly_matrix[i][j], mpoly_ctx));
        }
    }
    
    DET_PRINT("Matrix conversion complete\n");
}

// 函数2：在fq_nmod_mpoly_t矩阵上递归计算行列式
void compute_fq_nmod_mpoly_det_recursive(fq_nmod_mpoly_t det_result, 
                                        fq_nmod_mpoly_t **mpoly_matrix, 
                                        slong size, 
                                        fq_nmod_mpoly_ctx_t mpoly_ctx) {
    DET_PRINT("Computing %ld x %ld determinant on fq_nmod_mpoly matrix\n", size, size);
    
    if (size <= 0) {
        fq_nmod_mpoly_one(det_result, mpoly_ctx);
        return;
    }
    
    if (size == 1) {
        fq_nmod_mpoly_set(det_result, mpoly_matrix[0][0], mpoly_ctx);
        DET_PRINT("1x1 result: %ld terms\n", fq_nmod_mpoly_length(det_result, mpoly_ctx));
        return;
    }
    
    if (size == 2) {
        DET_PRINT("2x2: ad - bc using fq_nmod_mpoly operations\n");
        
        fq_nmod_mpoly_t ad, bc;
        fq_nmod_mpoly_init(ad, mpoly_ctx);
        fq_nmod_mpoly_init(bc, mpoly_ctx);
        
        // ad = a * d
        fq_nmod_mpoly_mul(ad, mpoly_matrix[0][0], mpoly_matrix[1][1], mpoly_ctx);
        // bc = b * c
        fq_nmod_mpoly_mul(bc, mpoly_matrix[0][1], mpoly_matrix[1][0], mpoly_ctx);
        // det = ad - bc
        fq_nmod_mpoly_sub(det_result, ad, bc, mpoly_ctx);
        
        DET_PRINT("2x2 result: %ld terms\n", fq_nmod_mpoly_length(det_result, mpoly_ctx));
        
        fq_nmod_mpoly_clear(ad, mpoly_ctx);
        fq_nmod_mpoly_clear(bc, mpoly_ctx);
        return;
    }
    
    // 3x3及更大：使用Laplace展开
    DET_PRINT("Laplace expansion along row 0 using fq_nmod_mpoly operations\n");
    
    fq_nmod_mpoly_zero(det_result, mpoly_ctx);
    
    fq_nmod_mpoly_t temp_result, cofactor, subdet;
    fq_nmod_mpoly_init(temp_result, mpoly_ctx);
    fq_nmod_mpoly_init(cofactor, mpoly_ctx);
    fq_nmod_mpoly_init(subdet, mpoly_ctx);
    
    for (slong col = 0; col < size; col++) {
        DET_PRINT("Column %ld expansion\n", col);
        
        // 跳过零元素
        if (fq_nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
            DET_PRINT("Skipping zero element\n");
            continue;
        }
        
        // 创建子矩阵
        fq_nmod_mpoly_t **submatrix = (fq_nmod_mpoly_t**) flint_malloc((size-1) * sizeof(fq_nmod_mpoly_t*));
        for (slong i = 0; i < size-1; i++) {
            submatrix[i] = (fq_nmod_mpoly_t*) flint_malloc((size-1) * sizeof(fq_nmod_mpoly_t));
            for (slong j = 0; j < size-1; j++) {
                fq_nmod_mpoly_init(submatrix[i][j], mpoly_ctx);
            }
        }
        
        // 填充子矩阵
        for (slong i = 1; i < size; i++) {
            slong sub_j = 0;
            for (slong j = 0; j < size; j++) {
                if (j != col) {
                    fq_nmod_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j], mpoly_ctx);
                    sub_j++;
                }
            }
        }
        
        // 递归计算子矩阵行列式
        compute_fq_nmod_mpoly_det_recursive(subdet, submatrix, size-1, mpoly_ctx);
        
        // 计算余子式：element * subdet
        fq_nmod_mpoly_mul(cofactor, mpoly_matrix[0][col], subdet, mpoly_ctx);
        
        // 应用符号并累加到结果
        if (col % 2 == 0) {
            fq_nmod_mpoly_add(temp_result, det_result, cofactor, mpoly_ctx);
        } else {
            fq_nmod_mpoly_sub(temp_result, det_result, cofactor, mpoly_ctx);
        }
        
        // 更新结果
        fq_nmod_mpoly_set(det_result, temp_result, mpoly_ctx);
        
        DET_PRINT("After column %ld: %ld terms\n", col, 
                  fq_nmod_mpoly_length(det_result, mpoly_ctx));
        
        // 清理子矩阵
        for (slong i = 0; i < size-1; i++) {
            for (slong j = 0; j < size-1; j++) {
                fq_nmod_mpoly_clear(submatrix[i][j], mpoly_ctx);
            }
            flint_free(submatrix[i]);
        }
        flint_free(submatrix);
    }
    
    DET_PRINT("Final result: %ld terms\n", fq_nmod_mpoly_length(det_result, mpoly_ctx));
    
    // 清理临时变量
    fq_nmod_mpoly_clear(temp_result, mpoly_ctx);
    fq_nmod_mpoly_clear(cofactor, mpoly_ctx);
    fq_nmod_mpoly_clear(subdet, mpoly_ctx);
}

void compute_fq_det_recursive_flint(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size) {
    if (size <= 0) {
        fq_mvpoly_init(result, matrix[0][0].nvars, matrix[0][0].npars, matrix[0][0].ctx);
        return;
    }
    
    DET_PRINT("Computing %ldx%ld determinant with FLINT (v2: batch conversion)\n", size, size);
    
    // 获取变量信息
    slong nvars = matrix[0][0].nvars;
    slong npars = matrix[0][0].npars;
    slong total_vars = nvars + npars;
    const fq_nmod_ctx_struct *ctx = matrix[0][0].ctx;
    
    fq_mvpoly_init(result, nvars, npars, ctx);
    
    // 创建fq_nmod_mpoly上下文
    fq_nmod_mpoly_ctx_t mpoly_ctx;
    fq_nmod_mpoly_ctx_init(mpoly_ctx, total_vars, ORD_LEX, ctx);
    
    // 分配fq_nmod_mpoly矩阵
    fq_nmod_mpoly_t **mpoly_matrix = (fq_nmod_mpoly_t**) flint_malloc(size * sizeof(fq_nmod_mpoly_t*));
    for (slong i = 0; i < size; i++) {
        mpoly_matrix[i] = (fq_nmod_mpoly_t*) flint_malloc(size * sizeof(fq_nmod_mpoly_t));
    }
    
    // 步骤1：批量转换整个矩阵
    clock_t conv_start = clock();
    fq_matrix_mvpoly_to_mpoly(mpoly_matrix, matrix, size, mpoly_ctx);
    clock_t conv_end = clock();
    double conv_time = (double)(conv_end - conv_start) / CLOCKS_PER_SEC;
    
    DET_PRINT("Matrix conversion time: %.6f seconds\n", conv_time);
    
    // 步骤2：在fq_nmod_mpoly矩阵上计算行列式
    fq_nmod_mpoly_t det_mpoly;
    fq_nmod_mpoly_init(det_mpoly, mpoly_ctx);
    
    clock_t det_start = clock();
    compute_fq_nmod_mpoly_det_recursive(det_mpoly, mpoly_matrix, size, mpoly_ctx);
    clock_t det_end = clock();
    double det_time = (double)(det_end - det_start) / CLOCKS_PER_SEC;
    
    DET_PRINT("Determinant computation time: %.6f seconds\n", det_time);
    
    // 步骤3：转换结果回fq_mvpoly格式
    clock_t result_start = clock();
    fq_nmod_mpoly_to_fq_mvpoly(result, det_mpoly, nvars, npars, mpoly_ctx, ctx);
    clock_t result_end = clock();
    double result_time = (double)(result_end - result_start) / CLOCKS_PER_SEC;
    
    DET_PRINT("Result conversion time: %.6f seconds\n", result_time);
    DET_PRINT("Total time: %.6f seconds\n", conv_time + det_time + result_time);
    DET_PRINT("Final result: %ld terms\n", result->nterms);
    
    // 清理fq_nmod_mpoly矩阵
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            fq_nmod_mpoly_clear(mpoly_matrix[i][j], mpoly_ctx);
        }
        flint_free(mpoly_matrix[i]);
    }
    flint_free(mpoly_matrix);
    
    // 清理其他临时对象
    fq_nmod_mpoly_clear(det_mpoly, mpoly_ctx);
    fq_nmod_mpoly_ctx_clear(mpoly_ctx);
}

// ============= 多项式矩阵行列式计算（使用fq_nmod_poly_mat_det_iter）=============

void compute_fq_det_polynomial_matrix_simple(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size) {
    if (size <= 0) {
        fq_mvpoly_init(result, matrix[0][0].nvars, matrix[0][0].npars, matrix[0][0].ctx);
        return;
    }
    
    DET_PRINT("Computing %ldx%ld determinant using simple polynomial matrix method\n", size, size);
    
    // Get context information
    slong nvars = matrix[0][0].nvars;
    slong npars = matrix[0][0].npars;
    const fq_nmod_ctx_struct *ctx = matrix[0][0].ctx;
    
    fq_mvpoly_init(result, nvars, npars, ctx);
    
    // Create fq_nmod_poly_mat
    fq_nmod_poly_mat_t poly_mat;
    fq_nmod_poly_mat_init(poly_mat, size, size, ctx);
    
    // Convert matrix: for univariate case, extract first variable as polynomial
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            fq_nmod_poly_struct *entry = fq_nmod_poly_mat_entry(poly_mat, i, j);
            fq_nmod_poly_zero(entry, ctx);
            
            // Convert mvpoly to polynomial (use first variable)
            for (slong k = 0; k < matrix[i][j].nterms; k++) {
                fq_monomial_t *term = &matrix[i][j].terms[k];
                slong degree = 0;
                if (term->var_exp && nvars > 0) {
                    degree = term->var_exp[0];
                }
                fq_nmod_poly_set_coeff(entry, degree, term->coeff, ctx);
            }
        }
    }
    
    // Compute determinant using polynomial matrix method
    fq_nmod_poly_t det_poly;
    fq_nmod_poly_init(det_poly, ctx);
    fq_nmod_poly_mat_det_iter(det_poly, poly_mat, ctx);
    
    // Convert result back to fq_mvpoly format
    slong degree = fq_nmod_poly_degree(det_poly, ctx);
    if (degree >= 0) {
        for (slong d = 0; d <= degree; d++) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ctx);
            fq_nmod_poly_get_coeff(coeff, det_poly, d, ctx);
            
            if (!fq_nmod_is_zero(coeff, ctx)) {
                slong *var_exp = NULL;
                slong *par_exp = NULL;
                
                if (nvars > 0) {
                    var_exp = (slong*) flint_calloc(nvars, sizeof(slong));
                    var_exp[0] = d;
                }
                if (npars > 0) {
                    par_exp = (slong*) flint_calloc(npars, sizeof(slong));
                }
                
                fq_mvpoly_add_term(result, var_exp, par_exp, coeff);
                
                if (var_exp) flint_free(var_exp);
                if (par_exp) flint_free(par_exp);
            }
            
            fq_nmod_clear(coeff, ctx);
        }
    }
    
    DET_PRINT("Polynomial matrix determinant completed: %ld terms\n", result->nterms);
    
    // Cleanup
    fq_nmod_poly_clear(det_poly, ctx);
    fq_nmod_poly_mat_clear(poly_mat, ctx);
}

// ============= 替换原有函数的接口 =============

// 替换原来的compute_fq_det_recursive函数
void compute_fq_det_recursive(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size) {
    DET_PRINT("Using FLINT-based determinant computation\n");
    
    // For univariate polynomial matrices, try the polynomial matrix method
    if (matrix[0][0].nvars == 1 && matrix[0][0].npars == 0 && size >= 3) {
        DET_PRINT("Trying polynomial matrix method for univariate case\n");
        compute_fq_det_polynomial_matrix_simple(result, matrix, size);
    } else {
        DET_PRINT("Using standard FLINT mpoly method\n");
        compute_fq_det_recursive_flint(result, matrix, size);
    }
}

#endif // FQ_DET_RECURSIVE_FLINT_H