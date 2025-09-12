#ifndef DIXON_FLINT_H
#define DIXON_FLINT_H
/* dixon_flint.h gcc -O3  -fsanitize=address -static-libasan  -march=native -o rescue rescue_attack_dixon.c -lflint -lmpfr -lgmp -lpthread -L/home/suohaohai02/mylinks -lflint -lstdc++ -lpml2 -fopenmp
 * Complete Dixon Resultant Implementation for Finite Extension Fields
 * 
 * Fixed compilation issues:
 * - Context pointer handling
 * - Function declaration order
 * - Assignment to array types
 * 
 * Compile with: gcc -O3 -march=native -o dixon_flint dixon_flint.c -lflint -lmpfr -lgmp -lpthread -L/home/suohaohai02/mylinks -lflint -lstdc++ -lpml2 -fopenmp -fsanitize=address -static-libasan 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include <unistd.h>
#include <ctype.h>

#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_poly.h>
#include <flint/fq_nmod_mat.h>
#include <flint/fq_nmod_mpoly.h>
#include <flint/nmod_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/profiler.h>

#include "fq_mvpoly.h"
#include "fq_mpoly_mat_det.h"
#include "fq_unified_interface.h"
#include "dixon_interface_flint.h"

/* Debug output control macro - set to 0 to disable all output */
#define DEBUG_OUTPUT_D 0

#if DEBUG_OUTPUT_D
    #define DEBUG_PRINT_D(...) printf(__VA_ARGS__)
#else
    #define DEBUG_PRINT_D(...) ((void)0)
#endif

// Method enumeration for determinant computation
typedef enum {
    DET_METHOD_RECURSIVE = 0,
    DET_METHOD_KRONECKER = 1,
    DET_METHOD_INTERPOLATION = 2,
    DET_METHOD_HUANG = 3
} det_method_t;

// ============ Matrix operations ============

// Build cancellation matrix in multivariate form
void build_fq_cancellation_matrix_mvpoly(fq_mvpoly_t ***M, fq_mvpoly_t *polys, 
                                        slong nvars, slong npars) {
    slong n = nvars + 1;
    
    // Allocate matrix space
    *M = (fq_mvpoly_t**) flint_malloc(n * sizeof(fq_mvpoly_t*));
    for (slong i = 0; i < n; i++) {
        (*M)[i] = (fq_mvpoly_t*) flint_malloc(n * sizeof(fq_mvpoly_t));
    }
    
    printf("Building %ld x %ld cancellation matrix (multivariate form)\n", n, n);
    
    // Build matrix entries
    for (slong i = 0; i < n; i++) {
        for (slong j = 0; j < n; j++) {
            fq_mvpoly_init(&(*M)[i][j], 2 * nvars, npars, polys[0].ctx);
            
            // Substitute variables according to row index
            for (slong t = 0; t < polys[j].nterms; t++) {
                slong *new_var_exp = (slong*) flint_calloc(2 * nvars, sizeof(slong));
                
                for (slong k = 0; k < nvars; k++) {
                    slong orig_exp = polys[j].terms[t].var_exp ? polys[j].terms[t].var_exp[k] : 0;
                    
                    if (k < i) {
                        // Use dual variable ~x_k
                        new_var_exp[nvars + k] = orig_exp;
                    } else {
                        // Use original variable x_k
                        new_var_exp[k] = orig_exp;
                    }
                }
                
                fq_mvpoly_add_term_fast(&(*M)[i][j], new_var_exp, polys[j].terms[t].par_exp, 
                                  polys[j].terms[t].coeff);
                flint_free(new_var_exp);
            }
        }
    }
}

void perform_fq_matrix_row_operations_mvpoly(fq_mvpoly_t ***new_matrix, fq_mvpoly_t ***original_matrix,
                                                   slong nvars, slong npars) {
    slong n = nvars + 1;
    
    // 分配新矩阵空间
    *new_matrix = (fq_mvpoly_t**) flint_malloc(n * sizeof(fq_mvpoly_t*));
    for (slong i = 0; i < n; i++) {
        (*new_matrix)[i] = (fq_mvpoly_t*) flint_malloc(n * sizeof(fq_mvpoly_t));
    }
    
    printf("Performing row operations on matrix:\n");
    
    // 第一行保持不变
    printf("Row 0: unchanged\n");
    for (slong j = 0; j < n; j++) {
        fq_mvpoly_copy(&(*new_matrix)[0][j], &(*original_matrix)[0][j]);
    }
    
    // 处理其余行
    for (slong i = 1; i < n; i++) {
        printf("Row %ld: (row[%ld] - row[%ld]) / (x%ld - ~x%ld)\n", i, i, i-1, i-1, i-1);
        
        for (slong j = 0; j < n; j++) {
            // 直接使用 fq_mvpoly_sub 计算差值
            fq_mvpoly_t diff;
            fq_mvpoly_sub(&diff, &(*original_matrix)[i][j], &(*original_matrix)[i-1][j]);
            
            // 执行除法
            fq_mvpoly_init(&(*new_matrix)[i][j], 2*nvars, npars, diff.ctx);
            
            if (diff.nterms > 0) {
                divide_by_fq_linear_factor_flint(&(*new_matrix)[i][j], &diff, 
                                               i-1, 2*nvars, npars);
            }
            
            fq_mvpoly_clear(&diff);
        }
    }
}


// Compute Dixon resultant degree bound
slong compute_fq_dixon_resultant_degree_bound(fq_mvpoly_t *polys, slong npolys, slong nvars, slong npars) {
    slong degree_product = 1;
    
    for (slong i = 0; i < npolys; i++) {
        slong max_total_deg = 0;
        
        // Find maximum total degree of polynomial i
        for (slong t = 0; t < polys[i].nterms; t++) {
            slong total_deg = 0;
            
            // Sum variable degrees
            if (polys[i].terms[t].var_exp) {
                for (slong j = 0; j < nvars; j++) {
                    total_deg += polys[i].terms[t].var_exp[j];
                }
            }
            
            // Sum parameter degrees
            if (polys[i].terms[t].par_exp && npars > 0) {
                for (slong j = 0; j < npars; j++) {
                    total_deg += polys[i].terms[t].par_exp[j];
                }
            }
            
            if (total_deg > max_total_deg) {
                max_total_deg = total_deg;
            }
        }
        
        degree_product *= max_total_deg;
    }
    
    return degree_product + 1;
}

void compute_fq_coefficient_matrix_det(fq_mvpoly_t *result, fq_mvpoly_t **coeff_matrix,
                                       slong size, slong npars, const fq_nmod_ctx_t ctx,
                                       det_method_t method, slong res_deg_bound) {
    if (size == 0) {
        fq_mvpoly_init(result, 0, npars, ctx);
        return;
    }
    
    fq_mvpoly_init(result, 0, npars, ctx);
    
    if (npars == 0) {
        // No parameters - scalar entries
        fq_nmod_mat_t scalar_mat;
        fq_nmod_mat_init(scalar_mat, size, size, ctx);
        
        for (slong i = 0; i < size; i++) {
            for (slong j = 0; j < size; j++) {
                if (coeff_matrix[i][j].nterms > 0) {
                    fq_nmod_set(fq_nmod_mat_entry(scalar_mat, i, j), 
                                coeff_matrix[i][j].terms[0].coeff, ctx);
                } else {
                    fq_nmod_zero(fq_nmod_mat_entry(scalar_mat, i, j), ctx);
                }
            }
        }
        
        printf("\nComputing Resultant\n");
        clock_t start = clock();
        
        fq_nmod_t det;
        fq_nmod_init(det, ctx);
        fq_nmod_mat_det(det, scalar_mat, ctx);
        
        clock_t end = clock();
        double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
        
        printf("End (%.3f seconds)\n", elapsed);
        
        if (!fq_nmod_is_zero(det, ctx)) {
            fq_mvpoly_add_term_fast(result, NULL, NULL, det);
        }
        
        fq_nmod_clear(det, ctx);
        fq_nmod_mat_clear(scalar_mat, ctx);
        
    } else if (npars == 1) {
        // One parameter - use Mulders-Storjohann algorithm
        printf("\nComputing Resultant using Mulders-Storjohann algorithm (univariate case)\n");
        clock_t start = clock();
        if (method == DET_METHOD_INTERPOLATION) {
            printf("\nMultiple parameters detected. Using fq_nmod interpolation method.\n");
            printf("  Parameters: %ld\n", npars);
            printf("  Matrix size: %ld x %ld\n", size, size);
            
            // Use the interpolation method
            fq_compute_det_by_interpolation(result, coeff_matrix, size,
                                           0, npars, ctx, res_deg_bound);
        }
        
        else {
            // Create polynomial matrix using our fq_nmod_poly_mat_t structure
            fq_nmod_poly_mat_t poly_mat;
            fq_nmod_poly_mat_init(poly_mat, size, size, ctx);
            
            // Convert mvpoly coefficients to polynomial matrix
            for (slong i = 0; i < size; i++) {
                for (slong j = 0; j < size; j++) {
                    fq_nmod_poly_struct *entry = fq_nmod_poly_mat_entry(poly_mat, i, j);
                    fq_nmod_poly_zero(entry, ctx);
                    
                    // Convert mvpoly to fq_nmod_poly (using parameter as polynomial variable)
                    for (slong t = 0; t < coeff_matrix[i][j].nterms; t++) {
                        slong deg = coeff_matrix[i][j].terms[t].par_exp ? 
                                   coeff_matrix[i][j].terms[t].par_exp[0] : 0;
                        fq_nmod_poly_set_coeff(entry, deg, 
                                              coeff_matrix[i][j].terms[t].coeff, ctx);
                    }
                }
            }
            
            // Compute determinant using Mulders-Storjohann algorithm
            fq_nmod_poly_t det_poly;
            fq_nmod_poly_init(det_poly, ctx);
            
            printf("  Matrix size: %ld x %ld\n", size, size);
            printf("  Using weak Popov form method...\n");
            
            fq_nmod_poly_mat_det_iter(det_poly, poly_mat, ctx);
            
            // Convert result back to fq_mvpoly
            slong det_deg = fq_nmod_poly_degree(det_poly, ctx);
            if (det_deg >= 0) {
                for (slong i = 0; i <= det_deg; i++) {
                    fq_nmod_t coeff;
                    fq_nmod_init(coeff, ctx);
                    fq_nmod_poly_get_coeff(coeff, det_poly, i, ctx);
                    if (!fq_nmod_is_zero(coeff, ctx)) {
                        slong par_exp[1] = {i};
                        fq_mvpoly_add_term_fast(result, NULL, par_exp, coeff);
                    }
                    fq_nmod_clear(coeff, ctx);
                }
            }
            
            printf("  Determinant degree: %ld\n", det_deg);
            printf("  Result terms: %ld\n", result->nterms);
            
            // Cleanup
            fq_nmod_poly_clear(det_poly, ctx);
            fq_nmod_poly_mat_clear(poly_mat, ctx);
        }

        
        clock_t end = clock();
        double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
        printf("End (%.3f seconds)\n", elapsed);
        
    } else {
        clock_t start = clock();
        switch (method) {
            case DET_METHOD_INTERPOLATION:
                printf("\nMultiple parameters detected. Using fq_nmod interpolation method.\n");
                printf("  Parameters: %ld\n", npars);
                printf("  Matrix size: %ld x %ld\n", size, size);
                
                // Use the interpolation method
                fq_compute_det_by_interpolation(result, coeff_matrix, size,
                                               0, npars, ctx, res_deg_bound);
                break;
                
            case DET_METHOD_RECURSIVE:
                printf("\nMultiple parameters detected. Using recursive expansion method.\n");
                printf("  Parameters: %ld\n", npars);
                printf("  Matrix size: %ld x %ld\n", size, size);
                
                // Use the recursive algorithm
                compute_fq_det_recursive(result, coeff_matrix, size);
                break;
                
            case DET_METHOD_KRONECKER:
                printf("\nMultiple parameters detected. Using Kronecker substitution method.\n");
                printf("  Parameters: %ld\n", npars);
                printf("  Matrix size: %ld x %ld\n", size, size);
                
                // Use the Kronecker algorithm
                compute_fq_det_kronecker(result, coeff_matrix, size);
                break;

            case DET_METHOD_HUANG:
                printf("\nMultiple parameters detected. Using Huang's interpolation method.\n");
                printf("  Parameters: %ld\n", npars);
                printf("  Matrix size: %ld x %ld\n", size, size);
                
                // Use the Kronecker algorithm
                compute_fq_det_huang_interpolation(result, coeff_matrix, size);
                break;
                
            default:
                printf("\nWarning: Unknown method %d, defaulting to interpolation.\n", method);
                fq_compute_det_by_interpolation(result, coeff_matrix, size,
                                               0, npars, ctx, res_deg_bound);
                break;
        }
        clock_t end = clock();
        double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
        printf("End (%.3f seconds)\n", elapsed);
    }
}



// ============ Submatrix finding functions ============


// Helper function to compute total degree of a polynomial
static slong compute_fq_polynomial_total_degree_old(fq_mvpoly_t *poly, slong npars) {
    if (poly == NULL || poly->nterms == 0) {
        return 0;
    }
    
    slong max_degree = -1;
    
    for (slong t = 0; t < poly->nterms; t++) {
        slong term_total_degree = 0;
        
        // Variable degrees
        if (poly->terms[t].var_exp) {
            for (slong v = 0; v < poly->nvars; v++) {
                term_total_degree += poly->terms[t].var_exp[v];
            }
        }
        
        // Parameter degrees
        if (poly->terms[t].par_exp && npars > 0) {
            for (slong p = 0; p < npars; p++) {
                term_total_degree += poly->terms[t].par_exp[p];
            }
        }
        
        if (max_degree == -1 || term_total_degree > max_degree) {
            max_degree = term_total_degree;
        }
    }
    
    return max_degree;
}
// Helper function to compute maximum degree (not total degree) of a polynomial
static slong compute_fq_polynomial_total_degree(fq_mvpoly_t *poly, slong npars) {
    if (poly == NULL || poly->nterms == 0) {
        return 0;
    }
    
    slong max_degree = 0;
    
    for (slong t = 0; t < poly->nterms; t++) {
        // Check variable degrees
        if (poly->terms[t].var_exp) {
            for (slong v = 0; v < poly->nvars; v++) {
                if (poly->terms[t].var_exp[v] > max_degree) {
                    max_degree = poly->terms[t].var_exp[v];
                }
            }
        }
        
        // Check parameter degrees
        if (poly->terms[t].par_exp && npars > 0) {
            for (slong p = 0; p < npars; p++) {
                if (poly->terms[t].par_exp[p] > max_degree) {
                    max_degree = poly->terms[t].par_exp[p];
                }
            }
        }
    }
    
    return max_degree;
}

// Compute row maximum total degree
static slong compute_fq_row_max_total_degree(fq_mvpoly_t **matrix_row, slong ncols, slong npars) {
    slong max_degree = -1;
    
    for (slong j = 0; j < ncols; j++) {
        slong poly_deg = compute_fq_polynomial_total_degree(matrix_row[j], npars);
        if (poly_deg > max_degree) {
            max_degree = poly_deg;
        }
    }
    
    return max_degree;
}

// Compute column maximum total degree
static slong compute_fq_col_max_total_degree(fq_mvpoly_t ***matrix, slong col_idx, slong nrows, slong npars) {
    slong max_degree = -1;
    
    for (slong i = 0; i < nrows; i++) {
        slong poly_deg = compute_fq_polynomial_total_degree(matrix[i][col_idx], npars);
        if (poly_deg > max_degree) {
            max_degree = poly_deg;
        }
    }
    
    return max_degree;
}


// 扩展的跟踪器结构，增加预分配的工作空间
typedef struct {
    field_elem_u *reduced_rows;    // 已约化的行向量
    slong *pivot_cols;             // 每行的主元列位置  
    slong *selected_indices;       // 已选择的原始行索引
    slong current_rank;            // 当前秩
    slong max_size;
    slong ncols;
    field_ctx_t *ctx;              // 统一字段上下文
    int initialized;               // 初始化标志
    
    // 预分配的工作空间，避免重复分配
    field_elem_u *work_row;        // 工作行
    field_elem_u *temp_vars;       // 临时变量池：[factor, temp, pivot_val, neg_temp]
    int workspace_initialized;     // 工作空间初始化标志
} unified_row_basis_tracker_t;

// 初始化优化的跟踪器
static void unified_row_basis_tracker_init(unified_row_basis_tracker_t *tracker, 
                                            slong max_size, slong ncols, 
                                            field_ctx_t *ctx) {
    tracker->max_size = max_size;
    tracker->ncols = ncols;
    tracker->ctx = ctx;
    tracker->current_rank = 0;
    tracker->initialized = 1;
    tracker->workspace_initialized = 0;
    
    void *ctx_ptr = (ctx->field_id == FIELD_ID_NMOD) ? 
                   (void*)&ctx->ctx.nmod_ctx : 
                   (void*)ctx->ctx.fq_ctx;
    
    // 分配主要存储
    tracker->reduced_rows = (field_elem_u*) flint_calloc(max_size * ncols, sizeof(field_elem_u));
    tracker->pivot_cols = (slong*) flint_calloc(max_size, sizeof(slong));
    tracker->selected_indices = (slong*) flint_calloc(max_size, sizeof(slong));
    
    // 预分配工作空间
    tracker->work_row = (field_elem_u*) flint_malloc(ncols * sizeof(field_elem_u));
    tracker->temp_vars = (field_elem_u*) flint_malloc(4 * sizeof(field_elem_u)); // factor, temp, pivot_val, neg_temp
    
    // 初始化所有字段元素
    for (slong i = 0; i < max_size * ncols; i++) {
        field_init_elem(&tracker->reduced_rows[i], ctx->field_id, ctx_ptr);
        field_set_zero(&tracker->reduced_rows[i], ctx->field_id, ctx_ptr);
    }
    
    // 初始化工作空间
    for (slong j = 0; j < ncols; j++) {
        field_init_elem(&tracker->work_row[j], ctx->field_id, ctx_ptr);
    }
    for (slong i = 0; i < 4; i++) {
        field_init_elem(&tracker->temp_vars[i], ctx->field_id, ctx_ptr);
    }
    tracker->workspace_initialized = 1;
    
    // 初始化主元列为-1
    for (slong i = 0; i < max_size; i++) {
        tracker->pivot_cols[i] = -1;
        tracker->selected_indices[i] = -1;
    }
}

// 清理优化的跟踪器
static void unified_row_basis_tracker_clear(unified_row_basis_tracker_t *tracker) {
    if (!tracker->initialized) return;
    
    void *ctx_ptr = (tracker->ctx->field_id == FIELD_ID_NMOD) ? 
                   (void*)&tracker->ctx->ctx.nmod_ctx : 
                   (void*)tracker->ctx->ctx.fq_ctx;
    
    // 清理主要存储
    if (tracker->reduced_rows) {
        for (slong i = 0; i < tracker->max_size * tracker->ncols; i++) {
            field_clear_elem(&tracker->reduced_rows[i], tracker->ctx->field_id, ctx_ptr);
        }
        flint_free(tracker->reduced_rows);
    }
    
    // 清理工作空间
    if (tracker->workspace_initialized) {
        for (slong j = 0; j < tracker->ncols; j++) {
            field_clear_elem(&tracker->work_row[j], tracker->ctx->field_id, ctx_ptr);
        }
        for (slong i = 0; i < 4; i++) {
            field_clear_elem(&tracker->temp_vars[i], tracker->ctx->field_id, ctx_ptr);
        }
        flint_free(tracker->work_row);
        flint_free(tracker->temp_vars);
    }
    
    if (tracker->pivot_cols) flint_free(tracker->pivot_cols);
    if (tracker->selected_indices) flint_free(tracker->selected_indices);
    
    tracker->initialized = 0;
    tracker->workspace_initialized = 0;
}

// 优化版本的添加行到基函数 - 数学逻辑完全不变
static int unified_try_add_row_to_basis(unified_row_basis_tracker_t *tracker, 
                                         const field_elem_u *unified_mat,
                                         slong new_row_idx, slong ncols) {
    if (!tracker->initialized || tracker->current_rank >= tracker->max_size) {
        return 0;
    }
    
    void *ctx_ptr = (tracker->ctx->field_id == FIELD_ID_NMOD) ? 
                   (void*)&tracker->ctx->ctx.nmod_ctx : 
                   (void*)tracker->ctx->ctx.fq_ctx;
    
    // 使用预分配的工作行，避免每次分配
    field_elem_u *work_row = tracker->work_row;
    
    // 复制输入行到工作行（重用已分配的空间）
    for (slong j = 0; j < ncols; j++) {
        field_set_elem(&work_row[j], &unified_mat[new_row_idx * ncols + j], 
                      tracker->ctx->field_id, ctx_ptr);
    }
    
    // 使用预分配的临时变量
    field_elem_u *factor = &tracker->temp_vars[0];
    field_elem_u *temp = &tracker->temp_vars[1];
    field_elem_u *pivot_val = &tracker->temp_vars[2];
    field_elem_u *neg_temp = &tracker->temp_vars[3];
    
    // 对已有的每个基向量进行消元（数学逻辑保持不变）
    for (slong i = 0; i < tracker->current_rank; i++) {
        slong pivot_col = tracker->pivot_cols[i];
        
        // 验证主元列索引
        if (pivot_col < 0 || pivot_col >= ncols) continue;
        
        // 如果工作行在主元位置非零，进行消元
        if (!field_is_zero(&work_row[pivot_col], tracker->ctx->field_id, ctx_ptr)) {
            // 获取基向量的主元值
            slong base_idx = i * ncols + pivot_col;
            field_set_elem(pivot_val, &tracker->reduced_rows[base_idx], 
                          tracker->ctx->field_id, ctx_ptr);
            
            // 计算消元因子 = work_row[pivot_col] / pivot_val
            field_inv(temp, pivot_val, tracker->ctx->field_id, ctx_ptr);
            field_mul(factor, &work_row[pivot_col], temp, tracker->ctx->field_id, ctx_ptr);
            
            // 执行消元：work_row -= factor * basis_row[i]
            // 优化：预先确定是否需要neg_temp，避免条件分支在内循环中
            int use_direct_add = (tracker->ctx->field_id >= FIELD_ID_GF28 && 
                                 tracker->ctx->field_id <= FIELD_ID_GF2128);
            
            for (slong j = 0; j < ncols; j++) {
                slong idx = i * ncols + j;
                field_mul(temp, factor, &tracker->reduced_rows[idx], 
                         tracker->ctx->field_id, ctx_ptr);
                
                if (use_direct_add) {
                    // 对于GF(2^n)，减法等于加法
                    field_add(&work_row[j], &work_row[j], temp, 
                             tracker->ctx->field_id, ctx_ptr);
                } else {
                    // 对于其他域，使用实际的减法
                    field_neg(neg_temp, temp, tracker->ctx->field_id, ctx_ptr);
                    field_add(&work_row[j], &work_row[j], neg_temp, 
                             tracker->ctx->field_id, ctx_ptr);
                }
            }
        }
    }
    
    // 寻找第一个非零位置（数学逻辑不变）
    slong first_nonzero = -1;
    for (slong j = 0; j < ncols; j++) {
        if (!field_is_zero(&work_row[j], tracker->ctx->field_id, ctx_ptr)) {
            first_nonzero = j;
            break;
        }
    }
    
    // 如果全零，则线性相关
    if (first_nonzero == -1) {
        return 0;  // 工作行不需要清理，因为是预分配的
    }
    
    // 标准化工作行（使第一个非零元素为1）
    field_inv(temp, &work_row[first_nonzero], tracker->ctx->field_id, ctx_ptr);
    
    // 存储标准化的行到基中
    slong base_row_start = tracker->current_rank * ncols;
    for (slong j = 0; j < ncols; j++) {
        slong idx = base_row_start + j;
        if (j < first_nonzero) {
            // 主元之前的位置应该是0
            field_set_zero(&tracker->reduced_rows[idx], tracker->ctx->field_id, ctx_ptr);
        } else {
            // 标准化并存储
            field_mul(&tracker->reduced_rows[idx], &work_row[j], temp, 
                     tracker->ctx->field_id, ctx_ptr);
        }
    }
    
    // 更新跟踪信息（数学逻辑不变）
    tracker->pivot_cols[tracker->current_rank] = first_nonzero;
    tracker->selected_indices[tracker->current_rank] = new_row_idx;
    tracker->current_rank++;
    
    return 1;
}

// 计算选定行子矩阵中某列的最大总度数
static slong compute_fq_selected_rows_col_max_total_degree(fq_mvpoly_t ***full_matrix, 
                                                           slong *selected_rows, 
                                                           slong num_selected_rows,
                                                           slong col_idx, 
                                                           slong npars) {
    slong max_degree = -1;
    
    for (slong i = 0; i < num_selected_rows; i++) {
        slong row_idx = selected_rows[i];
        slong poly_deg = compute_fq_polynomial_total_degree(full_matrix[row_idx][col_idx], npars);
        if (poly_deg > max_degree) {
            max_degree = poly_deg;
        }
    }
    
    return max_degree;
}

// 计算选定列子矩阵中某行的最大总度数
static slong compute_fq_selected_cols_row_max_total_degree(fq_mvpoly_t ***full_matrix, 
                                                           slong row_idx,
                                                           slong *selected_cols, 
                                                           slong num_selected_cols,
                                                           slong npars) {
    slong max_degree = -1;
    
    for (slong j = 0; j < num_selected_cols; j++) {
        slong col_idx = selected_cols[j];
        slong poly_deg = compute_fq_polynomial_total_degree(full_matrix[row_idx][col_idx], npars);
        if (poly_deg > max_degree) {
            max_degree = poly_deg;
        }
    }
    
    return max_degree;
}

// 检查两个索引数组是否相同
static int indices_equal(slong *indices1, slong *indices2, slong size) {
    for (slong i = 0; i < size; i++) {
        if (indices1[i] != indices2[i]) {
            return 0;
        }
    }
    return 1;
}
// 找到最大满秩子矩阵的主元行 - 确保线性无关
// 更高效的版本：使用增量秩检查
void find_pivot_rows_nmod_fixed(slong **selected_rows_out, slong *num_selected,
                                const nmod_mat_t mat) {
    
    slong nrows = nmod_mat_nrows(mat);
    slong ncols = nmod_mat_ncols(mat);
    nmod_t mod = mat->mod;
    
    if (nrows == 0 || ncols == 0) {
        *selected_rows_out = NULL;
        *num_selected = 0;
        return;
    }
    
    slong min_dim = FLINT_MIN(nrows, ncols);
    slong *selected_rows = (slong*) flint_malloc(min_dim * sizeof(slong));
    slong rank = 0;
    
    // Create working copy
    nmod_mat_t A;
    nmod_mat_init(A, nrows, ncols, mod.n);
    nmod_mat_set(A, mat);
    
    // Row permutation tracking
    slong *P = (slong*) flint_malloc(nrows * sizeof(slong));
    for (slong i = 0; i < nrows; i++) {
        P[i] = i;
    }
    
    // PROPER LU DECOMPOSITION using FLINT-style algorithm
    slong current_row = 0;
    for (slong col = 0; col < ncols && current_row < nrows; col++) {
        // printf("%d %d\n",col,ncols);
        // Find pivot in current column
        slong pivot_row = -1;
        for (slong i = current_row; i < nrows; i++) {
            if (nmod_mat_entry(A, i, col) != 0) {
                pivot_row = i;
                break;
            }
        }
        
        // No pivot found in this column
        if (pivot_row == -1) {
            continue;
        }
        
        // Record this pivot row
        selected_rows[rank] = P[pivot_row];
        rank++;
        
        // Swap rows if needed - use FLINT's efficient method
        if (pivot_row != current_row) {
            // Swap permutation
            slong temp_idx = P[current_row];
            P[current_row] = P[pivot_row];
            P[pivot_row] = temp_idx;
            
            // Swap matrix rows efficiently
            //nn_ptr temp_row = A->rows[current_row];
            //A->rows[current_row] = A->rows[pivot_row];
            //A->rows[pivot_row] = temp_row;
            // Swap matrix rows using compatible API
            for (slong j = 0; j < ncols; j++) {
                mp_limb_t temp_val = nmod_mat_entry(A, current_row, j);
                nmod_mat_entry(A, current_row, j) = nmod_mat_entry(A, pivot_row, j);
                nmod_mat_entry(A, pivot_row, j) = temp_val;
            }
        }
        
        // Eliminate below pivot using FLINT vectorized operations
        mp_limb_t pivot = nmod_mat_entry(A, current_row, col);
        mp_limb_t pivot_inv = n_invmod(pivot, mod.n);
        
        for (slong i = current_row + 1; i < nrows; i++) {
            mp_limb_t factor = nmod_mat_entry(A, i, col);
            if (factor == 0) continue;
            
            factor = n_mulmod2_preinv(factor, pivot_inv, mod.n, mod.ninv);
            mp_limb_t neg_factor = nmod_neg(factor, mod);
            
            // Use FLINT's vectorized subtraction
            /*
            slong remaining_cols = ncols - col;
            if (remaining_cols > 0) {
                _nmod_vec_scalar_addmul_nmod(A->rows[i] + col,
                                           A->rows[current_row] + col,
                                           remaining_cols, neg_factor, mod);
            }
            */
            mp_limb_t *row_i = &nmod_mat_entry(A, i, 0);
            mp_limb_t *row_current = &nmod_mat_entry(A, current_row, 0);
            if (ncols - col > 0) {
                _nmod_vec_scalar_addmul_nmod(row_i + col, 
                                           row_current + col, 
                                           ncols - col, 
                                           neg_factor, mod);
            }
        }
        
        current_row++;
    }
    
    // Set output
    *num_selected = rank;
    if (rank > 0) {
        *selected_rows_out = (slong*) flint_realloc(selected_rows, rank * sizeof(slong));
    } else {
        flint_free(selected_rows);
        *selected_rows_out = NULL;
    }
    
    // Cleanup
    nmod_mat_clear(A);
    flint_free(P);
}
// ============================================================================
// 优化的find_pivot_rows_simple - 直接调用nmod版本
// ============================================================================
void find_pivot_rows_simple(slong **selected_rows_out, slong *num_selected,
                                        const field_elem_u *unified_mat, 
                                        slong nrows, slong ncols,
                                        field_ctx_t *ctx) {
    
    // ============================================================================
    // 素数域快速路径：直接使用nmod_fixed实现（最优性能）
    // ============================================================================
    if (ctx->field_id == FIELD_ID_NMOD) {
        //printf("Using direct nmod_fixed implementation for prime field (optimal performance)\n");
        clock_t start = clock();
        
        // 转换为nmod_mat格式
        nmod_mat_t nmod_mat;
        nmod_mat_init(nmod_mat, nrows, ncols, ctx->ctx.nmod_ctx.n);
        
        // 高效拷贝数据（直接访问nmod字段）
        for (slong i = 0; i < nrows; i++) {
            for (slong j = 0; j < ncols; j++) {
                nmod_mat_entry(nmod_mat, i, j) = unified_mat[i * ncols + j].nmod;
            }
        }
        
        // 直接使用nmod_fixed版本（避免adaptive层的选择开销）
        find_pivot_rows_nmod_fixed(selected_rows_out, num_selected, nmod_mat);
        
        clock_t end = clock();
        double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
        
        nmod_mat_clear(nmod_mat);
        //printf("Prime field computation completed in %.4f seconds using direct nmod path\n", elapsed);
        return;
    }
    
    // ============================================================================
    // 非素数域：使用统一接口实现
    // ============================================================================
    //printf("Using unified interface for non-prime field computation\n");
    
    void *ctx_ptr = (ctx->field_id == FIELD_ID_FQ_ZECH) ? 
                   (void*)ctx->ctx.zech_ctx : 
                   (void*)ctx->ctx.fq_ctx;
    
    slong max_rank = FLINT_MIN(nrows, ncols);
    slong *selected_rows = (slong*) flint_malloc(max_rank * sizeof(slong));
    slong rank = 0;
    
    // 轻量级的线性无关检查器
    typedef struct {
        field_elem_u *rows;     
        slong *pivot_positions; 
        slong count;            
        slong ncols;
    } basis_tracker_t;
    
    basis_tracker_t tracker;
    tracker.rows = (field_elem_u*) flint_calloc(max_rank * ncols, sizeof(field_elem_u));
    tracker.pivot_positions = (slong*) flint_malloc(max_rank * sizeof(slong));
    tracker.count = 0;
    tracker.ncols = ncols;
    
    // 初始化所有field元素
    for (slong i = 0; i < max_rank * ncols; i++) {
        field_init_elem(&tracker.rows[i], ctx->field_id, ctx_ptr);
    }
    
    clock_t start = clock();
    
    // 逐行处理
    for (slong row = 0; row < nrows && rank < max_rank; row++) {
        // 创建测试向量
        field_elem_u *test_vec = (field_elem_u*) flint_malloc(ncols * sizeof(field_elem_u));
        for (slong j = 0; j < ncols; j++) {
            field_init_elem(&test_vec[j], ctx->field_id, ctx_ptr);
            field_set_elem(&test_vec[j], &unified_mat[row * ncols + j], ctx->field_id, ctx_ptr);
        }
        
        // 用现有基向量消元
        for (slong i = 0; i < tracker.count; i++) {
            slong pivot_col = tracker.pivot_positions[i];
            if (!field_is_zero(&test_vec[pivot_col], ctx->field_id, ctx_ptr)) {
                field_elem_u factor;
                field_init_elem(&factor, ctx->field_id, ctx_ptr);
                field_set_elem(&factor, &test_vec[pivot_col], ctx->field_id, ctx_ptr);
                
                for (slong j = 0; j < ncols; j++) {
                    field_elem_u temp;
                    field_init_elem(&temp, ctx->field_id, ctx_ptr);
                    field_mul(&temp, &factor, &tracker.rows[i * ncols + j], ctx->field_id, ctx_ptr);
                    field_sub(&test_vec[j], &test_vec[j], &temp, ctx->field_id, ctx_ptr);
                    field_clear_elem(&temp, ctx->field_id, ctx_ptr);
                }
                
                field_clear_elem(&factor, ctx->field_id, ctx_ptr);
            }
        }
        
        // 查找第一个非零位置
        slong pivot_pos = -1;
        for (slong j = 0; j < ncols; j++) {
            if (!field_is_zero(&test_vec[j], ctx->field_id, ctx_ptr)) {
                pivot_pos = j;
                break;
            }
        }
        
        if (pivot_pos >= 0) {
            // 线性无关，加入基
            selected_rows[rank] = row;
            tracker.pivot_positions[tracker.count] = pivot_pos;
            
            // 标准化并存储
            field_elem_u inv;
            field_init_elem(&inv, ctx->field_id, ctx_ptr);
            field_inv(&inv, &test_vec[pivot_pos], ctx->field_id, ctx_ptr);
            
            for (slong j = 0; j < ncols; j++) {
                field_mul(&tracker.rows[tracker.count * ncols + j], &test_vec[j], &inv, ctx->field_id, ctx_ptr);
            }
            
            field_clear_elem(&inv, ctx->field_id, ctx_ptr);
            tracker.count++;
            rank++;
        }
        
        // 清理测试向量
        for (slong j = 0; j < ncols; j++) {
            field_clear_elem(&test_vec[j], ctx->field_id, ctx_ptr);
        }
        flint_free(test_vec);
    }
    
    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    
    // 设置输出
    if (rank > 0) {
        *selected_rows_out = (slong*) flint_malloc(rank * sizeof(slong));
        memcpy(*selected_rows_out, selected_rows, rank * sizeof(slong));
    } else {
        *selected_rows_out = NULL;
    }
    *num_selected = rank;
    
    // 清理
    for (slong i = 0; i < max_rank * ncols; i++) {
        field_clear_elem(&tracker.rows[i], ctx->field_id, ctx_ptr);
    }
    flint_free(tracker.rows);
    flint_free(tracker.pivot_positions);
    flint_free(selected_rows);
    
    //printf("Non-prime field computation completed in %.4f seconds using unified interface\n", elapsed);
}
// 优化的 find_fq_optimal_maximal_rank_submatrix
void find_fq_optimal_maximal_rank_submatrix(fq_mvpoly_t ***full_matrix, 
                                           slong nrows, slong ncols,
                                           slong **row_indices_out, 
                                           slong **col_indices_out,
                                           slong *num_rows, slong *num_cols,
                                           slong npars) {
    printf("Finding maximal rank submatrix using iterative optimization method...\n");
    
    // 获取上下文
    const fq_nmod_ctx_struct *ctx = NULL;
    for (slong i = 0; i < nrows && !ctx; i++) {
        for (slong j = 0; j < ncols && !ctx; j++) {
            if (full_matrix[i][j] != NULL) {
                ctx = full_matrix[i][j]->ctx;
                break;
            }
        }
    }
    
    // 初始化统一字段上下文
    field_ctx_t unified_ctx;
    field_ctx_init(&unified_ctx, ctx);
    void *ctx_ptr = (unified_ctx.field_id == FIELD_ID_NMOD) ? 
                   (void*)&unified_ctx.ctx.nmod_ctx : 
                   (void*)unified_ctx.ctx.fq_ctx;
    
    // 评估矩阵并转换为统一格式
    printf("Evaluating matrix and converting to unified format...\n");
    clock_t conv_start = clock();
    
    // 生成评估参数
    fq_nmod_t *param_vals = (fq_nmod_t*) flint_malloc(npars * sizeof(fq_nmod_t));
    for (slong i = 0; i < npars; i++) {
        fq_nmod_init(param_vals[i], ctx);
        fq_nmod_set_si(param_vals[i], 2 + i, ctx);
    }
    
    // 分配统一格式矩阵
    field_elem_u *unified_mat = (field_elem_u*) flint_malloc(nrows * ncols * sizeof(field_elem_u));
    
    // 初始化并评估
    for (slong i = 0; i < nrows; i++) {
        for (slong j = 0; j < ncols; j++) {
            slong idx = i * ncols + j;
            field_init_elem(&unified_mat[idx], unified_ctx.field_id, ctx_ptr);
            
            // 评估多项式
            fq_nmod_t val;
            fq_nmod_init(val, ctx);
            // 安全评估：检查NULL指针
            if (full_matrix[i][j] == NULL) {
                fq_nmod_zero(val, ctx);  // NULL表示零多项式
            } else {
                evaluate_fq_mvpoly_at_params(val, full_matrix[i][j], param_vals);
            }            fq_nmod_to_field_elem(&unified_mat[idx], val, &unified_ctx);
            fq_nmod_clear(val, ctx);
        }
    }
    
    clock_t conv_end = clock();
    //printf("Conversion time: %.3f seconds\n", (double)(conv_end - conv_start) / CLOCKS_PER_SEC);
    
    // 初始化存储选择结果的数组
    slong *current_row_indices = NULL;
    slong *current_col_indices = NULL;
    slong *prev_row_indices = NULL;
    slong *prev_col_indices = NULL;
    slong current_size = 0;
    slong prev_size = 0;
    
    // 迭代优化
    const slong MAX_ITERATIONS = 10;
    slong iteration = 0;
    int converged = 0;
    
    clock_t iter_start = clock();
    
    while (iteration < MAX_ITERATIONS && !converged) {
        //printf("\n--- Iteration %ld ---\n", iteration + 1);
        
        // 保存上一次的结果
        if (iteration > 0) {
            prev_row_indices = (slong*) flint_malloc(current_size * sizeof(slong));
            prev_col_indices = (slong*) flint_malloc(current_size * sizeof(slong));
            memcpy(prev_row_indices, current_row_indices, current_size * sizeof(slong));
            memcpy(prev_col_indices, current_col_indices, current_size * sizeof(slong));
            prev_size = current_size;
        }
        
        slong row_rank_selected = 0;  // 记录行选择的秩
        

        if (iteration == 0) {
            // First iteration: direct pivot row selection using Gaussian elimination
            //printf("Initial row selection using Gaussian elimination...\n");
            
            // Use direct Gaussian elimination to find pivot rows
            slong *selected_rows = NULL;
            slong num_selected = 0;
            
            find_pivot_rows_simple(&selected_rows, &num_selected, 
                                  unified_mat, nrows, ncols, &unified_ctx);
            
            current_size = num_selected;
            row_rank_selected = num_selected;
            current_row_indices = (slong*) flint_malloc(current_size * sizeof(slong));
            memcpy(current_row_indices, selected_rows, current_size * sizeof(slong));
            
            //printf("Selected %ld rows\n", current_size);
            
            // 新增：如果选择的行数过多，直接转置求列主元
            if (current_size > 1000) {
                printf("Large matrix detected (%ld rows), using transpose method for column selection...\n", current_size);
                
                // 构建选定行的转置矩阵
                field_elem_u *transposed_mat = (field_elem_u*) flint_malloc(ncols * current_size * sizeof(field_elem_u));
                void *ctx_ptr = (unified_ctx.field_id == FIELD_ID_NMOD) ? 
                               (void*)&unified_ctx.ctx.nmod_ctx : 
                               (void*)unified_ctx.ctx.fq_ctx;
                
                // 初始化转置矩阵元素
                for (slong i = 0; i < ncols * current_size; i++) {
                    field_init_elem(&transposed_mat[i], unified_ctx.field_id, ctx_ptr);
                }
                
                // 填充转置矩阵：原矩阵的行变为转置矩阵的列
                for (slong i = 0; i < current_size; i++) {
                    slong orig_row = current_row_indices[i];
                    for (slong j = 0; j < ncols; j++) {
                        slong src_idx = orig_row * ncols + j;
                        slong dst_idx = j * current_size + i; // 转置索引
                        field_set_elem(&transposed_mat[dst_idx], &unified_mat[src_idx],
                                      unified_ctx.field_id, ctx_ptr);
                    }
                }
                
                // 对转置矩阵找主元行（即原矩阵的主元列）
                slong *selected_cols = NULL;
                slong num_selected_cols = 0;
                
                find_pivot_rows_simple(&selected_cols, &num_selected_cols, 
                                      transposed_mat, ncols, current_size, &unified_ctx);
                
                // 设置列索引
                current_col_indices = (slong*) flint_malloc(num_selected_cols * sizeof(slong));
                memcpy(current_col_indices, selected_cols, num_selected_cols * sizeof(slong));
                
                // 更新矩阵大小为较小的维度
                current_size = FLINT_MIN(current_size, num_selected_cols);
                
                //printf("Selected %ld columns via transpose method\n", num_selected_cols);
                //printf("Final submatrix size: %ld x %ld\n", current_size, current_size);
                
                // 清理转置矩阵
                for (slong i = 0; i < ncols * current_size; i++) {
                    field_clear_elem(&transposed_mat[i], unified_ctx.field_id, ctx_ptr);
                }
                flint_free(transposed_mat);
                flint_free(selected_rows);
                flint_free(selected_cols);
                
                // 直接跳出迭代循环，避免进一步的优化步骤
                converged = 1;
                //printf("Large matrix optimization completed, skipping further iterations.\n");
                break;
            } else {
                // 正常情况：释放临时数组，继续迭代过程
                flint_free(selected_rows);
            }
        } else {
            // 后续迭代：基于当前列选择重新选择行
            //printf("Re-selecting rows based on current columns...\n");
            
            // 基于当前选定的列计算所有行的度数
            fq_index_degree_pair *row_degrees = (fq_index_degree_pair*) flint_malloc(nrows * sizeof(fq_index_degree_pair));
            
            for (slong i = 0; i < nrows; i++) {
                row_degrees[i].index = i;
                row_degrees[i].degree = compute_fq_selected_cols_row_max_total_degree(
                    full_matrix, i, current_col_indices, current_size, npars);
            }
            
            qsort(row_degrees, nrows, sizeof(fq_index_degree_pair), compare_fq_degrees);
            
            // 构建当前列的转置子矩阵用于行选择
            field_elem_u *col_submat = (field_elem_u*) flint_malloc(nrows * current_size * sizeof(field_elem_u));
            for (slong i = 0; i < nrows * current_size; i++) {
                field_init_elem(&col_submat[i], unified_ctx.field_id, ctx_ptr);
            }
            
            for (slong i = 0; i < nrows; i++) {
                for (slong j = 0; j < current_size; j++) {
                    slong src_idx = i * ncols + current_col_indices[j];
                    slong dst_idx = i * current_size + j;
                    field_set_elem(&col_submat[dst_idx], &unified_mat[src_idx],
                                  unified_ctx.field_id, ctx_ptr);
                }
            }
            
            // 选择行
            unified_row_basis_tracker_t row_tracker;
            unified_row_basis_tracker_init(&row_tracker, current_size, current_size, &unified_ctx);
            
            for (slong i = 0; i < nrows && row_tracker.current_rank < current_size; i++) {
                slong row_idx = row_degrees[i].index;
                unified_try_add_row_to_basis(&row_tracker, col_submat, row_idx, current_size);
            }
            
            // 更新行索引
            flint_free(current_row_indices);
            current_row_indices = (slong*) flint_malloc(row_tracker.current_rank * sizeof(slong));
            memcpy(current_row_indices, row_tracker.selected_indices, row_tracker.current_rank * sizeof(slong));
            row_rank_selected = row_tracker.current_rank;
            current_size = row_tracker.current_rank;
            
            // 清理
            for (slong i = 0; i < nrows * current_size; i++) {
                field_clear_elem(&col_submat[i], unified_ctx.field_id, ctx_ptr);
            }
            flint_free(col_submat);
            unified_row_basis_tracker_clear(&row_tracker);
            flint_free(row_degrees);
            
            //printf("Selected %ld rows\n", current_size);
        }
        
        // 基于当前行选择列
        //printf("Selecting columns based on current rows...\n");
        
        // 基于选定的行计算列度数
        fq_index_degree_pair *col_degrees = (fq_index_degree_pair*) flint_malloc(ncols * sizeof(fq_index_degree_pair));
        
        for (slong j = 0; j < ncols; j++) {
            col_degrees[j].index = j;
            col_degrees[j].degree = compute_fq_selected_rows_col_max_total_degree(
                full_matrix, current_row_indices, current_size, j, npars);
        }
        
        qsort(col_degrees, ncols, sizeof(fq_index_degree_pair), compare_fq_degrees);
        
        // 构建选定行的子矩阵并转置
        field_elem_u *transposed = (field_elem_u*) flint_malloc(ncols * current_size * sizeof(field_elem_u));
        for (slong i = 0; i < ncols * current_size; i++) {
            field_init_elem(&transposed[i], unified_ctx.field_id, ctx_ptr);
        }
        
        for (slong i = 0; i < current_size; i++) {
            for (slong j = 0; j < ncols; j++) {
                slong src_idx = current_row_indices[i] * ncols + j;
                slong dst_idx = j * current_size + i;
                field_set_elem(&transposed[dst_idx], &unified_mat[src_idx],
                              unified_ctx.field_id, ctx_ptr);
            }
        }
        
        // 选择列基
        unified_row_basis_tracker_t col_tracker;
        unified_row_basis_tracker_init(&col_tracker, ncols, current_size, &unified_ctx);
        
        for (slong j = 0; j < ncols && col_tracker.current_rank < current_size; j++) {
            slong col_idx = col_degrees[j].index;
            unified_try_add_row_to_basis(&col_tracker, transposed, col_idx, current_size);
        }
        
        slong col_rank_selected = col_tracker.current_rank;
        
        // 更新列索引
        if (iteration == 0) {
            current_col_indices = (slong*) flint_malloc(col_tracker.current_rank * sizeof(slong));
        } else {
            flint_free(current_col_indices);
            current_col_indices = (slong*) flint_malloc(col_tracker.current_rank * sizeof(slong));
        }
        memcpy(current_col_indices, col_tracker.selected_indices, col_tracker.current_rank * sizeof(slong));
        current_size = FLINT_MIN(current_size, col_tracker.current_rank);
        
        //printf("Selected %ld columns\n", col_tracker.current_rank);
        
        // 清理
        for (slong i = 0; i < ncols * current_size; i++) {
            field_clear_elem(&transposed[i], unified_ctx.field_id, ctx_ptr);
        }
        flint_free(transposed);
        unified_row_basis_tracker_clear(&col_tracker);
        flint_free(col_degrees);
        
        // 检查收敛
        if (iteration > 0) {
            if (current_size == prev_size &&
                current_size == FLINT_MIN(col_rank_selected, row_rank_selected) &&
                indices_equal(current_row_indices, prev_row_indices, current_size) &&
                indices_equal(current_col_indices, prev_col_indices, current_size)) {
                converged = 1;
                //printf("Converged after %d iteration!\n", iteration);
            }
            
            flint_free(prev_row_indices);
            flint_free(prev_col_indices);
        }
        
        iteration++;
    }
    
    clock_t iter_end = clock();
    printf("Iterative optimization completed in %ld iterations (%.3f seconds)\n", 
           iteration, (double)(iter_end - iter_start) / CLOCKS_PER_SEC);
    
    // 设置输出
    *row_indices_out = (slong*) flint_malloc(current_size * sizeof(slong));
    *col_indices_out = (slong*) flint_malloc(current_size * sizeof(slong));
    
    memcpy(*row_indices_out, current_row_indices, current_size * sizeof(slong));
    memcpy(*col_indices_out, current_col_indices, current_size * sizeof(slong));
    
    *num_rows = current_size;
    *num_cols = current_size;
    
    // 验证结果
    printf("Verifying final %ld x %ld submatrix...\n", current_size, current_size);
    fq_nmod_mat_t final_mat;
    fq_nmod_mat_init(final_mat, current_size, current_size, ctx);
    
    for (slong i = 0; i < current_size; i++) {
        for (slong j = 0; j < current_size; j++) {
            fq_nmod_t temp;
            fq_nmod_init(temp, ctx);
            slong idx = (*row_indices_out)[i] * ncols + (*col_indices_out)[j];
            field_elem_to_fq_nmod(temp, &unified_mat[idx], &unified_ctx);
            fq_nmod_set(fq_nmod_mat_entry(final_mat, i, j), temp, ctx);
            fq_nmod_clear(temp, ctx);
        }
    }
    
    slong final_rank = fq_nmod_mat_rank(final_mat, ctx);
    printf("Final rank: %ld/%ld %s\n", final_rank, current_size, 
           (final_rank == current_size) ? "✓" : "✗");
    
    // 计算并显示最终子矩阵的最大度数
    slong max_degree = -1;
    for (slong i = 0; i < current_size; i++) {
        for (slong j = 0; j < current_size; j++) {
            slong deg = compute_fq_polynomial_total_degree(
                full_matrix[(*row_indices_out)[i]][(*col_indices_out)[j]], npars);
            if (deg > max_degree) {
                max_degree = deg;
            }
        }
    }
    printf("Maximum total degree in final submatrix: %ld\n", max_degree);
    
    // 打印选择的索引
    printf("Selected rows: ");
    for (slong i = 0; i < current_size; i++) {
        printf("%ld ", (*row_indices_out)[i]);
    }
    printf("\nSelected cols: ");
    for (slong i = 0; i < current_size; i++) {
        printf("%ld ", (*col_indices_out)[i]);
    }
    printf("\n");
    
    // 清理
    flint_free(current_row_indices);
    flint_free(current_col_indices);
    
    for (slong i = 0; i < npars; i++) {
        fq_nmod_clear(param_vals[i], ctx);
    }
    flint_free(param_vals);
    
    // 清理统一格式矩阵
    for (slong i = 0; i < nrows * ncols; i++) {
        field_clear_elem(&unified_mat[i], unified_ctx.field_id, ctx_ptr);
    }
    flint_free(unified_mat);
    
    fq_nmod_mat_clear(final_mat, ctx);
}

// Collect monomials
typedef struct {
    slong *exp;
    slong idx;
} monom_t;

// Simple hash table entry
typedef struct hash_entry {
    slong *exp;
    slong idx;
    struct hash_entry *next;
} hash_entry_t;

// Optimized monomial collection function - replaces the original O(n²) loop
void collect_unique_monomials(
    monom_t **x_monoms_out, slong *nx_monoms_out,
    monom_t **dual_monoms_out, slong *ndual_monoms_out,
    const fq_mvpoly_t *dixon_poly, 
    const slong *d0, const slong *d1, slong nvars) {
    
    if (dixon_poly->nterms == 0) {
        *x_monoms_out = NULL; *nx_monoms_out = 0;
        *dual_monoms_out = NULL; *ndual_monoms_out = 0;
        return;
    }
    
    // Simple hash table size - power of 2 for fast modulo
    slong hash_size = 1024;
    while (hash_size < dixon_poly->nterms) hash_size <<= 1;
    
    // Hash tables for x and dual monomials
    hash_entry_t **x_buckets = (hash_entry_t**) flint_calloc(hash_size, sizeof(hash_entry_t*));
    hash_entry_t **dual_buckets = (hash_entry_t**) flint_calloc(hash_size, sizeof(hash_entry_t*));
    
    // Pre-allocate storage
    monom_t *x_monoms = (monom_t*) flint_malloc(dixon_poly->nterms * sizeof(monom_t));
    monom_t *dual_monoms = (monom_t*) flint_malloc(dixon_poly->nterms * sizeof(monom_t));
    slong *x_exp_storage = (slong*) flint_calloc(dixon_poly->nterms * nvars, sizeof(slong));
    slong *dual_exp_storage = (slong*) flint_calloc(dixon_poly->nterms * nvars, sizeof(slong));
    
    slong nx_monoms = 0, ndual_monoms = 0;
    
    // Process each term
    for (slong i = 0; i < dixon_poly->nterms; i++) {
        const slong *var_exp = dixon_poly->terms[i].var_exp;
        if (!var_exp) continue;
        
        // Check degree bounds
        int valid = 1;
        for (slong k = 0; k < nvars && valid; k++) {
            if (var_exp[k] >= d0[k] || var_exp[nvars + k] >= d1[k]) {
                valid = 0;
            }
        }
        if (!valid) continue;
        
        // Process x-monomial (first nvars components)
        // Simple hash: sum of exponents * prime
        ulong x_hash = 0;
        for (slong k = 0; k < nvars; k++) {
            x_hash = x_hash * 31 + var_exp[k];
        }
        x_hash &= (hash_size - 1); // Fast modulo for power of 2
        
        // Check if x-monomial already exists
        hash_entry_t *entry = x_buckets[x_hash];
        int found = 0;
        while (entry && !found) {
            if (memcmp(entry->exp, var_exp, nvars * sizeof(slong)) == 0) {
                found = 1;
            } else {
                entry = entry->next;
            }
        }
        
        if (!found) {
            // Add new x-monomial
            slong *x_exp = &x_exp_storage[nx_monoms * nvars];
            memcpy(x_exp, var_exp, nvars * sizeof(slong));
            
            x_monoms[nx_monoms].exp = x_exp;
            x_monoms[nx_monoms].idx = nx_monoms;
            
            // Insert into hash table
            hash_entry_t *new_entry = (hash_entry_t*) flint_malloc(sizeof(hash_entry_t));
            new_entry->exp = x_exp;
            new_entry->idx = nx_monoms;
            new_entry->next = x_buckets[x_hash];
            x_buckets[x_hash] = new_entry;
            
            nx_monoms++;
        }
        
        // Process dual monomial (next nvars components)
        const slong *dual_exp_src = &var_exp[nvars];
        
        // Hash for dual monomial
        ulong dual_hash = 0;
        for (slong k = 0; k < nvars; k++) {
            dual_hash = dual_hash * 31 + dual_exp_src[k];
        }
        dual_hash &= (hash_size - 1);
        
        // Check if dual monomial already exists
        entry = dual_buckets[dual_hash];
        found = 0;
        while (entry && !found) {
            if (memcmp(entry->exp, dual_exp_src, nvars * sizeof(slong)) == 0) {
                found = 1;
            } else {
                entry = entry->next;
            }
        }
        
        if (!found) {
            // Add new dual monomial
            slong *dual_exp = &dual_exp_storage[ndual_monoms * nvars];
            memcpy(dual_exp, dual_exp_src, nvars * sizeof(slong));
            
            dual_monoms[ndual_monoms].exp = dual_exp;
            dual_monoms[ndual_monoms].idx = ndual_monoms;
            
            // Insert into hash table
            hash_entry_t *new_entry = (hash_entry_t*) flint_malloc(sizeof(hash_entry_t));
            new_entry->exp = dual_exp;
            new_entry->idx = ndual_monoms;
            new_entry->next = dual_buckets[dual_hash];
            dual_buckets[dual_hash] = new_entry;
            
            ndual_monoms++;
        }
    }
    
    // Resize to actual size
    if (nx_monoms > 0) {
        x_monoms = (monom_t*) flint_realloc(x_monoms, nx_monoms * sizeof(monom_t));
    } else {
        flint_free(x_monoms);
        x_monoms = NULL;
        flint_free(x_exp_storage);
    }
    
    if (ndual_monoms > 0) {
        dual_monoms = (monom_t*) flint_realloc(dual_monoms, ndual_monoms * sizeof(monom_t));
    } else {
        flint_free(dual_monoms);
        dual_monoms = NULL;
        flint_free(dual_exp_storage);
    }
    
    // Clean up hash tables
    for (slong i = 0; i < hash_size; i++) {
        hash_entry_t *entry = x_buckets[i];
        while (entry) {
            hash_entry_t *next = entry->next;
            flint_free(entry);
            entry = next;
        }
        
        entry = dual_buckets[i];
        while (entry) {
            hash_entry_t *next = entry->next;
            flint_free(entry);
            entry = next;
        }
    }
    flint_free(x_buckets);
    flint_free(dual_buckets);
    
    // Set outputs
    *x_monoms_out = x_monoms;
    *nx_monoms_out = nx_monoms;
    *dual_monoms_out = dual_monoms;
    *ndual_monoms_out = ndual_monoms;
    
    printf("Found %ld x-monomials and %ld ~x-monomials (after degree filtering)\n", 
           nx_monoms, ndual_monoms);
}
// 按需分配单个元素
fq_mvpoly_t* get_matrix_entry_lazy(fq_mvpoly_t ***matrix, slong i, slong j,
                                  slong npars, const fq_nmod_ctx_t ctx) {
    if (!matrix[i][j]) {
        matrix[i][j] = (fq_mvpoly_t*) flint_malloc(sizeof(fq_mvpoly_t));
        fq_mvpoly_init(matrix[i][j], 0, npars, ctx);
    }
    return matrix[i][j];
}
void fill_coefficient_matrix_optimized(fq_mvpoly_t ***full_matrix,
                                      monom_t *x_monoms, slong nx_monoms,
                                      monom_t *dual_monoms, slong ndual_monoms,
                                      const fq_mvpoly_t *dixon_poly,
                                      const slong *d0, const slong *d1, 
                                      slong nvars, slong npars) {
    
    printf("Filling coefficient matrix with lazy allocation...\n");
    
    // 为每个Dixon多项式项找到对应的矩阵位置
    for (slong t = 0; t < dixon_poly->nterms; t++) {
        if (!dixon_poly->terms[t].var_exp) continue;
        
        // 检查度数界限
        int valid = 1;
        for (slong k = 0; k < nvars; k++) {
            if (dixon_poly->terms[t].var_exp[k] >= d0[k] || 
                dixon_poly->terms[t].var_exp[nvars + k] >= d1[k]) {
                valid = 0;
                break;
            }
        }
        if (!valid) continue;
        
        // 找到x-单项式的行索引
        slong row = -1;
        for (slong i = 0; i < nx_monoms; i++) {
            if (memcmp(x_monoms[i].exp, dixon_poly->terms[t].var_exp, 
                      nvars * sizeof(slong)) == 0) {
                row = i;
                break;
            }
        }
        
        // 找到dual-单项式的列索引
        slong col = -1;
        for (slong j = 0; j < ndual_monoms; j++) {
            if (memcmp(dual_monoms[j].exp, &dixon_poly->terms[t].var_exp[nvars], 
                      nvars * sizeof(slong)) == 0) {
                col = j;
                break;
            }
        }
        
        // 只为实际需要的位置分配内存
        if (row >= 0 && col >= 0) {
            fq_mvpoly_t *entry = get_matrix_entry_lazy(full_matrix, row, col, 
                                                      npars, dixon_poly->ctx);
            fq_mvpoly_add_term_fast(entry, NULL, dixon_poly->terms[t].par_exp, 
                                   dixon_poly->terms[t].coeff);
        }
    }
}
// 优化版本的 find_fq_optimal_maximal_rank_submatrix
// ============ Extract coefficient matrix ============

void extract_fq_coefficient_matrix_from_dixon(fq_mvpoly_t ***coeff_matrix,
                                              slong *row_indices, slong *col_indices,
                                              slong *matrix_size,
                                              const fq_mvpoly_t *dixon_poly,
                                              slong nvars, slong npars) {
    printf("Extracting coefficient matrix from Dixon polynomial\n");
    
    // First compute degree bounds
    slong *d0 = (slong*) flint_calloc(nvars, sizeof(slong));
    slong *d1 = (slong*) flint_calloc(nvars, sizeof(slong));
    
    // Find maximum degree for each variable
    for (slong i = 0; i < dixon_poly->nterms; i++) {
        if (dixon_poly->terms[i].var_exp) {
            // Check x variables (first nvars)
            for (slong j = 0; j < nvars; j++) {
                if (dixon_poly->terms[i].var_exp[j] > d0[j]) {
                    d0[j] = dixon_poly->terms[i].var_exp[j];
                }
            }
            // Check ~x variables (next nvars)
            for (slong j = 0; j < nvars; j++) {
                if (dixon_poly->terms[i].var_exp[nvars + j] > d1[j]) {
                    d1[j] = dixon_poly->terms[i].var_exp[nvars + j];
                }
            }
        }
    }
    
    // Add 1 to each degree bound
    for (slong i = 0; i < nvars; i++) {
        d0[i]++;
        d1[i]++;
    }
    
    printf("Degree bounds - x vars: ");
    for (slong i = 0; i < nvars; i++) printf("%ld ", d0[i]);
    printf("\nDegree bounds - ~x vars: ");
    for (slong i = 0; i < nvars; i++) printf("%ld ", d1[i]);
    printf("\n");
    
    // Calculate expected matrix dimensions
    slong expected_rows = 1;
    slong expected_cols = 1;
    for (slong i = 0; i < nvars; i++) {
        expected_rows *= d0[i];
        expected_cols *= d1[i];
    }
    printf("Expected matrix size: %ld x %ld\n", expected_rows, expected_cols);
    
    monom_t *x_monoms = (monom_t*) flint_malloc(dixon_poly->nterms * sizeof(monom_t));
    monom_t *dual_monoms = (monom_t*) flint_malloc(dixon_poly->nterms * sizeof(monom_t));
    slong nx_monoms = 0, ndual_monoms = 0;
    
    collect_unique_monomials(&x_monoms, &nx_monoms,
                        &dual_monoms, &ndual_monoms,
                        dixon_poly, d0, d1, nvars);
    
    if (nx_monoms == 0 || ndual_monoms == 0) {
        printf("Warning: Empty coefficient matrix\n");
        *matrix_size = 0;
        flint_free(d0);
        flint_free(d1);
        flint_free(x_monoms);
        flint_free(dual_monoms);
        return;
    }
    // Build full coefficient matrix (entries are polynomials in parameters)
    fq_mvpoly_t ***full_matrix = (fq_mvpoly_t***) flint_malloc(nx_monoms * sizeof(fq_mvpoly_t**));
    for (slong i = 0; i < nx_monoms; i++) {
        full_matrix[i] = (fq_mvpoly_t**) flint_calloc(ndual_monoms, sizeof(fq_mvpoly_t*));
        // 注意：使用calloc初始化为NULL
    }
    // Fill the coefficient matrix
 fill_coefficient_matrix_optimized(full_matrix, x_monoms, nx_monoms, 
                                 dual_monoms, ndual_monoms, dixon_poly, 
                                 d0, d1, nvars, npars);
    // Find maximal rank submatrix
    slong *row_idx_array = NULL;
    slong *col_idx_array = NULL;
    slong num_rows, num_cols;
    
    if (npars == 0) {
        // Special case: no parameters, directly find rank
        fq_nmod_mat_t eval_mat;
        fq_nmod_mat_init(eval_mat, nx_monoms, ndual_monoms, dixon_poly->ctx);
        for (slong i = 0; i < nx_monoms; i++) {
            for (slong j = 0; j < ndual_monoms; j++) {
                // **CRITICAL FIX: Check for NULL pointer before dereferencing**
                if (full_matrix[i][j] != NULL && full_matrix[i][j]->nterms > 0) {
                    // Entry exists and has terms - use the coefficient
                    fq_nmod_set(fq_nmod_mat_entry(eval_mat, i, j), 
                               full_matrix[i][j]->terms[0].coeff, dixon_poly->ctx);
                } else {
                    // Entry is NULL (never allocated) or empty - set to zero
                    fq_nmod_zero(fq_nmod_mat_entry(eval_mat, i, j), dixon_poly->ctx);
                }
            }
        }
        slong rank = fq_nmod_mat_rank(eval_mat, dixon_poly->ctx);
        // Create identity selection for square case
        slong min_size = FLINT_MIN(nx_monoms, ndual_monoms);
        slong actual_size = FLINT_MIN(rank, min_size);
        
        row_idx_array = (slong*) flint_malloc(actual_size * sizeof(slong));
        col_idx_array = (slong*) flint_malloc(actual_size * sizeof(slong));
        
        for (slong i = 0; i < actual_size; i++) {
            row_idx_array[i] = i;
            col_idx_array[i] = i;
        }
        num_rows = actual_size;
        num_cols = actual_size;
        
        fq_nmod_mat_clear(eval_mat, dixon_poly->ctx);

    } else {
        // Check if matrix is small enough to use directly
        if (nx_monoms < 4 && ndual_monoms < 4) {
            printf("Small matrix detected (%ld x %ld), using original matrix directly\n", 
                   nx_monoms, ndual_monoms);
            
            // Create identity selection for the original matrix
            slong min_size = FLINT_MIN(nx_monoms, ndual_monoms);
            
            row_idx_array = (slong*) flint_malloc(min_size * sizeof(slong));
            col_idx_array = (slong*) flint_malloc(min_size * sizeof(slong));
            
            for (slong i = 0; i < min_size; i++) {
                row_idx_array[i] = i;
                col_idx_array[i] = i;
            }
            
            num_rows = min_size;
            num_cols = min_size;
        } else {
            // Call submatrix finding function for larger matrices
            find_fq_optimal_maximal_rank_submatrix(full_matrix, nx_monoms, ndual_monoms,
                                                  &row_idx_array, &col_idx_array, 
                                                  &num_rows, &num_cols,
                                                  npars);
        }
    }
    // Take square part if needed
    slong submat_rank = FLINT_MIN(num_rows, num_cols);
    
    if (submat_rank == 0) {
        printf("Warning: Matrix has rank 0\n");
        *matrix_size = 0;
        
        // Cleanup full_matrix with NULL checks (lazy allocation)
        if (full_matrix) {
            for (slong i = 0; i < nx_monoms; i++) {
                if (full_matrix[i]) {
                    for (slong j = 0; j < ndual_monoms; j++) {
                        if (full_matrix[i][j] != NULL) {  // **CRITICAL: Check for NULL**
                            fq_mvpoly_clear(full_matrix[i][j]);
                            flint_free(full_matrix[i][j]);
                        }
                    }
                    flint_free(full_matrix[i]);
                }
            }
            flint_free(full_matrix);
        }
        
        // Cleanup index arrays
        if (row_idx_array) flint_free(row_idx_array);
        if (col_idx_array) flint_free(col_idx_array);
        
        // Just free the monomial arrays themselves
        if (x_monoms) flint_free(x_monoms);
        if (dual_monoms) flint_free(dual_monoms);
        
        // Free degree bound arrays
        if (d0) flint_free(d0);
        if (d1) flint_free(d1);
        
        return;
    }
    

    // 安全的清理代码 - 修复NULL指针段错误
    
    printf("Extracted submatrix of size %ld x %ld\n", submat_rank, submat_rank);
    
    // Build the output submatrix
    *coeff_matrix = (fq_mvpoly_t**) flint_malloc(submat_rank * sizeof(fq_mvpoly_t*));
    for (slong i = 0; i < submat_rank; i++) {
        (*coeff_matrix)[i] = (fq_mvpoly_t*) flint_malloc(submat_rank * sizeof(fq_mvpoly_t));
        for (slong j = 0; j < submat_rank; j++) {
            // 安全复制：检查源矩阵元素是否为NULL
            fq_mvpoly_t *source = full_matrix[row_idx_array[i]][col_idx_array[j]];
            if (source != NULL) {
                fq_mvpoly_copy(&(*coeff_matrix)[i][j], source);
            } else {
                // NULL元素表示零多项式，直接初始化为零
                fq_mvpoly_init(&(*coeff_matrix)[i][j], 0, npars, dixon_poly->ctx);
            }
        }
    }
    
    // Copy indices
    for (slong i = 0; i < submat_rank; i++) {
        row_indices[i] = row_idx_array[i];
        col_indices[i] = col_idx_array[i];
    }
    *matrix_size = submat_rank;
    
    // 安全的cleanup - 添加NULL检查
    printf("Cleaning up matrix...\n");
    if (full_matrix) {
        for (slong i = 0; i < nx_monoms; i++) {
            if (full_matrix[i]) {
                for (slong j = 0; j < ndual_monoms; j++) {
                    if (full_matrix[i][j]) {  // 关键：检查NULL指针
                        fq_mvpoly_clear(full_matrix[i][j]);
                        flint_free(full_matrix[i][j]);
                    }
                }
                flint_free(full_matrix[i]);
            }
        }
        flint_free(full_matrix);
    }
    
    // 安全清理其他数组
    if (row_idx_array) flint_free(row_idx_array);
    if (col_idx_array) flint_free(col_idx_array);
    
    // 清理单项式数组 - 注意这里的内存管理
    if (x_monoms) {
        // x_monoms[i].exp 指向连续存储，不需要单独释放
        flint_free(x_monoms);
    }
    if (dual_monoms) {
        // dual_monoms[j].exp 指向连续存储，不需要单独释放  
        flint_free(dual_monoms);
    }
    
    // 清理度数界限数组
    if (d0) flint_free(d0);
    if (d1) flint_free(d1);
    
    printf("Cleanup completed safely\n");

}

// Compute determinant of cancellation matrix
void compute_fq_cancel_matrix_det(fq_mvpoly_t *result, fq_mvpoly_t **modified_M_mvpoly,
                                 slong nvars, slong npars, det_method_t method) {
    printf("Computing cancellation matrix determinant using recursive expansion\n");
    
    clock_t start = clock();
    compute_fq_det_recursive(result, modified_M_mvpoly, nvars + 1);
    clock_t end = clock();
    
    printf("End (%.3f seconds)\n", (double)(end - start) / CLOCKS_PER_SEC);
}

// ============ Compute determinant of coefficient matrix ============

// ============ Main Dixon resultant function ============
void fq_dixon_resultant(fq_mvpoly_t *result, fq_mvpoly_t *polys, 
                       slong nvars, slong npars) {

    printf("Dixon Resultant Computation over Finite Extension Fields\n");
    printf("========================================================\n");
    printf("\n=== Dixon Resultant over F_{p^d} ===\n");
    cleanup_unified_workspace();
    // Step 1: Build cancellation matrix
    printf("\nStep 1: Build Cancellation Matrix\n");
    fq_mvpoly_t **M_mvpoly;
    build_fq_cancellation_matrix_mvpoly(&M_mvpoly, polys, nvars, npars);
    
    // 显示原始矩阵的分析
    //analyze_fq_matrix_mvpoly(M_mvpoly, nvars + 1, nvars + 1, "Original Cancellation");
    
    // Step 2: Perform row operations
    printf("\nStep 2: Perform Matrix Row Operations\n");
    fq_mvpoly_t **modified_M_mvpoly;
    perform_fq_matrix_row_operations_mvpoly(&modified_M_mvpoly, &M_mvpoly, nvars, npars);
    
    // **新增：显示除法后的矩阵**
    //printf("\n=== Matrix After Row Operations ===\n");
    //analyze_fq_matrix_mvpoly(modified_M_mvpoly, nvars + 1, nvars + 1, "After Row Operations");
    /*
    // 如果矩阵不太大，显示详细内容
    if (nvars <= 3) {
        print_fq_matrix_mvpoly(modified_M_mvpoly, nvars + 1, nvars + 1, 
                               "After Row Operations (Detailed)", 1);
    } else {
        print_fq_matrix_mvpoly(modified_M_mvpoly, nvars + 1, nvars + 1, 
                               "After Row Operations (Summary)", 0);
    }
    
    // 比较原始矩阵和修改后矩阵的复杂度
    printf("\n=== Complexity Comparison ===\n");
    
    // 计算原始矩阵的总项数
    slong orig_total_terms = 0;
    for (slong i = 0; i <= nvars; i++) {
        for (slong j = 0; j <= nvars; j++) {
            orig_total_terms += M_mvpoly[i][j].nterms;
        }
    }
    
    // 计算修改后矩阵的总项数
    slong modified_total_terms = 0;
    for (slong i = 0; i <= nvars; i++) {
        for (slong j = 0; j <= nvars; j++) {
            modified_total_terms += modified_M_mvpoly[i][j].nterms;
        }
    }
    
    printf("Original matrix total terms: %ld\n", orig_total_terms);
    printf("Modified matrix total terms: %ld\n", modified_total_terms);
    printf("Complexity ratio: %.2f\n", (double)modified_total_terms / orig_total_terms);
    */
    // Step 3: Compute determinant of modified matrix
    printf("\nStep 3: Compute determinant of modified matrix\n");
    fq_mvpoly_t d_poly;
    compute_fq_cancel_matrix_det(&d_poly, modified_M_mvpoly, nvars, npars, DET_METHOD_RECURSIVE);
    
    printf("Dixon polynomial has %ld terms\n", d_poly.nterms);
    
    if (d_poly.nterms <= 100) {
        fq_mvpoly_print_expanded(&d_poly, "Dixon", 1);
    } else {
        printf("Dixon polynomial too large to display (%ld terms)\n", d_poly.nterms);
        /*
        // 显示部分项作为示例
        printf("First 5 terms of Dixon polynomial:\n");
        for (slong i = 0; i < FLINT_MIN(5, d_poly.nterms); i++) {
            printf("  Term %ld: ", i);
            fq_nmod_print_pretty(d_poly.terms[i].coeff, d_poly.ctx);
            
            // 显示变量指数
            if (d_poly.terms[i].var_exp) {
                for (slong j = 0; j < d_poly.nvars; j++) {
                    if (d_poly.terms[i].var_exp[j] > 0) {
                        if (j < d_poly.nvars / 2) {
                            printf(" * x%ld^%ld", j, d_poly.terms[i].var_exp[j]);
                        } else {
                            printf(" * ~x%ld^%ld", j - d_poly.nvars / 2, d_poly.terms[i].var_exp[j]);
                        }
                    }
                }
            }
            
            // 显示参数指数
            if (d_poly.terms[i].par_exp && d_poly.npars > 0) {
                for (slong j = 0; j < d_poly.npars; j++) {
                    if (d_poly.terms[i].par_exp[j] > 0) {
                        printf(" * p%ld^%ld", j, d_poly.terms[i].par_exp[j]);
                    }
                }
            }
            printf("\n");
        }
        */
    }

    for (slong i = 0; i <= nvars; i++) {
        for (slong j = 0; j <= nvars; j++) {
            fq_mvpoly_clear(&M_mvpoly[i][j]);
            fq_mvpoly_clear(&modified_M_mvpoly[i][j]);
        }
        flint_free(M_mvpoly[i]);
        flint_free(modified_M_mvpoly[i]);
    }
    flint_free(M_mvpoly);
    flint_free(modified_M_mvpoly);
    
    // Step 4: Extract coefficient matrix
    printf("\nStep 4: Extract coefficient matrix\n");
    
    fq_mvpoly_t **coeff_matrix;
    slong *row_indices = (slong*) flint_malloc(d_poly.nterms * sizeof(slong));
    slong *col_indices = (slong*) flint_malloc(d_poly.nterms * sizeof(slong));
    slong matrix_size;
    
    extract_fq_coefficient_matrix_from_dixon(&coeff_matrix, row_indices, col_indices,
                                            &matrix_size, &d_poly, nvars, npars);
    
    if (matrix_size > 0) {
        printf("\nStep 5: Compute determinant of coefficient matrix (fq_nmod)\n");
        
        slong res_deg_bound = compute_fq_dixon_resultant_degree_bound(polys, nvars+1, nvars, npars);
        printf("Resultant degree bound: %ld\n", res_deg_bound);
        
        ulong field_size = 1;
        for (slong i = 0; i < fq_nmod_ctx_degree(polys[0].ctx); i++) {
            field_size *= fq_nmod_ctx_prime(polys[0].ctx);
        }
        
        det_method_t coeff_method; // DET_METHOD_RECURSIVE DET_METHOD_KRONECKER DET_METHOD_INTERPOLATION DET_METHOD_HUANG
        if (field_size < (1UL << 16)) {
            coeff_method = DET_METHOD_INTERPOLATION;
        } else if (matrix_size < 9) {
            coeff_method = DET_METHOD_RECURSIVE;
        } else {
            coeff_method = DET_METHOD_KRONECKER;
        }
        //coeff_method = DET_METHOD_INTERPOLATION;
#ifdef DIXON_DET_METHOD
        coeff_method = DIXON_DET_METHOD;
        printf("Method overridden by macro: %d\n", DIXON_DET_METHOD);
#endif

        printf("Using coefficient matrix determinant method %d:", coeff_method);
        
        compute_fq_coefficient_matrix_det(result, coeff_matrix, matrix_size,
                                         npars, polys[0].ctx, coeff_method, res_deg_bound);
        
        printf("Resultant polynomial has %ld terms\n", result->nterms);
        if (result->nterms <= 100) {
            fq_mvpoly_print(result, "Final Resultant");
        } else {
            printf("Final resultant too large to display (%ld terms)\n", result->nterms);
            
            // 显示次数信息
            slong max_par_deg = 0;
            for (slong i = 0; i < result->nterms; i++) {
                if (result->terms[i].par_exp && result->npars > 0) {
                    for (slong j = 0; j < result->npars; j++) {
                        if (result->terms[i].par_exp[j] > max_par_deg) {
                            max_par_deg = result->terms[i].par_exp[j];
                        }
                    }
                }
            }
            printf("Maximum parameter degree in resultant: %ld\n", max_par_deg);
        }
        fq_mvpoly_make_monic(result);
        // Cleanup coefficient matrix
        for (slong i = 0; i < matrix_size; i++) {
            for (slong j = 0; j < matrix_size; j++) {
                fq_mvpoly_clear(&coeff_matrix[i][j]);
            }
            flint_free(coeff_matrix[i]);
        }
        flint_free(coeff_matrix);
    } else {
        fq_mvpoly_init(result, 0, npars, polys[0].ctx);
        printf("Warning: Empty coefficient matrix, resultant is 0\n");
    }
    
    // Cleanup
    flint_free(row_indices);
    flint_free(col_indices);
    
    fq_mvpoly_clear(&d_poly);
    
    printf("\n=== Dixon Resultant Computation Complete ===\n");
}

#endif