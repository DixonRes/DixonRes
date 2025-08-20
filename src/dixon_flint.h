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
                
                fq_mvpoly_add_term(&(*M)[i][j], new_var_exp, polys[j].terms[t].par_exp, 
                                  polys[j].terms[t].coeff);
                flint_free(new_var_exp);
            }
        }
    }
}

void perform_fq_matrix_row_operations_mvpoly(fq_mvpoly_t ***new_matrix, fq_mvpoly_t ***original_matrix,
                                            slong nvars, slong npars) {
    slong n = nvars + 1;
    
    // Allocate new matrix space
    *new_matrix = (fq_mvpoly_t**) flint_malloc(n * sizeof(fq_mvpoly_t*));
    for (slong i = 0; i < n; i++) {
        (*new_matrix)[i] = (fq_mvpoly_t*) flint_malloc(n * sizeof(fq_mvpoly_t));
    }
    
    printf("Performing row operations on matrix:\n");
    
    // First row remains unchanged
    printf("Row 0: unchanged\n");
    for (slong j = 0; j < n; j++) {
        fq_mvpoly_copy(&(*new_matrix)[0][j], &(*original_matrix)[0][j]);
    }
    
    // Process remaining rows
    for (slong i = 1; i < n; i++) {
        printf("Row %ld: (row[%ld] - row[%ld]) / (x%ld - ~x%ld)\n", i, i, i-1, i-1, i-1);
        
        for (slong j = 0; j < n; j++) {
            // Calculate difference
            fq_mvpoly_t diff;
            fq_mvpoly_init(&diff, 2 * nvars, npars, (*original_matrix)[0][0].ctx);
            
            // diff = M[i][j] - M[i-1][j]
            fq_mvpoly_copy(&diff, &(*original_matrix)[i][j]);
            for (slong t = 0; t < (*original_matrix)[i-1][j].nterms; t++) {
                fq_nmod_t neg_coeff;
                fq_nmod_init(neg_coeff, diff.ctx);
                fq_nmod_neg(neg_coeff, (*original_matrix)[i-1][j].terms[t].coeff, diff.ctx);
                
                fq_mvpoly_add_term(&diff, 
                                  (*original_matrix)[i-1][j].terms[t].var_exp,
                                  (*original_matrix)[i-1][j].terms[t].par_exp,
                                  neg_coeff);
                                  
                fq_nmod_clear(neg_coeff, diff.ctx);
            }
            
            // Perform division using improved method
            fq_mvpoly_init(&(*new_matrix)[i][j], 2 * nvars, npars, diff.ctx);
            
            if (diff.nterms > 0) {
                // Use the improved division method (could also use divide_by_fq_linear_factor_flint)
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
            fq_mvpoly_add_term(result, NULL, NULL, det);
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
                        fq_mvpoly_add_term(result, NULL, par_exp, coeff);
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
    if (poly->nterms == 0) {
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
    if (poly->nterms == 0) {
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

// 重写的行基跟踪器结构
typedef struct {
    field_elem_u *reduced_rows;    // 已约化的行向量（统一格式）
    slong *pivot_cols;             // 每行的主元列位置  
    slong *selected_indices;       // 已选择的原始行索引
    slong current_rank;            // 当前秩
    slong max_size;
    slong ncols;
    field_ctx_t *ctx;              // 统一字段上下文
    int initialized;               // 初始化标志
} unified_row_basis_tracker_t;

// 完整重写：初始化统一格式的行基跟踪器
static void unified_row_basis_tracker_init(unified_row_basis_tracker_t *tracker, 
                                          slong max_size, slong ncols, 
                                          field_ctx_t *ctx) {
    tracker->max_size = max_size;
    tracker->ncols = ncols;
    tracker->ctx = ctx;
    tracker->current_rank = 0;
    tracker->initialized = 1;
    
    void *ctx_ptr = (ctx->field_id == FIELD_ID_NMOD) ? 
                   (void*)&ctx->ctx.nmod_ctx : 
                   (void*)ctx->ctx.fq_ctx;
    
    // 分配内存
    tracker->reduced_rows = (field_elem_u*) flint_calloc(max_size * ncols, sizeof(field_elem_u));
    tracker->pivot_cols = (slong*) flint_calloc(max_size, sizeof(slong));
    tracker->selected_indices = (slong*) flint_calloc(max_size, sizeof(slong));
    
    // 初始化所有元素
    for (slong i = 0; i < max_size * ncols; i++) {
        field_init_elem(&tracker->reduced_rows[i], ctx->field_id, ctx_ptr);
        field_set_zero(&tracker->reduced_rows[i], ctx->field_id, ctx_ptr);
    }
    
    // 初始化主元列为-1
    for (slong i = 0; i < max_size; i++) {
        tracker->pivot_cols[i] = -1;
        tracker->selected_indices[i] = -1;
    }
}

// 完整重写：清理统一格式的行基跟踪器
static void unified_row_basis_tracker_clear(unified_row_basis_tracker_t *tracker) {
    if (!tracker->initialized) return;
    
    void *ctx_ptr = (tracker->ctx->field_id == FIELD_ID_NMOD) ? 
                   (void*)&tracker->ctx->ctx.nmod_ctx : 
                   (void*)tracker->ctx->ctx.fq_ctx;
    
    // 清理所有字段元素
    if (tracker->reduced_rows) {
        for (slong i = 0; i < tracker->max_size * tracker->ncols; i++) {
            field_clear_elem(&tracker->reduced_rows[i], tracker->ctx->field_id, ctx_ptr);
        }
        flint_free(tracker->reduced_rows);
    }
    
    if (tracker->pivot_cols) flint_free(tracker->pivot_cols);
    if (tracker->selected_indices) flint_free(tracker->selected_indices);
    
    tracker->initialized = 0;
}

static int unified_try_add_row_to_basis(unified_row_basis_tracker_t *tracker, 
                                       const field_elem_u *unified_mat,
                                       slong new_row_idx, slong ncols) {
    if (!tracker->initialized || tracker->current_rank >= tracker->max_size) {
        return 0;
    }
    
    void *ctx_ptr = (tracker->ctx->field_id == FIELD_ID_NMOD) ? 
                   (void*)&tracker->ctx->ctx.nmod_ctx : 
                   (void*)tracker->ctx->ctx.fq_ctx;
    
    // 分配工作行并初始化
    field_elem_u *work_row = (field_elem_u*) flint_malloc(ncols * sizeof(field_elem_u));
    for (slong j = 0; j < ncols; j++) {
        field_init_elem(&work_row[j], tracker->ctx->field_id, ctx_ptr);
    }
    
    // 复制输入行到工作行
    for (slong j = 0; j < ncols; j++) {
        field_set_elem(&work_row[j], &unified_mat[new_row_idx * ncols + j], 
                      tracker->ctx->field_id, ctx_ptr);
    }
    
    // 对已有的每个基向量进行消元
    for (slong i = 0; i < tracker->current_rank; i++) {
        slong pivot_col = tracker->pivot_cols[i];
        
        // 验证主元列索引
        if (pivot_col < 0 || pivot_col >= ncols) continue;
        
        // 如果工作行在主元位置非零，进行消元
        if (!field_is_zero(&work_row[pivot_col], tracker->ctx->field_id, ctx_ptr)) {
            // 分配临时变量
            field_elem_u factor, temp, pivot_val;
            field_init_elem(&factor, tracker->ctx->field_id, ctx_ptr);
            field_init_elem(&temp, tracker->ctx->field_id, ctx_ptr);
            field_init_elem(&pivot_val, tracker->ctx->field_id, ctx_ptr);
            
            // 获取基向量的主元值
            slong base_idx = i * ncols + pivot_col;
            field_set_elem(&pivot_val, &tracker->reduced_rows[base_idx], 
                          tracker->ctx->field_id, ctx_ptr);
            
            // 计算消元因子 = work_row[pivot_col] / pivot_val
            field_inv(&temp, &pivot_val, tracker->ctx->field_id, ctx_ptr);
            field_mul(&factor, &work_row[pivot_col], &temp, tracker->ctx->field_id, ctx_ptr);
            
            // 执行消元：work_row -= factor * basis_row[i]
            for (slong j = 0; j < ncols; j++) {
                slong idx = i * ncols + j;
                field_mul(&temp, &factor, &tracker->reduced_rows[idx], 
                         tracker->ctx->field_id, ctx_ptr);
                
                // 对于GF(2^n)，减法等于加法
                if (tracker->ctx->field_id >= FIELD_ID_GF28 && 
                    tracker->ctx->field_id <= FIELD_ID_GF2128) {
                    field_add(&work_row[j], &work_row[j], &temp, 
                             tracker->ctx->field_id, ctx_ptr);
                } else {
                    // 对于其他域，使用实际的减法
                    field_elem_u neg_temp;
                    field_init_elem(&neg_temp, tracker->ctx->field_id, ctx_ptr);
                    field_neg(&neg_temp, &temp, tracker->ctx->field_id, ctx_ptr);
                    field_add(&work_row[j], &work_row[j], &neg_temp, 
                             tracker->ctx->field_id, ctx_ptr);
                    field_clear_elem(&neg_temp, tracker->ctx->field_id, ctx_ptr);
                }
            }
            
            // 清理临时变量
            field_clear_elem(&factor, tracker->ctx->field_id, ctx_ptr);
            field_clear_elem(&temp, tracker->ctx->field_id, ctx_ptr);
            field_clear_elem(&pivot_val, tracker->ctx->field_id, ctx_ptr);
        }
    }
    
    // 寻找第一个非零位置
    slong first_nonzero = -1;
    for (slong j = 0; j < ncols; j++) {
        if (!field_is_zero(&work_row[j], tracker->ctx->field_id, ctx_ptr)) {
            first_nonzero = j;
            break;
        }
    }
    
    // 如果全零，则线性相关
    if (first_nonzero == -1) {
        // 清理工作行
        for (slong j = 0; j < ncols; j++) {
            field_clear_elem(&work_row[j], tracker->ctx->field_id, ctx_ptr);
        }
        flint_free(work_row);
        return 0;
    }
    
    // 标准化工作行（使第一个非零元素为1）
    field_elem_u inv_leading;
    field_init_elem(&inv_leading, tracker->ctx->field_id, ctx_ptr);
    field_inv(&inv_leading, &work_row[first_nonzero], tracker->ctx->field_id, ctx_ptr);
    
    // 存储标准化的行到基中
    slong base_row_start = tracker->current_rank * ncols;
    for (slong j = 0; j < ncols; j++) {
        slong idx = base_row_start + j;
        if (j < first_nonzero) {
            // 主元之前的位置应该是0
            field_set_zero(&tracker->reduced_rows[idx], tracker->ctx->field_id, ctx_ptr);
        } else {
            // 标准化并存储
            field_mul(&tracker->reduced_rows[idx], &work_row[j], &inv_leading, 
                     tracker->ctx->field_id, ctx_ptr);
        }
    }
    
    // 更新跟踪信息
    tracker->pivot_cols[tracker->current_rank] = first_nonzero;
    tracker->selected_indices[tracker->current_rank] = new_row_idx;
    tracker->current_rank++;
    
    // 清理
    field_clear_elem(&inv_leading, tracker->ctx->field_id, ctx_ptr);
    for (slong j = 0; j < ncols; j++) {
        field_clear_elem(&work_row[j], tracker->ctx->field_id, ctx_ptr);
    }
    flint_free(work_row);
    
    return 1;
}


// 优化的 find_fq_optimal_maximal_rank_submatrix
void find_fq_optimal_maximal_rank_submatrix(fq_mvpoly_t ***full_matrix, 
                                           slong nrows, slong ncols,
                                           slong **row_indices_out, 
                                           slong **col_indices_out,
                                           slong *num_rows, slong *num_cols,
                                           slong npars) {
    printf("Finding maximal rank submatrix using optimized incremental method...\n");
    
    // 获取上下文
    const fq_nmod_ctx_struct *ctx = full_matrix[0][0]->ctx;
    
    // 初始化统一字段上下文
    field_ctx_t unified_ctx;
    field_ctx_init(&unified_ctx, ctx);
    void *ctx_ptr = (unified_ctx.field_id == FIELD_ID_NMOD) ? 
                   (void*)&unified_ctx.ctx.nmod_ctx : 
                   (void*)unified_ctx.ctx.fq_ctx;
    
    // 1. 计算度数并排序
    fq_index_degree_pair *row_degrees = (fq_index_degree_pair*) flint_malloc(nrows * sizeof(fq_index_degree_pair));
    fq_index_degree_pair *col_degrees = (fq_index_degree_pair*) flint_malloc(ncols * sizeof(fq_index_degree_pair));
    
    for (slong i = 0; i < nrows; i++) {
        row_degrees[i].index = i;
        row_degrees[i].degree = compute_fq_row_max_total_degree(full_matrix[i], ncols, npars);
    }
    
    for (slong j = 0; j < ncols; j++) {
        col_degrees[j].index = j;
        col_degrees[j].degree = compute_fq_col_max_total_degree(full_matrix, j, nrows, npars);
    }
    
    
    qsort(row_degrees, nrows, sizeof(fq_index_degree_pair), compare_fq_degrees);
    qsort(col_degrees, ncols, sizeof(fq_index_degree_pair), compare_fq_degrees);
    
    // 2. 评估矩阵并转换为统一格式
    //printf("Evaluating matrix and converting to unified format...\n");
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
            evaluate_fq_mvpoly_at_params(val, full_matrix[i][j], param_vals);
            fq_nmod_to_field_elem(&unified_mat[idx], val, &unified_ctx);
            fq_nmod_clear(val, ctx);
        }
    }
    
    clock_t conv_end = clock();
    //printf("Conversion time: %.3f seconds\n", (double)(conv_end - conv_start) / CLOCKS_PER_SEC);
    
    // 3. 选择行
    //printf("Selecting rows using incremental Gaussian elimination...\n");
    clock_t row_start = clock();
    
    unified_row_basis_tracker_t row_tracker;
    unified_row_basis_tracker_init(&row_tracker, FLINT_MIN(nrows, ncols), ncols, &unified_ctx);
    
    for (slong i = 0; i < nrows && row_tracker.current_rank < FLINT_MIN(nrows, ncols); i++) {
        slong row_idx = row_degrees[i].index;
        unified_try_add_row_to_basis(&row_tracker, unified_mat, row_idx, ncols);
    }
    
    clock_t row_end = clock();
    //printf("Selected %ld linearly independent rows (time: %.3f seconds)\n", row_tracker.current_rank, (double)(row_end - row_start) / CLOCKS_PER_SEC);
    
    // 4. 选择列
    //printf("Selecting columns...\n");
    clock_t col_start = clock();
    
    // 构建选定行的子矩阵
    field_elem_u *selected_submat = (field_elem_u*) flint_malloc(row_tracker.current_rank * ncols * sizeof(field_elem_u));
    for (slong i = 0; i < row_tracker.current_rank * ncols; i++) {
        field_init_elem(&selected_submat[i], unified_ctx.field_id, ctx_ptr);
    }
    
    for (slong i = 0; i < row_tracker.current_rank; i++) {
        for (slong j = 0; j < ncols; j++) {
            slong src_idx = row_tracker.selected_indices[i] * ncols + j;
            slong dst_idx = i * ncols + j;
            field_set_elem(&selected_submat[dst_idx], &unified_mat[src_idx],
                          unified_ctx.field_id, ctx_ptr);
        }
    }
    
    // 转置以选择列
    field_elem_u *transposed = (field_elem_u*) flint_malloc(ncols * row_tracker.current_rank * sizeof(field_elem_u));
    for (slong i = 0; i < ncols * row_tracker.current_rank; i++) {
        field_init_elem(&transposed[i], unified_ctx.field_id, ctx_ptr);
    }
    
    for (slong i = 0; i < row_tracker.current_rank; i++) {
        for (slong j = 0; j < ncols; j++) {
            slong src_idx = i * ncols + j;
            slong dst_idx = j * row_tracker.current_rank + i;
            field_set_elem(&transposed[dst_idx], &selected_submat[src_idx],
                          unified_ctx.field_id, ctx_ptr);
        }
    }
    
    // 选择列基
    unified_row_basis_tracker_t col_tracker;
    unified_row_basis_tracker_init(&col_tracker, ncols, row_tracker.current_rank, &unified_ctx);
    
    for (slong j = 0; j < ncols && col_tracker.current_rank < row_tracker.current_rank; j++) {
        slong col_idx = col_degrees[j].index;
        unified_try_add_row_to_basis(&col_tracker, transposed, col_idx, row_tracker.current_rank);
    }
    
    clock_t col_end = clock();
    //printf("Selected %ld linearly independent columns (time: %.3f seconds)\n", col_tracker.current_rank, (double)(col_end - col_start) / CLOCKS_PER_SEC);
    
    // 5. 构建最终结果
    slong final_size = FLINT_MIN(row_tracker.current_rank, col_tracker.current_rank);
    
    *row_indices_out = (slong*) flint_malloc(final_size * sizeof(slong));
    *col_indices_out = (slong*) flint_malloc(final_size * sizeof(slong));
    
    for (slong i = 0; i < final_size; i++) {
        (*row_indices_out)[i] = row_tracker.selected_indices[i];
        (*col_indices_out)[i] = col_tracker.selected_indices[i];
    }
    
    *num_rows = final_size;
    *num_cols = final_size;
    
    // 6. 验证结果
    //printf("Verifying final %ld x %ld submatrix...\n", final_size, final_size);
    fq_nmod_mat_t final_mat;
    fq_nmod_mat_init(final_mat, final_size, final_size, ctx);
    
    for (slong i = 0; i < final_size; i++) {
        for (slong j = 0; j < final_size; j++) {
            fq_nmod_t temp;
            fq_nmod_init(temp, ctx);
            slong idx = (*row_indices_out)[i] * ncols + (*col_indices_out)[j];
            field_elem_to_fq_nmod(temp, &unified_mat[idx], &unified_ctx);
            fq_nmod_set(fq_nmod_mat_entry(final_mat, i, j), temp, ctx);
            fq_nmod_clear(temp, ctx);
        }
    }
    
    slong final_rank = fq_nmod_mat_rank(final_mat, ctx);
    printf("Final rank: %ld/%ld %s\n", final_rank, final_size, (final_rank == final_size) ? "✓" : "✗");
    /*
    // 打印选择的索引（用于调试）
    printf("Selected rows: ");
    for (slong i = 0; i < final_size; i++) {
        printf("%ld ", (*row_indices_out)[i]);
    }
    printf("\nSelected cols: ");
    for (slong i = 0; i < final_size; i++) {
        printf("%ld ", (*col_indices_out)[i]);
    }
    */
    // 清理
    unified_row_basis_tracker_clear(&row_tracker);
    unified_row_basis_tracker_clear(&col_tracker);
    flint_free(row_degrees);
    flint_free(col_degrees);
    
    for (slong i = 0; i < npars; i++) {
        fq_nmod_clear(param_vals[i], ctx);
    }
    flint_free(param_vals);
    
    // 清理统一格式矩阵
    for (slong i = 0; i < nrows * ncols; i++) {
        field_clear_elem(&unified_mat[i], unified_ctx.field_id, ctx_ptr);
    }
    flint_free(unified_mat);
    
    for (slong i = 0; i < row_tracker.current_rank * ncols; i++) {
        field_clear_elem(&selected_submat[i], unified_ctx.field_id, ctx_ptr);
    }
    flint_free(selected_submat);
    
    for (slong i = 0; i < ncols * row_tracker.current_rank; i++) {
        field_clear_elem(&transposed[i], unified_ctx.field_id, ctx_ptr);
    }
    flint_free(transposed);
    
    fq_nmod_mat_clear(final_mat, ctx);
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
    
    // Collect monomials
    typedef struct {
        slong *exp;
        slong idx;
    } monom_t;
    
    monom_t *x_monoms = (monom_t*) flint_malloc(dixon_poly->nterms * sizeof(monom_t));
    monom_t *dual_monoms = (monom_t*) flint_malloc(dixon_poly->nterms * sizeof(monom_t));
    slong nx_monoms = 0, ndual_monoms = 0;
    
    // Collect unique monomials that satisfy degree bounds
    for (slong i = 0; i < dixon_poly->nterms; i++) {
        // Check if this term satisfies degree bounds
        int valid = 1;
        if (dixon_poly->terms[i].var_exp) {
            for (slong k = 0; k < nvars; k++) {
                if (dixon_poly->terms[i].var_exp[k] >= d0[k] || 
                    dixon_poly->terms[i].var_exp[nvars + k] >= d1[k]) {
                    valid = 0;
                    break;
                }
            }
        }
        
        if (!valid) continue;
        
        // Check x-monomial
        int found = 0;
        for (slong j = 0; j < nx_monoms; j++) {
            int same = 1;
            for (slong k = 0; k < nvars; k++) {
                if (x_monoms[j].exp[k] != dixon_poly->terms[i].var_exp[k]) {
                    same = 0;
                    break;
                }
            }
            if (same) {
                found = 1;
                break;
            }
        }
        if (!found) {
            x_monoms[nx_monoms].exp = (slong*) flint_calloc(nvars, sizeof(slong));
            x_monoms[nx_monoms].idx = nx_monoms;
            for (slong k = 0; k < nvars; k++) {
                x_monoms[nx_monoms].exp[k] = dixon_poly->terms[i].var_exp[k];
            }
            nx_monoms++;
        }
        
        // Check ~x-monomial
        found = 0;
        for (slong j = 0; j < ndual_monoms; j++) {
            int same = 1;
            for (slong k = 0; k < nvars; k++) {
                if (dual_monoms[j].exp[k] != dixon_poly->terms[i].var_exp[nvars + k]) {
                    same = 0;
                    break;
                }
            }
            if (same) {
                found = 1;
                break;
            }
        }
        if (!found) {
            dual_monoms[ndual_monoms].exp = (slong*) flint_calloc(nvars, sizeof(slong));
            dual_monoms[ndual_monoms].idx = ndual_monoms;
            for (slong k = 0; k < nvars; k++) {
                dual_monoms[ndual_monoms].exp[k] = dixon_poly->terms[i].var_exp[nvars + k];
            }
            ndual_monoms++;
        }
    }
    
    printf("Found %ld x-monomials and %ld ~x-monomials (after degree filtering)\n", 
           nx_monoms, ndual_monoms);
    
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
        full_matrix[i] = (fq_mvpoly_t**) flint_malloc(ndual_monoms * sizeof(fq_mvpoly_t*));
        for (slong j = 0; j < ndual_monoms; j++) {
            full_matrix[i][j] = (fq_mvpoly_t*) flint_malloc(sizeof(fq_mvpoly_t));
            fq_mvpoly_init(full_matrix[i][j], 0, npars, dixon_poly->ctx);
        }
    }
    
    // Fill the coefficient matrix
    for (slong t = 0; t < dixon_poly->nterms; t++) {
        // Check degree bounds again
        int valid = 1;
        if (dixon_poly->terms[t].var_exp) {
            for (slong k = 0; k < nvars; k++) {
                if (dixon_poly->terms[t].var_exp[k] >= d0[k] || 
                    dixon_poly->terms[t].var_exp[nvars + k] >= d1[k]) {
                    valid = 0;
                    break;
                }
            }
        }
        
        if (!valid) continue;
        
        // Find row index (x-monomial)
        slong row = -1;
        for (slong i = 0; i < nx_monoms; i++) {
            int same = 1;
            for (slong k = 0; k < nvars; k++) {
                if (x_monoms[i].exp[k] != dixon_poly->terms[t].var_exp[k]) {
                    same = 0;
                    break;
                }
            }
            if (same) {
                row = i;
                break;
            }
        }
        
        // Find column index (~x-monomial)
        slong col = -1;
        for (slong j = 0; j < ndual_monoms; j++) {
            int same = 1;
            for (slong k = 0; k < nvars; k++) {
                if (dual_monoms[j].exp[k] != dixon_poly->terms[t].var_exp[nvars + k]) {
                    same = 0;
                    break;
                }
            }
            if (same) {
                col = j;
                break;
            }
        }
        
        // Add the parameter polynomial to the matrix entry
        if (row >= 0 && col >= 0) {
            fq_mvpoly_add_term(full_matrix[row][col], NULL, 
                              dixon_poly->terms[t].par_exp, dixon_poly->terms[t].coeff);
        }
    }
    
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
                if (full_matrix[i][j]->nterms > 0) {
                    fq_nmod_set(fq_nmod_mat_entry(eval_mat, i, j), 
                               full_matrix[i][j]->terms[0].coeff, dixon_poly->ctx);
                } else {
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
        // Call submatrix finding function
        find_fq_optimal_maximal_rank_submatrix(full_matrix, nx_monoms, ndual_monoms,
                                              &row_idx_array, &col_idx_array, 
                                              &num_rows, &num_cols,
                                              npars);
    }
    
    // Take square part if needed
    slong submat_rank = FLINT_MIN(num_rows, num_cols);
    
    if (submat_rank == 0) {
        printf("Warning: Matrix has rank 0\n");
        *matrix_size = 0;
        
        // Cleanup
        for (slong i = 0; i < nx_monoms; i++) {
            for (slong j = 0; j < ndual_monoms; j++) {
                fq_mvpoly_clear(full_matrix[i][j]);
                flint_free(full_matrix[i][j]);
            }
            flint_free(full_matrix[i]);
        }
        flint_free(full_matrix);
        
        if (row_idx_array) flint_free(row_idx_array);
        if (col_idx_array) flint_free(col_idx_array);
        
        for (slong i = 0; i < nx_monoms; i++) {
            flint_free(x_monoms[i].exp);
        }
        for (slong j = 0; j < ndual_monoms; j++) {
            flint_free(dual_monoms[j].exp);
        }
        flint_free(x_monoms);
        flint_free(dual_monoms);
        flint_free(d0);
        flint_free(d1);
        return;
    }
    
    printf("Extracted submatrix of size %ld x %ld\n", submat_rank, submat_rank);
    
    // Build the output submatrix
    *coeff_matrix = (fq_mvpoly_t**) flint_malloc(submat_rank * sizeof(fq_mvpoly_t*));
    for (slong i = 0; i < submat_rank; i++) {
        (*coeff_matrix)[i] = (fq_mvpoly_t*) flint_malloc(submat_rank * sizeof(fq_mvpoly_t));
        for (slong j = 0; j < submat_rank; j++) {
            fq_mvpoly_copy(&(*coeff_matrix)[i][j], full_matrix[row_idx_array[i]][col_idx_array[j]]);
        }
    }
    
    // Copy indices
    for (slong i = 0; i < submat_rank; i++) {
        row_indices[i] = row_idx_array[i];
        col_indices[i] = col_idx_array[i];
    }
    *matrix_size = submat_rank;
    
    // Cleanup
    for (slong i = 0; i < nx_monoms; i++) {
        for (slong j = 0; j < ndual_monoms; j++) {
            fq_mvpoly_clear(full_matrix[i][j]);
            flint_free(full_matrix[i][j]);
        }
        flint_free(full_matrix[i]);
    }
    flint_free(full_matrix);
    
    if (row_idx_array) flint_free(row_idx_array);
    if (col_idx_array) flint_free(col_idx_array);
    
    for (slong i = 0; i < nx_monoms; i++) {
        flint_free(x_monoms[i].exp);
    }
    for (slong j = 0; j < ndual_monoms; j++) {
        flint_free(dual_monoms[j].exp);
    }
    flint_free(x_monoms);
    flint_free(dual_monoms);
    flint_free(d0);
    flint_free(d1);
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
    }
    
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
        
        det_method_t coeff_method = DET_METHOD_INTERPOLATION; // DET_METHOD_RECURSIVE DET_METHOD_KRONECKER DET_METHOD_INTERPOLATION DET_METHOD_HUANG
        char *method_env = getenv("DIXON_DET_METHOD");
        if (method_env) {
            coeff_method = atoi(method_env);
        }

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
    fq_mvpoly_clear(&d_poly);
    
    printf("\n=== Dixon Resultant Computation Complete ===\n");
}

#endif