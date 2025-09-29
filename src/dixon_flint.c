/*
 * dixon_flint.c - Dixon Resultant Implementation for Finite Extension Fields
 *
 * This file contains the complete implementation of Dixon resultant computation
 * over finite fields using FLINT library for polynomial arithmetic and matrix operations.
 */

#include "dixon_flint.h"

// Global method selection variable definition
det_method_t dixon_global_method = -1;

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
            // Directly use fq_mvpoly_sub to compute difference
            fq_mvpoly_t diff;
            fq_mvpoly_sub(&diff, &(*original_matrix)[i][j], &(*original_matrix)[i-1][j]);
            
            // Perform division
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

// Extended tracker structure with pre-allocated workspace
// Initialize optimized tracker
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
    
    // Allocate main storage
    tracker->reduced_rows = (field_elem_u*) flint_calloc(max_size * ncols, sizeof(field_elem_u));
    tracker->pivot_cols = (slong*) flint_calloc(max_size, sizeof(slong));
    tracker->selected_indices = (slong*) flint_calloc(max_size, sizeof(slong));
    
    // Pre-allocate workspace
    tracker->work_row = (field_elem_u*) flint_malloc(ncols * sizeof(field_elem_u));
    tracker->temp_vars = (field_elem_u*) flint_malloc(4 * sizeof(field_elem_u)); // factor, temp, pivot_val, neg_temp
    
    // Initialize all field elements
    for (slong i = 0; i < max_size * ncols; i++) {
        field_init_elem(&tracker->reduced_rows[i], ctx->field_id, ctx_ptr);
        field_set_zero(&tracker->reduced_rows[i], ctx->field_id, ctx_ptr);
    }
    
    // Initialize workspace
    for (slong j = 0; j < ncols; j++) {
        field_init_elem(&tracker->work_row[j], ctx->field_id, ctx_ptr);
    }
    for (slong i = 0; i < 4; i++) {
        field_init_elem(&tracker->temp_vars[i], ctx->field_id, ctx_ptr);
    }
    tracker->workspace_initialized = 1;
    
    // Initialize pivot columns to -1
    for (slong i = 0; i < max_size; i++) {
        tracker->pivot_cols[i] = -1;
        tracker->selected_indices[i] = -1;
    }
}

// Clear optimized tracker
static void unified_row_basis_tracker_clear(unified_row_basis_tracker_t *tracker) {
    if (!tracker->initialized) return;
    
    void *ctx_ptr = (tracker->ctx->field_id == FIELD_ID_NMOD) ? 
                   (void*)&tracker->ctx->ctx.nmod_ctx : 
                   (void*)tracker->ctx->ctx.fq_ctx;
    
    // Clear main storage
    if (tracker->reduced_rows) {
        for (slong i = 0; i < tracker->max_size * tracker->ncols; i++) {
            field_clear_elem(&tracker->reduced_rows[i], tracker->ctx->field_id, ctx_ptr);
        }
        flint_free(tracker->reduced_rows);
    }
    
    // Clear workspace
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

// Optimized version of adding row to basis - mathematical logic completely unchanged
static int unified_try_add_row_to_basis(unified_row_basis_tracker_t *tracker, 
                                         const field_elem_u *unified_mat,
                                         slong new_row_idx, slong ncols) {
    if (!tracker->initialized || tracker->current_rank >= tracker->max_size) {
        return 0;
    }
    
    void *ctx_ptr = (tracker->ctx->field_id == FIELD_ID_NMOD) ? 
                   (void*)&tracker->ctx->ctx.nmod_ctx : 
                   (void*)tracker->ctx->ctx.fq_ctx;
    
    // Use pre-allocated work row to avoid allocation each time
    field_elem_u *work_row = tracker->work_row;
    
    // Copy input row to work row (reuse allocated space)
    for (slong j = 0; j < ncols; j++) {
        field_set_elem(&work_row[j], &unified_mat[new_row_idx * ncols + j], 
                      tracker->ctx->field_id, ctx_ptr);
    }
    
    // Use pre-allocated temporary variables
    field_elem_u *factor = &tracker->temp_vars[0];
    field_elem_u *temp = &tracker->temp_vars[1];
    field_elem_u *pivot_val = &tracker->temp_vars[2];
    field_elem_u *neg_temp = &tracker->temp_vars[3];
    
    // Perform elimination for each existing basis vector (mathematical logic unchanged)
    for (slong i = 0; i < tracker->current_rank; i++) {
        slong pivot_col = tracker->pivot_cols[i];
        
        // Validate pivot column index
        if (pivot_col < 0 || pivot_col >= ncols) continue;
        
        // If work row is non-zero at pivot position, perform elimination
        if (!field_is_zero(&work_row[pivot_col], tracker->ctx->field_id, ctx_ptr)) {
            // Get pivot value of basis vector
            slong base_idx = i * ncols + pivot_col;
            field_set_elem(pivot_val, &tracker->reduced_rows[base_idx], 
                          tracker->ctx->field_id, ctx_ptr);
            
            // Calculate elimination factor = work_row[pivot_col] / pivot_val
            field_inv(temp, pivot_val, tracker->ctx->field_id, ctx_ptr);
            field_mul(factor, &work_row[pivot_col], temp, tracker->ctx->field_id, ctx_ptr);
            
            // Perform elimination: work_row -= factor * basis_row[i]
            // Optimization: pre-determine if neg_temp is needed to avoid conditional branches in inner loop
            int use_direct_add = (tracker->ctx->field_id >= FIELD_ID_GF28 && 
                                 tracker->ctx->field_id <= FIELD_ID_GF2128);
            
            for (slong j = 0; j < ncols; j++) {
                slong idx = i * ncols + j;
                field_mul(temp, factor, &tracker->reduced_rows[idx], 
                         tracker->ctx->field_id, ctx_ptr);
                
                if (use_direct_add) {
                    // For GF(2^n), subtraction equals addition
                    field_add(&work_row[j], &work_row[j], temp, 
                             tracker->ctx->field_id, ctx_ptr);
                } else {
                    // For other fields, use actual subtraction
                    field_neg(neg_temp, temp, tracker->ctx->field_id, ctx_ptr);
                    field_add(&work_row[j], &work_row[j], neg_temp, 
                             tracker->ctx->field_id, ctx_ptr);
                }
            }
        }
    }
    
    // Find first non-zero position (mathematical logic unchanged)
    slong first_nonzero = -1;
    for (slong j = 0; j < ncols; j++) {
        if (!field_is_zero(&work_row[j], tracker->ctx->field_id, ctx_ptr)) {
            first_nonzero = j;
            break;
        }
    }
    
    // If all zeros, then linearly dependent
    if (first_nonzero == -1) {
        return 0;  // Work row doesn't need cleanup because it's pre-allocated
    }
    
    // Normalize work row (make first non-zero element 1)
    field_inv(temp, &work_row[first_nonzero], tracker->ctx->field_id, ctx_ptr);
    
    // Store normalized row to basis
    slong base_row_start = tracker->current_rank * ncols;
    for (slong j = 0; j < ncols; j++) {
        slong idx = base_row_start + j;
        if (j < first_nonzero) {
            // Positions before pivot should be 0
            field_set_zero(&tracker->reduced_rows[idx], tracker->ctx->field_id, ctx_ptr);
        } else {
            // Normalize and store
            field_mul(&tracker->reduced_rows[idx], &work_row[j], temp, 
                     tracker->ctx->field_id, ctx_ptr);
        }
    }
    
    // Update tracking information (mathematical logic unchanged)
    tracker->pivot_cols[tracker->current_rank] = first_nonzero;
    tracker->selected_indices[tracker->current_rank] = new_row_idx;
    tracker->current_rank++;
    
    return 1;
}

// Compute maximum total degree of a column in selected rows submatrix
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

// Compute maximum total degree of a row in selected columns submatrix
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

// Check if two index arrays are identical
static int indices_equal(slong *indices1, slong *indices2, slong size) {
    for (slong i = 0; i < size; i++) {
        if (indices1[i] != indices2[i]) {
            return 0;
        }
    }
    return 1;
}
// Find pivot rows for maximal rank submatrix - ensure linear independence
// More efficient version: using incremental rank checking
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
// Optimized find_pivot_rows_simple - directly calls nmod version
// ============================================================================
void find_pivot_rows_simple(slong **selected_rows_out, slong *num_selected,
                                        const field_elem_u *unified_mat, 
                                        slong nrows, slong ncols,
                                        field_ctx_t *ctx) {
    
    // ============================================================================
    // Prime field fast path: directly use nmod_fixed implementation (optimal performance)
    // ============================================================================
    if (ctx->field_id == FIELD_ID_NMOD) {
        //printf("Using direct nmod_fixed implementation for prime field (optimal performance)\n");
        clock_t start = clock();
        
        // Convert to nmod_mat format
        nmod_mat_t nmod_mat;
        nmod_mat_init(nmod_mat, nrows, ncols, ctx->ctx.nmod_ctx.n);
        
        // Efficient data copy (direct access to nmod field)
        for (slong i = 0; i < nrows; i++) {
            for (slong j = 0; j < ncols; j++) {
                nmod_mat_entry(nmod_mat, i, j) = unified_mat[i * ncols + j].nmod;
            }
        }
        
        // Directly use nmod_fixed version (avoid adaptive layer selection overhead)
        find_pivot_rows_nmod_fixed(selected_rows_out, num_selected, nmod_mat);
        
        clock_t end = clock();
        double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
        
        nmod_mat_clear(nmod_mat);
        //printf("Prime field computation completed in %.4f seconds using direct nmod path\n", elapsed);
        return;
    }
    
    // ============================================================================
    // Non-prime field: use unified interface implementation
    // ============================================================================
    //printf("Using unified interface for non-prime field computation\n");
    
    void *ctx_ptr = (ctx->field_id == FIELD_ID_FQ_ZECH) ? 
                   (void*)ctx->ctx.zech_ctx : 
                   (void*)ctx->ctx.fq_ctx;
    
    slong max_rank = FLINT_MIN(nrows, ncols);
    slong *selected_rows = (slong*) flint_malloc(max_rank * sizeof(slong));
    slong rank = 0;
    
    // Lightweight linear independence checker
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
    
    // Initialize all field elements
    for (slong i = 0; i < max_rank * ncols; i++) {
        field_init_elem(&tracker.rows[i], ctx->field_id, ctx_ptr);
    }
    
    clock_t start = clock();
    
    // Process row by row
    for (slong row = 0; row < nrows && rank < max_rank; row++) {
        // Create test vector
        field_elem_u *test_vec = (field_elem_u*) flint_malloc(ncols * sizeof(field_elem_u));
        for (slong j = 0; j < ncols; j++) {
            field_init_elem(&test_vec[j], ctx->field_id, ctx_ptr);
            field_set_elem(&test_vec[j], &unified_mat[row * ncols + j], ctx->field_id, ctx_ptr);
        }
        
        // Eliminate using existing basis vectors
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
        
        // Find first non-zero position
        slong pivot_pos = -1;
        for (slong j = 0; j < ncols; j++) {
            if (!field_is_zero(&test_vec[j], ctx->field_id, ctx_ptr)) {
                pivot_pos = j;
                break;
            }
        }
        
        if (pivot_pos >= 0) {
            // Linearly independent, add to basis
            selected_rows[rank] = row;
            tracker.pivot_positions[tracker.count] = pivot_pos;
            
            // Normalize and store
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
        
        // Clean up test vector
        for (slong j = 0; j < ncols; j++) {
            field_clear_elem(&test_vec[j], ctx->field_id, ctx_ptr);
        }
        flint_free(test_vec);
    }
    
    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    
    // Set output
    if (rank > 0) {
        *selected_rows_out = (slong*) flint_malloc(rank * sizeof(slong));
        memcpy(*selected_rows_out, selected_rows, rank * sizeof(slong));
    } else {
        *selected_rows_out = NULL;
    }
    *num_selected = rank;
    
    // Cleanup
    for (slong i = 0; i < max_rank * ncols; i++) {
        field_clear_elem(&tracker.rows[i], ctx->field_id, ctx_ptr);
    }
    flint_free(tracker.rows);
    flint_free(tracker.pivot_positions);
    flint_free(selected_rows);
    
    //printf("Non-prime field computation completed in %.4f seconds using unified interface\n", elapsed);
}
// Optimized find_fq_optimal_maximal_rank_submatrix
void find_fq_optimal_maximal_rank_submatrix(fq_mvpoly_t ***full_matrix, 
                                           slong nrows, slong ncols,
                                           slong **row_indices_out, 
                                           slong **col_indices_out,
                                           slong *num_rows, slong *num_cols,
                                           slong npars) {
    printf("Finding maximal rank submatrix using iterative optimization method...\n");
    
    // Get context
    const fq_nmod_ctx_struct *ctx = NULL;
    for (slong i = 0; i < nrows && !ctx; i++) {
        for (slong j = 0; j < ncols && !ctx; j++) {
            if (full_matrix[i][j] != NULL) {
                ctx = full_matrix[i][j]->ctx;
                break;
            }
        }
    }
    
    // Initialize unified field context
    field_ctx_t unified_ctx;
    field_ctx_init(&unified_ctx, ctx);
    void *ctx_ptr = (unified_ctx.field_id == FIELD_ID_NMOD) ? 
                   (void*)&unified_ctx.ctx.nmod_ctx : 
                   (void*)unified_ctx.ctx.fq_ctx;
    
    // Evaluate matrix and convert to unified format
    printf("Evaluating matrix and converting to unified format...\n");
    clock_t conv_start = clock();
    /*
    // Generate evaluation parameters
    fq_nmod_t *param_vals = (fq_nmod_t*) flint_malloc(npars * sizeof(fq_nmod_t));
    for (slong i = 0; i < npars; i++) {
        fq_nmod_init(param_vals[i], ctx);
        fq_nmod_set_si(param_vals[i], 2 + i, ctx);
    }
    */
        
    // Generate evaluation parameters using t^(i+1) + c*(i+1) where t is the multiplicative group generator
    fq_nmod_t *param_vals = (fq_nmod_t*) flint_malloc(npars * sizeof(fq_nmod_t));
    
    // Get field characteristic
    mp_limb_t p = fq_nmod_ctx_prime(ctx);
    slong d = fq_nmod_ctx_degree(ctx);
    //printf("Field: F_%lu^%ld, characteristic p = %lu\n", p, d, p);
    
    // Get the multiplicative group generator (primitive element)
    fq_nmod_t gen;
    fq_nmod_init(gen, ctx);
    fq_nmod_gen(gen, ctx);  // This gives us the generator t of F_{p^d}
    
    //printf("Field generator t = ");
    //fq_nmod_print_pretty(gen, ctx);
    //printf("\n");
    
    //printf("Using evaluation points t^(i+1) + c*(i+1) to avoid using 1:\n");
    
    // Choose a multiplier c that gives good distribution
    slong multiplier = (p > 3) ? ((p + 1) / 3) : 2;
    
    // Generate evaluation points using i+1 instead of i to avoid t^0 = 1
    for (slong i = 0; i < npars; i++) {
        fq_nmod_init(param_vals[i], ctx);
        
        fq_nmod_t power_val, addend;
        fq_nmod_init(power_val, ctx);
        fq_nmod_init(addend, ctx);
        
        slong exp = i + 1;  // Use i+1 so we start from t^1, not t^0
        
        // Use t^(i+1)
        fq_nmod_pow_ui(power_val, gen, exp, ctx);
        
        // Use multiplier*(i+1) mod p
        slong coeff = (multiplier * exp) % p;
        if (coeff == 0) coeff = exp % p;  // Fallback to exp if multiplier gives 0
        if (coeff == 0) coeff = 2;       // Final fallback
        
        fq_nmod_set_si(addend, coeff, ctx);
        fq_nmod_add(param_vals[i], power_val, addend, ctx);  // t^(i+1) + c*(i+1)
        
        // Output the evaluation point
        //printf("  Parameter %ld: ", i);
        //fq_nmod_print_pretty(param_vals[i], ctx);
        //printf(" (t^%ld + %ld)\n", exp, coeff);
        
        fq_nmod_clear(power_val, ctx);
        fq_nmod_clear(addend, ctx);
    }
    
    fq_nmod_clear(gen, ctx);
        
    // Allocate unified format matrix
    field_elem_u *unified_mat = (field_elem_u*) flint_malloc(nrows * ncols * sizeof(field_elem_u));
    
    // Initialize and evaluate
    for (slong i = 0; i < nrows; i++) {
        for (slong j = 0; j < ncols; j++) {
            slong idx = i * ncols + j;
            field_init_elem(&unified_mat[idx], unified_ctx.field_id, ctx_ptr);
            
            // Evaluate polynomial
            fq_nmod_t val;
            fq_nmod_init(val, ctx);
            // Safe evaluation: check for NULL pointer
            if (full_matrix[i][j] == NULL) {
                fq_nmod_zero(val, ctx);  // NULL represents zero polynomial
            } else {
                evaluate_fq_mvpoly_at_params(val, full_matrix[i][j], param_vals);
            }            fq_nmod_to_field_elem(&unified_mat[idx], val, &unified_ctx);
            fq_nmod_clear(val, ctx);
        }
    }
    
    clock_t conv_end = clock();
    //printf("Conversion time: %.3f seconds\n", (double)(conv_end - conv_start) / CLOCKS_PER_SEC);
    
    // Initialize arrays to store selection results
    slong *current_row_indices = NULL;
    slong *current_col_indices = NULL;
    slong *prev_row_indices = NULL;
    slong *prev_col_indices = NULL;
    slong current_size = 0;
    slong prev_size = 0;
    
    // Iterative optimization
    const slong MAX_ITERATIONS = 10;
    slong iteration = 0;
    int converged = 0;
    
    clock_t iter_start = clock();
    
    while (iteration < MAX_ITERATIONS && !converged) {
        //printf("\n--- Iteration %ld ---\n", iteration + 1);
        
        // Save previous results
        if (iteration > 0) {
            prev_row_indices = (slong*) flint_malloc(current_size * sizeof(slong));
            prev_col_indices = (slong*) flint_malloc(current_size * sizeof(slong));
            memcpy(prev_row_indices, current_row_indices, current_size * sizeof(slong));
            memcpy(prev_col_indices, current_col_indices, current_size * sizeof(slong));
            prev_size = current_size;
        }
        
        slong row_rank_selected = 0;  // Record row selection rank
        

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
            
            // New addition: if too many rows selected, directly use transpose method for column selection
            if (current_size > 1000) {
                printf("Large matrix detected (%ld rows), using transpose method for column selection...\n", current_size);
                
                // Build transpose matrix of selected rows
                field_elem_u *transposed_mat = (field_elem_u*) flint_malloc(ncols * current_size * sizeof(field_elem_u));
                void *ctx_ptr = (unified_ctx.field_id == FIELD_ID_NMOD) ? 
                               (void*)&unified_ctx.ctx.nmod_ctx : 
                               (void*)unified_ctx.ctx.fq_ctx;
                
                // Initialize transpose matrix elements
                for (slong i = 0; i < ncols * current_size; i++) {
                    field_init_elem(&transposed_mat[i], unified_ctx.field_id, ctx_ptr);
                }
                
                // Fill transpose matrix: original matrix rows become transpose matrix columns
                for (slong i = 0; i < current_size; i++) {
                    slong orig_row = current_row_indices[i];
                    for (slong j = 0; j < ncols; j++) {
                        slong src_idx = orig_row * ncols + j;
                        slong dst_idx = j * current_size + i; // Transpose index
                        field_set_elem(&transposed_mat[dst_idx], &unified_mat[src_idx],
                                      unified_ctx.field_id, ctx_ptr);
                    }
                }
                
                // Find pivot rows for transpose matrix (i.e., pivot columns of original matrix)
                slong *selected_cols = NULL;
                slong num_selected_cols = 0;
                
                find_pivot_rows_simple(&selected_cols, &num_selected_cols, 
                                      transposed_mat, ncols, current_size, &unified_ctx);
                
                // Set column indices
                current_col_indices = (slong*) flint_malloc(num_selected_cols * sizeof(slong));
                memcpy(current_col_indices, selected_cols, num_selected_cols * sizeof(slong));
                
                // Update matrix size to smaller dimension
                current_size = FLINT_MIN(current_size, num_selected_cols);
                
                //printf("Selected %ld columns via transpose method\n", num_selected_cols);
                //printf("Final submatrix size: %ld x %ld\n", current_size, current_size);
                
                // Clean up transpose matrix
                for (slong i = 0; i < ncols * current_size; i++) {
                    field_clear_elem(&transposed_mat[i], unified_ctx.field_id, ctx_ptr);
                }
                flint_free(transposed_mat);
                flint_free(selected_rows);
                flint_free(selected_cols);
                
                // Directly exit iteration loop, avoid further optimization steps
                converged = 1;
                //printf("Large matrix optimization completed, skipping further iterations.\n");
                break;
            } else {
                // Normal case: free temporary arrays, continue iteration process
                flint_free(selected_rows);
            }
        } else {
            // Subsequent iterations: re-select rows based on current columns
            //printf("Re-selecting rows based on current columns...\n");
            
            // Calculate degrees of all rows based on currently selected columns
            fq_index_degree_pair *row_degrees = (fq_index_degree_pair*) flint_malloc(nrows * sizeof(fq_index_degree_pair));
            
            for (slong i = 0; i < nrows; i++) {
                row_degrees[i].index = i;
                row_degrees[i].degree = compute_fq_selected_cols_row_max_total_degree(
                    full_matrix, i, current_col_indices, current_size, npars);
            }
            
            qsort(row_degrees, nrows, sizeof(fq_index_degree_pair), compare_fq_degrees);
            
            // Build transpose submatrix of current columns for row selection
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
            
            // Select rows
            unified_row_basis_tracker_t row_tracker;
            unified_row_basis_tracker_init(&row_tracker, current_size, current_size, &unified_ctx);
            
            for (slong i = 0; i < nrows && row_tracker.current_rank < current_size; i++) {
                slong row_idx = row_degrees[i].index;
                unified_try_add_row_to_basis(&row_tracker, col_submat, row_idx, current_size);
            }
            
            // Update row indices
            flint_free(current_row_indices);
            current_row_indices = (slong*) flint_malloc(row_tracker.current_rank * sizeof(slong));
            memcpy(current_row_indices, row_tracker.selected_indices, row_tracker.current_rank * sizeof(slong));
            row_rank_selected = row_tracker.current_rank;
            current_size = row_tracker.current_rank;
            
            // Cleanup
            for (slong i = 0; i < nrows * current_size; i++) {
                field_clear_elem(&col_submat[i], unified_ctx.field_id, ctx_ptr);
            }
            flint_free(col_submat);
            unified_row_basis_tracker_clear(&row_tracker);
            flint_free(row_degrees);
            
            //printf("Selected %ld rows\n", current_size);
        }
        
        // Select columns based on current rows
        //printf("Selecting columns based on current rows...\n");
        
        // Calculate column degrees based on selected rows
        fq_index_degree_pair *col_degrees = (fq_index_degree_pair*) flint_malloc(ncols * sizeof(fq_index_degree_pair));
        
        for (slong j = 0; j < ncols; j++) {
            col_degrees[j].index = j;
            col_degrees[j].degree = compute_fq_selected_rows_col_max_total_degree(
                full_matrix, current_row_indices, current_size, j, npars);
        }
        
        qsort(col_degrees, ncols, sizeof(fq_index_degree_pair), compare_fq_degrees);
        
        // Build submatrix of selected rows and transpose
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
        
        // Select column basis
        unified_row_basis_tracker_t col_tracker;
        unified_row_basis_tracker_init(&col_tracker, ncols, current_size, &unified_ctx);
        
        for (slong j = 0; j < ncols && col_tracker.current_rank < current_size; j++) {
            slong col_idx = col_degrees[j].index;
            unified_try_add_row_to_basis(&col_tracker, transposed, col_idx, current_size);
        }
        
        slong col_rank_selected = col_tracker.current_rank;
        
        // Update column indices
        if (iteration == 0) {
            current_col_indices = (slong*) flint_malloc(col_tracker.current_rank * sizeof(slong));
        } else {
            flint_free(current_col_indices);
            current_col_indices = (slong*) flint_malloc(col_tracker.current_rank * sizeof(slong));
        }
        memcpy(current_col_indices, col_tracker.selected_indices, col_tracker.current_rank * sizeof(slong));
        current_size = FLINT_MIN(current_size, col_tracker.current_rank);
        
        //printf("Selected %ld columns\n", col_tracker.current_rank);
        
        // Cleanup
        for (slong i = 0; i < ncols * current_size; i++) {
            field_clear_elem(&transposed[i], unified_ctx.field_id, ctx_ptr);
        }
        flint_free(transposed);
        unified_row_basis_tracker_clear(&col_tracker);
        flint_free(col_degrees);
        
        // Check convergence
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
    
    // Set output
    *row_indices_out = (slong*) flint_malloc(current_size * sizeof(slong));
    *col_indices_out = (slong*) flint_malloc(current_size * sizeof(slong));
    
    memcpy(*row_indices_out, current_row_indices, current_size * sizeof(slong));
    memcpy(*col_indices_out, current_col_indices, current_size * sizeof(slong));
    
    *num_rows = current_size;
    *num_cols = current_size;
    
    // Verify results
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
    
    // Calculate and display maximum degree of final submatrix
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
    
    // Print selected indices
    printf("Selected rows: ");
    for (slong i = 0; i < current_size; i++) {
        printf("%ld ", (*row_indices_out)[i]);
    }
    printf("\nSelected cols: ");
    for (slong i = 0; i < current_size; i++) {
        printf("%ld ", (*col_indices_out)[i]);
    }
    printf("\n");
    
    // Cleanup
    flint_free(current_row_indices);
    flint_free(current_col_indices);
    
    for (slong i = 0; i < npars; i++) {
        fq_nmod_clear(param_vals[i], ctx);
    }
    flint_free(param_vals);
    
    // Clean up unified format matrix
    for (slong i = 0; i < nrows * ncols; i++) {
        field_clear_elem(&unified_mat[i], unified_ctx.field_id, ctx_ptr);
    }
    flint_free(unified_mat);
    
    fq_nmod_mat_clear(final_mat, ctx);
}

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
// Allocate single element on demand
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
    
    // Find corresponding matrix positions for each Dixon polynomial term
    for (slong t = 0; t < dixon_poly->nterms; t++) {
        if (!dixon_poly->terms[t].var_exp) continue;
        
        // Check degree bounds
        int valid = 1;
        for (slong k = 0; k < nvars; k++) {
            if (dixon_poly->terms[t].var_exp[k] >= d0[k] || 
                dixon_poly->terms[t].var_exp[nvars + k] >= d1[k]) {
                valid = 0;
                break;
            }
        }
        if (!valid) continue;
        
        // Find x-monomial row index
        slong row = -1;
        for (slong i = 0; i < nx_monoms; i++) {
            if (memcmp(x_monoms[i].exp, dixon_poly->terms[t].var_exp, 
                      nvars * sizeof(slong)) == 0) {
                row = i;
                break;
            }
        }
        
        // Find dual-monomial column index
        slong col = -1;
        for (slong j = 0; j < ndual_monoms; j++) {
            if (memcmp(dual_monoms[j].exp, &dixon_poly->terms[t].var_exp[nvars], 
                      nvars * sizeof(slong)) == 0) {
                col = j;
                break;
            }
        }
        
        // Only allocate memory for positions that are actually needed
        if (row >= 0 && col >= 0) {
            fq_mvpoly_t *entry = get_matrix_entry_lazy(full_matrix, row, col, 
                                                      npars, dixon_poly->ctx);
            fq_mvpoly_add_term_fast(entry, NULL, dixon_poly->terms[t].par_exp, 
                                   dixon_poly->terms[t].coeff);
        }
    }
}
// Optimized version of find_fq_optimal_maximal_rank_submatrix
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
        // Note: use calloc to initialize to NULL
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
        if (nx_monoms < 4 && ndual_monoms < 4 && expected_rows < 4 && expected_cols < 4) {
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
    

    // Safe cleanup code - fix NULL pointer segfault
    
    printf("Extracted submatrix of size %ld x %ld\n", submat_rank, submat_rank);
    
    // Build the output submatrix
    *coeff_matrix = (fq_mvpoly_t**) flint_malloc(submat_rank * sizeof(fq_mvpoly_t*));
    for (slong i = 0; i < submat_rank; i++) {
        (*coeff_matrix)[i] = (fq_mvpoly_t*) flint_malloc(submat_rank * sizeof(fq_mvpoly_t));
        for (slong j = 0; j < submat_rank; j++) {
            // Safe copy: check if source matrix element is NULL
            fq_mvpoly_t *source = full_matrix[row_idx_array[i]][col_idx_array[j]];
            if (source != NULL) {
                fq_mvpoly_copy(&(*coeff_matrix)[i][j], source);
            } else {
                // NULL element represents zero polynomial, directly initialize to zero
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
    
    // Safe cleanup - add NULL checks
    printf("Cleaning up matrix...\n");
    if (full_matrix) {
        for (slong i = 0; i < nx_monoms; i++) {
            if (full_matrix[i]) {
                for (slong j = 0; j < ndual_monoms; j++) {
                    if (full_matrix[i][j]) {  // Critical: check NULL pointer
                        fq_mvpoly_clear(full_matrix[i][j]);
                        flint_free(full_matrix[i][j]);
                    }
                }
                flint_free(full_matrix[i]);
            }
        }
        flint_free(full_matrix);
    }
    
    // Safe cleanup of other arrays
    if (row_idx_array) flint_free(row_idx_array);
    if (col_idx_array) flint_free(col_idx_array);
    
    // Clean up monomial arrays - note memory management here
    if (x_monoms) {
        // x_monoms[i].exp points to contiguous storage, no need to free individually
        flint_free(x_monoms);
    }
    if (dual_monoms) {
        // dual_monoms[j].exp points to contiguous storage, no need to free individually  
        flint_free(dual_monoms);
    }
    
    // Clean up degree bound arrays
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
    
    // Display analysis of original matrix
    //analyze_fq_matrix_mvpoly(M_mvpoly, nvars + 1, nvars + 1, "Original Cancellation");
    
    // Step 2: Perform row operations
    printf("\nStep 2: Perform Matrix Row Operations\n");
    fq_mvpoly_t **modified_M_mvpoly;
    perform_fq_matrix_row_operations_mvpoly(&modified_M_mvpoly, &M_mvpoly, nvars, npars);
    
    // Step 3: Compute determinant of modified matrix
    printf("\nStep 3: Compute determinant of modified matrix\n");
    fq_mvpoly_t d_poly;
    compute_fq_cancel_matrix_det(&d_poly, modified_M_mvpoly, nvars, npars, DET_METHOD_RECURSIVE);
    
    printf("Dixon polynomial has %ld terms\n", d_poly.nterms);
    
    if (d_poly.nterms <= 100) {
        fq_mvpoly_print_expanded(&d_poly, "Dixon", 1);
    } else {
        printf("Dixon polynomial too large to display (%ld terms)\n", d_poly.nterms);
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
        #ifdef _OPENMP
        if (npars > 1) {
            coeff_method = DET_METHOD_INTERPOLATION;
        } else 
        #endif
        if (matrix_size < 9) {
            coeff_method = DET_METHOD_RECURSIVE;
        } else {
            coeff_method = DET_METHOD_KRONECKER;
        }
        //coeff_method = DET_METHOD_INTERPOLATION;
        if (dixon_global_method != -1) {
            coeff_method = dixon_global_method;
            printf("Method overridden by global variable: %d\n", dixon_global_method);
        }
        printf("Using coefficient matrix determinant method %d:", coeff_method);
        
        compute_fq_coefficient_matrix_det(result, coeff_matrix, matrix_size,
                                         npars, polys[0].ctx, coeff_method, res_deg_bound);
        
        printf("Resultant polynomial has %ld terms\n", result->nterms);
        if (result->nterms <= 100) {
            fq_mvpoly_print(result, "Final Resultant");
        } else {
            printf("Final resultant too large to display (%ld terms)\n", result->nterms);
            
            // Display degree information
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

void fq_dixon_resultant_with_names(fq_mvpoly_t *result, fq_mvpoly_t *polys, 
                                  slong nvars, slong npars,
                                  char **var_names, char **par_names, 
                                  const char *gen_name) {

    printf("Dixon Resultant Computation over Finite Extension Fields\n");
    printf("========================================================\n");
    printf("\n=== Dixon Resultant over F_{p^d} ===\n");
    cleanup_unified_workspace();
    
    // Step 1: Build cancellation matrix
    printf("\nStep 1: Build Cancellation Matrix\n");
    fq_mvpoly_t **M_mvpoly;
    build_fq_cancellation_matrix_mvpoly(&M_mvpoly, polys, nvars, npars);
    
    // Step 2: Perform row operations
    printf("\nStep 2: Perform Matrix Row Operations\n");
    fq_mvpoly_t **modified_M_mvpoly;
    perform_fq_matrix_row_operations_mvpoly(&modified_M_mvpoly, &M_mvpoly, nvars, npars);
    //print_fq_matrix_mvpoly(modified_M_mvpoly, nvars + 1, nvars + 1, "After Row Operations (Detailed)", 1);
    // Step 3: Compute determinant of modified matrix
    printf("\nStep 3: Compute determinant of modified matrix\n");
    fq_mvpoly_t d_poly;
    compute_fq_cancel_matrix_det(&d_poly, modified_M_mvpoly, nvars, npars, DET_METHOD_RECURSIVE);
    
    printf("Dixon polynomial has %ld terms\n", d_poly.nterms);
    
    if (d_poly.nterms <= 100) {
        // Use the enhanced print function with proper names
        fq_mvpoly_print_with_names(&d_poly, "Dixon", var_names, par_names, gen_name, 1);
    } else {
        printf("Dixon polynomial too large to display (%ld terms)\n", d_poly.nterms);
    }

    // Continue with existing logic...
    
    // Step 4: Extract coefficient matrix
    printf("\nStep 4: Extract coefficient matrix\n");
    
    fq_mvpoly_t **coeff_matrix;
    slong *row_indices = (slong*) flint_malloc(d_poly.nterms * sizeof(slong));
    slong *col_indices = (slong*) flint_malloc(d_poly.nterms * sizeof(slong));
    slong matrix_size;
    
    extract_fq_coefficient_matrix_from_dixon(&coeff_matrix, row_indices, col_indices,
                                            &matrix_size, &d_poly, nvars, npars);
    //print_fq_matrix_mvpoly(coeff_matrix, matrix_size, matrix_size, "Dixon Matrix", 1);
    if (matrix_size > 0) {
        printf("\nStep 5: Compute determinant of coefficient matrix (fq_nmod)\n");
        
        slong res_deg_bound = compute_fq_dixon_resultant_degree_bound(polys, nvars+1, nvars, npars);
        printf("Resultant degree bound: %ld\n", res_deg_bound);
        
        det_method_t coeff_method;
        #ifdef _OPENMP
        if (npars > 1) {
            coeff_method = DET_METHOD_INTERPOLATION;
        } else 
        #endif
        if (matrix_size < 9) {
            coeff_method = DET_METHOD_RECURSIVE;
        } else {
            coeff_method = DET_METHOD_KRONECKER;
        }
        
        if (dixon_global_method != -1) {
            coeff_method = dixon_global_method;
            printf("Method overridden by global variable: %d\n", dixon_global_method);
        }
        printf("Using coefficient matrix determinant method %d:", coeff_method);
        
        compute_fq_coefficient_matrix_det(result, coeff_matrix, matrix_size,
                                         npars, polys[0].ctx, coeff_method, res_deg_bound);
        
        printf("Resultant polynomial has %ld terms\n", result->nterms);
        if (result->nterms <= 100) {
            // Use enhanced print with original parameter names
            fq_mvpoly_print_with_names(result, "Final Resultant", NULL, par_names, gen_name, 0);
        } else {
            printf("Final resultant too large to display (%ld terms)\n", result->nterms);
            
            // Display degree information
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
    
    // Cleanup matrices
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
    
    printf("\n=== Dixon Resultant Computation Complete ===\n");
}

