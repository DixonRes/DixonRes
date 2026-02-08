#ifndef DIXON_COMPLEXITY_H
#define DIXON_COMPLEXITY_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include "dixon_interface_flint.h"

// Define omega parameter for complexity calculation
#define DIXON_OMEGA 2.3

// Helper structure to collect polynomial information
typedef struct {
    char **all_vars;           // All variables found in polynomials
    slong num_all_vars;        // Total number of unique variables
    slong max_vars;            // Allocated space for variables
    long *degrees;             // Degree of each polynomial
    slong num_polys;           // Number of polynomials
    const fq_nmod_ctx_struct *ctx;
} poly_analysis_t;

// Variable hash table entry for fast variable lookups
typedef struct var_entry {
    char *name;
    slong index;
    struct var_entry *next;
} var_entry_t;

// Hash table for variable management
typedef struct {
    var_entry_t **buckets;
    slong bucket_count;
    slong count;
} var_hash_table_t;

// Lightweight parser state for degree calculation only
typedef struct {
    const char *input;
    size_t pos;
    size_t len;
    var_hash_table_t var_table;
    long max_degree_found;
    const fq_nmod_ctx_struct *ctx;
    const char *generator_name;
} lightweight_parser_t;

// Function declarations

// Basic utility functions
int compare_desc(const void *a, const void *b);
int all_equal(const long *arr, int len);

// Core Dixon complexity calculations
void solve_inequality_system(fmpz_t result, const long *t_values, int n);
void fuss_catalan(fmpz_t result, long n, long d);
void dixon_size(fmpz_t result, const long *a_values, int len, int show_details);
double dixon_complexity(const long *a_values, int len, int n, double omega);

// Polynomial analysis functions
static void poly_analysis_init(poly_analysis_t *analysis, slong num_polys, const fq_nmod_ctx_t ctx);
static void poly_analysis_clear(poly_analysis_t *analysis);
static slong find_variable(poly_analysis_t *analysis, const char *var_name);
static void add_variable(poly_analysis_t *analysis, const char *var_name);

// Variable hash table functions
static slong hash_string(const char *str, slong bucket_count);
static void var_hash_init(var_hash_table_t *table, slong initial_buckets);
static void var_hash_clear(var_hash_table_t *table);
static slong var_hash_find(var_hash_table_t *table, const char *name);
static slong var_hash_add(var_hash_table_t *table, const char *name);

// Optimized polynomial analysis functions
static slong find_variable_optimized(poly_analysis_t *analysis, const char *var_name);
static int add_variable_optimized(poly_analysis_t *analysis, const char *var_name);
static int parse_and_extract_degree(lightweight_parser_t *parser);
static void analyze_single_polynomial(poly_analysis_t *analysis, slong poly_idx, const char *poly_str);
static void analyze_single_polynomial_old(poly_analysis_t *analysis, slong poly_idx, const char *poly_str);
static int is_elimination_var(const char *var_name, const char **elim_vars, slong num_elim_vars);

// Main Dixon complexity analysis functions
char* dixon_complexity_auto(const char **poly_strings, slong num_polys,
                           const char **elim_vars, slong num_elim_vars,
                           const fq_nmod_ctx_t ctx);
char* dixon_complexity_auto_str(const char *poly_string, const char *vars_string, const fq_nmod_ctx_t ctx);

// Complexity extraction functions
static long extract_constant_term(const char *poly_str);
double extract_max_complexity(const char **poly_strings, slong num_polys);
double extract_max_complexity_str(const char *poly_string);

// Test function
int test_dixon_complexity(void);

#endif // DIXON_COMPLEXITY_H