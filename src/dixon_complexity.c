#include "dixon_complexity.h"

// Comparison function for descending order sorting
int compare_desc(const void *a, const void *b) {
    long long_a = *(const long*)a;
    long long_b = *(const long*)b;
    return (long_b > long_a) - (long_b < long_a);
}

// Solve inequality system to count non-negative integer solutions
void solve_inequality_system(fmpz_t result, const long *t_values, int n) {
    fmpz_set_ui(result, 0);
    
    if (n == 0) {
        return;
    }
    
    // Check for negative values
    for (int i = 0; i < n; i++) {
        if (t_values[i] < 0) {
            return;
        }
    }
    
    // Find maximum value for array sizing
    long max_t = t_values[0];
    for (int i = 1; i < n; i++) {
        if (t_values[i] > max_t) {
            max_t = t_values[i];
        }
    }
    
    // Allocate dynamic programming arrays
    fmpz *prev_dp = (fmpz*)malloc((max_t + 1) * sizeof(fmpz));
    fmpz *curr_dp = (fmpz*)malloc((max_t + 1) * sizeof(fmpz));
    
    if (!prev_dp || !curr_dp) {
        printf("Memory allocation failed!\n");
        if (prev_dp) free(prev_dp);
        if (curr_dp) free(curr_dp);
        return;
    }
    
    // Initialize arrays
    for (long i = 0; i <= max_t; i++) {
        fmpz_init(prev_dp + i);
        fmpz_init(curr_dp + i);
        fmpz_set_ui(prev_dp + i, 0);
        fmpz_set_ui(curr_dp + i, 0);
    }
    
    // Initialize first layer: S₀ = x₀, range is 0 to t₀
    for (long s = 0; s <= t_values[0]; s++) {
        fmpz_set_ui(prev_dp + s, 1);
    }
    
    // Calculate layer by layer
    for (int i = 1; i < n; i++) {
        long t_i = t_values[i];
        long t_prev = t_values[i-1];
        
        // Clear current layer
        for (long s = 0; s <= t_i; s++) {
            fmpz_set_ui(curr_dp + s, 0);
        }
        
        // Optimization: precompute cumulative sum
        fmpz_t cumsum;
        fmpz_init(cumsum);
        fmpz_set_ui(cumsum, 0);
        
        long limit = (t_prev + 1 < t_i + 1) ? t_prev + 1 : t_i + 1;
        for (long prev_s = 0; prev_s < limit; prev_s++) {
            fmpz_add(cumsum, cumsum, prev_dp + prev_s);
            if (prev_s <= t_i) {
                fmpz_set(curr_dp + prev_s, cumsum);
            }
        }
        
        // Handle remaining states
        for (long s = limit; s <= t_i; s++) {
            long inner_limit = (s + 1 < t_prev + 1) ? s + 1 : t_prev + 1;
            for (long prev_s = 0; prev_s < inner_limit; prev_s++) {
                fmpz_add(curr_dp + s, curr_dp + s, prev_dp + prev_s);
            }
        }
        
        // Swap arrays
        fmpz *temp = prev_dp;
        prev_dp = curr_dp;
        curr_dp = temp;
        
        fmpz_clear(cumsum);
    }
    
    // Calculate final result
    fmpz_set_ui(result, 0);
    for (long i = 0; i <= t_values[n-1]; i++) {
        fmpz_add(result, result, prev_dp + i);
    }
    
    // Cleanup memory
    for (long i = 0; i <= max_t; i++) {
        fmpz_clear(prev_dp + i);
        fmpz_clear(curr_dp + i);
    }
    free(prev_dp);
    free(curr_dp);
}

// Calculate Fuss-Catalan number: (1/(n(d-1)+1)) * C(nd, n)
void fuss_catalan(fmpz_t result, long n, long d) {
    if (n == 0 || d == 0) {
        fmpz_set_ui(result, 1);
        return;
    }
    
    fmpz_t numerator, denominator;
    fmpz_init(numerator);
    fmpz_init(denominator);
    
    // Calculate binomial coefficient C(nd, n)
    fmpz_bin_uiui(numerator, n * d, n);
    
    // Calculate denominator n(d-1)+1
    long denom_val = n * (d - 1) + 1;
    fmpz_set_ui(denominator, denom_val);
    
    // Calculate result = numerator / denominator
    fmpz_fdiv_q(result, numerator, denominator);
    
    fmpz_clear(numerator);
    fmpz_clear(denominator);
}

// Check if all elements in array are equal
int all_equal(const long *arr, int len) {
    if (len <= 1) return 1;
    for (int i = 1; i < len; i++) {
        if (arr[i] != arr[0]) return 0;
    }
    return 1;
}

// Calculate Dixon matrix size
void dixon_size(fmpz_t result, const long *a_values, int len, int show_details) {
    fmpz_set_ui(result, 0);
    
    if (len == 0) {
        return;
    }
    
    // Step 1: Copy all elements, sort in descending order, then remove the minimum value
    long *sorted_values = (long*)malloc(len * sizeof(long));
    if (!sorted_values) {
        printf("Memory allocation failed!\n");
        return;
    }
    
    for (int i = 0; i < len; i++) {
        sorted_values[i] = a_values[i];
    }
    qsort(sorted_values, len, sizeof(long), compare_desc);
    
    // Remove the minimum value (last element), keep the first len-1 elements
    int new_len = len - 1;
    if (new_len <= 0) {
        free(sorted_values);
        return;
    }
    
    long *b_values = (long*)malloc(new_len * sizeof(long));
    if (!b_values) {
        printf("Memory allocation failed!\n");
        free(sorted_values);
        return;
    }
    
    for (int i = 0; i < new_len; i++) {
        b_values[i] = sorted_values[i];
    }
    
    if (show_details) {
        printf("Input sequence: ");
        for (int i = 0; i < len; i++) {
            printf("%ld ", a_values[i]);
        }
        printf("\nSorted sequence (before removing minimum): ");
        for (int i = 0; i < len; i++) {
            printf("%ld ", sorted_values[i]);
        }
        printf("\nSorted sequence (after removing minimum): ");
        for (int i = 0; i < new_len; i++) {
            printf("%ld ", b_values[i]);
        }
        printf("\n");
    }
    
    // If all elements are equal, use Fuss-Catalan numbers
    if (all_equal(b_values, new_len)) {
        fuss_catalan(result, len, b_values[0]);
        if (show_details) {
            printf("All elements are equal, using Fuss-Catalan number: ");
            fmpz_print(result);
            printf("\n");
        }
        free(sorted_values);
        free(b_values);
        return;
    }
    
    // Step 2: Generate cumulative sum sequence
    long *t_values = (long*)malloc(new_len * sizeof(long));
    if (!t_values) {
        printf("Memory allocation failed!\n");
        free(sorted_values);
        free(b_values);
        return;
    }
    
    long cumsum = 0;
    for (int i = 0; i < new_len; i++) {
        cumsum += (b_values[i] - 1);  // Subtract 1 from each element
        t_values[i] = cumsum;
    }
    
    if (show_details) {
        printf("Cumulative sum t_values: ");
        for (int i = 0; i < new_len; i++) {
            printf("%ld ", t_values[i]);
        }
        printf("\n\n");
    }
    
    // Step 3: Solve the inequality system
    solve_inequality_system(result, t_values, new_len);
    
    if (show_details) {
        printf("Inequality system:\n");
        for (int i = 0; i < new_len; i++) {
            printf("  ");
            for (int j = 0; j <= i; j++) {
                if (j > 0) printf(" + ");
                printf("x_%d", j);
            }
            printf(" <= %ld\n", t_values[i]);
        }
        printf("\nNumber of solutions: ");
        fmpz_print(result);
        printf("\n");
    }
    
    free(sorted_values);
    free(b_values);
    free(t_values);
}

// Calculate Dixon complexity
double dixon_complexity(const long *a_values, int len, int n, double omega) {
    fmpz_t size;
    fmpz_init(size);
    
    dixon_size(size, a_values, len, 0);
    
    // Check if size is zero to avoid computation issues
    if (fmpz_is_zero(size)) {
        fmpz_clear(size);
        return 0.0;
    }
    
    long d = 0;
    int m = len;
    
    if (m == n + 1) {
        d = 1;
    } else if (m == n) {
        for (int i = 0; i < len; i++) {
            d += a_values[i];
        }
    } else if (m < n) {
        for (int i = 0; i < len; i++) {
            d += a_values[i];
        }
        // Calculate (d+1)^(n-m+1)
        long exponent = n - m + 1;
        long base = d + 1;
        
        // Use double to avoid integer overflow
        double d_double = pow((double)base, (double)exponent);
        
        if (!isfinite(d_double) || d_double <= 0) {
            fmpz_clear(size);
            return INFINITY;
        }
        
        // Calculate log₂(d * size^omega)
        double size_double = fmpz_get_d(size);
        
        if (size_double <= 0 || !isfinite(size_double)) {
            fmpz_clear(size);
            return 0.0;
        }
        
        double powered_size = pow(size_double, omega);
        
        if (!isfinite(powered_size) || powered_size <= 0) {
            fmpz_clear(size);
            return INFINITY;
        }
        
        double product = d_double * powered_size;
        
        if (!isfinite(product) || product <= 0) {
            fmpz_clear(size);
            return INFINITY;
        }
        
        double result = log2(product);
        fmpz_clear(size);
        return result;
        
    } else {
        fmpz_clear(size);
        return 0.0;
    }
    
    // Handle m >= n cases
    if (d <= 0) {
        fmpz_clear(size);
        return 0.0;
    }
    
    double size_double = fmpz_get_d(size);
    
    if (size_double <= 0 || !isfinite(size_double)) {
        fmpz_clear(size);
        return 0.0;
    }
    
    double powered_size = pow(size_double, omega);
    
    if (!isfinite(powered_size) || powered_size <= 0) {
        fmpz_clear(size);
        return INFINITY;
    }
    
    double product = (double)d * powered_size;
    
    if (!isfinite(product) || product <= 0) {
        fmpz_clear(size);
        return INFINITY;
    }
    
    double result = log2(product);
    
    fmpz_clear(size);
    return result;
}

// Initialize polynomial analysis structure
static void poly_analysis_init(poly_analysis_t *analysis, slong num_polys, const fq_nmod_ctx_t ctx) {
    analysis->max_vars = 16;
    analysis->all_vars = (char**) malloc(analysis->max_vars * sizeof(char*));
    analysis->num_all_vars = 0;
    analysis->degrees = (long*) calloc(num_polys, sizeof(long));
    analysis->num_polys = num_polys;
    analysis->ctx = ctx;
}

// Clear polynomial analysis structure
static void poly_analysis_clear(poly_analysis_t *analysis) {
    for (slong i = 0; i < analysis->num_all_vars; i++) {
        free(analysis->all_vars[i]);
    }
    free(analysis->all_vars);
    free(analysis->degrees);
}

// Check if variable already exists in the list
static slong find_variable(poly_analysis_t *analysis, const char *var_name) {
    for (slong i = 0; i < analysis->num_all_vars; i++) {
        if (strcmp(analysis->all_vars[i], var_name) == 0) {
            return i;
        }
    }
    return -1;
}

// Add variable to the list if not already present
static void add_variable(poly_analysis_t *analysis, const char *var_name) {
    if (find_variable(analysis, var_name) >= 0) {
        return; // Already exists
    }
    
    // Expand array if needed
    if (analysis->num_all_vars >= analysis->max_vars) {
        analysis->max_vars *= 2;
        analysis->all_vars = (char**) realloc(analysis->all_vars, 
                                             analysis->max_vars * sizeof(char*));
    }
    
    analysis->all_vars[analysis->num_all_vars] = strdup(var_name);
    analysis->num_all_vars++;
}

// Simple hash function for variable names
static slong hash_string(const char *str, slong bucket_count) {
    slong hash = 5381;
    while (*str) {
        hash = ((hash << 5) + hash) + *str++;
    }
    return hash < 0 ? (-hash) % bucket_count : hash % bucket_count;
}

// Initialize variable hash table
static void var_hash_init(var_hash_table_t *table, slong initial_buckets) {
    table->bucket_count = initial_buckets;
    table->buckets = (var_entry_t**) calloc(initial_buckets, sizeof(var_entry_t*));
    table->count = 0;
}

// Clear variable hash table
static void var_hash_clear(var_hash_table_t *table) {
    for (slong i = 0; i < table->bucket_count; i++) {
        var_entry_t *entry = table->buckets[i];
        while (entry) {
            var_entry_t *next = entry->next;
            free(entry->name);
            free(entry);
            entry = next;
        }
    }
    free(table->buckets);
    table->buckets = NULL;
    table->bucket_count = 0;
    table->count = 0;
}

// Find variable in hash table (returns index or -1 if not found)
static slong var_hash_find(var_hash_table_t *table, const char *name) {
    if (!table->buckets) return -1;
    
    slong bucket = hash_string(name, table->bucket_count);
    var_entry_t *entry = table->buckets[bucket];
    
    while (entry) {
        if (strcmp(entry->name, name) == 0) {
            return entry->index;
        }
        entry = entry->next;
    }
    return -1;
}

// Add variable to hash table (returns index)
static slong var_hash_add(var_hash_table_t *table, const char *name) {
    // Check if already exists
    slong existing = var_hash_find(table, name);
    if (existing >= 0) return existing;
    
    // Add new entry
    slong bucket = hash_string(name, table->bucket_count);
    var_entry_t *entry = (var_entry_t*) malloc(sizeof(var_entry_t));
    if (!entry) return -1;
    
    entry->name = strdup(name);
    if (!entry->name) {
        free(entry);
        return -1;
    }
    entry->index = table->count++;
    entry->next = table->buckets[bucket];
    table->buckets[bucket] = entry;
    
    return entry->index;
}

// Optimized version of find_variable using hash table
static slong find_variable_optimized(poly_analysis_t *analysis, const char *var_name) {
    // For small numbers of variables, linear search is still fast
    if (analysis->num_all_vars < 16) {
        for (slong i = 0; i < analysis->num_all_vars; i++) {
            if (strcmp(analysis->all_vars[i], var_name) == 0) {
                return i;
            }
        }
        return -1;
    }
    
    // For larger numbers, we'd use a hash table (implementation depends on context)
    // For now, keep linear search but with early termination optimizations
    for (slong i = 0; i < analysis->num_all_vars; i++) {
        const char *existing = analysis->all_vars[i];
        // Quick length check first
        if (strlen(existing) == strlen(var_name) && 
            strcmp(existing, var_name) == 0) {
            return i;
        }
    }
    return -1;
}

// Optimized version of add_variable with better memory management
static int add_variable_optimized(poly_analysis_t *analysis, const char *var_name) {
    // Check if already exists using optimized search
    if (find_variable_optimized(analysis, var_name) >= 0) {
        return 1; // Already exists
    }
    
    // Expand array if needed (grow by 50% instead of doubling for better memory usage)
    if (analysis->num_all_vars >= analysis->max_vars) {
        slong new_max = analysis->max_vars + (analysis->max_vars >> 1) + 8;
        char **new_vars = (char**) realloc(analysis->all_vars, 
                                           new_max * sizeof(char*));
        if (!new_vars) {
            return 0; // Memory allocation failed
        }
        analysis->all_vars = new_vars;
        analysis->max_vars = new_max;
    }
    
    // Add the variable
    analysis->all_vars[analysis->num_all_vars] = strdup(var_name);
    if (!analysis->all_vars[analysis->num_all_vars]) {
        return 0; // strdup failed
    }
    
    analysis->num_all_vars++;
    return 1;
}

// Fast degree-only parsing without full polynomial construction
static int parse_and_extract_degree(lightweight_parser_t *parser) {
    parser->max_degree_found = 0;
    var_hash_init(&parser->var_table, 16);
    
    // Simple tokenizer focused on finding variables and their powers
    while (parser->pos < parser->len) {
        char c = parser->input[parser->pos];
        
        // Skip whitespace and operators
        if (isspace(c) || c == '+' || c == '-' || c == '*' || c == '(' || c == ')') {
            parser->pos++;
            continue;
        }
        
        // Handle numbers
        if (isdigit(c)) {
            while (parser->pos < parser->len && isdigit(parser->input[parser->pos])) {
                parser->pos++;
            }
            continue;
        }
        
        // Handle identifiers (variables)
        if (isalpha(c) || c == '_') {
            size_t start = parser->pos;
            while (parser->pos < parser->len && 
                   (isalnum(parser->input[parser->pos]) || parser->input[parser->pos] == '_')) {
                parser->pos++;
            }
            
            // Extract variable name
            size_t len = parser->pos - start;
            char *var_name = (char*) malloc(len + 1);
            if (!var_name) return 0;
            
            strncpy(var_name, parser->input + start, len);
            var_name[len] = '\0';
            
            // Skip if it's the generator
            if (parser->generator_name && strcmp(var_name, parser->generator_name) == 0) {
                free(var_name);
                continue;
            }
            
            // Add to variable table
            slong var_index = var_hash_add(&parser->var_table, var_name);
            free(var_name);
            
            if (var_index < 0) return 0; // Error
            
            // Check for power operator
            long power = 1;
            if (parser->pos < parser->len && parser->input[parser->pos] == '^') {
                parser->pos++; // skip '^'
                power = 0;
                while (parser->pos < parser->len && isdigit(parser->input[parser->pos])) {
                    power = power * 10 + (parser->input[parser->pos] - '0');
                    parser->pos++;
                }
                if (power == 0) power = 1; // Handle edge case
            }
            
            // Update max degree for this term (simple heuristic)
            if (power > parser->max_degree_found) {
                parser->max_degree_found = power;
            }
        } else {
            parser->pos++; // Skip unknown characters
        }
    }
    
    return 1;
}

// Main optimized function
static void analyze_single_polynomial(poly_analysis_t *analysis, slong poly_idx, 
                                               const char *poly_str) {
    // Input validation
    if (!analysis || !poly_str || poly_idx >= analysis->num_polys) {
        if (analysis && poly_idx < analysis->num_polys) {
            analysis->degrees[poly_idx] = 0;
        }
        return;
    }
    
    // Early exit for empty polynomial
    size_t poly_len = strlen(poly_str);
    if (poly_len == 0) {
        analysis->degrees[poly_idx] = 0;
        return;
    }
    
    // Get generator name (cached if possible)
    char *gen_name = get_generator_name(analysis->ctx);
    if (!gen_name) {
        analysis->degrees[poly_idx] = 0;
        return;
    }
    
    // Use lightweight parser for degree calculation
    lightweight_parser_t light_parser;
    light_parser.input = poly_str;
    light_parser.pos = 0;
    light_parser.len = poly_len;
    light_parser.ctx = analysis->ctx;
    light_parser.generator_name = gen_name;
    
    int parse_success = parse_and_extract_degree(&light_parser);
    
    if (parse_success) {
        // Extract variables from the hash table and add to global analysis
        for (slong i = 0; i < light_parser.var_table.bucket_count; i++) {
            var_entry_t *entry = light_parser.var_table.buckets[i];
            while (entry) {
                if (!add_variable_optimized(analysis, entry->name)) {
                    // Handle allocation failure gracefully
                    break;
                }
                entry = entry->next;
            }
        }
        
        analysis->degrees[poly_idx] = light_parser.max_degree_found;
    } else {
        // Fallback to original method if lightweight parsing fails
        printf("  Lightweight parsing failed, using full parser for polynomial %ld\n", poly_idx);
        
        // Initialize full parser state
        parser_state_t state;
        state.var_names = NULL;
        state.nvars = 0;
        state.npars = 0;
        state.max_pars = 16;
        state.par_names = (char**) malloc(state.max_pars * sizeof(char*));
        if (!state.par_names) {
            analysis->degrees[poly_idx] = 0;
            goto cleanup;
        }
        
        state.ctx = analysis->ctx;
        state.current.str = NULL;
        fq_nmod_init(state.current.value, analysis->ctx);
        state.generator_name = strdup(gen_name);
        
        // Parse polynomial to discover all variables/parameters
        fq_mvpoly_t temp_poly;
        fq_mvpoly_init(&temp_poly, 0, state.max_pars, analysis->ctx);
        
        state.input = poly_str;
        state.pos = 0;
        state.len = poly_len;
        next_token(&state);
        
        parse_expression(&state, &temp_poly);
        
        // Add discovered parameters as variables
        for (slong i = 0; i < state.npars; i++) {
            add_variable_optimized(analysis, state.par_names[i]);
        }
        
        // Calculate total degree efficiently
        long max_degree = 0;
        for (slong i = 0; i < temp_poly.nterms; i++) {
            long term_degree = 0;
            
            // Sum degrees from parameters (variables in our context)
            if (temp_poly.terms[i].par_exp) {
                for (slong j = 0; j < state.npars; j++) {
                    term_degree += temp_poly.terms[i].par_exp[j];
                }
            }
            
            // Update max degree
            if (term_degree > max_degree) {
                max_degree = term_degree;
            }
        }
        
        analysis->degrees[poly_idx] = max_degree;
        
        // Cleanup full parser state
        fq_mvpoly_clear(&temp_poly);
        for (slong i = 0; i < state.npars; i++) {
            free(state.par_names[i]);
        }
        free(state.par_names);
        if (state.generator_name) free(state.generator_name);
        fq_nmod_clear(state.current.value, analysis->ctx);
        if (state.current.str) free(state.current.str);
    }
    
cleanup:
    // Cleanup lightweight parser
    var_hash_clear(&light_parser.var_table);
    free(gen_name);
}

static void analyze_single_polynomial_old(poly_analysis_t *analysis, slong poly_idx, 
                                     const char *poly_str) {
    // Get generator name
    char *gen_name = get_generator_name(analysis->ctx);
    
    // Initialize parser state for variable discovery
    parser_state_t state;
    state.var_names = NULL;
    state.nvars = 0;
    state.npars = 0;
    state.max_pars = 16;
    state.par_names = (char**) malloc(state.max_pars * sizeof(char*));
    state.ctx = analysis->ctx;
    state.current.str = NULL;
    fq_nmod_init(state.current.value, analysis->ctx);
    state.generator_name = strdup(gen_name);
    
    // Parse polynomial to discover all variables/parameters
    fq_mvpoly_t temp_poly;
    fq_mvpoly_init(&temp_poly, 0, state.max_pars, analysis->ctx);
    
    state.input = poly_str;
    state.pos = 0;
    state.len = strlen(poly_str);
    next_token(&state);
    
    parse_expression(&state, &temp_poly);
    
    // Add all discovered parameters as variables to our global list
    for (slong i = 0; i < state.npars; i++) {
        add_variable(analysis, state.par_names[i]);
    }
    
    // Calculate total degree of this polynomial
    long max_degree = 0;
    for (slong i = 0; i < temp_poly.nterms; i++) {
        long term_degree = 0;
        
        // Sum degrees from parameters (variables in our context)
        if (temp_poly.terms[i].par_exp) {
            for (slong j = 0; j < state.npars; j++) {
                term_degree += temp_poly.terms[i].par_exp[j];
            }
        }
        
        if (term_degree > max_degree) {
            max_degree = term_degree;
        }
    }
    
    analysis->degrees[poly_idx] = max_degree;
    
    // Cleanup
    fq_mvpoly_clear(&temp_poly);
    for (slong i = 0; i < state.npars; i++) {
        free(state.par_names[i]);
    }
    free(state.par_names);
    if (state.generator_name) free(state.generator_name);
    fq_nmod_clear(state.current.value, analysis->ctx);
    if (state.current.str) free(state.current.str);
    free(gen_name);
}

// Check if a variable is in the elimination list
static int is_elimination_var(const char *var_name, const char **elim_vars, slong num_elim_vars) {
    for (slong i = 0; i < num_elim_vars; i++) {
        if (strcmp(var_name, elim_vars[i]) == 0) {
            return 1;
        }
    }
    return 0;
}

// Main function: dixon_complexity_auto with complexity encoding
char* dixon_complexity_auto(const char **poly_strings, slong num_polys,
                           const char **elim_vars, slong num_elim_vars,
                           const fq_nmod_ctx_t ctx) {
    
    printf("\n=== Dixon Complexity Analysis ===\n");
    
    // Initialize analysis structure
    poly_analysis_t analysis;
    poly_analysis_init(&analysis, num_polys, ctx);
    
    // Analyze each polynomial
    printf("Analyzing %ld polynomials...\n", num_polys);
    for (slong i = 0; i < num_polys; i++) {
        printf("  Polynomial %ld: ", i + 1);
        analyze_single_polynomial(&analysis, i, poly_strings[i]);
        printf("degree = %ld\n", analysis.degrees[i]);
    }
    
    // Print all discovered variables
    printf("\nDiscovered variables (%ld): ", analysis.num_all_vars);
    for (slong i = 0; i < analysis.num_all_vars; i++) {
        if (i > 0) printf(", ");
        printf("%s", analysis.all_vars[i]);
    }
    printf("\n");
    
    printf("Elimination variables (%ld): ", num_elim_vars);
    for (slong i = 0; i < num_elim_vars; i++) {
        if (i > 0) printf(", ");
        printf("%s", elim_vars[i]);
    }
    printf("\n");
    
    // Identify remaining variables (not being eliminated)
    char **remaining_vars = (char**) malloc(analysis.num_all_vars * sizeof(char*));
    slong num_remaining = 0;
    
    for (slong i = 0; i < analysis.num_all_vars; i++) {
        if (!is_elimination_var(analysis.all_vars[i], elim_vars, num_elim_vars)) {
            remaining_vars[num_remaining] = strdup(analysis.all_vars[i]);
            num_remaining++;
        }
    }
    
    printf("Remaining variables (%ld): ", num_remaining);
    for (slong i = 0; i < num_remaining; i++) {
        if (i > 0) printf(", ");
        printf("%s", remaining_vars[i]);
    }
    printf("\n");
    
    // Calculate total degree product (td)
    long td = 1;
    printf("\nDegree calculation: ");
    for (slong i = 0; i < num_polys; i++) {
        if (i > 0) printf(" * ");
        printf("%ld", analysis.degrees[i]);
        td *= analysis.degrees[i];
    }
    printf(" = %ld\n", td);
    
    // Compute Dixon matrix size and complexity
    printf("\n=== Dixon Matrix Analysis ===\n");
    
    fmpz_t matrix_size;
    fmpz_init(matrix_size);
    dixon_size(matrix_size, analysis.degrees, num_polys, 1);
    
    printf("Dixon matrix size: ");
    fmpz_print(matrix_size);
    printf("\n");
    
    // Calculate complexity
    double complexity = dixon_complexity(analysis.degrees, num_polys, 
                                       analysis.num_all_vars, DIXON_OMEGA);
    printf("Dixon complexity (log_2): %.6f\n", complexity);
    printf("Omega parameter: %.1f\n", DIXON_OMEGA);
    
    // Encode complexity as integer (complexity * 1000)
    long complexity_encoded = (long)(complexity * 1000.0 + 0.5); // Round to nearest integer
    
    // Generate evaluation polynomial with encoded complexity
    printf("\n=== Generating Evaluation Polynomial ===\n");
    
    size_t eval_poly_size = 1024; // Initial buffer size
    char *eval_poly = (char*) malloc(eval_poly_size);
    eval_poly[0] = '\0';
    
    if (num_remaining == 0) {
        // No remaining variables, just output the encoded complexity as constant
        sprintf(eval_poly, "%ld", complexity_encoded);
        printf("No remaining variables, evaluation polynomial: %s\n", eval_poly);
    } else {
        // Build polynomial: x₁^td + x₂^td + ... + xₙ^td + complexity_encoded
        for (slong i = 0; i < num_remaining; i++) {
            if (i > 0) {
                strcat(eval_poly, " + ");
            }
            
            // Check if we need more space
            size_t needed = strlen(eval_poly) + strlen(remaining_vars[i]) + 50;
            if (needed >= eval_poly_size) {
                eval_poly_size = needed * 2;
                eval_poly = (char*) realloc(eval_poly, eval_poly_size);
            }
            
            strcat(eval_poly, remaining_vars[i]);
            if (td > 1) {
                char exp_str[32];
                sprintf(exp_str, "^%ld", td);
                strcat(eval_poly, exp_str);
            }
        }
        
        // Add encoded complexity as constant term
        char complexity_str[64];
        sprintf(complexity_str, " + %ld", complexity_encoded);
        
        // Ensure buffer is large enough
        size_t needed = strlen(eval_poly) + strlen(complexity_str) + 1;
        if (needed >= eval_poly_size) {
            eval_poly_size = needed * 2;
            eval_poly = (char*) realloc(eval_poly, eval_poly_size);
        }
        
        strcat(eval_poly, complexity_str);
        printf("Evaluation polynomial: %s\n", eval_poly);
    }
    
    // Cleanup remaining variables
    for (slong i = 0; i < num_remaining; i++) {
        free(remaining_vars[i]);
    }
    free(remaining_vars);
    
    // Cleanup analysis
    poly_analysis_clear(&analysis);
    fmpz_clear(matrix_size);
    
    printf("=== Analysis Complete ===\n\n");
    
    return eval_poly;
}

// String interface version (matching dixon_str signature)
char* dixon_complexity_auto_str(const char *poly_string,    // comma-separated polynomials
                                const char *vars_string,     // comma-separated elimination variables
                                const fq_nmod_ctx_t ctx) {
    
    // Split input strings
    slong num_polys, num_vars;
    char **poly_array = split_string(poly_string, &num_polys);
    char **vars_array = split_string(vars_string, &num_vars);
    
    // Convert to const char**
    const char **poly_strings = (const char**) malloc(num_polys * sizeof(char*));
    const char **elim_vars = (const char**) malloc(num_vars * sizeof(char*));
    
    for (slong i = 0; i < num_polys; i++) {
        poly_strings[i] = poly_array[i];
    }
    for (slong i = 0; i < num_vars; i++) {
        elim_vars[i] = vars_array[i];
    }
    
    // Call main function
    char *result = dixon_complexity_auto(poly_strings, num_polys, elim_vars, num_vars, ctx);
    
    // Cleanup
    free_split_strings(poly_array, num_polys);
    free_split_strings(vars_array, num_vars);
    free(poly_strings);
    free(elim_vars);
    
    return result;
}

// Extract constant term from a polynomial string
static long extract_constant_term(const char *poly_str) {
    if (!poly_str || strlen(poly_str) == 0) {
        return 0;
    }
    
    // Simple case: if the entire string is just a number
    char *endptr;
    long value = strtol(poly_str, &endptr, 10);
    if (*endptr == '\0') {
        // Entire string is a number
        return value;
    }
    
    // More complex case: find constant terms in the polynomial
    long max_constant = 0;
    const char *ptr = poly_str;
    
    while (*ptr) {
        // Skip whitespace
        while (*ptr && isspace(*ptr)) ptr++;
        if (!*ptr) break;
        
        // Skip '+' sign
        if (*ptr == '+') {
            ptr++;
            while (*ptr && isspace(*ptr)) ptr++;
        }
        
        // Check if this is a number (constant term)
        if (isdigit(*ptr) || (*ptr == '-' && isdigit(*(ptr+1)))) {
            char *term_end;
            long term_value = strtol(ptr, &term_end, 10);
            
            // Check if this is truly a constant (no variables following)
            const char *check_ptr = term_end;
            while (*check_ptr && isspace(*check_ptr)) check_ptr++;
            
            if (!*check_ptr || *check_ptr == '+' || *check_ptr == '-') {
                // This is a constant term
                if (term_value > max_constant) {
                    max_constant = term_value;
                }
            }
            
            ptr = term_end;
        } else {
            // Skip to next term
            while (*ptr && *ptr != '+' && *ptr != '-') ptr++;
            if (*ptr == '-') ptr++; // Don't skip the minus sign as it might be part of the next term
        }
    }
    
    return max_constant;
}

// Extract maximum complexity from multiple polynomial strings
double extract_max_complexity(const char **poly_strings, slong num_polys) {
    if (!poly_strings || num_polys == 0) {
        return 0.0;
    }
    
    long max_encoded = 0;
    
    printf("\n=== Extracting Complexity Values ===\n");
    
    for (slong i = 0; i < num_polys; i++) {
        long constant = extract_constant_term(poly_strings[i]);
        printf("Polynomial %ld: constant term = %ld\n", i + 1, constant);
        
        if (constant > max_encoded) {
            max_encoded = constant;
        }
    }
    
    double complexity = (double)max_encoded / 1000.0;
    printf("Maximum encoded value: %ld\n", max_encoded);
    printf("Decoded complexity: %.6f\n", complexity);
    printf("=== Extraction Complete ===\n\n");
    
    return complexity;
}

// String interface for complexity extraction
double extract_max_complexity_str(const char *poly_string) {
    if (!poly_string || strlen(poly_string) == 0) {
        return 0.0;
    }
    
    // Split input string by commas
    slong num_polys;
    char **poly_array = split_string(poly_string, &num_polys);
    
    if (!poly_array || num_polys == 0) {
        return 0.0;
    }
    
    // Convert to const char**
    const char **poly_strings = (const char**) malloc(num_polys * sizeof(char*));
    for (slong i = 0; i < num_polys; i++) {
        poly_strings[i] = poly_array[i];
    }
    
    // Extract complexity
    double complexity = extract_max_complexity(poly_strings, num_polys);
    
    // Cleanup
    free_split_strings(poly_array, num_polys);
    free(poly_strings);
    
    return complexity;
}

int test_dixon_complexity() {
    // Test data
    long a1[] = {1000, 1000, 1000, 1001, 1002, 1003};
    int len = sizeof(a1) / sizeof(a1[0]);
    double omega = 2.3;
    
    printf("Dixon Complexity Results:\n");
    for (int n = 0; n < 10; n++) {
        double complexity = dixon_complexity(a1, len, n, omega);
        printf("n=%d: %.6f\n", n, complexity);
    }
    
    // Test dixon_size function
    printf("\nTesting dixon_size with details:\n");
    fmpz_t test_result;
    fmpz_init(test_result);
    dixon_size(test_result, a1, len, 1);
    fmpz_clear(test_result);
    
    return 0;
}
