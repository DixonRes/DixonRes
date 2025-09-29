// gcc -O3 -march=native -o dixon main.c -lflint -lmpfr -lgmp -lpthread -L/home/suohaohai02/mylinks -lflint -lstdc++ -lpml2 -fopenmp
//x86_64-w64-mingw32-gcc -O3 -march=native -fopenmp -static -o dixon.exe dixon.c -I../include -L../lib -lflint -lmpfr -lgmp -lm -lpthread -lstdc++
//gcc -O3 -march=native -fopenmp -DHAVE_FLINT -DHAVE_PML -fcommon -I/home/ideal/Attack/Dixon/pml-main/flint-extras/include -I/home/ideal/Attack/Dixon/flint-main/include -o dixon *.c -L/home/ideal/Attack/Dixon/flint-main/lib -L/home/ideal/Attack/Dixon/flint-main -L/home/ideal/Attack/Dixon/pml-main/flint-extras/lib -L/home/ideal/Attack/Dixon/pml-main/flint-extras -lflint -lpml -lmpfr -lgmp -lm -lpthread -lstdc++ -fopenmp
//x86_64-w64-mingw32-gcc -O3 -march=native -fopenmp -static -o dixon.exe *.c -mwindows -I../include -L../lib -lflint -lmpfr -lgmp -lm -lpthread -lstdc++ -lcomctl32 -lcomdlg32 -luser32 -lgdi32

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/fmpz.h>
#include <flint/fmpz_factor.h>
#include <flint/fq_nmod.h>
#include <ctype.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <time.h>

#include "dixon_flint.h"
#include "dixon_interface_flint.h"
#include "fq_mvpoly.h"
#include "fq_unified_interface.h"
#include "unified_mpoly_resultant.h"
#include "dixon_with_ideal_reduction.h"
#include "dixon_complexity.h"
#include "polynomial_system_solver.h"

void test_dixon();

// Print usage instructions
static void print_usage(const char *prog_name) {
    printf("Usage:\n");
    printf("  Command line - Basic Dixon:\n");
    printf("    %s \"polynomials\" \"eliminate_vars\" field_size\n", prog_name);
    printf("\n");
    printf("  Command line - Dixon with ideal reduction:\n");
    printf("    %s \"polynomials\" \"eliminate_vars\" \"ideal_generators\" \"all_variables\" field_size\n", 
           prog_name);
    printf("\n");
    printf("  Command line - Polynomial System Solver:\n");
    printf("    %s --solve \"polynomials\" field_size\n", prog_name);
    printf("    (Solves n equations in n variables - equation count MUST equal variable count)\n");
    printf("\n");
    printf("  File input:\n");
    printf("    %s input_file\n", prog_name);
    printf("    (Result will be saved to input_file with '_solution' suffix)\n");
    printf("\n");
    printf("  File input - Polynomial System Solver:\n");
    printf("    %s --solve input_file\n", prog_name);
    printf("    (File format: field_size on first line, polynomials on remaining lines)\n");
    printf("\n");
    printf("  Silent mode (no console output, only result file):\n");
    printf("    %s --silent \"polynomials\" \"eliminate_vars\" field_size\n", prog_name);
    printf("    %s --silent --solve \"polynomials\" field_size\n", prog_name);
    printf("    %s --silent input_file\n", prog_name);
    printf("\n");
    printf("  Test:\n");
    printf("    %s --test\n", prog_name);
    printf("    %s --test-solver\n", prog_name);
    printf("\n");
    printf("File format (basic Dixon - multiline):\n");
    printf("  Line 1: field size (prime < 2^63, prime power, or p^k format)\n");
    printf("  Lines 2 to n-1: polynomials (can span multiple lines, comma separated)\n");
    printf("  Last line: variables TO ELIMINATE (comma separated) - NOT all variables!\n");
    printf("\n");
    printf("File format (polynomial solver - multiline):\n");
    printf("  Line 1: field size (prime < 2^63, prime power, or p^k format)\n");
    printf("  Lines 2 to n: polynomials (one per line or comma separated)\n");
    printf("  (All variables will be automatically detected and solved)\n");
    printf("\n");
    printf("IMPORTANT NOTES:\n");
    printf("  Dixon mode:\n");
    printf("    - The variable list contains ONLY the variables to ELIMINATE, not all variables\n");
    printf("    - Number of variables to eliminate MUST be (number of equations - 1)\n");
    printf("  Solver mode:\n");
    printf("    - Number of equations MUST equal number of variables (n×n system)\n");
    printf("    - All variables are automatically detected from the polynomials\n");
    printf("    - System is solved completely for all variables\n");
    printf("  - For field extensions (e.g., F_256 = F_2^8), the generator is 't' by default\n");
    printf("  - Field size can be input as: 257, 256, 2^8, 3^5, etc.\n");
    printf("  - Prime fields must be less than 2^63 (max: 9223372036854775783)\n");
    printf("\n");
    printf("Examples:\n");
    printf("  Basic Dixon (3 equations, eliminate 2 variables):\n");
    printf("    %s \"x + y + z, x*y + y*z + z*x, x*y*z + 1\" \"x, y\" 257\n", prog_name);
    printf("    (eliminates x and y, keeps z in the resultant)\n");
    printf("\n");
    printf("  Polynomial System Solver (2*2 system):\n");
    printf("    %s --solve \"x + y - 3, 2*x - y\" 257\n", prog_name);
    printf("    (solves for both x and y)\n");
    printf("\n");
    printf("  Polynomial System Solver (3*3 system):\n");
    printf("    %s --solve \"x^2 + y^2 + z^2 - 5, x + y + z - 3, x*y*z - 1\" 257\n", prog_name);
    printf("    (solves for x, y, and z)\n");
    printf("\n");
    printf("  Extension field (2 equations, eliminate 1 variable):\n");
    printf("    %s \"x + y^2 + t, x*y + t*y + 1\" \"x\" 2^8\n", prog_name);
    printf("    (F_256 = F_2^8, 't' is the field extension generator)\n");
    printf("\n");
    printf("  Dixon with ideal reduction:\n");
    printf("    %s \"a1^2 + a2^2 + a3^2 + a4^2 - 10, a4^3 - a1 - a2*a3 - 5\" \"a4\" \"a2^3 = 2*a1 + 1, a3^3 = a1*a2 + 3, a4^3 = a1 + a2*a3 + 5\" 257\n", prog_name);
    printf("    (eliminates a4 and reduces the result using the triangular ideal)\n");
    printf("\n");
    printf("  Silent polynomial solver:\n");
    printf("    %s --silent --solve \"x^2 + y^2 - 1, x + y - 1\" 7\n", prog_name);
    printf("    -> Creates solution.dat with result, shows only timing\n");
    printf("\n");
    printf("  File input:\n");
    printf("    %s example.dat\n", prog_name);
    printf("    -> Output saved to: example_solution.dat\n");
}

// Trim whitespace from both ends of a string
static char* trim(char *str) {
    char *end;
    
    // Trim leading space
    while(isspace((unsigned char)*str)) str++;
    
    if(*str == 0)  // All spaces
        return str;
    
    // Trim trailing space
    end = str + strlen(str) - 1;
    while(end > str && isspace((unsigned char)*end)) end--;
    
    // Write new null terminator
    end[1] = '\0';
    
    return str;
}

// Check if a number is a prime power, return prime and power
static int check_prime_power(const fmpz_t n, fmpz_t prime, ulong *power) {
    if (fmpz_cmp_ui(n, 1) <= 0) return 0;
    
    // Check if it's directly a prime
    if (fmpz_is_probabprime(n)) {
        fmpz_set(prime, n);
        *power = 1;
        return 1;
    }
    
    // Use FLINT's factorization functionality
    fmpz_factor_t factors;
    fmpz_factor_init(factors);
    fmpz_factor(factors, n);
    
    // Check if there's only one prime factor (i.e., it's a prime power)
    if (factors->num == 1) {
        fmpz_set(prime, factors->p + 0);  // First (and only) prime factor
        *power = factors->exp[0];         // Corresponding power
        fmpz_factor_clear(factors);
        return 1;
    }
    
    fmpz_factor_clear(factors);
    return 0;
}

// Parse field size string, supporting formats like "2^8", "256", etc.
static int parse_field_size(const char *field_str, fmpz_t prime, ulong *power) {
    if (!field_str || strlen(field_str) == 0) {
        return 0;
    }
    
    // Check if it contains '^' (power notation)
    const char *caret = strchr(field_str, '^');
    if (caret) {
        // Format: p^k
        char *prime_str = malloc(caret - field_str + 1);
        strncpy(prime_str, field_str, caret - field_str);
        prime_str[caret - field_str] = '\0';
        
        // Parse prime part
        fmpz_t p;
        fmpz_init(p);
        int success = fmpz_set_str(p, prime_str, 10);
        if (success != 0) {
            fmpz_clear(p);
            free(prime_str);
            return 0;
        }
        
        // Parse power part
        char *endptr;
        ulong k = strtoul(caret + 1, &endptr, 10);
        if (*endptr != '\0' || k == 0) {
            fmpz_clear(p);
            free(prime_str);
            return 0;
        }
        
        // Check if p is prime
        if (!fmpz_is_probabprime(p)) {
            fmpz_clear(p);
            free(prime_str);
            return 0;
        }
        
        fmpz_set(prime, p);
        *power = k;
        fmpz_clear(p);
        free(prime_str);
        return 1;
    } else {
        // Format: direct number (could be p or p^k)
        fmpz_t field_size;
        fmpz_init(field_size);
        
        int success = fmpz_set_str(field_size, field_str, 10);
        if (success != 0) {
            fmpz_clear(field_size);
            return 0;
        }
        
        // Check if it's a prime power
        int result = check_prime_power(field_size, prime, power);
        fmpz_clear(field_size);
        return result;
    }
}

// Generate output filename from input filename
// example.dat -> example_solution.dat
// p1.txt -> p1_solution.txt
static char* generate_output_filename(const char *input_filename) {
    if (!input_filename) return NULL;
    
    // Find the last dot for file extension
    const char *dot = strrchr(input_filename, '.');
    
    if (dot) {
        // Has extension
        size_t base_len = dot - input_filename;
        size_t ext_len = strlen(dot);
        char *output = malloc(base_len + 9 + ext_len + 1);  // 9 for "_solution"
        
        strncpy(output, input_filename, base_len);
        output[base_len] = '\0';
        strcat(output, "_solution");
        strcat(output, dot);
        
        return output;
    } else {
        // No extension
        size_t len = strlen(input_filename);
        char *output = malloc(len + 9 + 1);
        strcpy(output, input_filename);
        strcat(output, "_solution");
        return output;
    }
}

static char* read_entire_line(FILE *fp) {
    if (!fp) return NULL;
    
    size_t capacity = 4096;
    size_t length = 0;
    char *line = malloc(capacity);
    if (!line) return NULL;
    
    int c;
    while ((c = fgetc(fp)) != EOF && c != '\n' && c != '\r') {
        if (length + 1 >= capacity) {
            capacity *= 2;
            char *new_line = realloc(line, capacity);
            if (!new_line) {
                free(line);
                return NULL;
            }
            line = new_line;
        }
        line[length++] = (char)c;
    }
    
    // Handle Windows-style \r\n
    if (c == '\r') {
        int next = fgetc(fp);
        if (next != '\n' && next != EOF) {
            ungetc(next, fp);
        }
    }
    
    if (length == 0 && c == EOF) {
        free(line);
        return NULL;
    }
    
    line[length] = '\0';
    return line;
}

static int read_multiline_file(FILE *fp, char **field_str, char **polys_str, 
                              char **vars_str, char **ideal_str, char **allvars_str) {
    char **lines = NULL;
    int line_count = 0;
    int line_capacity = 10;
    
    // Allocate initial array for lines
    lines = malloc(line_capacity * sizeof(char*));
    if (!lines) return 0;
    
    char *line;
    while ((line = read_entire_line(fp)) != NULL) {

        char *trimmed = trim(line);
        
        // Skip empty lines
        if (strlen(trimmed) == 0) {
            free(line);
            continue;
        }
        
        // Skip comment lines
        if (trimmed[0] == '#') {
            free(line);
            continue;
        }
        
        // Expand array if needed
        if (line_count >= line_capacity) {
            line_capacity *= 2;
            char **new_lines = realloc(lines, line_capacity * sizeof(char*));
            if (!new_lines) {
                for (int i = 0; i < line_count; i++) {
                    free(lines[i]);
                }
                free(lines);
                free(line);
                return 0;
            }
            lines = new_lines;
        }
        
        // Store the line
        lines[line_count++] = strdup(trimmed);
        free(line); 
    }
    // Check minimum line count
    if (line_count < 3) {
        fprintf(stderr, "Error: File must contain at least 3 non-empty lines\n");
        fprintf(stderr, "  Line 1: field size\n");
        fprintf(stderr, "  Lines 2 to n-1: polynomials\n");
        fprintf(stderr, "  Line n: variables to ELIMINATE (must be #equations - 1)\n");
        for (int i = 0; i < line_count; i++) {
            free(lines[i]);
        }
        free(lines);
        return 0;
    }
    
    // First line is field size
    *field_str = lines[0];
    
    // Last line is variables to eliminate
    *vars_str = lines[line_count - 1];
    
    // Lines 2 to n-1 are polynomials (join them)
    size_t total_len = 0;
    for (int i = 1; i < line_count - 1; i++) {
        total_len += strlen(lines[i]) + 2;  // +2 for space or comma
    }
    
    // Use dynamic allocation to avoid truncation
    char *poly_buffer = malloc(total_len + 1);
    if (!poly_buffer) {
        for (int i = 0; i < line_count; i++) {
            free(lines[i]);
        }
        free(lines);
        return 0;
    }
    poly_buffer[0] = '\0';
    
    for (int i = 1; i < line_count - 1; i++) {
        if (i > 1) {
            // Add space between lines
            strcat(poly_buffer, " ");
        }
        strcat(poly_buffer, lines[i]);
    }
    
    *polys_str = poly_buffer;
    
    // Free the lines array (but not the strings we're returning)
    for (int i = 1; i < line_count - 1; i++) {
        free(lines[i]);
    }
    free(lines);
    
    // For now, we don't support ideal reduction in multiline format
    *ideal_str = NULL;
    *allvars_str = NULL;
    
    return 1;
}

// Read file for polynomial solver (simpler format: field size, then polynomials)
static int read_solver_file(FILE *fp, char **field_str, char **polys_str) {
    char **lines = NULL;
    int line_count = 0;
    int line_capacity = 10;
    
    // Allocate initial array for lines
    lines = malloc(line_capacity * sizeof(char*));
    if (!lines) return 0;
    
    // Read all lines
    char *line;
    while ((line = read_entire_line(fp)) != NULL) {
        char *trimmed = trim(line);
        
        // Skip empty lines and comments
        if (strlen(trimmed) == 0 || trimmed[0] == '#') {
            free(line);
            continue;
        }
        
        // Expand array if needed
        if (line_count >= line_capacity) {
            line_capacity *= 2;
            char **new_lines = realloc(lines, line_capacity * sizeof(char*));
            if (!new_lines) {
                for (int i = 0; i < line_count; i++) {
                    free(lines[i]);
                }
                free(lines);
                free(line);
                return 0;
            }
            lines = new_lines;
        }
        
        lines[line_count++] = strdup(trimmed);
        free(line);
    }
    
    // Check minimum line count (field size + at least one polynomial)
    if (line_count < 2) {
        fprintf(stderr, "Error: Solver file must contain at least 2 non-empty lines\n");
        fprintf(stderr, "  Line 1: field size\n");
        fprintf(stderr, "  Lines 2+: polynomials\n");
        for (int i = 0; i < line_count; i++) {
            free(lines[i]);
        }
        free(lines);
        return 0;
    }
    
    // First line is field size
    *field_str = lines[0];
    
    // Rest are polynomials - join them intelligently
    size_t total_len = 0;
    for (int i = 1; i < line_count; i++) {
        total_len += strlen(lines[i]) + 3;  // +3 for possible ", " separator
    }
    
    char *poly_buffer = malloc(total_len + 1);
    if (!poly_buffer) {
        for (int i = 0; i < line_count; i++) {
            free(lines[i]);
        }
        free(lines);
        return 0;
    }
    poly_buffer[0] = '\0';
    
    for (int i = 1; i < line_count; i++) {
        if (i > 1) {
            // Check if previous line already ends with comma, or current line starts with comma
            size_t prev_len = strlen(poly_buffer);
            int prev_ends_with_comma = (prev_len > 0 && 
                (poly_buffer[prev_len-1] == ',' || 
                 (prev_len > 1 && poly_buffer[prev_len-2] == ',' && poly_buffer[prev_len-1] == ' ')));
            int curr_starts_with_comma = (lines[i][0] == ',');
            
            if (!prev_ends_with_comma && !curr_starts_with_comma) {
                // Need to add comma separator
                strcat(poly_buffer, ", ");
            } else if (!prev_ends_with_comma && curr_starts_with_comma) {
                // Current starts with comma, just add space
                strcat(poly_buffer, " ");
            }
            // If previous ends with comma, don't add anything extra
        }
        strcat(poly_buffer, lines[i]);
    }
    
    *polys_str = poly_buffer;
    
    // Free the lines array (except the first line which we're returning)
    for (int i = 1; i < line_count; i++) {
        free(lines[i]);
    }
    free(lines);
    
    return 1;
}

// Save polynomial solver result to file
static void save_solver_result_to_file(const char *filename, const char *polys_str, 
                                       mp_limb_t prime, ulong power, 
                                       const polynomial_solutions_t *sols,
                                       double computation_time) {
    FILE *out_fp = fopen(filename, "w");
    if (!out_fp) {
        fprintf(stderr, "Warning: Could not create output file '%s'\n", filename);
        return;
    }
    
    // Calculate actual field size for display
    mp_limb_t field_size = 1;
    for (ulong i = 0; i < power; i++) {
        field_size *= prime;
    }
    
    // Write header information
    fprintf(out_fp, "Polynomial System Solver\n");
    fprintf(out_fp, "========================\n");
    fprintf(out_fp, "Field: F_%lu", prime);
    if (power > 1) {
        fprintf(out_fp, "^%lu (size %lu)", power, field_size);
        fprintf(out_fp, "\nField extension generator: t");
    }
    fprintf(out_fp, "\n");
    
    fprintf(out_fp, "Polynomials: %s\n", polys_str);
    fprintf(out_fp, "Computation time: %.3f seconds\n", computation_time);
    fprintf(out_fp, "\nSolutions:\n");
    fprintf(out_fp, "==========\n");
    
    // Write solutions directly to file
    if (!sols) {
        fprintf(out_fp, "Solution structure is null\n");
        fclose(out_fp);
        return;
    }
    
    fprintf(out_fp, "\n=== Polynomial System Solutions ===\n");
    
    if (!sols->is_valid) {
        fprintf(out_fp, "Solving failed");
        if (sols->error_message) {
            fprintf(out_fp, ": %s", sols->error_message);
        }
        fprintf(out_fp, "\n");
        fclose(out_fp);
        return;
    }
    
    if (sols->has_no_solutions) {
        fprintf(out_fp, "System has no solutions over the finite field\n");
        fclose(out_fp);
        return;
    }
    
    if (sols->num_variables == 0) {
        fprintf(out_fp, "No variables\n");
        fclose(out_fp);
        return;
    }
    
    if (sols->num_solution_sets == 0) {
        fprintf(out_fp, "No solutions found\n");
        fclose(out_fp);
        return;
    }
    
    fprintf(out_fp, "Found %ld complete solution set(s):\n", sols->num_solution_sets);
    
    for (slong set = 0; set < sols->num_solution_sets; set++) {
        fprintf(out_fp, "\nSolution set %ld:\n", set + 1);
        
        for (slong var = 0; var < sols->num_variables; var++) {
            fprintf(out_fp, "  %s = ", sols->variable_names[var]);
            
            slong num_sols = sols->solutions_per_var[set * sols->num_variables + var];
            if (num_sols == 0) {
                fprintf(out_fp, "no solution");
            } else if (num_sols == 1) {
                // Convert solution to string for file output
                char *sol_str = fq_nmod_get_str_pretty(sols->solution_sets[set][var][0], sols->ctx);
                fprintf(out_fp, "%s", sol_str);
                free(sol_str);
            } else {
                fprintf(out_fp, "{");
                for (slong sol = 0; sol < num_sols; sol++) {
                    if (sol > 0) fprintf(out_fp, ", ");
                    char *sol_str = fq_nmod_get_str_pretty(sols->solution_sets[set][var][sol], sols->ctx);
                    fprintf(out_fp, "%s", sol_str);
                    free(sol_str);
                }
                fprintf(out_fp, "}");
            }
            fprintf(out_fp, "\n");
        }
    }
    
    // Also show the compatibility view
    fprintf(out_fp, "\n=== Compatibility View ===\n");
    for (slong var = 0; var < sols->num_variables; var++) {
        fprintf(out_fp, "%s = {", sols->variable_names[var]);
        
        slong total_printed = 0;
        for (slong set = 0; set < sols->num_solution_sets; set++) {
            slong num_sols = sols->solutions_per_var[set * sols->num_variables + var];
            for (slong sol = 0; sol < num_sols; sol++) {
                if (total_printed > 0) fprintf(out_fp, ", ");
                char *sol_str = fq_nmod_get_str_pretty(sols->solution_sets[set][var][sol], sols->ctx);
                fprintf(out_fp, "%s", sol_str);
                free(sol_str);
                total_printed++;
            }
        }
        fprintf(out_fp, "}");
        if (total_printed > 1) {
            fprintf(out_fp, " (%ld solutions)", total_printed);
        } else if (total_printed == 0) {
            fprintf(out_fp, " (no solutions)");
        }
        fprintf(out_fp, "\n");
    }
    
    fprintf(out_fp, "=== Solution Complete ===\n\n");
    
    fclose(out_fp);
}

// Save result to file
static void save_result_to_file(const char *filename, const char *polys_str, 
                               const char *vars_str, const char *ideal_str, 
                               const char *allvars_str, mp_limb_t prime, 
                               ulong power, const char *result, 
                               double computation_time) {
    FILE *out_fp = fopen(filename, "w");
    if (!out_fp) {
        fprintf(stderr, "Warning: Could not create output file '%s'\n", filename);
        return;
    }
    
    // Calculate actual field size for display
    mp_limb_t field_size = 1;
    for (ulong i = 0; i < power; i++) {
        field_size *= prime;
    }
    
    // Write header information
    fprintf(out_fp, "Dixon Resultant Computation\n");
    fprintf(out_fp, "==========================\n");
    fprintf(out_fp, "Field: F_%lu", prime);
    if (power > 1) {
        fprintf(out_fp, "^%lu (size %lu)", power, field_size);
        fprintf(out_fp, "\nField extension generator: t");
    }
    fprintf(out_fp, "\n");
    
    if (ideal_str && allvars_str) {
        fprintf(out_fp, "Mode: Dixon with ideal reduction\n");
        fprintf(out_fp, "Ideal generators: %s\n", ideal_str);
        fprintf(out_fp, "All variables: %s\n", allvars_str);
    } else {
        fprintf(out_fp, "Mode: Basic Dixon resultant\n");
    }
    
    fprintf(out_fp, "Variables eliminated: %s\n", vars_str);
    fprintf(out_fp, "Polynomials: %s\n", polys_str);
    fprintf(out_fp, "Computation time: %.3f seconds\n", computation_time);
    fprintf(out_fp, "\nResultant:\n");
    fprintf(out_fp, "%s\n", result);
    
    fclose(out_fp);
}

// Count number of comma-separated items in a string
static int count_comma_separated_items(const char *str) {
    if (!str || strlen(str) == 0) return 0;
    
    int count = 1;
    for (const char *p = str; *p; p++) {
        if (*p == ',') count++;
    }
    return count;
}

int main(int argc, char *argv[]) {
    clock_t start_time = clock();
    
    if (argc == 1) {
        print_usage(argv[0]);
        return 0;
    }
    
    // Check for flags
    int silent_mode = 0;
    int solve_mode = 0;
    int arg_offset = 0;
    
    // Parse flags
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--silent") == 0) {
            silent_mode = 1;
            arg_offset++;
        } else if (strcmp(argv[i], "--solve") == 0) {
            solve_mode = 1;
            arg_offset++;
        } else {
            break; // Stop at first non-flag argument
        }
    }
    
    // Check if help is needed
    if ((argc == 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)) ||
        (argc - arg_offset == 1 && (strcmp(argv[1+arg_offset], "--help") == 0 || 
                                    strcmp(argv[1+arg_offset], "-h") == 0))) {
        if (!silent_mode) print_usage(argv[0]);
        return 0;
    }

    // Check for test modes
    if ((argc == 2 && strcmp(argv[1], "--test") == 0) ||
        (argc - arg_offset == 1 && strcmp(argv[1+arg_offset], "--test") == 0)) {
        if (!silent_mode) test_dixon();
        return 0;
    }
    
    if ((argc == 2 && strcmp(argv[1], "--test-solver") == 0) ||
        (argc - arg_offset == 1 && strcmp(argv[1+arg_offset], "--test-solver") == 0)) {
        if (!silent_mode) test_polynomial_solver();
        return 0;
    }
    
    char *polys_str = NULL;
    char *vars_str = NULL;
    char *ideal_str = NULL;
    char *allvars_str = NULL;
    char *field_str = NULL;
    int need_free = 0;
    char *input_filename = NULL;
    char *output_filename = NULL;
    mp_limb_t prime = 0;
    ulong power = 0;
    int is_file_input = 0;
    
    //fq_interpolation_use_half_threads();
    
    // Adjust argc and argv for flags
    int effective_argc = argc - arg_offset;
    char **effective_argv = argv + arg_offset;
    
    // Determine input mode
    if (solve_mode) {
        // Polynomial system solver mode
        if (effective_argc == 2) {
            // File input for solver
            FILE *fp = fopen(effective_argv[1], "r");
            if (!fp) {
                if (!silent_mode) {
                    fprintf(stderr, "Error: Cannot open file '%s'\n", effective_argv[1]);
                }
                return 1;
            }
            
            if (!silent_mode) {
                printf("Reading polynomial system from file: %s\n", effective_argv[1]);
            }
            
            input_filename = strdup(effective_argv[1]);
            output_filename = generate_output_filename(input_filename);
            is_file_input = 1;
            
            if (!read_solver_file(fp, &field_str, &polys_str)) {
                fclose(fp);
                return 1;
            }
            
            fclose(fp);
            need_free = 1;
            
        } else if (effective_argc == 3) {
            // Command line solver mode
            polys_str = effective_argv[1];
            field_str = effective_argv[2];
            output_filename = strdup("solution.dat");
            
        } else {
            if (!silent_mode) {
                fprintf(stderr, "Error: Solver mode requires either:\n");
                fprintf(stderr, "  %s --solve \"polynomials\" field_size\n", argv[0]);
                fprintf(stderr, "  %s --solve input_file\n", argv[0]);
            }
            return 1;
        }
        
    } else {
        // Original Dixon modes
        if (effective_argc == 2) {
            // Try to read as file
            FILE *fp = fopen(effective_argv[1], "r");
            if (!fp) {
                if (!silent_mode) {
                    fprintf(stderr, "Error: Cannot open file '%s'\n", effective_argv[1]);
                }
                return 1;
            }
            
            if (!silent_mode) {
                printf("Reading from file: %s\n", effective_argv[1]);
            }
            
            // Save input filename for output file generation
            input_filename = strdup(effective_argv[1]);
            output_filename = generate_output_filename(input_filename);
            is_file_input = 1;
            
            // Use new multiline reading function
            if (!read_multiline_file(fp, &field_str, &polys_str, &vars_str, &ideal_str, &allvars_str)) {
                fclose(fp);
                return 1;
            }
            
            fclose(fp);
            need_free = 1;
            
        } else if (effective_argc == 4 || effective_argc == 5) {
            // Command line arguments
            polys_str = effective_argv[1];
            vars_str = effective_argv[2];
            
            if (effective_argc == 5) {
                ideal_str = effective_argv[3];
                field_str = effective_argv[4];
            } else {
                field_str = effective_argv[3];
            }
            
            // For command line input, always use solution.dat
            output_filename = strdup("solution.dat");
            
        } else {
            if (!silent_mode) {
                print_usage(argv[0]);
            }
            return 1;
        }
    }
    
    // Parse finite field size using new function
    fmpz_t p_fmpz;
    fmpz_init(p_fmpz);
    
    if (!parse_field_size(field_str, p_fmpz, &power)) {
        if (!silent_mode) {
            fprintf(stderr, "Error: Invalid field size '%s'\n", field_str);
            fprintf(stderr, "Field size must be:\n");
            fprintf(stderr, "  - A prime number (e.g., 257)\n");
            fprintf(stderr, "  - A prime power (e.g., 256 = 2^8)\n");
            fprintf(stderr, "  - In p^k format (e.g., 2^8, 3^5)\n");
        }
        if (need_free) {
            free(field_str);
            free(polys_str);
            if (vars_str) free(vars_str);
            if (ideal_str) free(ideal_str);
        }
        if (input_filename) free(input_filename);
        if (output_filename) free(output_filename);
        fmpz_clear(p_fmpz);
        return 1;
    }
    
    prime = fmpz_get_ui(p_fmpz);
    
    // Calculate actual field size for display
    mp_limb_t field_size = 1;
    for (ulong i = 0; i < power; i++) {
        field_size *= prime;
    }
    
    if (!silent_mode) {
        printf("\n=== %s ===\n", solve_mode ? "Polynomial System Solver" : "Dixon Resultant Computation");
        printf("Field: F_%lu", prime);
        if (power > 1) {
            printf("^%lu (size %lu)", power, field_size);
            printf("\nField extension generator: t");
        }
        printf("\n");
    }
    
    // Initialize finite field
    fq_nmod_ctx_t ctx;
    fq_nmod_ctx_init(ctx, p_fmpz, power, "t");
    
    char *result = NULL;
    polynomial_solutions_t *solutions = NULL;
    
    if (solve_mode) {
        // Polynomial system solver mode
        if (!silent_mode) {
            printf("\nMode: Polynomial System Solver\n");
            
            // Count polynomials
            int poly_count = count_comma_separated_items(polys_str);
            printf("Number of polynomials: %d\n", poly_count);
            printf("Polynomials: %s\n", polys_str);
            printf("Note: Number of variables must equal number of equations\n");
            printf("--------------------------------\n");
        }
        
        // Redirect output if silent
        int original_stdout = -1;
        int original_stderr = -1;
        int dev_null = -1;
        
        if (silent_mode) {
            fflush(stdout);
            fflush(stderr);
            
            original_stdout = dup(STDOUT_FILENO);
            original_stderr = dup(STDERR_FILENO);
            
            dev_null = open("/dev/null", O_WRONLY);
            if (dev_null != -1) {
                dup2(dev_null, STDOUT_FILENO);
                dup2(dev_null, STDERR_FILENO);
                close(dev_null);
            }
        }
        
        // Solve the polynomial system
        solutions = solve_polynomial_system_string(polys_str, ctx);
        
        if (silent_mode && original_stdout != -1) {
            fflush(stdout);
            fflush(stderr);
            
            dup2(original_stdout, STDOUT_FILENO);
            dup2(original_stderr, STDERR_FILENO);
            close(original_stdout);
            close(original_stderr);
        }
        
    } else if (ideal_str) {
        // Dixon with ideal reduction
        if (!silent_mode) {
            printf("\nMode: Dixon with ideal reduction\n");
            printf("Polynomials: %s\n", polys_str);
            printf("Eliminate: %s\n", vars_str);
            printf("Ideal generators: %s\n", ideal_str);
            printf("--------------------------------\n");
        }
        
        result = dixon_with_ideal_reduction_str(polys_str, vars_str, 
                                               ideal_str, ctx);
    } else {
        // Basic Dixon computation
        if (!silent_mode) {
            printf("\nMode: Basic Dixon resultant\n");
            
            // Count polynomials and variables
            int poly_count = count_comma_separated_items(polys_str);
            int var_count = count_comma_separated_items(vars_str);
            
            printf("Number of polynomials: %d\n", poly_count);
            printf("Variables to ELIMINATE: %s (count: %d)\n", vars_str, var_count);
            if (var_count != poly_count - 1) {
                printf("WARNING: Dixon method requires eliminating exactly %d variables for %d equations!\n", 
                       poly_count - 1, poly_count);
            }
            printf("--------------------------------\n");
        }
                
        int original_stdout = -1;
        int original_stderr = -1;
        int dev_null = -1;
        
        if (silent_mode) {
            // Flush all buffers
            fflush(stdout);
            fflush(stderr);
            
            // Save original file descriptors
            original_stdout = dup(STDOUT_FILENO);
            original_stderr = dup(STDERR_FILENO);
            
            // Open /dev/null
            dev_null = open("/dev/null", O_WRONLY);
            if (dev_null != -1) {
                dup2(dev_null, STDOUT_FILENO);
                dup2(dev_null, STDERR_FILENO);
                close(dev_null);
            }
        }
        
        result = dixon_str(polys_str, vars_str, ctx);
        
        if (silent_mode && original_stdout != -1) {
            // Restore original output
            fflush(stdout);
            fflush(stderr);
            
            dup2(original_stdout, STDOUT_FILENO);
            dup2(original_stderr, STDERR_FILENO);
            close(original_stdout);
            close(original_stderr);
        }
    }
    
    // Calculate computation time
    clock_t end_time = clock();
    double computation_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    
    // Output result
    if (solve_mode) {
        // Handle polynomial solver results
        if (solutions) {
            if (silent_mode || !silent_mode) {
                printf("\n=== POLYNOMIAL SYSTEM SOLUTIONS ===\n");
                print_polynomial_solutions(solutions);
                printf("====================================\n");
            }
            
            // Save result to output file
            if (output_filename) {
                save_solver_result_to_file(output_filename, polys_str, prime, power, 
                                          solutions, computation_time);
                
                if (!silent_mode) {
                    printf("\nResult saved to: %s\n", output_filename);
                }
            }
            
            polynomial_solutions_clear(solutions);
            free(solutions);
        } else {
            if (!silent_mode) {
                fprintf(stderr, "\nError: Polynomial system solving failed\n");
            }
        }
        
    } else {
        // Handle Dixon resultant results
        if (result) {
            if (0 && !silent_mode) {
                printf("\n=== RESULT ===\n");
                printf("%s\n", result);
                printf("==============\n");
            }
            
            // Save result to output file
            if (output_filename) {
                save_result_to_file(output_filename, polys_str, vars_str, 
                                   ideal_str, allvars_str, prime, power, 
                                   result, computation_time);
                
                if (!silent_mode) {
                    printf("\nResult saved to: %s\n", output_filename);
                }
            }
            
            free(result);
        } else {
            if (!silent_mode) {
                fprintf(stderr, "\nError: Computation failed\n");
            }
        }
    }
    
    // Always print computation time
    printf("Total computation time: %.3f seconds\n", computation_time);
    
    // Cleanup
    fq_nmod_ctx_clear(ctx);
    fmpz_clear(p_fmpz);
    
    if (need_free) {
        free(field_str);
        free(polys_str);
        if (vars_str) free(vars_str);
        if (ideal_str) free(ideal_str);
        if (allvars_str) free(allvars_str);
    }
    
    if (input_filename) {
        free(input_filename);
    }
    
    if (output_filename) {
        free(output_filename);
    }
    
    cleanup_unified_workspace();
    flint_cleanup();
    return 0;
}