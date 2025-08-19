// gcc -O3 -march=native -o dixon main.c -lflint -lmpfr -lgmp -lpthread -L/home/suohaohai02/mylinks -lflint -lstdc++ -lpml2 -fopenmp
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/fmpz.h>
#include <flint/fq_nmod.h>
#include <ctype.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <time.h>
#include "dixon_test.h"

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
    printf("  File input:\n");
    printf("    %s input_file\n", prog_name);
    printf("    (Result will be saved to input_file with '_solution' suffix)\n");
    printf("\n");
    printf("  Silent mode (no console output, only result file):\n");
    printf("    %s --silent \"polynomials\" \"eliminate_vars\" field_size\n", prog_name);
    printf("    %s --silent input_file\n", prog_name);
    printf("\n");
    printf("  Test:\n");
    printf("    %s --test\n", prog_name);
    printf("\n");
    printf("File format (basic Dixon - multiline):\n");
    printf("  Line 1: field size (prime, prime power, or p^k format)\n");
    printf("  Lines 2 to n-1: polynomials (can span multiple lines, comma separated)\n");
    printf("  Last line: variables TO ELIMINATE (comma separated) - NOT all variables!\n");
    printf("\n");
    printf("IMPORTANT NOTES:\n");
    printf("  - The variable list contains ONLY the variables to ELIMINATE, not all variables\n");
    printf("  - Number of variables to eliminate MUST be (number of equations - 1)\n");
    printf("  - For field extensions (e.g., F_256 = F_2^8), the generator is 't' by default\n");
    printf("  - Field size can be input as: 257, 2^8, 3^5, etc.\n");
    printf("  - Command line inputs will create 'solution.dat' file\n");
    printf("  - --silent mode suppresses all output except timing and creates result file only\n");
    printf("\n");
    printf("Examples:\n");
    printf("  Basic field (3 equations, eliminate 2 variables):\n");
    printf("    %s \"x + y + z, x*y + y*z + z*x, x*y*z + 1\" \"x, y\" 257\n", prog_name);
    printf("    (eliminates x and y, keeps z in the resultant)\n");
    printf("\n");
    printf("  Extension field (2 equations, eliminate 1 variable):\n");
    printf("    %s \"x + y^2 + t, x*y + t*y + 1\" \"x\" 2^8\n", prog_name);
    printf("    (F_256 = F_2^8, 't' is the field extension generator)\n");
    printf("\n");
    printf("  Silent mode:\n");
    printf("    %s --silent \"x + y, x*y + 1\" \"x\" 257\n", prog_name);
    printf("    → Creates solution.dat with result, shows only timing\n");
    printf("\n");
    printf("  File input:\n");
    printf("    %s example.dat\n", prog_name);
    printf("    → Output saved to: example_solution.dat\n");
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

// 检查是否是素数幂，返回素数和幂次
static int check_prime_power(mp_limb_t n, mp_limb_t *prime, ulong *power) {
    if (n <= 1) return 0;
    
    // 直接是素数的情况
    if (n_is_prime(n)) {
        *prime = n;
        *power = 1;
        return 1;
    }
    
    // 检查是否是素数的幂
    mp_limb_t root;
    int is_perfect = n_is_perfect_power(&root, n);
    
    if (is_perfect && n_is_prime(root)) {
        *prime = root;
        // 计算幂次
        *power = 1;
        mp_limb_t temp = root;
        while (temp < n) {
            temp *= root;
            (*power)++;
        }
        return 1;
    }
    
    // 使用 n_factor_t 结构进行因式分解
    n_factor_t factors;
    n_factor_init(&factors);  // 初始化
    n_factor(&factors, n, 1);  // 1 表示需要证明素性
    
    // 检查是否只有一个素因子（即是素数的幂）
    if (factors.num == 1) {
        *prime = factors.p[0];
        *power = factors.exp[0];
        return 1;
    }
    
    // n_factor_t 不需要清理，因为它不分配动态内存
    return 0;
}

// Parse field size string, supporting formats like "2^8", "256", etc.
static int parse_field_size(const char *field_str, mp_limb_t *prime, ulong *power) {
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
        
        char *endptr;
        mp_limb_t p = strtoul(prime_str, &endptr, 10);
        if (*endptr != '\0' || p == 0) {
            free(prime_str);
            return 0;
        }
        
        ulong k = strtoul(caret + 1, &endptr, 10);
        if (*endptr != '\0' || k == 0) {
            free(prime_str);
            return 0;
        }
        
        // Check if p is prime
        if (!n_is_prime(p)) {
            free(prime_str);
            return 0;
        }
        
        *prime = p;
        *power = k;
        free(prime_str);
        return 1;
    } else {
        // Format: direct number (could be p or p^k)
        char *endptr;
        mp_limb_t field_size = strtoul(field_str, &endptr, 10);
        if (*endptr != '\0' || field_size == 0) {
            return 0;
        }
        
        // Check if it's a prime power
        return check_prime_power(field_size, prime, power);
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
    
    // 处理 Windows 风格的 \r\n
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
    
    // 使用新的读行函数，没有长度限制
    char *line;
    while ((line = read_entire_line(fp)) != NULL) {
        // 去除首尾空白
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
                // 清理已分配的内存
                for (int i = 0; i < line_count; i++) {
                    free(lines[i]);
                }
                free(lines);
                free(line);
                return 0;
            }
            lines = new_lines;
        }
        
        // Store the line（使用strdup复制trimmed部分）
        lines[line_count++] = strdup(trimmed);
        free(line);  // 释放原始行
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
    
    // 使用动态分配，避免截断
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

int main(int argc, char *argv[]) {
    clock_t start_time = clock();
    
    // Check for silent mode
    int silent_mode = 0;
    int arg_offset = 0;
    if (argc > 1 && strcmp(argv[1], "--silent") == 0) {
        silent_mode = 1;
        arg_offset = 1;
    }
    
    // Check if help is needed
    if ((argc == 2 && strcmp(argv[1], "--help") == 0) || 
        (argc == 3 && strcmp(argv[1], "--silent") == 0 && strcmp(argv[2], "--help") == 0)) {
        if (!silent_mode) print_usage(argv[0]);
        return 0;
    }

    if ((argc == 2 && strcmp(argv[1], "--test") == 0) ||
        (argc == 3 && strcmp(argv[1], "--silent") == 0 && strcmp(argv[2], "--test") == 0)) {
        if (!silent_mode) test_dixon();
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
    
    fq_interpolation_use_half_threads();
    
    // Adjust argc and argv for silent mode
    int effective_argc = argc - arg_offset;
    char **effective_argv = argv + arg_offset;
    
    // Determine if file input or command line input
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
        
    } else if (effective_argc == 4 || effective_argc == 6) {
        // Command line arguments
        polys_str = effective_argv[1];
        vars_str = effective_argv[2];
        
        if (effective_argc == 6) {
            ideal_str = effective_argv[3];
            allvars_str = effective_argv[4];
            field_str = effective_argv[5];
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
    
    // Parse finite field size using new function
    if (!parse_field_size(field_str, &prime, &power)) {
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
            free(vars_str);
            if (ideal_str) free(ideal_str);
            if (allvars_str) free(allvars_str);
        }
        if (input_filename) free(input_filename);
        if (output_filename) free(output_filename);
        return 1;
    }
    
    // Calculate actual field size for display
    mp_limb_t field_size = 1;
    for (ulong i = 0; i < power; i++) {
        field_size *= prime;
    }
    
    if (!silent_mode) {
        printf("\n=== Dixon Resultant Computation ===\n");
        printf("Field: F_%lu", prime);
        if (power > 1) {
            printf("^%lu (size %lu)", power, field_size);
            printf("\nField extension generator: t");
        }
        printf("\n");
    }
    
    // Initialize finite field
    fmpz_t p;
    fmpz_init_set_ui(p, prime);
    
    fq_nmod_ctx_t ctx;
    fq_nmod_ctx_init(ctx, p, power, "t");
    
    char *result = NULL;
    
    if (ideal_str && allvars_str) {
        // Dixon with ideal reduction
        if (!silent_mode) {
            printf("\nMode: Dixon with ideal reduction\n");
            printf("Polynomials: %s\n", polys_str);
            printf("Eliminate: %s\n", vars_str);
            printf("Ideal generators: %s\n", ideal_str);
            printf("All variables: %s\n", allvars_str);
            printf("--------------------------------\n");
        }
        
        result = dixon_with_ideal_reduction_str(polys_str, vars_str, 
                                               ideal_str, allvars_str, ctx);
    } else {
        // Basic Dixon computation
        if (!silent_mode) {
            printf("\nMode: Basic Dixon resultant\n");
            
            // Count polynomials and variables
            int poly_count = 1;
            for (const char *p = polys_str; *p; p++) {
                if (*p == ',') poly_count++;
            }
            int var_count = 1;
            for (const char *p = vars_str; *p; p++) {
                if (*p == ',') var_count++;
            }
            
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
            // 刷新所有缓冲区
            fflush(stdout);
            fflush(stderr);
            
            // 保存原始文件描述符
            original_stdout = dup(STDOUT_FILENO);
            original_stderr = dup(STDERR_FILENO);
            
            // 打开 /dev/null
            dev_null = open("/dev/null", O_WRONLY);
            if (dev_null != -1) {
                dup2(dev_null, STDOUT_FILENO);
                dup2(dev_null, STDERR_FILENO);
                close(dev_null);
            }
        }
        
        result = dixon_str(polys_str, vars_str, ctx);
        
        if (silent_mode && original_stdout != -1) {
            // 恢复原始输出
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
    
    // Always print computation time
    printf("Total computation time: %.3f seconds\n", computation_time);
    
    // Cleanup
    fq_nmod_ctx_clear(ctx);
    fmpz_clear(p);
    
    if (need_free) {
        free(field_str);
        free(polys_str);
        free(vars_str);
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