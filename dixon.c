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
#include "dixon_test.h"

#define PROGRAM_VERSION "0.0.2"

/* =========================================================================
 * Print usage
 * ========================================================================= */
static void print_usage(const char *prog_name)
{
    printf("===============================================\n");
    printf("DixonRes v%s\n", PROGRAM_VERSION);
    printf("FLINT version: %s (Recommended: 3.4.0)\n", FLINT_VERSION);
#ifdef HAVE_PML
    printf("PML support: ENABLED\n");
#else
    printf("PML support: DISABLED\n");
#endif
    printf("===============================================\n");

    printf("USAGE:\n");
    printf("  Basic Dixon:\n");
    printf("    %s \"polynomials\" \"eliminate_vars\" field_size\n", prog_name);
    printf("    -> Default output file: solution+timestamp.dat\n");

    printf("  Polynomial system solver:\n");
    printf("    %s --solve \"polynomials\" field_size\n", prog_name);
    printf("    -> Writes all solutions to solution+timestamp.dat\n");

    printf("  Complexity analysis:\n");
    printf("    %s --comp \"polynomials\" \"eliminate_vars\" field_size\n", prog_name);
    printf("    %s -c    \"polynomials\" \"eliminate_vars\" field_size\n", prog_name);
    printf("    %s --comp input_file\n", prog_name);
    printf("    -> Prints complexity info; saves to comp+timestamp.dat\n");
    printf("    Add --omega <value> (or -w <value>) to set omega (default: %.4g)\n",
           DIXON_OMEGA);
//     printf("  Dixon with ideal reduction:\n");
//     printf("    %s \"polynomials\" \"eliminate_vars\" \"ideal_generators\" \"all_variables\" field_size\n", prog_name);
//     printf("    -> Writes reduced resultant to solution+timestamp.dat\n");
    printf("  File input:\n");
    printf("    %s input_file\n", prog_name);
    printf("    %s --solve input_file\n", prog_name);
    printf("    -> Output saved as input_file_solution+timestamp.dat\n");

    printf("  Silent mode:\n");
    printf("    %s --silent [--solve|--comp|-c] <args>\n", prog_name);
    printf("    -> No console output; solution file is still generated\n");

    printf("FILE FORMAT (Basic Dixon / Complexity, multiline):\n");
    printf("  Line 1 : field size (prime or p^k; generator defaults to 't')\n");
    printf("  Line 2+: polynomials (comma-separated, may span multiple lines)\n");
    printf("  Last   : variables TO ELIMINATE (comma-separated)\n");

    printf("FILE FORMAT (Polynomial Solver, multiline):\n");
    printf("  Line 1 : field size\n");
    printf("  Line 2+: polynomials (one per line or comma-separated)\n");

    printf("EXAMPLES:\n");
    printf("  %s \"x+y+z, x*y+y*z+z*x, x*y*z+1\" \"x,y\" 257\n", prog_name);
    printf("  %s --comp \"x^2+y^2+1, x*y+z, x+y+z^2\" \"x,y\" 257\n", prog_name);
    printf("  %s --solve \"x^2+y^2+z^2-6, x+y+z-4, x*y*z-x-1\" 257\n", prog_name);
    printf("  %s --solve \"x^2 + t*y, x*y + t^2\" \"2^8: t^8 + t^4 + t^3 + t + 1\"\n", prog_name);
    printf("  (AES polynomial for F_256, 't' is the field extension generator)\n");
    printf("  %s --silent \"x+y^2+t, x*y+t*y+1\" \"x\" 2^8\n", prog_name);
    printf("  %s example.dat\n", prog_name);
}

/* =========================================================================
 * Utility helpers
 * ========================================================================= */
static char *trim(char *str)
{
    char *end;
    while (isspace((unsigned char)*str)) str++;
    if (*str == 0) return str;
    end = str + strlen(str) - 1;
    while (end > str && isspace((unsigned char)*end)) end--;
    end[1] = '\0';
    return str;
}

static int check_prime_power(const fmpz_t n, fmpz_t prime, ulong *power)
{
    if (fmpz_cmp_ui(n, 1) <= 0) return 0;
    if (fmpz_is_probabprime(n)) {
        fmpz_set(prime, n);
        *power = 1;
        return 1;
    }
    fmpz_factor_t factors;
    fmpz_factor_init(factors);
    fmpz_factor(factors, n);
    if (factors->num == 1) {
        fmpz_set(prime, factors->p + 0);
        *power = factors->exp[0];
        fmpz_factor_clear(factors);
        return 1;
    }
    fmpz_factor_clear(factors);
    return 0;
}

static int parse_field_polynomial(nmod_poly_t modulus, const char *poly_str,
                                  mp_limb_t prime, const char *var_name)
{
    if (!poly_str || !var_name) return 0;
    nmod_poly_zero(modulus);
    char *work_str = strdup(poly_str);

    char *dst = work_str;
    for (char *src = work_str; *src; src++)
        if (!isspace(*src)) *dst++ = *src;
    *dst = '\0';

    char *token = work_str;
    while (*token) {
        if (*token == '+' || *token == '-') { token++; continue; }

        char *term_end = token;
        while (*term_end && *term_end != '+' && *term_end != '-') term_end++;

        char term_char = *term_end;
        *term_end = '\0';

        mp_limb_t coeff  = 1;
        ulong      degree = 0;
        char      *var_pos = strstr(token, var_name);

        if (!var_pos) {
            if (strlen(token) > 0) coeff = strtoul(token, NULL, 10);
            degree = 0;
        } else {
            if (var_pos > token) {
                size_t coeff_len = var_pos - token;
                char  *coeff_str = malloc(coeff_len + 1);
                strncpy(coeff_str, token, coeff_len);
                coeff_str[coeff_len] = '\0';
                size_t len = strlen(coeff_str);
                if (len > 0 && coeff_str[len - 1] == '*') coeff_str[len - 1] = '\0';
                if (strlen(coeff_str) > 0) coeff = strtoul(coeff_str, NULL, 10);
                else                        coeff = 1;
                free(coeff_str);
            }
            char *degree_pos = var_pos + strlen(var_name);
            if (*degree_pos == '^') degree = strtoul(degree_pos + 1, NULL, 10);
            else                    degree = 1;
        }

        mp_limb_t existing = nmod_poly_get_coeff_ui(modulus, degree);
        nmod_poly_set_coeff_ui(modulus, degree, (existing + coeff) % prime);

        *term_end = term_char;
        token = term_end;
    }
    free(work_str);
    return 1;
}

static int parse_field_size(const char *field_str, fmpz_t prime, ulong *power,
                            char **field_poly, char **gen_var)
{
    if (!field_str || strlen(field_str) == 0) return 0;
    if (field_poly) *field_poly = NULL;
    if (gen_var)    *gen_var    = NULL;

    const char *colon = strchr(field_str, ':');
    if (colon) {
        size_t  size_len   = colon - field_str;
        char   *size_part  = malloc(size_len + 1);
        strncpy(size_part, field_str, size_len);
        size_part[size_len] = '\0';
        char *trimmed_size = trim(size_part);

        if (field_poly) {
            const char *poly_start = colon + 1;
            while (isspace(*poly_start)) poly_start++;
            *field_poly = strdup(poly_start);
            if (gen_var && *field_poly) {
                const char *p = *field_poly;
                while (*p && !isalpha(*p)) p++;
                if (*p && isalpha(*p)) {
                    char var_name[2] = { *p, '\0' };
                    *gen_var = strdup(var_name);
                }
            }
        }
        int result = parse_field_size(trimmed_size, prime, power, NULL, NULL);
        free(size_part);
        return result;
    }

    const char *caret = strchr(field_str, '^');
    if (caret) {
        char *prime_str = malloc(caret - field_str + 1);
        strncpy(prime_str, field_str, caret - field_str);
        prime_str[caret - field_str] = '\0';

        fmpz_t p;
        fmpz_init(p);
        if (fmpz_set_str(p, prime_str, 10) != 0) {
            fmpz_clear(p); free(prime_str); return 0;
        }
        char *endptr;
        ulong k = strtoul(caret + 1, &endptr, 10);
        if (*endptr != '\0' || k == 0) {
            fmpz_clear(p); free(prime_str); return 0;
        }
        if (!fmpz_is_probabprime(p)) {
            fmpz_clear(p); free(prime_str); return 0;
        }
        fmpz_set(prime, p);
        *power = k;
        fmpz_clear(p);
        free(prime_str);
        return 1;
    } else {
        fmpz_t field_size;
        fmpz_init(field_size);
        int success = fmpz_set_str(field_size, field_str, 10);
        if (success != 0) { fmpz_clear(field_size); return 0; }
        int result = check_prime_power(field_size, prime, power);
        fmpz_clear(field_size);
        return result;
    }
}

static char *generate_output_filename(const char *input_filename)
{
    if (!input_filename) return NULL;
    const char *dot = strrchr(input_filename, '.');
    if (dot) {
        size_t base_len = dot - input_filename;
        size_t ext_len  = strlen(dot);
        char  *output   = malloc(base_len + 9 + ext_len + 1);
        strncpy(output, input_filename, base_len);
        output[base_len] = '\0';
        strcat(output, "_solution");
        strcat(output, dot);
        return output;
    } else {
        size_t len    = strlen(input_filename);
        char  *output = malloc(len + 9 + 1);
        strcpy(output, input_filename);
        strcat(output, "_solution");
        return output;
    }
}

/* Generate _comp variant of the output filename */
static char *generate_comp_filename(const char *input_filename)
{
    if (!input_filename) return NULL;
    const char *dot = strrchr(input_filename, '.');
    if (dot) {
        size_t base_len = dot - input_filename;
        size_t ext_len  = strlen(dot);
        char  *output   = malloc(base_len + 6 + ext_len + 1); /* 6 = "_comp" */
        strncpy(output, input_filename, base_len);
        output[base_len] = '\0';
        strcat(output, "_comp");
        strcat(output, dot);
        return output;
    } else {
        size_t len    = strlen(input_filename);
        char  *output = malloc(len + 6 + 1);
        strcpy(output, input_filename);
        strcat(output, "_comp");
        return output;
    }
}

/* =========================================================================
 * Complexity analysis: save to file
 * ========================================================================= */
static void save_comp_result_to_file(
        const char   *filename,
        const char   *polys_str,
        const char   *vars_str,
        mp_limb_t     prime,
        ulong         power,
        slong         num_polys,
        slong         num_all_vars,
        char        **all_vars,
        char        **elim_var_list,
        slong         num_elim,
        const long   *degrees,
        const fmpz_t  matrix_size,
        long          bezout_bound,
        double        complexity,
        double        omega,
        double        comp_time)
{
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Warning: Cannot create output file '%s'\n", filename);
        return;
    }

    mp_limb_t field_size = 1;
    for (ulong i = 0; i < power; i++) field_size *= prime;

    fprintf(fp, "Dixon Complexity Analysis\n");
    fprintf(fp, "=========================\n");
    fprintf(fp, "Field: F_%lu", prime);
    if (power > 1) fprintf(fp, "^%lu (size %lu)", power, field_size);
    fprintf(fp, "\n");
    fprintf(fp, "Polynomials: %s\n", polys_str);
    fprintf(fp, "Eliminate:   %s\n", vars_str);
    fprintf(fp, "Computation time: %.3f seconds\n\n", comp_time);

    fprintf(fp, "--- Input Summary ---\n");
    fprintf(fp, "Equations : %ld\n", num_polys);
    fprintf(fp, "Variables (%ld): ", num_all_vars);
    for (slong i = 0; i < num_all_vars; i++) {
        if (i > 0) fprintf(fp, ", ");
        fprintf(fp, "%s", all_vars[i]);
    }
    fprintf(fp, "\n");
    fprintf(fp, "Elim vars (%ld): ", num_elim);
    for (slong i = 0; i < num_elim; i++) {
        if (i > 0) fprintf(fp, ", ");
        fprintf(fp, "%s", elim_var_list[i]);
    }
    fprintf(fp, "\n");

    /* Remaining vars */
    fprintf(fp, "Remaining vars: ");
    int first = 1;
    for (slong i = 0; i < num_all_vars; i++) {
        int is_elim = 0;
        for (slong j = 0; j < num_elim; j++)
            if (strcmp(all_vars[i], elim_var_list[j]) == 0) { is_elim = 1; break; }
        if (!is_elim) {
            if (!first) fprintf(fp, ", ");
            fprintf(fp, "%s", all_vars[i]);
            first = 0;
        }
    }
    if (first) fprintf(fp, "(none)");
    fprintf(fp, "\n\n");

    fprintf(fp, "--- Degree & Size ---\n");
    fprintf(fp, "Degree sequence: [");
    for (slong i = 0; i < num_polys; i++) {
        if (i > 0) fprintf(fp, ", ");
        fprintf(fp, "%ld", degrees[i]);
    }
    fprintf(fp, "]\n");
    fprintf(fp, "Bezout bound (degree product): %ld\n", bezout_bound);
    fprintf(fp, "Dixon matrix size: ");
    fmpz_fprint(fp, matrix_size);
    fprintf(fp, "\n");
    fprintf(fp, "Resultant degree estimate (Bezout): %ld\n", bezout_bound);
    fprintf(fp, "\n--- Complexity ---\n");
    fprintf(fp, "Complexity (log2, omega=%.4g): %.6f\n", omega, complexity);

    fclose(fp);
}

/* =========================================================================
 * Complexity analysis: run and display
 * ========================================================================= */
static void run_complexity_analysis(
        const char      *polys_str,
        const char      *vars_str,
        mp_limb_t        prime,
        ulong            power,
        const fq_nmod_ctx_t ctx,
        const char      *output_filename,
        int              silent_mode,
        double           comp_time,
        double           omega)
{
    /* ---- split polynomials and elimination variables ---- */
    slong   num_polys;
    char  **poly_arr = split_string(polys_str, &num_polys);

    slong   num_elim;
    char  **elim_arr = split_string(vars_str, &num_elim);

    /* ---- get generator name ---- */
    char *gen_name = get_generator_name(ctx);  /* may return NULL for prime fields */

    /* ---- collect all variables ---- */
    char  **all_vars;
    slong   num_all_vars;
    collect_variables((const char **)poly_arr, num_polys,
                      gen_name, &all_vars, &num_all_vars);

    /* ---- compute degree of each polynomial ---- */
    if (num_polys <= 0) {
        if (!silent_mode) fprintf(stderr, "Error: no polynomials to analyze\n");
        free_split_strings(poly_arr, num_polys);
        free_split_strings(elim_arr, num_elim);
        if (gen_name) free(gen_name);
        for (slong i = 0; i < num_all_vars; i++) free(all_vars[i]);
        free(all_vars);
        return;
    }
    long *degrees = calloc((size_t)num_polys, sizeof(long));
    for (slong i = 0; i < num_polys; i++)
        degrees[i] = get_poly_total_degree(poly_arr[i], gen_name);

    /* ---- Bezout bound = product of degrees ---- */
    long bezout = 1;
    for (slong i = 0; i < num_polys; i++) bezout *= degrees[i];

    /* ---- Dixon matrix size via Hessenberg recurrence ---- */
    fmpz_t matrix_size;
    fmpz_init(matrix_size);
    /* suppress internal prints from dixon_size */
    {
        int orig_stdout = dup(STDOUT_FILENO);
        int devnull     = open("/dev/null", O_WRONLY);
        if (devnull != -1) { dup2(devnull, STDOUT_FILENO); close(devnull); }
        dixon_size(matrix_size, degrees, (int)num_polys, 0);
        fflush(stdout);
        dup2(orig_stdout, STDOUT_FILENO);
        close(orig_stdout);
    }

    /* ---- Dixon complexity ---- */
    double complexity = dixon_complexity(degrees, (int)num_polys,
                                        (int)num_all_vars, omega);

    /* ---- console output (concise) ---- */
    if (!silent_mode) {
        printf("\n=== Complexity Analysis ===\n");
        printf("Equations: %ld  |  Variables: %ld  |  Eliminate: %ld\n",
               num_polys, num_all_vars, num_elim);

        printf("All vars : ");
        for (slong i = 0; i < num_all_vars; i++) {
            if (i > 0) printf(", ");
            printf("%s", all_vars[i]);
        }
        printf("\n");

        printf("Degrees  : [");
        for (slong i = 0; i < num_polys; i++) {
            if (i > 0) printf(", ");
            printf("%ld", degrees[i]);
        }
        printf("]\n");

        printf("Bezout bound      : %ld\n", bezout);

        printf("Dixon matrix size : ");
        fmpz_print(matrix_size);
        printf("\n");

        printf("Resultant deg est : %ld  (Bezout bound)\n", bezout);

        if (isfinite(complexity))
            printf("Complexity (log2) : %.4f  (omega=%.4g)\n",
                   complexity, omega);
        else
            printf("Complexity (log2) : inf / undefined\n");

        if (output_filename)
            printf("Report saved to   : %s\n", output_filename);
        printf("===========================\n");
    }

    /* ---- save to file ---- */
    if (output_filename) {
        save_comp_result_to_file(
            output_filename, polys_str, vars_str,
            prime, power,
            num_polys, num_all_vars, all_vars,
            elim_arr, num_elim,
            degrees, matrix_size, bezout,
            complexity, omega, comp_time);
    }

    /* ---- cleanup ---- */
    fmpz_clear(matrix_size);
    free(degrees);
    for (slong i = 0; i < num_all_vars; i++) free(all_vars[i]);
    free(all_vars);
    if (gen_name) free(gen_name);
    free_split_strings(poly_arr, num_polys);
    free_split_strings(elim_arr, num_elim);
}

/* =========================================================================
 * File reading helpers (unchanged from original)
 * ========================================================================= */
static char *read_entire_line(FILE *fp)
{
    if (!fp) return NULL;
    size_t capacity = 4096, length = 0;
    char  *line = malloc(capacity);
    if (!line) return NULL;
    int c;
    while ((c = fgetc(fp)) != EOF && c != '\n' && c != '\r') {
        if (length + 1 >= capacity) {
            capacity *= 2;
            char *nl = realloc(line, capacity);
            if (!nl) { free(line); return NULL; }
            line = nl;
        }
        line[length++] = (char)c;
    }
    if (c == '\r') {
        int next = fgetc(fp);
        if (next != '\n' && next != EOF) ungetc(next, fp);
    }
    if (length == 0 && c == EOF) { free(line); return NULL; }
    line[length] = '\0';
    return line;
}

static int read_multiline_file(FILE *fp, char **field_str, char **polys_str,
                               char **vars_str, char **ideal_str,
                               char **allvars_str)
{
    char **lines   = NULL;
    int   line_count = 0, line_capacity = 10;
    lines = malloc(line_capacity * sizeof(char *));
    if (!lines) return 0;

    char *line;
    while ((line = read_entire_line(fp)) != NULL) {
        char *trimmed = trim(line);
        if (strlen(trimmed) == 0 || trimmed[0] == '#') { free(line); continue; }
        if (line_count >= line_capacity) {
            line_capacity *= 2;
            char **nl = realloc(lines, line_capacity * sizeof(char *));
            if (!nl) {
                for (int i = 0; i < line_count; i++) free(lines[i]);
                free(lines); free(line); return 0;
            }
            lines = nl;
        }
        lines[line_count++] = strdup(trimmed);
        free(line);
    }

    if (line_count < 3) {
        fprintf(stderr, "Error: File must contain at least 3 non-empty lines\n");
        fprintf(stderr, "  Line 1: field size\n");
        fprintf(stderr, "  Lines 2 to n-1: polynomials\n");
        fprintf(stderr, "  Line n: variables to ELIMINATE\n");
        for (int i = 0; i < line_count; i++) free(lines[i]);
        free(lines);
        return 0;
    }

    *field_str = lines[0];
    *vars_str  = lines[line_count - 1];

    size_t total_len = 0;
    for (int i = 1; i < line_count - 1; i++)
        total_len += strlen(lines[i]) + 2;

    char *poly_buffer = malloc(total_len + 1);
    if (!poly_buffer) {
        for (int i = 0; i < line_count; i++) free(lines[i]);
        free(lines);
        return 0;
    }
    poly_buffer[0] = '\0';
    for (int i = 1; i < line_count - 1; i++) {
        if (i > 1) strcat(poly_buffer, " ");
        strcat(poly_buffer, lines[i]);
    }
    *polys_str = poly_buffer;

    for (int i = 1; i < line_count - 1; i++) free(lines[i]);
    free(lines);

    *ideal_str    = NULL;
    *allvars_str  = NULL;
    return 1;
}

static int read_solver_file(FILE *fp, char **field_str, char **polys_str)
{
    char **lines   = NULL;
    int   line_count = 0, line_capacity = 10;
    lines = malloc(line_capacity * sizeof(char *));
    if (!lines) return 0;

    char *line;
    while ((line = read_entire_line(fp)) != NULL) {
        char *trimmed = trim(line);
        if (strlen(trimmed) == 0 || trimmed[0] == '#') { free(line); continue; }
        if (line_count >= line_capacity) {
            line_capacity *= 2;
            char **nl = realloc(lines, line_capacity * sizeof(char *));
            if (!nl) {
                for (int i = 0; i < line_count; i++) free(lines[i]);
                free(lines); free(line); return 0;
            }
            lines = nl;
        }
        lines[line_count++] = strdup(trimmed);
        free(line);
    }

    if (line_count < 2) {
        fprintf(stderr, "Error: Solver file must contain at least 2 non-empty lines\n");
        for (int i = 0; i < line_count; i++) free(lines[i]);
        free(lines);
        return 0;
    }

    *field_str = lines[0];

    size_t total_len = 0;
    for (int i = 1; i < line_count; i++)
        total_len += strlen(lines[i]) + 3;

    char *poly_buffer = malloc(total_len + 1);
    if (!poly_buffer) {
        for (int i = 0; i < line_count; i++) free(lines[i]);
        free(lines);
        return 0;
    }
    poly_buffer[0] = '\0';
    for (int i = 1; i < line_count; i++) {
        if (i > 1) {
            size_t prev_len = strlen(poly_buffer);
            int prev_comma = (prev_len > 0 &&
                (poly_buffer[prev_len - 1] == ',' ||
                 (prev_len > 1 && poly_buffer[prev_len - 2] == ',' &&
                  poly_buffer[prev_len - 1] == ' ')));
            int curr_comma = (lines[i][0] == ',');
            if (!prev_comma && !curr_comma)      strcat(poly_buffer, ", ");
            else if (!prev_comma && curr_comma)  strcat(poly_buffer, " ");
        }
        strcat(poly_buffer, lines[i]);
    }
    *polys_str = poly_buffer;

    for (int i = 1; i < line_count; i++) free(lines[i]);
    free(lines);
    return 1;
}

/* =========================================================================
 * Result saving helpers (unchanged from original)
 * ========================================================================= */
static void save_solver_result_to_file(const char *filename,
                                       const char *polys_str,
                                       mp_limb_t prime, ulong power,
                                       const polynomial_solutions_t *sols,
                                       double computation_time)
{
    FILE *out_fp = fopen(filename, "w");
    if (!out_fp) {
        fprintf(stderr, "Warning: Could not create output file '%s'\n", filename);
        return;
    }
    mp_limb_t field_size = 1;
    for (ulong i = 0; i < power; i++) field_size *= prime;

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
    fprintf(out_fp, "\nSolutions:\n==========\n");

    if (!sols) { fprintf(out_fp, "Solution structure is null\n"); fclose(out_fp); return; }

    fprintf(out_fp, "\n=== Polynomial System Solutions ===\n");
    if (!sols->is_valid) {
        fprintf(out_fp, "Solving failed");
        if (sols->error_message) fprintf(out_fp, ": %s", sols->error_message);
        fprintf(out_fp, "\n");
        fclose(out_fp); return;
    }
    if (sols->has_no_solutions) {
        fprintf(out_fp, "System has no solutions over the finite field\n");
        fclose(out_fp); return;
    }
    if (sols->num_variables == 0) {
        fprintf(out_fp, "No variables\n"); fclose(out_fp); return;
    }
    if (sols->num_solution_sets == 0) {
        fprintf(out_fp, "No solutions found\n"); fclose(out_fp); return;
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
                char *sol_str = fq_nmod_get_str_pretty(sols->solution_sets[set][var][0], sols->ctx);
                fprintf(out_fp, "%s", sol_str); free(sol_str);
            } else {
                fprintf(out_fp, "{");
                for (slong sol = 0; sol < num_sols; sol++) {
                    if (sol > 0) fprintf(out_fp, ", ");
                    char *sol_str = fq_nmod_get_str_pretty(sols->solution_sets[set][var][sol], sols->ctx);
                    fprintf(out_fp, "%s", sol_str); free(sol_str);
                }
                fprintf(out_fp, "}");
            }
            fprintf(out_fp, "\n");
        }
    }

    fprintf(out_fp, "\n=== Compatibility View ===\n");
    for (slong var = 0; var < sols->num_variables; var++) {
        fprintf(out_fp, "%s = {", sols->variable_names[var]);
        slong total_printed = 0;
        for (slong set = 0; set < sols->num_solution_sets; set++) {
            slong num_sols = sols->solutions_per_var[set * sols->num_variables + var];
            for (slong sol = 0; sol < num_sols; sol++) {
                if (total_printed > 0) fprintf(out_fp, ", ");
                char *sol_str = fq_nmod_get_str_pretty(sols->solution_sets[set][var][sol], sols->ctx);
                fprintf(out_fp, "%s", sol_str); free(sol_str);
                total_printed++;
            }
        }
        fprintf(out_fp, "}");
        if (total_printed > 1)      fprintf(out_fp, " (%ld solutions)", total_printed);
        else if (total_printed == 0) fprintf(out_fp, " (no solutions)");
        fprintf(out_fp, "\n");
    }
    fprintf(out_fp, "=== Solution Complete ===\n\n");
    fclose(out_fp);
}

static void save_result_to_file(const char *filename,
                                const char *polys_str,
                                const char *vars_str,
                                const char *ideal_str,
                                const char *allvars_str,
                                mp_limb_t prime, ulong power,
                                const char *result,
                                double computation_time)
{
    FILE *out_fp = fopen(filename, "w");
    if (!out_fp) {
        fprintf(stderr, "Warning: Could not create output file '%s'\n", filename);
        return;
    }
    mp_limb_t field_size = 1;
    for (ulong i = 0; i < power; i++) field_size *= prime;

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
    fprintf(out_fp, "\nResultant:\n%s\n", result);
    fclose(out_fp);
}

static int count_comma_separated_items(const char *str)
{
    if (!str || strlen(str) == 0) return 0;
    int count = 1;
    for (const char *p = str; *p; p++)
        if (*p == ',') count++;
    return count;
}

/* =========================================================================
 * main
 * ========================================================================= */
int main(int argc, char *argv[])
{
    clock_t start_time = clock();

    if (argc == 1) { print_usage(argv[0]); return 0; }

    /* ---- parse leading flags ---- */
    int    silent_mode = 0;
    int    solve_mode  = 0;
    int    comp_mode   = 0;
    int    arg_offset  = 0;
    double omega       = DIXON_OMEGA;   /* default, overridden by --omega */

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--silent") == 0) {
            silent_mode = 1; arg_offset++;
        } else if (strcmp(argv[i], "--solve") == 0) {
            solve_mode = 1; arg_offset++;
        } else if (strcmp(argv[i], "--comp") == 0 ||
                   strcmp(argv[i], "-c")     == 0) {
            comp_mode = 1; arg_offset++;
        } else if ((strcmp(argv[i], "--omega") == 0 ||
                    strcmp(argv[i], "-w")      == 0) && i + 1 < argc) {
            char *endptr = NULL;
            double val = strtod(argv[i + 1], &endptr);
            if (endptr && *endptr == '\0' && val > 0.0) {
                omega = val;
            } else {
                fprintf(stderr, "Warning: invalid --omega value '%s', "
                                "using default %.4g\n", argv[i + 1], omega);
            }
            arg_offset += 2;
            i++;          /* skip the value token */
        } else {
            break;
        }
    }

    int    effective_argc = argc - arg_offset;
    char **effective_argv = argv + arg_offset;

    /* ---- help ---- */
    if (effective_argc >= 2 &&
        (strcmp(effective_argv[1], "--help") == 0 ||
         strcmp(effective_argv[1], "-h")     == 0)) {
        if (!silent_mode) print_usage(argv[0]);
        return 0;
    }

    /* ---- version banner ---- */
    if (!silent_mode) {
        printf("=================================================\n");
        printf("Dixon Resultant & Polynomial Solver v%s\n", PROGRAM_VERSION);
        printf("FLINT version: %s\n", FLINT_VERSION);
#ifdef HAVE_PML
        printf("PML support: ENABLED\n");
#else
        printf("PML support: DISABLED\n");
#endif
        printf("=================================================\n\n");
    }

    /* ---- test modes ---- */
    if (argc >= 2 && strcmp(argv[1], "--test") == 0) {
        if (!silent_mode) {
            if (argc >= 3) test_dixon(atoi(argv[2]));
            else           test_dixon(0);
        }
        return 0;
    }
    if (effective_argc >= 2 &&
        strcmp(effective_argv[1], "--test-solver") == 0) {
        if (!silent_mode) test_polynomial_solver();
        return 0;
    }

    /* ---- input variables ---- */
    char *polys_str     = NULL;
    char *vars_str      = NULL;
    char *ideal_str     = NULL;
    char *allvars_str   = NULL;
    char *field_str     = NULL;
    int   need_free     = 0;
    char *input_filename  = NULL;
    char *output_filename = NULL;
    mp_limb_t prime = 0;
    ulong     power = 0;

    /* ---- determine input mode ---- */
    if (comp_mode) {
        /* Complexity analysis mode: same argument format as basic Dixon */
        if (effective_argc == 2) {
            /* file input */
            FILE *fp = fopen(effective_argv[1], "r");
            if (!fp) {
                if (!silent_mode)
                    fprintf(stderr, "Error: Cannot open file '%s'\n",
                            effective_argv[1]);
                return 1;
            }
            if (!silent_mode)
                printf("Reading from file: %s\n", effective_argv[1]);

            input_filename  = strdup(effective_argv[1]);
            output_filename = generate_comp_filename(input_filename);

            if (!read_multiline_file(fp, &field_str, &polys_str,
                                     &vars_str, &ideal_str, &allvars_str)) {
                fclose(fp); return 1;
            }
            fclose(fp);
            need_free = 1;

        } else if (effective_argc == 4) {
            /* command line: polys vars field_size */
            polys_str = effective_argv[1];
            vars_str  = effective_argv[2];
            field_str = effective_argv[3];

            char buffer[128];
            time_t now = time(NULL);
            struct tm *t = localtime(&now);
            if (t) strftime(buffer, sizeof(buffer), "comp_%Y%m%d_%H%M%S.dat", t);
            else   strcpy(buffer, "comp.dat");
            output_filename = strdup(buffer);

        } else {
            if (!silent_mode) {
                fprintf(stderr, "Error: Complexity mode requires:\n");
                fprintf(stderr, "  %s --comp \"polynomials\" \"eliminate_vars\" field_size\n",
                        argv[0]);
                fprintf(stderr, "  %s --comp input_file\n", argv[0]);
            }
            return 1;
        }

    } else if (solve_mode) {
        if (effective_argc == 2) {
            FILE *fp = fopen(effective_argv[1], "r");
            if (!fp) {
                if (!silent_mode)
                    fprintf(stderr, "Error: Cannot open file '%s'\n",
                            effective_argv[1]);
                return 1;
            }
            if (!silent_mode)
                printf("Reading polynomial system from file: %s\n",
                       effective_argv[1]);

            input_filename  = strdup(effective_argv[1]);
            output_filename = generate_output_filename(input_filename);

            if (!read_solver_file(fp, &field_str, &polys_str)) {
                fclose(fp); return 1;
            }
            fclose(fp);
            need_free = 1;

        } else if (effective_argc == 3) {
            polys_str = effective_argv[1];
            field_str = effective_argv[2];
            char buffer[128];
            time_t now = time(NULL);
            struct tm *t = localtime(&now);
            if (t) strftime(buffer, sizeof(buffer), "solution_%Y%m%d_%H%M%S.dat", t);
            else   strcpy(buffer, "solution.dat");
            output_filename = strdup(buffer);

        } else {
            if (!silent_mode) {
                fprintf(stderr, "Error: Solver mode requires either:\n");
                fprintf(stderr, "  %s --solve \"polynomials\" field_size\n", argv[0]);
                fprintf(stderr, "  %s --solve input_file\n", argv[0]);
            }
            return 1;
        }

    } else {
        /* Basic Dixon */
        if (effective_argc == 2) {
            FILE *fp = fopen(effective_argv[1], "r");
            if (!fp) {
                if (!silent_mode)
                    fprintf(stderr, "Error: Cannot open file '%s'\n",
                            effective_argv[1]);
                return 1;
            }
            if (!silent_mode)
                printf("Reading from file: %s\n", effective_argv[1]);

            input_filename  = strdup(effective_argv[1]);
            output_filename = generate_output_filename(input_filename);

            if (!read_multiline_file(fp, &field_str, &polys_str,
                                     &vars_str, &ideal_str, &allvars_str)) {
                fclose(fp); return 1;
            }
            fclose(fp);
            need_free = 1;

        } else if (effective_argc == 4 || effective_argc == 5) {
            polys_str = effective_argv[1];
            vars_str  = effective_argv[2];
            if (effective_argc == 5) {
                ideal_str = effective_argv[3];
                field_str = effective_argv[4];
            } else {
                field_str = effective_argv[3];
            }
            char buffer[128];
            time_t now = time(NULL);
            struct tm *t = localtime(&now);
            if (t) strftime(buffer, sizeof(buffer), "solution_%Y%m%d_%H%M%S.dat", t);
            else   strcpy(buffer, "solution.dat");
            output_filename = strdup(buffer);

        } else {
            if (!silent_mode) print_usage(argv[0]);
            return 1;
        }
    }

    /* ---- parse field size ---- */
    fmpz_t p_fmpz;
    fmpz_init(p_fmpz);
    char *field_poly_str = NULL;
    char *gen_var_name   = NULL;

    if (!parse_field_size(field_str, p_fmpz, &power,
                          &field_poly_str, &gen_var_name)) {
        if (!silent_mode) {
            fprintf(stderr, "Error: Invalid field size '%s'\n", field_str);
            fprintf(stderr, "Field size must be a prime, prime power (e.g. 256), or p^k (e.g. 2^8)\n");
        }
        if (need_free) {
            free(field_str); free(polys_str);
            if (vars_str)    free(vars_str);
            if (ideal_str)   free(ideal_str);
        }
        if (input_filename)  free(input_filename);
        if (output_filename) free(output_filename);
        fmpz_clear(p_fmpz);
        if (field_poly_str) free(field_poly_str);
        if (gen_var_name)   free(gen_var_name);
        return 1;
    }

    prime = fmpz_get_ui(p_fmpz);

    mp_limb_t field_size = 1;
    for (ulong i = 0; i < power; i++) field_size *= prime;

    if (!silent_mode) {
        printf("=== %s ===\n",
               comp_mode  ? "Complexity Analysis" :
               solve_mode ? "Polynomial System Solver" :
                            "Dixon Resultant Computation");
        printf("Field: F_%lu", prime);
        if (power > 1) {
            printf("^%lu (size %lu)", power, field_size);
            printf("\nField extension generator: t");
        }
        printf("\n");
    }

    /* ---- initialize finite field ---- */
    fq_nmod_ctx_t ctx;

    if (power > 1 && field_poly_str) {
        const char *var_name = gen_var_name ? gen_var_name : "t";
        if (!silent_mode) {
            printf("Using custom field polynomial: %s\n", field_poly_str);
            printf("Generator variable: %s\n", var_name);
        }
        nmod_poly_t modulus;
        nmod_poly_init(modulus, prime);
        if (!parse_field_polynomial(modulus, field_poly_str, prime, var_name)) {
            fprintf(stderr, "Error: Failed to parse field polynomial\n");
            nmod_poly_clear(modulus); fmpz_clear(p_fmpz);
            if (field_poly_str) free(field_poly_str);
            if (gen_var_name)   free(gen_var_name);
            return 1;
        }
        if (nmod_poly_degree(modulus) != (slong)power) {
            fprintf(stderr, "Error: Polynomial degree %ld doesn't match power %lu\n",
                    nmod_poly_degree(modulus), power);
            nmod_poly_clear(modulus); fmpz_clear(p_fmpz);
            if (field_poly_str) free(field_poly_str);
            if (gen_var_name)   free(gen_var_name);
            return 1;
        }
        fq_nmod_ctx_init_modulus(ctx, modulus, var_name);
        nmod_poly_clear(modulus);

    } else {
        fq_nmod_ctx_init(ctx, p_fmpz, power, "t");

        if (!silent_mode && power > 1) {
            printf("Using FLINT's default irreducible polynomial:\n  ");
            const nmod_poly_struct *modulus = ctx->modulus;
            int first_term = 1;
            for (slong i = nmod_poly_degree(modulus); i >= 0; i--) {
                mp_limb_t coeff = nmod_poly_get_coeff_ui(modulus, i);
                if (coeff != 0) {
                    if (!first_term) printf(" + ");
                    if      (i == 0) printf("%lu", coeff);
                    else if (i == 1) { if (coeff == 1) printf("t"); else printf("%lu*t", coeff); }
                    else             { if (coeff == 1) printf("t^%ld", i); else printf("%lu*t^%ld", coeff, i); }
                    first_term = 0;
                }
            }
            printf("\n");
        }
    }

    /* ======================================================
     * EXECUTE requested mode
     * ====================================================== */

    char *result = NULL;
    polynomial_solutions_t *solutions = NULL;

    if (comp_mode) {
        /* ---- Complexity analysis ---- */
        if (!silent_mode) {
            int poly_count = count_comma_separated_items(polys_str);
            int var_count  = count_comma_separated_items(vars_str);
            printf("\nPolynomials (%d): %s\n", poly_count, polys_str);
            printf("Eliminate   (%d): %s\n", var_count,  vars_str);
            printf("Omega            : %.4g\n", omega);
            printf("--------------------------------\n");
        }

        clock_t end_time  = clock();
        double comp_time  = (double)(end_time - start_time) / CLOCKS_PER_SEC;

        run_complexity_analysis(polys_str, vars_str,
                                prime, power, ctx,
                                output_filename, silent_mode,
                                comp_time, omega);

    } else if (solve_mode) {
        /* ---- Polynomial system solver ---- */
        if (!silent_mode) {
            printf("\nMode: Polynomial System Solver\n");
            int poly_count = count_comma_separated_items(polys_str);
            printf("Number of polynomials: %d\n", poly_count);
            printf("Polynomials: %s\n", polys_str);
            printf("Note: Number of variables must equal number of equations\n");
            printf("--------------------------------\n");
        }

        int orig_stdout = -1, orig_stderr = -1, devnull = -1;
        if (silent_mode) {
            fflush(stdout); fflush(stderr);
            orig_stdout = dup(STDOUT_FILENO);
            orig_stderr = dup(STDERR_FILENO);
            devnull = open("/dev/null", O_WRONLY);
            if (devnull != -1) {
                dup2(devnull, STDOUT_FILENO);
                dup2(devnull, STDERR_FILENO);
                close(devnull);
            }
        }

        solutions = solve_polynomial_system_string(polys_str, ctx);

        if (silent_mode && orig_stdout != -1) {
            fflush(stdout); fflush(stderr);
            dup2(orig_stdout, STDOUT_FILENO);
            dup2(orig_stderr, STDERR_FILENO);
            close(orig_stdout); close(orig_stderr);
        }

    } else if (ideal_str) {
        /* ---- Dixon with ideal reduction ---- */
        if (!silent_mode) {
            printf("\nMode: Dixon with ideal reduction\n");
            printf("Polynomials: %s\n", polys_str);
            printf("Eliminate: %s\n", vars_str);
            printf("Ideal generators: %s\n", ideal_str);
            printf("--------------------------------\n");
        }
        result = dixon_with_ideal_reduction_str(polys_str, vars_str, ideal_str, ctx);

    } else {
        /* ---- Basic Dixon resultant ---- */
        if (!silent_mode) {
            printf("\nMode: Basic Dixon resultant\n");
            int poly_count = count_comma_separated_items(polys_str);
            int var_count  = count_comma_separated_items(vars_str);
            printf("Number of polynomials: %d\n", poly_count);
            printf("Variables to ELIMINATE: %s (count: %d)\n", vars_str, var_count);
            if (var_count != poly_count - 1)
                printf("WARNING: Dixon method requires eliminating exactly %d variables "
                       "for %d equations!\n", poly_count - 1, poly_count);
            printf("--------------------------------\n");
        }

        int orig_stdout = -1, orig_stderr = -1, devnull = -1;
        if (silent_mode) {
            fflush(stdout); fflush(stderr);
            orig_stdout = dup(STDOUT_FILENO);
            orig_stderr = dup(STDERR_FILENO);
            devnull = open("/dev/null", O_WRONLY);
            if (devnull != -1) {
                dup2(devnull, STDOUT_FILENO);
                dup2(devnull, STDERR_FILENO);
                close(devnull);
            }
        }

        result = dixon_str(polys_str, vars_str, ctx);

        if (silent_mode && orig_stdout != -1) {
            fflush(stdout); fflush(stderr);
            dup2(orig_stdout, STDOUT_FILENO);
            dup2(orig_stderr, STDERR_FILENO);
            close(orig_stdout); close(orig_stderr);
        }
    }

    /* ---- compute total time ---- */
    clock_t end_time       = clock();
    double  computation_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    /* ---- output results ---- */
    if (comp_mode) {
        /* already printed inside run_complexity_analysis */
        ;
    } else if (solve_mode) {
        if (solutions) {
            printf("\n=== POLYNOMIAL SYSTEM SOLUTIONS ===\n");
            print_polynomial_solutions(solutions);
            printf("====================================\n");

            if (output_filename) {
                save_solver_result_to_file(output_filename, polys_str,
                                           prime, power, solutions,
                                           computation_time);
                if (!silent_mode)
                    printf("\nResult saved to: %s\n", output_filename);
            }
            polynomial_solutions_clear(solutions);
            free(solutions);
        } else {
            if (!silent_mode)
                fprintf(stderr, "\nError: Polynomial system solving failed\n");
        }
    } else {
        if (result) {
            if (output_filename) {
                save_result_to_file(output_filename, polys_str, vars_str,
                                    ideal_str, allvars_str, prime, power,
                                    result, computation_time);
                if (!silent_mode)
                    printf("\nResult saved to: %s\n", output_filename);
            }
            free(result);
        } else {
            if (!silent_mode)
                fprintf(stderr, "\nError: Computation failed\n");
        }
    }

    printf("Total computation time: %.3f seconds\n", computation_time);

    /* ---- cleanup ---- */
    fq_nmod_ctx_clear(ctx);
    if (field_poly_str) free(field_poly_str);
    if (gen_var_name)   free(gen_var_name);
    fmpz_clear(p_fmpz);

    if (need_free) {
        free(field_str);
        free(polys_str);
        if (vars_str)   free(vars_str);
        if (ideal_str)  free(ideal_str);
        if (allvars_str) free(allvars_str);
    }
    if (input_filename)  free(input_filename);
    if (output_filename) free(output_filename);

    cleanup_unified_workspace();
    flint_cleanup();
    return 0;
}
