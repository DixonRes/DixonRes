// dixon_magma_implementation.c
// Pure Dixon resultant implementation for the Magma script computation
// gcc -O3 -march=native -fopenmp -DHAVE_PML -o vision vision.c -L/home/suohaohai02/mylinks -lflint -lmpfr -lgmp -lpthread -lstdc++ -L/home/suohaohai02/mylinks -lpml -fopenmp
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <flint/flint.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_poly.h>
#include <flint/fq_nmod_mat.h>
#include <flint/fmpz.h>

// Include the Dixon interface headers
#include "dixon_flint.h"
#include "dixon_interface_flint.h"

// Helper function to setup GF(2^8) with the given irreducible polynomial
void setup_gf28_context(fq_nmod_ctx_t ctx) {
    // Irreducible polynomial: z8^8 + z8^4 + z8^3 + z8 + 1
    // In binary: 100011011 = 0x11B
    nmod_poly_t modulus;
    nmod_poly_init(modulus, 2);
    
    // Set coefficients for z8^8 + z8^4 + z8^3 + z8 + 1
    nmod_poly_set_coeff_ui(modulus, 0, 1);  // constant term
    nmod_poly_set_coeff_ui(modulus, 2, 1);  // z8
    nmod_poly_set_coeff_ui(modulus, 3, 1);  // z8^3
    nmod_poly_set_coeff_ui(modulus, 4, 1);  // z8^4
    nmod_poly_set_coeff_ui(modulus, 8, 1);  // z8^8
    
    fq_nmod_ctx_init_modulus(ctx, modulus, "z8");
    nmod_poly_clear(modulus);
}

// Define all polynomials as global constants
const char* POLY_P0 = "(z8^7 + z8^5 + z8^3)*x1^4*x3^4 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^3 + z8^2)*x2^4*x3^4 + (z8^6 + z8^3 + z8^2 + z8)*x1^4*x3^2 + (z8^6 + z8^5 + z8^3 + 1)*x2^4*x3^2 + (z8^6 + z8^3 + z8^2 + z8)*x1^2*x3^4 + (z8^6 + z8^5 + z8^3 + 1)*x2^2*x3^4 + (z8^7 + z8^6 + z8^5 + 1)*x1^4*x3 + (z8^4 + z8^3 + z8^2 + z8 + 1)*x2^4*x3 + (z8^7 + z8^6 + z8^5 + 1)*x1*x3^4 + (z8^4 + z8^3 + z8^2 + z8 + 1)*x2*x3^4 + (z8^6 + z8^3 + z8 + 1)*x1^2*x3^2 + (z8^7 + z8^6 + z8^5)*x2^2*x3^2 + (z8^7 + z8^6 + z8^4 + z8^3 + z8^2 + 1)*x3^4 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^3 + z8 + 1)*x1^2*x3 + (z8^3)*x2^2*x3 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^3 + z8 + 1)*x1*x3^2 + (z8^3)*x2*x3^2 + (z8^6 + z8^3 + z8^2 + z8)*x1*x3 + (z8^6 + z8^5 + z8^3 + 1)*x2*x3 + (z8^6 + z8^5 + z8^4 + 1)*x3^2 + (z8^7 + z8^5 + z8^2)*x3 + 1";

const char* POLY_P1 = "(z8^7 + z8^6 + z8^5 + z8^2 + 1)*x1^4*x4^4 + (z8^7 + z8^5 + z8^4 + 1)*x2^4*x4^4 + (z8^7 + z8^6 + z8^4 + z8)*x1^4*x4^2 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^2 + 1)*x2^4*x4^2 + (z8^7 + z8^6 + z8^4 + z8)*x1^2*x4^4 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^2 + 1)*x2^2*x4^4 + (z8^5 + z8^4 + z8^3 + z8^2 + z8)*x1^4*x4 + (z8^7 + z8^6)*x2^4*x4 + (z8^5 + z8^4 + z8^3 + z8^2 + z8)*x1*x4^4 + (z8^7 + z8^6)*x2*x4^4 + (z8^7 + z8^6 + z8^4 + z8^3 + z8^2 + 1)*x1^2*x4^2 + (z8^6 + z8^5 + z8^4 + z8^2 + z8)*x2^2*x4^2 + (z8^7 + z8^5 + z8^3)*x4^4 + (z8^4)*x1^2*x4 + (z8^7 + z8^6 + z8^5 + z8 + 1)*x2^2*x4 + (z8^4)*x1*x4^2 + (z8^7 + z8^6 + z8^5 + z8 + 1)*x2*x4^2 + (z8^7 + z8^6 + z8^4 + z8)*x1*x4 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^2 + 1)*x2*x4 + (z8^6 + z8^3 + z8^2 + z8)*x4^2 + (z8^7 + z8^6 + z8^5 + 1)*x4 + 1";

const char* POLY_P2 = "(z8^7 + z8^5 + z8^3)*x5^4*x7^4 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^3 + z8^2)*x6^4*x7^4 + (z8^6 + z8^3 + z8^2 + z8)*x5^4*x7^2 + (z8^6 + z8^5 + z8^3 + 1)*x6^4*x7^2 + (z8^6 + z8^3 + z8^2 + z8)*x5^2*x7^4 + (z8^6 + z8^5 + z8^3 + 1)*x6^2*x7^4 + (z8^7 + z8^6 + z8^5 + 1)*x5^4*x7 + (z8^4 + z8^3 + z8^2 + z8 + 1)*x6^4*x7 + (z8^7 + z8^6 + z8^5 + 1)*x5*x7^4 + (z8^4 + z8^3 + z8^2 + z8 + 1)*x6*x7^4 + (z8^6 + z8^3 + z8 + 1)*x5^2*x7^2 + (z8^7 + z8^6 + z8^5)*x6^2*x7^2 + (z8^6 + z8^4 + z8^3 + z8^2 + z8 + 1)*x7^4 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^3 + z8 + 1)*x5^2*x7 + (z8^3)*x6^2*x7 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^3 + z8 + 1)*x5*x7^2 + (z8^3)*x6*x7^2 + (z8^6 + z8^3 + z8^2 + z8)*x5*x7 + (z8^6 + z8^5 + z8^3 + 1)*x6*x7 + (z8^7 + z8^5 + z8)*x7^2 + (z8^5 + z8^4 + z8^3 + z8)*x7 + 1";

const char* POLY_P3 = "(z8^7 + z8^6 + z8^5 + z8^2 + 1)*x5^4*x8^4 + (z8^7 + z8^5 + z8^4 + 1)*x6^4*x8^4 + (z8^7 + z8^6 + z8^4 + z8)*x5^4*x8^2 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^2 + 1)*x6^4*x8^2 + (z8^7 + z8^6 + z8^4 + z8)*x5^2*x8^4 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^2 + 1)*x6^2*x8^4 + (z8^5 + z8^4 + z8^3 + z8^2 + z8)*x5^4*x8 + (z8^7 + z8^6)*x6^4*x8 + (z8^5 + z8^4 + z8^3 + z8^2 + z8)*x5*x8^4 + (z8^7 + z8^6)*x6*x8^4 + (z8^7 + z8^6 + z8^4 + z8^3 + z8^2 + 1)*x5^2*x8^2 + (z8^6 + z8^5 + z8^4 + z8^2 + z8)*x6^2*x8^2 + (z8^6 + z8^5 + z8^2 + z8)*x8^4 + (z8^4)*x5^2*x8 + (z8^7 + z8^6 + z8^5 + z8 + 1)*x6^2*x8 + (z8^4)*x5*x8^2 + (z8^7 + z8^6 + z8^5 + z8 + 1)*x6*x8^2 + (z8^7 + z8^6 + z8^4 + z8)*x5*x8 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^2 + 1)*x6*x8 + (z8^6 + z8^3 + z8^2)*x8^2 + (z8^7 + z8^5 + z8^4 + z8^3 + z8^2)*x8 + 1";

const char* POLY_P4 = "z8*x0*x1 + (z8^7 + z8^6 + z8^4 + z8^2 + z8)*x1 + 1";

const char* POLY_P5 = "(z8^2 + z8)*x0*x2 + (z8^6 + z8^4 + z8)*x2 + 1";

const char* POLY_P6 = "z8*x3*x5 + (z8 + 1)*x4*x5 + (z8^6 + z8^2 + z8)*x5 + 1";

const char* POLY_P7 = "(z8^2 + z8)*x3*x6 + (z8^2 + z8 + 1)*x4*x6 + (z8^7 + z8^3 + z8^2 + 1)*x6 + 1";

const char* POLY_P8 = "z8*x7*x9 + (z8 + 1)*x8*x9 + (z8^5 + z8^4 + z8)*x9 + 1";

const char* POLY_P9 = "(z8^2 + z8)*x7*x10 + (z8^2 + z8 + 1)*x8*x10 + (z8^5 + z8^4 + z8^3 + z8 + 1)*x10 + 1";

const char* POLY_P10 = "(z8^4 + z8^3 + z8^2)*x9^4 + (z8^4 + z8)*x10^4 + (z8^7 + z8^6 + z8^4 + z8^3 + z8^2 + z8 + 1)*x9^2 + (z8^5 + z8^4 + z8^3 + z8^2 + z8)*x10^2 + (z8^6 + z8^3 + z8^2 + 1)*x9 + (z8^7 + z8^6 + z8^5 + z8^2 + 1)*x10 + (z8^7 + z8^6 + z8^5 + z8^3 + z8^2 + z8)";

// Main implementation following the Magma script
int main() {
    printf("=== Pure Dixon Implementation of Magma Script ===\n\n");
    
    // Setup GF(2^8) context
    fq_nmod_ctx_t ctx;
    setup_gf28_context(ctx);
    printf("Field: GF(2^8) with modulus z8^8 + z8^4 + z8^3 + z8 + 1\n\n");
    
    // Allocate buffer for polynomial strings
    size_t buffer_size = 1024 * 1024 * 10; // 10MB buffer for large polynomials
    char* poly_string = (char*) malloc(buffer_size);
    if (!poly_string) {
        printf("Memory allocation failed!\n");
        return 1;
    }
    
    // Step 1: r1 = Resultant(p4, p5, x0)
    printf("Step 1: Computing r1 = Resultant(p4, p5, x0)\n");
    clock_t start = clock();
    
    snprintf(poly_string, buffer_size, "%s,%s", POLY_P4, POLY_P5);
    char* r1 = dixon_str(poly_string, "x0", ctx);
    
    clock_t end = clock();
    printf("Time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("r1 length: %zu characters\n\n", strlen(r1));
    
    // Step 2: r2 = dixon([r1, p0, p1], [x1, x2])
    printf("Step 2: Computing r2 = dixon([r1, p0, p1], [x1, x2])\n");
    start = clock();
    
    snprintf(poly_string, buffer_size, "%s,%s,%s", r1, POLY_P0, POLY_P1);
    char* r2 = dixon_str(poly_string, "x1,x2", ctx);
    
    end = clock();
    printf("Time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("r2 length: %zu characters\n\n", strlen(r2));
    
    // Step 3: r3 = dixon([r2, p6, p7], [x3, x4])
    printf("Step 3: Computing r3 = dixon([r2, p6, p7], [x3, x4])\n");
    start = clock();
    
    snprintf(poly_string, buffer_size, "%s,%s,%s", r2, POLY_P6, POLY_P7);
    char* r3 = dixon_str(poly_string, "x3,x4", ctx);
    
    end = clock();
    printf("Time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("r3 length: %zu characters\n\n", strlen(r3));
    
    // Step 4: s1 = Resultant(p9, p10, x10)
    printf("Step 4: Computing s1 = Resultant(p9, p10, x10)\n");
    start = clock();
    
    snprintf(poly_string, buffer_size, "%s,%s", POLY_P9, POLY_P10);
    char* s1 = dixon_str(poly_string, "x10", ctx);
    
    end = clock();
    printf("Time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("s1 length: %zu characters\n\n", strlen(s1));
    
    // Step 5: s2 = Resultant(s1, p8, x9)
    printf("Step 5: Computing s2 = Resultant(s1, p8, x9)\n");
    start = clock();
    
    snprintf(poly_string, buffer_size, "%s,%s", s1, POLY_P8);
    char* s2 = dixon_str(poly_string, "x9", ctx);
    
    end = clock();
    printf("Time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("s2 length: %zu characters\n\n", strlen(s2));
    
    // Step 6: s3 = Resultant(s2, p3, x8)
    printf("Step 6: Computing s3 = Resultant(s2, p3, x8)\n");
    start = clock();
    
    snprintf(poly_string, buffer_size, "%s,%s", s2, POLY_P3);
    char* s3 = dixon_str(poly_string, "x8", ctx);
    
    end = clock();
    printf("Time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("s3 length: %zu characters\n\n", strlen(s3));
    
    // Step 7: s4 = Resultant(s3, p2, x7)
    printf("Step 7: Computing s4 = Resultant(s3, p2, x7)\n");
    start = clock();
    
    snprintf(poly_string, buffer_size, "%s,%s", s3, POLY_P2);
    char* s4 = dixon_str(poly_string, "x7", ctx);
    
    end = clock();
    printf("Time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("s4 length: %zu characters\n\n", strlen(s4));
    
    // Step 8: d = Resultant(r3, s4, x5)
    printf("Step 8: Computing d = Resultant(r3, s4, x5)\n");
    start = clock();
    
    snprintf(poly_string, buffer_size, "%s,%s", r3, s4);
    char* d = dixon_ori(poly_string, "x5", ctx); // Eliminate x5 to get univariate polynomial in x6
    
    end = clock();
    printf("Time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("Final resultant d length: %zu characters\n\n", strlen(d));
    
    // The result d should be either 0 or 1
    printf("Final result: %s\n", d);
    
    // If the final result is a univariate polynomial, find its roots
    if (strcmp(d, "0") != 0 && strcmp(d, "1") != 0) {
        printf("\nThe final resultant is not a constant.\n");
        printf("This might indicate remaining variables or an unexpected result.\n");
    }
    
    // Cleanup
    free(poly_string);
    free(r1);
    free(r2);
    free(r3);
    free(s1);
    free(s2);
    free(s3);
    free(s4);
    free(d);
    
    fq_nmod_ctx_clear(ctx);
    
    printf("\n=== Computation Complete ===\n");
    
    return 0;
}