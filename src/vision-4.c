// dixon_solver.c
// Complete implementation for solving the polynomial system from Magma
// Compile: gcc -O3 -march=native -fopenmp -o dixon_solver dixon_solver.c -lflint -lmpfr -lgmp -lpthread
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

// Include the Dixon interface headers (assumed to be available)
#include "dixon_flint.h"
#include "dixon_interface_flint.h"

// Helper function to setup GF(2^8) with the irreducible polynomial from Magma
void setup_gf28_context(fq_nmod_ctx_t ctx) {
    // Irreducible polynomial: z8^8 + z8^4 + z8^3 + z8^2 + 1
    nmod_poly_t modulus;
    nmod_poly_init(modulus, 2);
    
    // Set coefficients for z8^8 + z8^4 + z8^3 + z8^2 + 1
    nmod_poly_set_coeff_ui(modulus, 0, 1);  // constant term
    nmod_poly_set_coeff_ui(modulus, 2, 1);  // z8^2
    nmod_poly_set_coeff_ui(modulus, 3, 1);  // z8^3
    nmod_poly_set_coeff_ui(modulus, 4, 1);  // z8^4
    nmod_poly_set_coeff_ui(modulus, 8, 1);  // z8^8
    
    fq_nmod_ctx_init_modulus(ctx, modulus, "z8");
    nmod_poly_clear(modulus);
}

// Define all 15 polynomials from the Magma script
// These are the exact polynomials from ps[1] through ps[15]

const char* POLY_P1 = "(z8^4 + z8 + 1)*x1^4*x3^4 + (z8^7 + z8^4 + z8^2)*x2^4*x3^4 + "
    "(z8^7 + z8^5 + z8^3 + z8^2 + z8)*x1^4*x3^2 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^3 + 1)*x2^4*x3^2 + "
    "(z8^7 + z8^5 + z8^3 + z8^2 + z8)*x1^2*x3^4 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^3 + 1)*x2^2*x3^4 + "
    "(z8^6 + z8^5 + z8^3 + z8^2 + z8)*x1^4*x3 + (z8^6 + z8^4 + z8^3 + 1)*x2^4*x3 + "
    "(z8^6 + z8^5 + z8^3 + z8^2 + z8)*x1*x3^4 + (z8^6 + z8^4 + z8^3 + 1)*x2*x3^4 + "
    "(z8^7 + z8^6 + z8^3 + z8^2 + z8)*x1^2*x3^2 + (z8^7 + z8^5 + z8^3 + 1)*x2^2*x3^2 + "
    "(z8^7 + z8^6 + z8^4 + z8^3 + 1)*x3^4 + (z8^6 + z8^5 + z8^3)*x1^2*x3 + "
    "(z8^6 + z8^4 + z8^3 + z8^2)*x2^2*x3 + (z8^6 + z8^5 + z8^3)*x1*x3^2 + "
    "(z8^6 + z8^4 + z8^3 + z8^2)*x2*x3^2 + (z8^7 + z8^3)*x1*x3 + "
    "(z8^7 + z8^6 + z8^3 + z8^2)*x2*x3 + (z8^7 + z8^4 + z8^3)*x3^2 + "
    "(z8^7 + z8^2 + z8 + 1)*x3 + 1";

const char* POLY_P2 = "(z8^5 + z8^4 + z8^2 + 1)*x1^4*x4^4 + (z8^7 + z8^5 + z8^4 + z8)*x2^4*x4^4 + "
    "(z8^7 + z8^6 + z8^5 + z8^3 + z8^2 + z8 + 1)*x1^4*x4^2 + (z8^7 + z8^5 + z8^4 + z8^3)*x2^4*x4^2 + "
    "(z8^7 + z8^6 + z8^5 + z8^3 + z8^2 + z8 + 1)*x1^2*x4^4 + (z8^7 + z8^5 + z8^4 + z8^3)*x2^2*x4^4 + "
    "(z8^7 + z8^5 + z8^4 + z8)*x1^4*x4 + (z8^7 + z8^2 + 1)*x2^4*x4 + "
    "(z8^7 + z8^5 + z8^4 + z8)*x1*x4^4 + (z8^7 + z8^2 + 1)*x2*x4^4 + "
    "(z8^6 + z8^3 + z8^2 + z8 + 1)*x1^2*x4^2 + (z8^5 + z8^3)*x2^2*x4^2 + "
    "(z8^7 + z8^6 + z8^5 + z8^3)*x4^4 + (z8^7 + z8^5 + z8^4 + z8^3)*x1^2*x4 + "
    "(z8^7 + z8^3 + z8^2)*x2^2*x4 + (z8^7 + z8^5 + z8^4 + z8^3)*x1*x4^2 + "
    "(z8^7 + z8^3 + z8^2)*x2*x4^2 + (z8^7 + z8^2 + 1)*x1*x4 + "
    "(z8^7 + z8^6 + 1)*x2*x4 + (z8^7 + z8^6 + z8^4 + z8^3 + z8^2)*x4^2 + "
    "(z8^7 + z8^6 + z8^2)*x4 + 1";

const char* POLY_P3 = "(z8^4 + z8 + 1)*x5^4*x7^4 + (z8^7 + z8^4 + z8^2)*x6^4*x7^4 + "
    "(z8^7 + z8^5 + z8^3 + z8^2 + z8)*x5^4*x7^2 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^3 + 1)*x6^4*x7^2 + "
    "(z8^7 + z8^5 + z8^3 + z8^2 + z8)*x5^2*x7^4 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^3 + 1)*x6^2*x7^4 + "
    "(z8^6 + z8^5 + z8^3 + z8^2 + z8)*x5^4*x7 + (z8^6 + z8^4 + z8^3 + 1)*x6^4*x7 + "
    "(z8^6 + z8^5 + z8^3 + z8^2 + z8)*x5*x7^4 + (z8^6 + z8^4 + z8^3 + 1)*x6*x7^4 + "
    "(z8^7 + z8^6 + z8^3 + z8^2 + z8)*x5^2*x7^2 + (z8^7 + z8^5 + z8^3 + 1)*x6^2*x7^2 + "
    "(z8^7 + z8^6 + z8^5 + z8^4 + 1)*x7^4 + (z8^6 + z8^5 + z8^3)*x5^2*x7 + "
    "(z8^6 + z8^4 + z8^3 + z8^2)*x6^2*x7 + (z8^6 + z8^5 + z8^3)*x5*x7^2 + "
    "(z8^6 + z8^4 + z8^3 + z8^2)*x6*x7^2 + (z8^7 + z8^3)*x5*x7 + "
    "(z8^7 + z8^6 + z8^3 + z8^2)*x6*x7 + (z8^7 + z8^6 + z8^5 + z8)*x7^2 + "
    "(z8^7 + z8^5 + z8 + 1)*x7 + 1";

const char* POLY_P4 = "(z8^5 + z8^4 + z8^2 + 1)*x5^4*x8^4 + (z8^7 + z8^5 + z8^4 + z8)*x6^4*x8^4 + "
    "(z8^7 + z8^6 + z8^5 + z8^3 + z8^2 + z8 + 1)*x5^4*x8^2 + (z8^7 + z8^5 + z8^4 + z8^3)*x6^4*x8^2 + "
    "(z8^7 + z8^6 + z8^5 + z8^3 + z8^2 + z8 + 1)*x5^2*x8^4 + (z8^7 + z8^5 + z8^4 + z8^3)*x6^2*x8^4 + "
    "(z8^7 + z8^5 + z8^4 + z8)*x5^4*x8 + (z8^7 + z8^2 + 1)*x6^4*x8 + "
    "(z8^7 + z8^5 + z8^4 + z8)*x5*x8^4 + (z8^7 + z8^2 + 1)*x6*x8^4 + "
    "(z8^6 + z8^3 + z8^2 + z8 + 1)*x5^2*x8^2 + (z8^5 + z8^3)*x6^2*x8^2 + "
    "(z8^6 + z8^5 + z8^4 + z8^3 + z8^2 + z8)*x8^4 + (z8^7 + z8^5 + z8^4 + z8^3)*x5^2*x8 + "
    "(z8^7 + z8^3 + z8^2)*x6^2*x8 + (z8^7 + z8^5 + z8^4 + z8^3)*x5*x8^2 + "
    "(z8^7 + z8^3 + z8^2)*x6*x8^2 + (z8^7 + z8^2 + 1)*x5*x8 + "
    "(z8^7 + z8^6 + 1)*x6*x8 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^3 + z8^2 + z8)*x8^2 + "
    "(z8^6 + z8^5 + z8^3 + z8 + 1)*x8 + 1";

const char* POLY_P5 = "(z8^4 + z8 + 1)*x9^4*x11^4 + (z8^7 + z8^4 + z8^2)*x10^4*x11^4 + "
    "(z8^7 + z8^5 + z8^3 + z8^2 + z8)*x9^4*x11^2 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^3 + 1)*x10^4*x11^2 + "
    "(z8^7 + z8^5 + z8^3 + z8^2 + z8)*x9^2*x11^4 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^3 + 1)*x10^2*x11^4 + "
    "(z8^6 + z8^5 + z8^3 + z8^2 + z8)*x9^4*x11 + (z8^6 + z8^4 + z8^3 + 1)*x10^4*x11 + "
    "(z8^6 + z8^5 + z8^3 + z8^2 + z8)*x9*x11^4 + (z8^6 + z8^4 + z8^3 + 1)*x10*x11^4 + "
    "(z8^7 + z8^6 + z8^3 + z8^2 + z8)*x9^2*x11^2 + (z8^7 + z8^5 + z8^3 + 1)*x10^2*x11^2 + "
    "(z8^7 + z8^5 + z8 + 1)*x11^4 + (z8^6 + z8^5 + z8^3)*x9^2*x11 + "
    "(z8^6 + z8^4 + z8^3 + z8^2)*x10^2*x11 + (z8^6 + z8^5 + z8^3)*x9*x11^2 + "
    "(z8^6 + z8^4 + z8^3 + z8^2)*x10*x11^2 + (z8^7 + z8^3)*x9*x11 + "
    "(z8^7 + z8^6 + z8^3 + z8^2)*x10*x11 + (z8^7 + z8^6 + z8^3 + z8^2 + 1)*x11^2 + "
    "(z8^4 + z8^3 + z8^2 + 1)*x11 + 1";

const char* POLY_P6 = "(z8^5 + z8^4 + z8^2 + 1)*x9^4*x12^4 + (z8^7 + z8^5 + z8^4 + z8)*x10^4*x12^4 + "
    "(z8^7 + z8^6 + z8^5 + z8^3 + z8^2 + z8 + 1)*x9^4*x12^2 + (z8^7 + z8^5 + z8^4 + z8^3)*x10^4*x12^2 + "
    "(z8^7 + z8^6 + z8^5 + z8^3 + z8^2 + z8 + 1)*x9^2*x12^4 + (z8^7 + z8^5 + z8^4 + z8^3)*x10^2*x12^4 + "
    "(z8^7 + z8^5 + z8^4 + z8)*x9^4*x12 + (z8^7 + z8^2 + 1)*x10^4*x12 + "
    "(z8^7 + z8^5 + z8^4 + z8)*x9*x12^4 + (z8^7 + z8^2 + 1)*x10*x12^4 + "
    "(z8^6 + z8^3 + z8^2 + z8 + 1)*x9^2*x12^2 + (z8^5 + z8^3)*x10^2*x12^2 + "
    "(z8^7 + z8^6 + z8^4 + z8^2 + 1)*x12^4 + (z8^7 + z8^5 + z8^4 + z8^3)*x9^2*x12 + "
    "(z8^7 + z8^3 + z8^2)*x10^2*x12 + (z8^7 + z8^5 + z8^4 + z8^3)*x9*x12^2 + "
    "(z8^7 + z8^3 + z8^2)*x10*x12^2 + (z8^7 + z8^2 + 1)*x9*x12 + "
    "(z8^7 + z8^6 + 1)*x10*x12 + (z8^6 + z8^5 + z8^4 + z8^3)*x12^2 + "
    "(z8^7 + z8^3 + 1)*x12 + 1";

const char* POLY_P7 = "z8*x0*x1 + (z8^7 + z8^6 + z8^4 + z8^3 + z8 + 1)*x1 + 1";

const char* POLY_P8 = "(z8^2 + z8)*x0*x2 + (z8^7 + z8^6 + z8^4 + z8^3 + z8 + 1)*x2 + 1";

const char* POLY_P9 = "z8*x3*x5 + (z8 + 1)*x4*x5 + (z8^7 + z8^6 + z8^4 + z8^2 + 1)*x5 + 1";

const char* POLY_P10 = "(z8^2 + z8)*x3*x6 + (z8^2 + z8 + 1)*x4*x6 + (z8^7 + z8^3 + z8^2 + z8 + 1)*x6 + 1";

const char* POLY_P11 = "z8*x7*x9 + (z8 + 1)*x8*x9 + (z8^7 + z8^2 + 1)*x9 + 1";

const char* POLY_P12 = "(z8^2 + z8)*x7*x10 + (z8^2 + z8 + 1)*x8*x10 + (z8^7 + z8^5 + z8^3 + z8^2)*x10 + 1";

const char* POLY_P13 = "z8*x11*x13 + (z8 + 1)*x12*x13 + (z8^6 + z8^5 + z8^3)*x13 + 1";

const char* POLY_P14 = "(z8^2 + z8)*x11*x14 + (z8^2 + z8 + 1)*x12*x14 + (z8^7 + z8^6 + z8^3)*x14 + 1";

const char* POLY_P15 = "(z8^7 + z8^5 + z8^3 + 1)*x13^4 + (z8^6 + z8^5 + z8^4 + z8 + 1)*x14^4 + "
    "(z8^6 + z8^4 + z8^3 + z8^2 + 1)*x13^2 + (z8^7 + z8^6 + z8^5 + z8^4 + z8^3 + z8^2 + 1)*x14^2 + "
    "(z8^4 + z8^2)*x13 + (z8^4 + z8^3 + z8^2 + z8)*x14 + z8^4";

// Main implementation following the Magma script structure
int main() {
    printf("=== Dixon Resultant Solver for Magma Polynomial System ===\n\n");
    
    // Setup GF(2^8) context
    fq_nmod_ctx_t ctx;
    setup_gf28_context(ctx);
    printf("Field: GF(2^8) with modulus z8^8 + z8^4 + z8^3 + z8^2 + 1\n");
    printf("Computing elimination sequence from Magma script...\n\n");
    
    // Allocate large buffer for intermediate results
    size_t buffer_size = 1024 * 1024 * 50; // 50MB for very large polynomials
    char* poly_string = (char*) malloc(buffer_size);
    if (!poly_string) {
        printf("Memory allocation failed!\n");
        return 1;
    }
    
    clock_t total_start = clock();
    
    // ===== First Chain: Eliminate x0, x1, x2, x3, x4 =====
    printf("====== FIRST CHAIN ======\n");
    
    // Step 1: r1 = Resultant(p7, p8, x0)
    printf("Step 1: Computing r1 = Resultant(p7, p8, x0)\n");
    clock_t start = clock();
    
    snprintf(poly_string, buffer_size, "%s,%s", POLY_P7, POLY_P8);
    char* r1 = dixon_str(poly_string, "x0", ctx);
    
    clock_t end = clock();
    printf("  Time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("  Result length: %zu characters\n\n", strlen(r1));
    
    // Step 2: r2 = dixon([r1, p1, p2], [x1, x2])
    printf("Step 2: Computing r2 = dixon([r1, p1, p2], [x1, x2])\n");
    start = clock();
    
    snprintf(poly_string, buffer_size, "%s,%s,%s", r1, POLY_P1, POLY_P2);
    char* r2 = dixon_str(poly_string, "x1,x2", ctx);
    
    end = clock();
    printf("  Time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("  Result length: %zu characters\n\n", strlen(r2));
    
    // Step 3: r3 = dixon([r2, p9, p10], [x3, x4])
    printf("Step 3: Computing r3 = dixon([r2, p9, p10], [x3, x4])\n");
    start = clock();
    
    snprintf(poly_string, buffer_size, "%s,%s,%s", r2, POLY_P9, POLY_P10);
    char* r3 = dixon_str(poly_string, "x3,x4", ctx);
    
    end = clock();
    printf("  Time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("  Result length: %zu characters\n\n", strlen(r3));
    
    // Step 4: r4 = dixon([r3, p3, p4], [x5, x6])
    printf("Step 4: Computing r4 = dixon([r3, p3, p4], [x5, x6])\n");
    start = clock();
    
    snprintf(poly_string, buffer_size, "%s,%s,%s", r3, POLY_P3, POLY_P4);
    char* r4 = dixon_str(poly_string, "x5,x6", ctx);
    
    end = clock();
    printf("  Time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("  Result length: %zu characters\n\n", strlen(r4));
    
    // ===== Second Chain: Eliminate x14, x13, x12, x11, x10, x9 =====
    printf("====== SECOND CHAIN ======\n");
    
    // Step 5: s1 = Resultant(p15, p14, x14)
    printf("Step 5: Computing s1 = Resultant(p15, p14, x14)\n");
    start = clock();
    
    snprintf(poly_string, buffer_size, "%s,%s", POLY_P15, POLY_P14);
    char* s1 = dixon_str(poly_string, "x14", ctx);
    
    end = clock();
    printf("  Time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("  Result length: %zu characters\n\n", strlen(s1));
    
    // Step 6: s2 = Resultant(s1, p13, x13)
    printf("Step 6: Computing s2 = Resultant(s1, p13, x13)\n");
    start = clock();
    
    snprintf(poly_string, buffer_size, "%s,%s", s1, POLY_P13);
    char* s2 = dixon_str(poly_string, "x13", ctx);
    
    end = clock();
    printf("  Time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("  Result length: %zu characters\n\n", strlen(s2));
    
    // Step 7: s3 = Resultant(s2, p6, x12)
    printf("Step 7: Computing s3 = Resultant(s2, p6, x12)\n");
    start = clock();
    
    snprintf(poly_string, buffer_size, "%s,%s", s2, POLY_P6);
    char* s3 = dixon_str(poly_string, "x12", ctx);
    
    end = clock();
    printf("  Time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("  Result length: %zu characters\n\n", strlen(s3));
    
    // Step 8: s4 = Resultant(s3, p5, x11)
    printf("Step 8: Computing s4 = Resultant(s3, p5, x11)\n");
    start = clock();
    
    snprintf(poly_string, buffer_size, "%s,%s", s3, POLY_P5);
    char* s4 = dixon_str(poly_string, "x11", ctx);
    
    end = clock();
    printf("  Time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("  Result length: %zu characters\n\n", strlen(s4));
    
    // Step 9: s5 = Resultant(s4, p12, x10)
    printf("Step 9: Computing s5 = Resultant(s4, p12, x10)\n");
    start = clock();
    
    snprintf(poly_string, buffer_size, "%s,%s", s4, POLY_P12);
    char* s5 = dixon_str(poly_string, "x10", ctx);
    
    end = clock();
    printf("  Time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("  Result length: %zu characters\n\n", strlen(s5));
    
    // Step 10: s6 = Resultant(s5, p11, x9)
    printf("Step 10: Computing s6 = Resultant(s5, p11, x9)\n");
    start = clock();
    
    concat_polynomials(s5, POLY_P11, &poly_string, &buffer_size);
    // snprintf(poly_string, buffer_size, "%s,%s", s5, POLY_P11);
    char* s6 = dixon_str(poly_string, "x9", ctx);
    
    end = clock();
    printf("  Time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("  Result length: %zu characters\n\n", strlen(s6));
    
    // ===== Final Step: Combine the two chains =====
    printf("====== FINAL STEP ======\n");
    
    // Step 11: d = Resultant(r4, s6, x7)
    printf("Step 11: Computing d = Resultant(r4, s6, x7)\n");
    printf("  This eliminates x7 and should give us a polynomial in x8 only.\n");
    start = clock();

    concat_polynomials(r4, s6, &poly_string, &buffer_size);
    // snprintf(poly_string, buffer_size, "%s,%s", r4, s6);
    char* d = dixon_ori(poly_string, "x7", ctx);
    
    end = clock();
    printf("  Time: %.3f seconds\n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("  Final resultant length: %zu characters\n\n", strlen(d));
    
    clock_t total_end = clock();
    
    // Check if the result is a constant
    if (strcmp(d, "0") == 0) {
        printf("=== RESULT ===\n");
        printf("The final resultant is 0.\n");
        printf("This means the system has a solution in the algebraic closure.\n");
    } else if (strcmp(d, "1") == 0) {
        printf("=== RESULT ===\n");
        printf("The final resultant is 1.\n");
        printf("This means the system has no solutions.\n");
    } else {
        printf("=== RESULT ===\n");
        printf("The final resultant is a polynomial in x8.\n");
        printf("The system has solutions if and only if this polynomial has roots.\n");
        
        // The roots should be found by the dixon_str function
        // which calls find_and_print_roots_of_univariate_resultant
    }
    
    printf("\n=== COMPUTATION SUMMARY ===\n");
    printf("Total computation time: %.3f seconds\n", (double)(total_end - total_start) / CLOCKS_PER_SEC);
    printf("Field: GF(2^8)\n");
    printf("Variables eliminated: x0, x1, x2, x3, x4, x5, x6, x7, x9, x10, x11, x12, x13, x14\n");
    printf("Remaining variable: x8\n");
    
    // Cleanup
    free(poly_string);
    free(r1);
    free(r2);
    free(r3);
    free(r4);
    free(s1);
    free(s2);
    free(s3);
    free(s4);
    free(s5);
    free(s6);
    free(d);
    
    fq_nmod_ctx_clear(ctx);
    
    printf("\n=== Computation Complete ===\n");
    
    return 0;
}