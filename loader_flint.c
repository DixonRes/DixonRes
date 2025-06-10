#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/perm.h>
#include <flint/ulong_extras.h>
#include <nmod_poly_mat_extra.h>
#include <nmod_poly_mat_utils.h>
//gcc flint.c -L/home/suohaohai02/mylinks -lflint -lstdc++ -lpml2 -o loader2
int main() {
    FILE* fin = fopen("dixon_matrix", "r");
    if (!fin) {
        fprintf(stderr, "Error opening file\n");
        return 1;
    }
    
    char magic[16];
    fscanf(fin, "%15s", magic);
    
    long p, m, d;
    fscanf(fin, "%ld %ld %ld", &p, &m, &d);
    
    nmod_poly_mat_t mat;
    nmod_poly_mat_init(mat, m, m, p);
    
    for (long i = 0; i < m; i++) {
        for (long j = 0; j < m; j++) {
            nmod_poly_init(nmod_poly_mat_entry(mat, i, j), p);
        }
    }
    
    char line[1024];
    while (fgets(line, sizeof(line), fin)) {
        if (line[0] == '\0' || line[0] == '\n') continue;
        
        if (line[0] == 'P') {
            long i, j, num_terms;
            sscanf(line, "P %ld %ld %ld", &i, &j, &num_terms);
            
            nmod_poly_struct* poly = nmod_poly_mat_entry(mat, i, j);
            
            for (long k = 0; k < num_terms; k++) {
                if (!fgets(line, sizeof(line), fin)) break;
                
                long exp, coeff;
                if (sscanf(line, "T %ld %ld", &exp, &coeff) == 2) {
                    nmod_poly_set_coeff_ui(poly, exp, coeff % p);
                }
            }
        }
    }
    fclose(fin);
    
    nmod_poly_t det;
    nmod_poly_init(det, p);
    
    clock_t start = clock();
    printf("determinant_via_iter\n");
    
    nmod_poly_mat_det_iter(det, mat);
    
    printf("1\n");
    printf("Determinant: ");
    //nmod_poly_print_pretty(det, "x");
    printf("\n");
    clock_t end = clock();
    double duration = (double)(end - start) * 1000.0 / CLOCKS_PER_SEC;
    printf("time: %.0f ms\n", duration);

    nmod_poly_mat_clear(mat);
    nmod_poly_clear(det);
    
    return 0;
}
