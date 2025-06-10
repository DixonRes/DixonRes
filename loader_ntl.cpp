#include <fstream>
#include <sstream>
#include <NTL/lzz_pXFactoring.h> 
#include <mat_lzz_pX_determinant.h>
#include <chrono>
#include <exception>
using namespace NTL;

int main() {
    std::ifstream fin("dixon_matrix");
    std::string magic;
    fin >> magic;
    
    long p, m, d;
    fin.ignore(1024, '\n');
    fin >> p >> m >> d;
    
    zz_p::init(p);
    Mat<zz_pX> A;
    A.SetDims(m, m);
    
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        
        if (line.find("P") == 0) {
            std::istringstream iss(line);
            std::string tag;
            long i, j, num_terms;
            iss >> tag >> i >> j >> num_terms;
            
            for (long k = 0; k < num_terms; k++) {
                std::getline(fin, line);
                long exp, coeff;
                sscanf(line.c_str(), "T %ld %ld", &exp, &coeff);
                SetCoeff(A[i][j], exp, coeff);
            }
        }
    }


    
    
    zz_pX det;
    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "determinant_via_linsolve" << std::endl;
    //bool is_true = PML::determinant_generic_knowing_degree(det, A, d);
    //std::cout << is_true << std::endl;
    PML::determinant_via_linsolve(det, A);
    //determinant_via_linsolve //determinant_generic_knowing_degree //determinant_via_evaluation_FFT
    //determinant_via_evaluation_general //determinant_via_evaluation_geometric
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "time: " << duration.count() << " ms" << std::endl;

    //std::cout << det << std::endl;
/*
    Vec<zz_p> roots;
    //roots = FindRoots(det);
    roots.SetLength(0);
    for (long i = 0; i < p; i++) {
        zz_p x;
        conv(x, i); 
        if (eval(det, x) == 0) {
            roots.append(x);
        }
    }
    
    std::cout << "Roots (" << roots.length() << "):" << std::endl;
    for (long i = 0; i < roots.length(); i++) {
        std::cout << "Root #" << i << ": " << roots[i] << std::endl;
    }
*/
    return 0;
}
//g++ loader.cpp -L/home/suohaohai02/mylinks -lpml -lntl -o loader