#include <omp.h>
int main(void){ return omp_get_max_threads() > 0 ? 0 : 1; }
