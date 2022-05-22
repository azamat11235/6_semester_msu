#include <stdio.h>
#include "debug.h"
#include "qr.h"
#include "qr_omp.h"


#define LINE "--------------------------------\n"

int main() {
    printf(LINE);
    if (test(qr_omp) == 0)
        printf("test (qr_omp): OK\n");
    if (test(qr) == 0)
        printf("test (qr):     OK\n");
    printf(LINE);
    
    printf("> qr_omp:\n");   
    for (int n = 256; n <= 2048; n *= 2)
        printf("    size = %4d, time = %f\n", n, compute_time(qr_omp, n));
    printf(LINE);
    printf("> qr:\n");
    for (int n = 256; n <= 2048; n *= 2)
        printf("    size = %4d, time = %f\n", n, compute_time(qr, n));
    printf(LINE); 

    return 0;
}
