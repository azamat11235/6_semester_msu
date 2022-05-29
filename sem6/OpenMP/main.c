#include <stdio.h>
#include "debug.h"
#include "qr.h"
#include "qr_omp.h"


int main() {
    if (test(qr_omp) == 0)
        printf("test (qr_omp): OK\n");
    printf("\n");
    printf("%-5s\t%-10s\t%-9s\n", "size", "proc_count", "time (s.)");
    printf("---------------------------------\n");
    for (int size = 256; size <= 2048; size *= 2) {
    	for (int proc_count = 1; proc_count <= 8; proc_count *= 2) {
		printf("%-5d\t%10d\t%-.6f\n", size, proc_count, compute_time(qr_omp, size, proc_count));
	}
	printf("\n");
    }

    return 0;
}
