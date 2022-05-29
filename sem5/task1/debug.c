#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <cblas.h>
#include "debug.h"
#include "routine.h"
#include "qr.h"
#include "qr_cblas.h"
#include "qr_batch.h"

typedef void qr_func_type(double*, double*, int);
typedef void restore_q_func_type(double*, double*, int);


double compute_time(qr_func_type* qr_func, int size) {
    const int num_of_iters = 5;
    double t = 0;
    double* matrix;
    double* q;
    allocMatrix(&q, size);
    allocMatrix(&matrix, size);
    for (int i = 0; i < num_of_iters; ++i) {
        clock_t start;
        clock_t end;
        fillMatrix(matrix, size);
        start = clock();
        (*qr_func)(matrix, q, size);
        end = clock();
        t += end - start;
    }
    free(matrix);
    free(q);
    return t / num_of_iters / CLOCKS_PER_SEC;
}

int test(qr_func_type* qr_func, restore_q_func_type* restore_q_func)  {
    int flag = 0;
    const int size = 256;
    const double EPS = 1e-10;
    double* a;
    double* q_;
    double* r;
    double* q;
    double* qr;
    double* qt;
    double* qqt;
    double* eye;
    allocMatrix(&a, size);
    allocMatrix(&q_, size);
    allocMatrix(&r, size);
    allocMatrix(&q, size);
    allocMatrix(&qr, size);
    allocMatrix(&qt, size);
    allocMatrix(&qqt, size);
    allocMatrix(&eye, size);
    
    get_eye(eye, size);
    fillMatrix(a, size);
    memcpy(r, a, sizeof(double)*size*size);
    (*qr_func)(r, q_, size);
    restore_q_func(q_, q, size);
    for (int icol = 0; icol < size && flag == 0; ++icol)
        for (int irow = icol + 1; irow < size && flag == 0; ++irow)
            if (abs(r[size * irow + icol]) > 1e-5) {
                printf("Warning: R is not an upper triangular.\n");
                flag = 1;
            }
    mat_mul(q, r, qr, size);
    if (!are_equal(a, qr, size, EPS)) {
        printf("Warning: A != QR.\n");
        flag = 1;
    }
    memcpy(qt, q, sizeof(double) * size * size);
    matrix_transpose(qt, size);
    mat_mul(qt, q, qqt, size);
    if (!are_equal(qqt, eye, size, EPS)) {
        printf("Warning: Q is not an orthogonal matrix.\n");
        flag = 1;
    }

    free(a);
    free(q_);
    free(r);
    free(q);
    free(qr);
    free(qt);
    free(qqt);
    free(eye);

    return flag;
}
