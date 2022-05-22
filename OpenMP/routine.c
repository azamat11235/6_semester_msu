#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "parameters.h"
#include "routine.h"


void bcache(double* a, int na, double* cache, int i, int j, int k) {
    for (int ii = 0; ii < _b; ++ii) {
            int row_abs = i + ii;
            for (int jj = 0; jj < _b; ++jj) {
                int col_abs = j + jj;
                cache[k*_b*_b + ii*_b + jj] = a[row_abs*na + col_abs];
            }
        }
}

void bflush(double* a, int na, double* cache, int i, int j, int k) {
    for (int ii = 0; ii < _b; ++ii) {
            int row_abs = i + ii;
            for (int jj = 0; jj < _b; ++jj) {
                int col_abs = j + jj;
                a[row_abs*na + col_abs] = cache[k*_b*_b + ii*_b + jj];
            }
        }
}

void allocMatrix(double** matrix, int n) {
    *matrix = malloc(sizeof(double)*n*n);
}

void fillMatrix(double* matrix, int n) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            matrix[n * i + j] = (double)(rand() % 1000) / 100;
}

int are_equal(double* a, double* b, int n, double eps) {
    double norm_fro = 0;
    for (int i = 0; i < n * n; ++i)
        norm_fro += pow((a[i] - b[i]), 2);
    norm_fro = sqrt(norm_fro);
    if (norm_fro > eps)
        return 0;
    else
        return 1;
}

double* get_eye(double* matrix, int size) {
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            matrix[size * i + j] = (i == j);
    return matrix;
}

double* matrix_transpose(double* matrix, int n) {
    double d;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            d = matrix[n * i + j];
            matrix[n * i + j] = matrix[n * j + i];
            matrix[n * j + i] = d;
        }
    }
    return matrix;
}

void mat_mul(double* a, double* b, double* c, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            c[n * i + j] = 0;
            for (int k = 0; k < n; ++k)
                c[n * i + j] += a[n * i + k] * b[n * k + j];
        }
    }
}
