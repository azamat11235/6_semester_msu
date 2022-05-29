#include <cblas.h>
#include <math.h>
#include "qr.h"

    
void qr(double* a, double* q, int n) {
    for (int j = 0; j < n - 1; ++j) {
        for (int i = j + 1; i < n; ++i) {
            double ajj =  a[j*n + j];
            double aij =  a[i*n + j];
            double c = ajj / sqrt(ajj*ajj + aij*aij);
            double s = -aij / sqrt(ajj*ajj + aij*aij);
            q[i*n + j] = c;
            q[j*n + i] = s;
            for (int k = j; k < n; ++k) {
                double ajk = a[j*n + k];
                double aik = a[i*n + k];
                a[j*n + k] = ajk*c - aik*s;
                a[i*n + k] = ajk*s + aik*c;
            }
        }
    }
}

void restore_q(double* q, double* restored_q, int n) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            restored_q[i*n + j] = (double)(i == j);
    for (int j = n-2; j >= 0; --j) {
        for (int i = n-1; i > j; --i) {
            double c = q[i*n + j];
            double s = -q[j*n + i];
            for (int k = j; k < n; ++k) {
                double qjk = restored_q[j*n + k];
                double qik = restored_q[i*n + k];
                restored_q[j*n + k] = qjk*c - qik*s;
                restored_q[i*n + k] = qjk*s + qik*c;
            }
        }
    }
}
