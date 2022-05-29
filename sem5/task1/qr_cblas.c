#include <cblas.h>
#include "qr_cblas.h"


void qr_cblas(double* a, double* q, int n) {
    for (int j = 0; j < n - 1; ++j) {
      for (int i = j + 1; i < n; ++i) {
          double c;
          double s;
          cblas_drotg(&a[n * j + j], &a[n * i + j], &c, &s);
          a[n * i + j] = 0;
          cblas_drot(n-j-1, &a[j*n + j+1], 1, &a[i*n + j+1], 1, c, s);
          q[n*i + j] = c;
          q[n*j + i] = -s;
      }
    }
}
