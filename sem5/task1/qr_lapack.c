#include <lapacke.h>
#include "qr_lapack.h"


void qr_lapack(double* a, double* tau, int n) {
  LAPACKE_dgeqrfp(LAPACK_ROW_MAJOR, n, n, a, n, tau);
}
