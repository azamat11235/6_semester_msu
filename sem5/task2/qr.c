#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "qr.h"

void qr(double* a, int ma, int na, double* q, int nq, double* r, int mr, int nr) {
    double* r_ = malloc(sizeof(double)*ma*na);
    memcpy(r_, a, sizeof(double)*ma*na);
    for (int icol = 0; icol < na; ++icol) {
        for (int irow = icol + 1; irow < ma; ++irow) {
            double c;
            double s;
            double rii = r_[na * icol + icol];
            double rji = r_[na * irow + icol];
            c = rii / sqrt(rii * rii + rji * rji);
            s = -rji / sqrt(rii * rii + rji * rji);
            *q++ = c;
            *q++ = -s;
            for (int k = icol; k < na; ++k) {
                double xi = r_[na * icol + k];
                double xj = r_[na * irow + k];
                r_[na * icol + k] = xi * c - xj * s;
                r_[na * irow + k] = xi * s + xj * c;
            }
        }
    }
    for (int i = 0; i < mr; ++i)
        for (int j = 0; j < nr; ++j)
           r[i*nr + j] = r_[i*nr + j];
    free(r_)
}
