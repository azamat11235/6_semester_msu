#include <math.h>
#include <omp.h>
#include "routine.h"
#include "parameters.h"
#include "qr_omp.h"


void compute_params2(double aii, double aji, double* c, double* s) {
    *c = aii / sqrt(aii * aii + aji * aji);
    *s = -aji / sqrt(aii * aii + aji * aji);
}

void rotate2(double* xi, double* xj, double c, double s) {
    double xi_ = (*xi) * c - (*xj) * s;
    double xj_ = (*xi) * s + (*xj) * c;
    *xi = xi_;
    *xj = xj_;
}

void qr_omp(double* a, double* q, int n) {
    double cache[4*_b*_b] = {0};
    for (int jb = 0; jb < n; jb += _b) {
        bcache(a, n, cache, jb, jb, 2); // кешируем диаг. блок
        // вращаем диаг. блок
        for (int j = 0; j < _b-1; ++j) {
            for (int i = j+1; i < _b; ++i) {
                double ajj = cache[2*_b*_b + j*_b + j];
                double aij = cache[2*_b*_b + i*_b + j];
                double c;
                double s;
                compute_params2(ajj, aij, &c, &s);
                cache[i*_b + j] = c;
                cache[j*_b + i] = s;
                for (int k = j; k < _b; ++k) {
                    int jk = 2*_b*_b + j*_b + k;
                    int ik = 2*_b*_b + i*_b + k;
                    rotate2(&cache[jk], &cache[ik], c, s);
                }

            }
        }
        bflush(q, n, cache, jb, jb, 0); // sin, cos
        // вращаем столбец (поддиаг. блоки)
        for (int ib = jb+_b; ib < n; ib += _b) {
            bcache(a, n, cache, ib, jb, 3); // поддиаг. блок
            for (int j = 0; j < _b; ++j) {
                for (int i = 0; i < _b; ++i) {
                    double ajj = cache[2*_b*_b + j*_b + j];
                    double aij = cache[3*_b*_b + i*_b + j];
                    double c;
                    double s;
                    compute_params2(ajj, aij, &c, &s);
                    cache[i*_b + j] = c;
                    cache[_b*_b + j*_b + i] = s;
                    for (int k = j; k < _b; ++k) {
                       int jk = 2*_b*_b + j*_b + k;
                       int ik = 3*_b*_b + i*_b + k;
                       rotate2(&cache[jk], &cache[ik], c, s);
                    }
                }
            }
            bflush(q, n, cache, ib, jb, 0); // cos
            bflush(q, n, cache, jb, ib, 1); // sin
            bflush(a, n, cache, ib, jb, 3); // поддиаг. блок
        }
        bflush(a, n, cache, jb, jb, 2); // диаг. блок

        #pragma omp parallel num_threads(4)
        {
            // обновляем строку (блоки справа от диаг.)
            #pragma omp for private(cache)
            for (int jb2 = jb+_b; jb2 < n; jb2 += _b) {
                bcache(a, n, cache, jb, jb2, 2); // внедиаг. блок
                bcache(q, n, cache, jb, jb, 0);  // cos, sin диаг. блока
                // вращения диаг. блока
                for (int j = 0; j < _b-1; ++j) {
                    for (int i = j+1; i < _b; ++i) {
                        double c;
                        double s;
                        c = cache[i*_b + j];
                        s = cache[j*_b + i];
                        for (int k = 0; k < _b; ++k) {
                            int jk = 2*_b*_b + j*_b + k;
                            int ik = 2*_b*_b + i*_b + k;
                            rotate2(&cache[jk], &cache[ik], c, s);
                        }
                    }
                }
                // вращения поддиаг. блоков
                for (int ib = jb+_b; ib < n; ib += _b) {
                    bcache(a, n, cache, ib, jb2, 3); // внедиаг. блок (нижний)
                    bcache(q, n, cache, ib, jb, 0);  // cos поддиаг. блока
                    bcache(q, n, cache, jb, ib, 1);  // sin поддиаг. блока
                    for (int j = 0; j < _b; ++j) {
                        for (int i = 0; i < _b; ++i) {
                            double c;
                            double s;
                            c = cache[i*_b + j];
                            s = cache[_b*_b + j*_b + i];
                            for (int k = 0; k < _b; ++k) {
                                int jk = 2*_b*_b + j*_b + k;
                                int ik = 3*_b*_b + i*_b + k;
                                rotate2(&cache[jk], &cache[ik], c, s);
                            }
                        }
                    }
                    bflush(a, n, cache, ib, jb2, 3); // внедиаг. блок (нижний)
                }
                bflush(a, n, cache, jb, jb2, 2); // недиаг. блок
            }
        }
    }
}
