#ifndef _ROUTINE_
#define _ROUTINE_

void bcache(double* a, int na, double* cache, int i, int j, int k); //
void bflush(double* a, int na, double* cache, int i, int j, int k); //
void allocMatrix(double** matrix, int n);
void fillMatrix(double* matrix, int n);
int are_equal(double* a, double* b, int n, double eps);
double* get_eye(double* matrix, int size);
double* matrix_transpose(double* matrix, int n);
void mat_mul(double* a, double* b, double* c, int n);

#endif
