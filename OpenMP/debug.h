#ifndef _DEBUG_
#define _DEBUG_

typedef void qr_func_type(double*, double*, int, int);
double compute_time(qr_func_type* qr_func, int size, int proc_count);
int test(qr_func_type* f);

#endif
