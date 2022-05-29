#ifndef _QRBATCH_
#define _QRBATCH_

void qr_batch(double* a, double* q, int n);
void restore_q_batch(double* q, double* restored_q, int n);

#endif
