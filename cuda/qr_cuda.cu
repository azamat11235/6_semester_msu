#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <stdio.h>
#include <cuda_runtime.h>
#include <cuda.h>
#include <ctime>

#define BLOCK_SIZE 64

double* matrix_transpose(double* matrix, int n);
void mat_mul(double* a, double* b, double* c, int n);


void restore_q(double* q, double* restored_q, int n) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            restored_q[n * i + j] = (double)(i == j);
     for (int jb = n - 1; jb >= 0; jb -= BLOCK_SIZE) {
        for (int ib = n - 1; ib >= 0; ib -= BLOCK_SIZE) {
            for (int j = 0; j < BLOCK_SIZE; ++j) {
                for (int i = 0; i < BLOCK_SIZE; ++i) {
                    int row_abs = ib - i;
                    int col_abs = jb - j;
                    if (row_abs > col_abs) {
                        double c = q[row_abs*n + col_abs];
                        double s = -q[col_abs*n + row_abs];
                        for (int k = 0; k < n; ++k) {
                            double q0k =  restored_q[col_abs*n + k] * c + restored_q[row_abs*n + k] * s;
                            double q1k = -restored_q[col_abs*n + k] * s + restored_q[row_abs*n + k] * c;
                            restored_q[col_abs*n + k] = q0k;
                            restored_q[row_abs*n + k] = q1k;
                        }
                    }
                }
            }
        }
    }
}

bool check_result(double* A, double *q, double *R, int size) {
    const double eps = 1e-10;
    double *Q = new double[size*size];
    restore_q(q, Q, size);
    bool OK = true;

    for (int i = 0; i < size && OK; ++i) {
        for (int j = 0; j < i && OK; ++j) {
            if (R[i*size + j] > eps) {
                // printf("R is not upper triangular\n");
                OK = false;
            }
        }
    }

    double *Qt = new double[size*size];
    memcpy(Qt, Q, sizeof(double)*size*size);
    matrix_transpose(Qt, size);

    double *QtQ = new double[size*size];
    mat_mul(Qt, Q, QtQ, size);

    for (int i = 0; i < size && OK; ++i) {
        for (int j = 0; j < i && OK; ++j) {
            if (std::abs((QtQ[i*size + j] - (i==j))) > eps) {
                // printf("Q^T*Q != I\n")
                OK = false;
            }
        }
    }

    double *QR = new double[size*size];
    mat_mul(Q, R, QR, size);

    for (int i = 0; i < size && OK; ++i) {
        for (int j = 0; j < i && OK; ++j) {
            if (std::abs(QR[i*size + j] - A[i*size + j]) > eps) {
                // printf("Q^T*Q != I\n")
                OK = false;
            }
        }
    }

    delete[] Q;
    delete[] Qt; 
    delete[] QtQ;
    delete[] QR;

    return !OK;
}

void fillMatrix(double* matrix, int n) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            matrix[n * i + j] = (double)(rand() % 1000) / 100;
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


int sign(double x) {
    return (x > 0) - (x < 0);
}

bool isUpperTriangular(double *a, int size) {
    double eps = 1e-10;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < i; ++j) {
            if (a[i*size + j] > eps)
                return false;
        }
    }
    return true;
}

void pm(double *a, int m=8, int n=8) {
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%8.4f ", a[i*n + j]);
        }
        std::cout << "\n";
    }
    std::cout << "-------------------\n";
}

__global__
void to_d_buf(double *d_a, double *d_buf, int i0, int j0, int size_d_a) {
    int i = i0 + BLOCK_SIZE*blockIdx.y  + threadIdx.y;
    int j = j0 + threadIdx.x;
    d_buf[(i-i0)*BLOCK_SIZE + (j-j0)] = d_a[i*size_d_a + j];
}


__global__
void upd(double *d_a, int size_a, double *d_cos, double *d_sin, int j) {
    int my_col = (j+BLOCK_SIZE) + blockIdx.x*BLOCK_SIZE + threadIdx.x;
    for (int i = j; i < size_a; i += BLOCK_SIZE) {
        for (int jj = 0; jj < BLOCK_SIZE; ++jj) {
            for (int ii = 0; ii < BLOCK_SIZE; ++ii) {
                int col = j + jj;
                int row = i + ii;
                if (row > col) {
                    double c =  d_cos[(row-j)*BLOCK_SIZE + jj];
                    double s = -d_sin[(row-j)*BLOCK_SIZE + jj];
                    // std::cout << c << ' ' << s << ' ' << row << ' ' << col << "\n";
                    double a0k =  d_a[col*size_a + my_col] * c + d_a[row*size_a + my_col] * s;
                    double a1k = -d_a[col*size_a + my_col] * s + d_a[row*size_a + my_col] * c;
                    d_a[col*size_a + my_col] = a0k;
                    d_a[row*size_a + my_col] = a1k;
                }
            }
        }
    }
}

void qr_cuda(double *a, double *q, int size) {
    double *buf = new double[size*BLOCK_SIZE * 3]; // a, q_cos, q_sin
    double *d_a;
    double *d_buf;
    double *d_cos;
    double *d_sin;
    cudaMalloc(&d_a, sizeof(double)*size*size);
    cudaMalloc(&d_cos, sizeof(double)*size*BLOCK_SIZE);
    cudaMalloc(&d_sin, sizeof(double)*size*BLOCK_SIZE);
    cudaMalloc(&d_buf, sizeof(double)*size*BLOCK_SIZE);
    cudaMemcpy(d_a, a, sizeof(double)*size*size, cudaMemcpyHostToDevice);

    cudaEvent_t syncEvent;
    cudaEventCreate(&syncEvent);

    for (int j = 0; j < size; j += BLOCK_SIZE) {
        for (int i = j; i < size; i += BLOCK_SIZE) {
            for (int jj = 0; jj < BLOCK_SIZE; ++jj) {
                int col = j + jj;
                for (int ii = 0; ii < BLOCK_SIZE; ++ii) {
                    int row = i + ii;
                    if (row > col) {
                        double a0 = a[col*size + col];
                        double a1 = a[row*size + col];
                        double sigma = sign(a1);
                        if (std::abs(a0) > std::abs(a1))
                            sigma = sign(a0);
                        double den = sigma * std::sqrt(a0*a0 + a1*a1);
                        double c = a0 / den;
                        double s = a1 / den;
                        // q[row*size + col] =  c;
                        // q[col*size + row] = -s;
                        int c_ind = (size*BLOCK_SIZE) + (row-j)*BLOCK_SIZE + (col-j);
                        int s_ind = c_ind + size*BLOCK_SIZE;
                        buf[c_ind] =  c;
                        buf[s_ind] = -s;
                        for (int k = col; k < j+BLOCK_SIZE; ++k) {
                            double a0k =  a[col*size + k] * c + a[row*size + k] * s;
                            double a1k = -a[col*size + k] * s + a[row*size + k] * c;
                            a[col*size + k] = a0k;
                            a[row*size + k] = a1k;
                        }
                    }
                }
            }
        }
        // cos sin to device
        cudaMemcpy(d_cos, &buf[size*BLOCK_SIZE], sizeof(double)*size*BLOCK_SIZE, cudaMemcpyHostToDevice);
        cudaMemcpy(d_sin, &buf[2*size*BLOCK_SIZE], sizeof(double)*size*BLOCK_SIZE, cudaMemcpyHostToDevice);

        dim3 gridSize = dim3((size-j-BLOCK_SIZE)/BLOCK_SIZE, 1, 1);
        dim3 blockSize = dim3(BLOCK_SIZE, 1, 1);
        cudaEventRecord(syncEvent, 0);
        upd<<<gridSize, blockSize>>>(d_a, size, d_cos, d_sin, j);

        for (int i = j; i < size; i += BLOCK_SIZE) { // cos&sin из buf в q
            for (int jj = 0; jj < BLOCK_SIZE; ++jj) {
                int col = j + jj;
                for (int ii = 0; ii < BLOCK_SIZE; ++ii) {
                    int row = i + ii;
                    if (row > col) {
                        int c_ind = (size*BLOCK_SIZE) + (row-j)*BLOCK_SIZE + (col-j);
                        int s_ind = c_ind + size*BLOCK_SIZE;
                        double c = buf[c_ind];
                        double s = buf[s_ind];
                        q[row*size + col] = c;
                        q[col*size + row] = s;
                    }
                }
            }
        }

        cudaEventSynchronize(syncEvent);

        // один столбцовый блок из device в host
        gridSize = dim3(1, size/BLOCK_SIZE, 1);
        blockSize = dim3(BLOCK_SIZE, BLOCK_SIZE, 1);
        cudaEventRecord(syncEvent, 0);
        to_d_buf<<<gridSize, blockSize>>>(d_a, d_buf, 0, j+BLOCK_SIZE, size);
        cudaEventSynchronize(syncEvent);
        cudaMemcpy(buf, d_buf, sizeof(double)*size*BLOCK_SIZE, cudaMemcpyDeviceToHost);
        for (int ii = 0; ii < size; ++ii) {
            for (int jj = j+BLOCK_SIZE; jj < j + 2*BLOCK_SIZE; ++jj) {
                a[ii*size + jj] = buf[ii*BLOCK_SIZE + (jj-j-BLOCK_SIZE)];
            }
        }
    }
    delete[] buf;
    cudaFree(d_a);
    cudaFree(d_cos);
    cudaFree(d_sin);
    cudaFree(d_buf);
}

int main() {
    printf("%-5s\t%-10s\t%12s\n", "size", "time (s.)", "check_result");
    printf("---------------------------------\n");
    for (int size = 256; size <= 2048; size *= 2) {
        double *a = new double[size*size];
        double *q = new double[size*size];
        fillMatrix(a, size);
        double *r = new double[size*size];
        memcpy(r, a, sizeof(double)*size*size);

        qr_cuda(r, q, size);

        printf("%-5s\t%-10s\t", size, "time (s.)");
        if (!check_result(a, q, r, size)) {
            printf("%12s", "OK\n");
        }
        else {
            printf("%12s", "Error!\n");
        }
        

    }
    
    std::cout << isUpperTriangular(a, size) << "\n";


    return 0;
}
