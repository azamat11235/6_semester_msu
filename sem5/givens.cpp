#include <iostream>
#include <ctime>
#include <cstring>
#include <cmath>

using std::cout;
using std::endl;

void print_matrix(const double* matrix, int size) {
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j)
            cout << matrix[size * i + j] << ' ';
        cout << endl;
    }
}

void fill_matrix(double* matrix, int size) {
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            matrix[size * i + j] = (double)(rand() % 1000) / 100;
}

double* matrix_transpose(double* matrix, int size) {
    double* matrix_T = new double[size * size];
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            matrix_T[size * i + j] = matrix[size * j + i];
        }
    }
    return matrix_T;
}

double* mat_mul(double* a, double* b, int size) {
    double* res = new double[size * size];
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            res[size * i + j] = 0;
            for (int k = 0; k < size; ++k)
                res[size * i + j] += a[size * i + k] * b[size * k + j];
        }
    }
    return res;
}

double* mat_sub(double* a, double* b, int size) {
    double* c = new double[size * size];
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            c[size * i + j] = a[size * i + j] - b[size * i + j];
    return c;
}

double* get_eye(int size) {
    double* eye = new double[size * size];
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            eye[size * i + j] = (i == j);
    return eye;
}

bool are_equal(double* a, double* b, int size) {
    const double EPS = 1e-10;
    double norm_fro = 0;
    for (int i = 0; i < size * size; ++i)
        norm_fro += std::pow((a[i] - b[i]), 2);
    norm_fro = std::sqrt(norm_fro);
    if (norm_fro > EPS)
        return false;
    else
        return true;
}

void compute_params(double aii, double aji, double* cos, double* sin) {
    *cos = aii / std::sqrt(aii * aii + aji * aji);
    *sin = -aji / std::sqrt(aii * aii + aji * aji);
}

void rotate(double* xi, double* xj, double cos, double sin) {
    double xi_ = (*xi) * cos - (*xj) * sin;
    double xj_ = (*xi) * sin + (*xj) * cos;
    *xi = xi_;
    *xj = xj_;
}

double* qr(double* a, int size) {
    int q_size = size * (size - 1);
    double* q = new double[q_size];
    double* iq = q;
    for (int icolumn = 0; icolumn < size - 1; ++icolumn) {
        for (int irow = icolumn + 1; irow < size; ++irow) {
            double cos;
            double sin;
            compute_params(a[size * icolumn + icolumn], a[size * irow + icolumn], &cos, &sin);
            *iq++ = cos;
            *iq++ = -sin;
            for (int k = icolumn; k < size; ++k)
                rotate(&a[size * icolumn + k], &a[size * irow + k], cos, sin);
        }
    }
    return q;
}

double* qx_mul(double* q, double* x, int size) {
    double* qx = new double[size * size];
    memcpy(qx, x, sizeof(x) * size * size);
    int q_size = size * (size - 1);
    q += q_size - 1;
    for (int icolumn = size - 1; icolumn >= 0; --icolumn) {
        for (int irow = size - 1; irow >= icolumn + 1; --irow) {
            double sin = *q--;
            double cos = *q--;
            for (int k = 0; k < size; ++k)
                rotate(&qx[size * icolumn + k], &qx[size * irow + k], cos, sin);
        }
    }
    return qx;
}

bool test(int size = 100, int num_of_tests = 5) {
    double* matrix = new double[size * size];
    double* r = new double[size * size];
    for (int i = 0; i < num_of_tests; ++i) {
        fill_matrix(matrix, size);
        memcpy(r, matrix, sizeof(matrix) * size * size);
        double* q = qr(r, size);
        double* qr_ = qx_mul(q, r, size);
        if (!are_equal(matrix, qr_, size)) {
            cout << "Warning: A != QR." << endl;
            return false;
        }
        double* eye = get_eye(size);
        double* Q = qx_mul(q, eye, size);
        double* Q_T = matrix_transpose(Q, size);
        double* QQ_T = mat_mul(Q, Q_T, size);
        if (!are_equal(QQ_T, eye, size)) {
            cout << "Warning: Q is not an orthogonal matrix." << endl;
            return false;
        }
        for (int icolumn = 0; icolumn < size; ++icolumn)
            for (int irow = icolumn + 1; irow < size; ++irow)
                if (r[size * irow + icolumn] > 1e-10) {
                    cout << "Warning: R is not an upper triangular." << endl;
                    return false;
                }
        delete[] q;
        delete[] qr_;
        delete[] eye;
        delete[] Q;
        delete[] Q_T;
        delete[] QQ_T;
    }
    delete[] matrix;
    delete[] r;
    return true;
}

double get_time(int size, int num_of_iters = 10) {
    double t = 0;
    double* matrix = new double[size * size];
    double* tmp = new double[size * size];
    fill_matrix(matrix, size);
    for (int i = 0; i < num_of_iters; ++i) {
        clock_t start;
        clock_t end;
        memcpy(tmp, matrix, sizeof(matrix) * size * size);
        start = clock();
        double* q = qr(tmp, size);
        end = clock();
        t += end - start;
        delete[] q;
    }
    delete[] matrix;
    delete[] tmp;
    return t / num_of_iters / CLOCKS_PER_SEC;
}

int main() {
    srand((unsigned)time(NULL));

    if (test())
        cout << "All tests passed!" << endl;
    cout << get_time(256) << endl;
    cout << get_time(512) << endl;
    cout << get_time(1024) << endl;
    cout << get_time(2048) << endl;

    return 0;
}