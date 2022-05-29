#include <iostream>
#include <mpi.h>
#include <cblas.h>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "parameters.h"

void save_res(int mpirank, double* A, double* Q, int m, int n) {
    std::string num;
    if (mpirank < 10)
        num = "0" + std::to_string(mpirank);
    else
        num = std::to_string(mpirank);
    std::string fname_a("data/a_part_" + num + ".dat");
    std::string fname_q("data/q_part_" + num + ".dat");
    std::ofstream out_a;
    std::ofstream out_q;
    out_a.open(fname_a);
    out_q.open(fname_q);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            out_a << A[i*n + j] << ' ';
            out_q << Q[i*n + j] << ' ';
        }
        out_a << '\n';
        out_q << '\n';
    }
}

int main(int argc, char **argv) {
    const int size = atoi(argv[1]);
    const int proc_num = atoi(argv[2]);

    int block_count_glob = std::ceil((float)size / _BLOCK_SIZE);
    const std::string fname("data/data.dat");

    MPI_Init(&argc, &argv);

    int mpirank;
    int mpisize;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
    bool mpiroot = (mpirank == 0);

    int m = size;
    int n = size / proc_num;
    int block_count_loc = std::ceil((float)n / _BLOCK_SIZE);
    double *A = new double[m*n];
    double *Q = new double[m*n];
    
    /* read and distr matr */
    double *buf = new double[m*_BLOCK_SIZE];
    if (mpiroot) {
        double *A_glob = new double[size*size];
        std::ifstream file(fname);
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                double d;
                file >> A_glob[i*size + j];
            }
        }

        for (int j_block_glob = 0; j_block_glob < block_count_glob; ++j_block_glob) {
            int cur_proc = j_block_glob % proc_num;
            int j_block_loc = j_block_glob / proc_num;
            for (int i = 0; i < m; ++i) {
                for (int j = 0; j < _BLOCK_SIZE; ++j) {
                    int j_glob = j_block_glob*_BLOCK_SIZE + j;
                    int j_loc = j_block_loc*_BLOCK_SIZE + j;
                    if (cur_proc == 0)
                        A[i*n + j_loc] = A_glob[i*size + j_glob];
                    else
                        buf[i*_BLOCK_SIZE + j] = A_glob[i*size + j_glob];
                }
            }
            if (cur_proc != 0) {
                MPI_Send(buf, m*_BLOCK_SIZE, MPI_DOUBLE, cur_proc, 666, MPI_COMM_WORLD);       
            }
        }

        delete A_glob;
    }
    else {
        for (int j_block = 0; j_block < block_count_loc; ++j_block) {
            MPI_Recv(buf, m*_BLOCK_SIZE, MPI_DOUBLE, mpiroot, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 0; i < m; ++i) {
                for (int j = 0; j < _BLOCK_SIZE; ++j) {
                    int j2 = j_block*_BLOCK_SIZE + j;
                    A[i*n + j2] = buf[i*_BLOCK_SIZE + j];
                }
            }
        }
    }
    delete buf;


    double start = MPI_Wtime();
    double *buf_c = new double[_BLOCK_SIZE * _BLOCK_SIZE];
    double *buf_s = new double[_BLOCK_SIZE * _BLOCK_SIZE];
    for (int j_block_glob = 0; j_block_glob < block_count_glob; ++j_block_glob) {

        int root = j_block_glob % proc_num;
        int j_block_loc = j_block_glob / proc_num;
        for (int i_block = j_block_glob; i_block < block_count_glob; ++i_block) {
            // вращение блока
            if (mpirank == root) {
                for (int j = 0; j < _BLOCK_SIZE; ++j) {
                    int col_loc = j_block_loc*_BLOCK_SIZE + j;
                    int col_glob = j_block_glob*_BLOCK_SIZE + j;
                    for (int i = 0; i < _BLOCK_SIZE; ++i) {
                        int row = i_block*_BLOCK_SIZE + i;
                        if (row <= col_glob)
                            continue;
                        int jj = col_glob*n + col_loc;
                        int ij = row*n + col_loc;
                        double ajj = A[jj];
                        double aij = A[ij];
                        double c;
                        double s;
                        cblas_drotg(&ajj, &aij, &c, &s);
                        cblas_drot(_BLOCK_SIZE-j, &A[jj], 1, &A[ij], 1, c, s);
                        buf_c[i*_BLOCK_SIZE + j] = c;
                        buf_s[i*_BLOCK_SIZE + j] = -s;
                    }
                }

                for (int i = 0; i < _BLOCK_SIZE; ++i) {
                    int row = i_block*_BLOCK_SIZE + i;
                    for (int j = 0; j < _BLOCK_SIZE; ++j) {
                        int col = j_block_loc*_BLOCK_SIZE + j;
                        if (j_block_glob == i_block && i < j)
                            Q[row*n + col] = buf_s[j*_BLOCK_SIZE + i];
                        else
                            Q[row*n + col] = buf_c[i*_BLOCK_SIZE + j];
                    }
                }
            }

            MPI_Bcast(buf_c, _BLOCK_SIZE * _BLOCK_SIZE,  MPI_DOUBLE, root, MPI_COMM_WORLD);
            MPI_Bcast(buf_s, _BLOCK_SIZE * _BLOCK_SIZE,  MPI_DOUBLE, root, MPI_COMM_WORLD);

            if (i_block != j_block_glob && i_block % proc_num == mpirank) {
                int j_block_loc2 = i_block / proc_num;
                for (int i = 0; i < _BLOCK_SIZE; ++i) {
                    int row = j_block_glob*_BLOCK_SIZE + i;
                    for (int j = 0; j < _BLOCK_SIZE; ++j) {
                        int col = j_block_loc2*_BLOCK_SIZE + j;
                        Q[row*n + col] = buf_s[j*_BLOCK_SIZE + i];
                    }
                }
            }

            // обновление строки
            int j0 = j_block_loc + (mpirank <= root);
            for (int j_block_loc2 = j0; j_block_loc2 < block_count_loc; ++j_block_loc2) {
                for (int j = 0; j < _BLOCK_SIZE; ++j) {
                    for (int i = 0; i < _BLOCK_SIZE; ++i) {
                        if (j_block_glob == i_block && i <= j)
                            continue;
                        int col = j_block_loc2*_BLOCK_SIZE;
                        int row1 = (j_block_glob*_BLOCK_SIZE + j)*n + col;
                        int row2 = (i_block*_BLOCK_SIZE + i)*n + col;
                        double c =  buf_c[i*_BLOCK_SIZE + j];
                        double s = -buf_s[i*_BLOCK_SIZE + j];
                        cblas_drot(_BLOCK_SIZE, &A[row1], 1, &A[row2], 1, c, s);
                    }
                }
            }
        }
    }
    delete buf_c;
    delete buf_s;
    double end = MPI_Wtime();
    if (mpiroot)
        std::cout << end-start;

    save_res(mpirank, A, Q, m, n);

    MPI_Finalize();

    return 0;
}