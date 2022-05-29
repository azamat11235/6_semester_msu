import numpy as np
import scipy.linalg
import os
import sys
import glob

def restore_q_batch(q, block_size):
    m, n = q.shape
    res = np.zeros(q.shape)
    for i in range(min(q.shape)):
        res[i, i] = 1
    for jb in range(n-1, -1, -block_size):
        for ib in range(m-1, -1, -block_size):
            for j in range(block_size):
                for i in range(block_size):
                    row = ib - i
                    col = jb - j
                    if row > col:
                        c = q[row, col]
                        s = q[col, row]
                        scipy.linalg.blas.drot(res[col], res[row], c, s, overwrite_x=1, overwrite_y=1)
    return res

size = 128
if len(sys.argv) > 1:
    size = int(sys.argv[1])
    proc_num = int(sys.argv[2])
# os.system(f'cda; mpic++ qr_mpi.cpp parameters.h -lcblas -lblas; mpirun -n 2 ./a.out {size};')

a = []
for line in open('data/data.dat'):
    a.append(np.asarray(list(map(float, line.split()))))
a = np.asarray(a)

for line in open('parameters.h'):
    if '_BLOCK_SIZE' in line:
        block_size = int(line.split()[-1])
        break
# print(block_size, proc_num)

a_files = sorted(glob.glob('data/a_part_*.dat'))
q_files = sorted(glob.glob('data/q_part_*.dat'))
r = np.empty(a.shape)
q = np.empty(a.shape)
for fname in a_files:
    mpirank = int(os.path.splitext(fname)[0][-2:])
    ai = []
    for line in open(fname):
        ai.append(np.asarray(list(map(float, line.split()))))
    ai = np.asarray(ai)
    m, n = ai.shape
    if n % block_size != 0:
        print(r"WARNING: n % block_size != 0")
    for j_block in range(n//block_size):
        j0 = (mpirank + j_block*proc_num) * block_size
        j1 = j_block*block_size
        r[:, j0 : j0+block_size] = ai[:, j1 : j1+block_size]

for fname in q_files:
    mpirank = int(os.path.splitext(fname)[0][-2:])
    qi = []
    for line in open(fname):
        qi.append(np.asarray(list(map(float, line.split()))))
    qi = np.asarray(qi)
    m, n = qi.shape
    if n % block_size != 0:
        print(r"WARNING: n % block_size != 0")
    for j_block in range(n//block_size):
        j0 = (mpirank + j_block*proc_num) * block_size
        j1 = j_block*block_size
        q[:, j0 : j0+block_size] = qi[:, j1 : j1+block_size]

q = restore_q_batch(q, block_size)

print("||A - QR||    =", np.linalg.norm(a- q@r))
print("||Q.T*Q - I|| =", np.linalg.norm(q.T @ q - np.eye(*q.shape)))