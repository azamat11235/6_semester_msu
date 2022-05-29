import numpy as np
import scipy.linalg
import os
import sys
import glob


def restore_q(q, proc_num):
    m, n = q.shape
    res = np.zeros(q.shape)
    for i in range(min(q.shape)):
        res[i, i] = 1
    for j in range(n-2, -1, -1):
        for i in range(proc_num-1, 0, -1):
            ii = j + i
            if ii >= m: continue
            c = q[ii, j]
            s = q[j, ii]
            scipy.linalg.blas.drot(res[j], res[ii], c, s, offx=j, offy=j, overwrite_x=1, overwrite_y=1)

        for i in range(proc_num):
            i0 = j + i
            i1 = (m - i0 - 1) // proc_num * proc_num + i0
            if i0 >= m or i1 >= m: break
            for ii in range(i1, i0, -proc_num):
                c = q[ii, j]
                s = q[j, ii]
                scipy.linalg.blas.drot(res[i0], res[ii], c, s, offx=j, offy=j, overwrite_x=1, overwrite_y=1)
    return res

a = np.load('data/a.npy')

r_files = sorted(glob.glob('data/r_part_*.npy'))
q_s_files = sorted(glob.glob('data/q_s_part_*.npy'))
r = np.empty(a.shape)
q = np.empty(a.shape)
proc_num = len(r_files)

for fname in r_files:
    mpirank = int(os.path.splitext(fname)[0][-2:])
    ri = np.load(fname)
    r[mpirank::proc_num, :] = ri

for fname in q_s_files:
    mpirank = int(os.path.splitext(fname)[0][-2:])
    q_s = np.load(fname)
    q_c = np.load(fname.replace('q_s', 'q_c'))
    for i in range(mpirank, q.shape[0], proc_num):
        for j in range(min(i, q.shape[1])):
            i_loc = (i - mpirank) // proc_num
            q[i, j] = q_c[i_loc, j]
            q[j, i] = q_s[i_loc, j]        

q0 = q.copy()
q = restore_q(q, proc_num)

for i in range(r.shape[0]):
    for j in range(r.shape[1]):
        if i > j and abs(r[i,j]) > 1e-15:
            print('!!!', i, j, r[i,j])
            break
print("||A - QR||_F    =", np.linalg.norm(a - q @ r))
print("||Q.T*Q - I||_F =", np.linalg.norm(q.T @ q - np.eye(q.shape[0])))
