import numpy as np
import scipy.linalg
import sys
import math
from mpi4py import MPI

def save_res(R, Q_c, Q_s):
    if my_rank < 10:
        num = '0' + str(my_rank)
    else:
        num = str(my_rank)
    fname_r = 'data/r_part_' + num + '.npy'
    fname_q = 'data/q_%s_part_' + num + '.npy'
    
    with open(fname_r, 'wb') as f:
        np.save(f, R)
    with open(fname_q % 'c', 'wb') as f:
        np.save(f, Q_c)
    with open(fname_q % 's', 'wb') as f:
        np.save(f, Q_s)


size = int(sys.argv[1])
proc_num = int(sys.argv[2])
comm = MPI.COMM_WORLD
my_rank = comm.Get_rank()
proc_name = MPI.Get_processor_name()
comm_size = comm.Get_size()
mpiroot = (my_rank == 0)

m = size // proc_num
n = size
Q_c = np.empty((m, n))
Q_s = np.empty((m, n))

'''if mpiroot:
    # print('matrix distribution...')
    fname = 'data/a.npy'
    A_glob = np.load(fname)
    for i in range(1, proc_num):
        comm.Send(A_glob[i::proc_num, :].flatten(), dest=i)
    A = A_glob[::proc_num, :]
else:
    A = np.empty(m*n)
    comm.Recv(A, source=0)
    A = A.reshape(m, n)
comm.Barrier()'''
A = np.random.rand(m, n)

buf = np.empty(n)
# if mpiroot: print('qr...')
t1 = MPI.Wtime()
for j in range(n-1):
    i0 = (j - my_rank) / proc_num
    main_proc = j % proc_num
    i0 = 0 if i0 < 0 else math.ceil(i0)
    for i in range(i0 + 1, A.shape[0]):
        c, s = scipy.linalg.blas.drotg(A[i0, j], A[i, j])
        scipy.linalg.blas.drot(A[i0], A[i], c, s, offx=j, offy=j, overwrite_x=1, overwrite_y=1)
        Q_c[i, j], Q_s[i, j] = c, -s
    if my_rank == main_proc:
        proc = (my_rank + 1) % proc_num
        for i in range(proc_num - 1):
            if j + i + 1 >= size: break
            comm.Send(A[i0, j:], dest=proc)
            comm.Recv(A[i0, j:], source=proc)
            proc = (proc + 1) % proc_num
    elif i0 < m:
            comm.Recv(buf[j:], source=main_proc)
            c, s = scipy.linalg.blas.drotg(buf[j], A[i0, j])
            scipy.linalg.blas.drot(buf, A[i0], c, s, offx=j, offy=j, overwrite_x=1, overwrite_y=1)
            Q_c[i0, j], Q_s[i0, j] = c, -s
            comm.Send(buf[j:], dest=main_proc)
comm.Barrier()
t2 = MPI.Wtime()

# save_res(A, Q_c, Q_s)
if mpiroot:
    print('%-5d\t%5d\t%-.6f\n' % (size, proc_num, t2-t1)) # print(t2 - t1)

comm.Barrier()
print('Node: %s, rank: %d/%d' % (proc_name, my_rank, comm_size))
