import numpy as np
import sys

fname = 'data/a.npy'

m = 256
n = 256
if len(sys.argv) > 1:
    m = int(sys.argv[1])
    n = int(sys.argv[2])

a = np.random.rand(m, n)

with open(fname, 'wb') as f:
    np.save(f, a)

# print('Matrix(%d, %d)' % (m, n))
