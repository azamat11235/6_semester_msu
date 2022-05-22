import numpy as np
import sys

fname = 'data.dat'
f = open(fname, 'w')

m = 256
n = 256
if len(sys.argv) > 1:
    m = int(sys.argv[1])
    n = int(sys.argv[2])

a = np.random.rand(m, n)

for row in a:
    for elem in row:
        f.write('%f ' % elem)
    f.write('\n')

f.close()
print('Matrix(%d, %d)' % (m, n))