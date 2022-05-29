import os

dname = 'output/'
out = 'RESULTS_nodes.txt'
out = open(out, 'w')

print('%-5s\t%5s\t%-9s' % ('size', 'nodes', 'time (s.)'), file=out)
print('-'*33, file=out)
for fname in sorted(os.listdir(dname)):
    with open(dname+fname) as f:
        line = f.readline()
        size, n, t = line.split()
        print('%-5s\t%5s\t%-9s' % (size, n, t), file=out)
    if n == '4' and size != '2048':
        print(file=out)

out.close()

