#!/bin/bash

size=${1}
n=${2}

mpiexec -n $n python3 qr.py $size $n
