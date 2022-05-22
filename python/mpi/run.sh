#!/bin/bash

printf '%-5s\t%-8s\t%-9s\n' size proc_num 'time (s.)'
echo '-----------------------------'
for ((size=256; size <= 2048; size*=2)) do
    rm data/*
    python3 gen_data.py $size $size
    for ((n=1; n <= 64; n*=2)) do
        t=`mpiexec -n $n python3 qr.py $size $n`
        printf '%-5s\t%8s\t%-.6f\n' $size $n $t
        sleep 16
    done
    echo
done
