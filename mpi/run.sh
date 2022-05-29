#!/bin/bash

printf '%-5s\t%-10s\t%-9s\n' size proc_count 'time (s.)'
echo '---------------------------------'
for ((size=256; size <= 2048; size*=2)) do
    rm -f data/data.dat
    rm -f data/a_part_*
    rm -f data/q_part_*
    python3 gen_data.py $size $size
    for ((proc_num=1; proc_num <= 4; proc_num*=2)) do
        mpic++ qr_mpi.cpp parameters.h -lcblas -lblas
        t=`mpirun -n ${proc_num} ./a.out ${size} ${proc_num}`
        printf '%-5s\t%10s\t%-.6f\n' $size $proc_num $t
        sleep 10
    done
    echo
done
rm a.out
