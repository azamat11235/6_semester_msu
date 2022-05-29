#!/bin/bash

printf '%-5s\t%-10s\t%-9s\n' size proc_count 'time (s.)'
echo '---------------------------------'
for ((size=256; size <= 2048; size*=2)) do
    if (test -f data/a.npy)
    then rm data/*
    fi
    python3 gen_data.py $size $size
    for ((n=1; n <= 64; n*=2)) do
        t=`mpiexec -n $n --oversubscribe python3 qr.py $size $n`
        printf '%-5s\t%10s\t%-.6f\n' $size $n $t
        #sleep 16
    done
    echo
    rm data/*
done
