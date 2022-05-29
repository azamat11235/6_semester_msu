#!/bin/bash

# printf '%-5s\t%5s\t%-9s\n' size nodes 'time (s.)'
# echo '---------------------------------'
touch output/a
rm output/*
for ((size=256; size <= 2048; size*=2)) do
    for ((n=1; n <= 4; n*=2)) do
    	if [ $size -lt 1024 ]
	then sbatch -p new -o output/slurm_0${size}_${n}.out --nodes=$n --job-name=0${size}_${n} run_nodes.sh $size $n > /dev/null
        else sbatch -p new -o output/slurm_${size}_${n}.out --nodes=$n --job-name=${size}_${n} run_nodes.sh $size $n > /dev/null
    	fi
    done
done
echo Submitted batch jobs
