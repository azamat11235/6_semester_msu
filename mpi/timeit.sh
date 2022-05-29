SIZE=${1}
PROC_NUM=${2}

rm -f data/*.dat

python3 gen_data.py ${SIZE} ${SIZE}

printf '%-5s\t%-8s\t%-9s\n' size proc_num 'time (s.)'

mpic++ qr_mpi.cpp parameters.h -lcblas -lblas
t=`mpirun -n ${PROC_NUM} ./a.out ${SIZE} ${PROC_NUM}`
printf '%-5s\t%8s\t%-.6f\n' $SIZE $PROC_NUM $t

rm -f data/*.dat
rm a.out