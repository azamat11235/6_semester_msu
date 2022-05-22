SIZE=128
PROC_NUM=4

rm -f data.dat
rm -f a_part_*
rm -f q_part_*

python3 gen_data.py ${SIZE} ${SIZE}

mpic++ qr_mpi.cpp parameters.h -lcblas -lblas
mpirun -n ${PROC_NUM} ./a.out ${SIZE} ${PROC_NUM}

python3 test.py ${SIZE} ${PROC_NUM}

rm -f data.dat
rm -f a_part_*
rm -f q_part_*
rm a.out
