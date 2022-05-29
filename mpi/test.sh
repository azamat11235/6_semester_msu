SIZE=256
PROC_NUM=4

rm -f data/*.dat

python3 gen_data.py ${SIZE} ${SIZE}
mpic++ qr_mpi.cpp parameters.h -lcblas -lblas
mpirun -n ${PROC_NUM} ./a.out ${SIZE} ${PROC_NUM} > /dev/null

python3 test.py ${SIZE} ${PROC_NUM}

rm -f data/*.dat
rm a.out
