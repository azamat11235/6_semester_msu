SIZE=${1}
PROC_NUM=${2}

rm -f data.dat
rm -f a_part_*
rm -f q_part_*

python3 gen_data.py ${SIZE} ${SIZE}

mpic++ qr_mpi.cpp parameters.h -lcblas -lblas
mpirun -n ${PROC_NUM} ./a.out ${SIZE} ${PROC_NUM}
echo "--------------------"

rm -f data.dat
rm -f a_part_*
rm -f q_part_*
rm a.out