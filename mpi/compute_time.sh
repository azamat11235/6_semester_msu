PROC_NUM=${1}

rm -f data.dat
rm -f a_part_*
rm -f q_part_*

echo "PROC_NUM = ${PROC_NUM}"
echo "--------------------"
SIZE=256
python3 gen_data.py ${SIZE} ${SIZE}
mpic++ qr_mpi.cpp parameters.h -lcblas -lblas
mpirun -n ${PROC_NUM} ./a.out ${SIZE} ${PROC_NUM}
echo "--------------------"

SIZE=512
python3 gen_data.py ${SIZE} ${SIZE}
mpic++ qr_mpi.cpp parameters.h -lcblas -lblas
mpirun -n ${PROC_NUM} ./a.out ${SIZE} ${PROC_NUM}
echo "--------------------"

SIZE=1024
python3 gen_data.py ${SIZE} ${SIZE}
mpic++ qr_mpi.cpp parameters.h -lcblas -lblas
mpirun -n ${PROC_NUM} ./a.out ${SIZE} ${PROC_NUM}
echo "--------------------"

SIZE=2048
python3 gen_data.py ${SIZE} ${SIZE}
mpic++ qr_mpi.cpp parameters.h -lcblas -lblas
mpirun -n ${PROC_NUM} ./a.out ${SIZE} ${PROC_NUM}
echo "--------------------"


PROC_NUM=1
echo "PROC_NUM = 1"
echo "--------------------"
SIZE=256
python3 gen_data.py ${SIZE} ${SIZE}
mpic++ qr_mpi.cpp parameters.h -lcblas -lblas
mpirun -n ${PROC_NUM} ./a.out ${SIZE} ${PROC_NUM}
echo "--------------------"

SIZE=512
python3 gen_data.py ${SIZE} ${SIZE}
mpic++ qr_mpi.cpp parameters.h -lcblas -lblas
mpirun -n ${PROC_NUM} ./a.out ${SIZE} ${PROC_NUM}
echo "--------------------"

SIZE=1024
python3 gen_data.py ${SIZE} ${SIZE}
mpic++ qr_mpi.cpp parameters.h -lcblas -lblas
mpirun -n ${PROC_NUM} ./a.out ${SIZE} ${PROC_NUM}
echo "--------------------"

SIZE=2048
python3 gen_data.py ${SIZE} ${SIZE}
mpic++ qr_mpi.cpp parameters.h -lcblas -lblas
mpirun -n ${PROC_NUM} ./a.out ${SIZE} ${PROC_NUM}
echo "--------------------"

rm -f data.dat
rm -f a_part_*
rm -f q_part_*
rm a.out