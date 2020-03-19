#!/bin/bash

export OMP_NUM_THREADS=8

ITER=1000

JACOBI="J_jacobi_up jacobi_omp jacobi_mod jacobi_mod_omp"

echo ""
for prog in $JACOBI;
do
    echo "start ${prog}"
    ./${prog} $ITER
    echo "end ${prog}"
    echo ""
done
