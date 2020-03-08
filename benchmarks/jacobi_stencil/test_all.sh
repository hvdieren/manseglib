#!/bin/bash

export OMP_NUM_THREADS=8

ITER=1000

JACOBI="jacobi_mod jacobi_mod_omp"

echo ""
for prog in $JACOBI;
do
    echo "start ${prog}"
    ./${prog} $ITER
    echo "end ${prog}"
    echo ""
done

diff --speed-large-files --strip-trailing-cr --suppress-common-lines \
    ./jacobi_mod_values.txt ./jacobi_mod_omp_values.txt | grep '>' -c
