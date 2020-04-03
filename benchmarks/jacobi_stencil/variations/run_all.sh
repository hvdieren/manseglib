#!/bin/bash

VAR="nb16_b256" # nb32_b128 nb32_b512 nb64_b64 nb64_b256 nb128_b32 nb128_b128 nb256_b16 nb256_b64 nb512_b32"

JS="J_jacobi_up jacobi_omp jacobi_mod jacobi_mod_omp"

ITER=1000 #0

for var in $VAR;
do
	echo "runs for ${var}"
	date
	for j in $JS;
	do
		echo "start ${j}"
		date
		./${var}/${j} $ITER
		date
		echo "end ${j}"
		echo ""
	done
	date
	echo "finished runs for ${var}"
	echo "----------------------------------------"
done