#!/bin/bash

# MTX="bcsstk01"
# MTX="bcsstk06"
MTX=("bcsstk01" "bcsstm27" "nasa2910" "Pres_Poisson" "obstclae" "Andrews") #"BenElechi1" "tmt_sym" "G3_Circuit")
TOL=1e-7
CG_ITER=2000
IR_ITER=10
STEP_CHECK=("10" "10" "10" "100" "100" "100") #"1000" "1000" "1000")

INDICES=(0 1 2 3 4 5) # 6 7 8)

for i in ${INDICES[@]}
do
	echo "ref start ${MTX[$i]}"
	date
	./mpir/sparsesolve ./data/${MTX[$i]}.mtx $IR_ITER $TOL $CG_ITER $TOL ${STEP_CHECK[$i]} > ./results/ref_${MTX[$i]}.out
	date
	echo "ref end ${MTX[$i]}"

	echo "manseg start ${MTX[i]}"
	date
	./mpir_class_manseg/sparsesolve ./data/${MTX[$i]}.mtx $IR_ITER $TOL $CG_ITER $TOL ${STEP_CHECK[$i]} > ./results/manseg_${MTX[$i]}.out
	date
	echo "manseg end ${MTX[$i]}"
done