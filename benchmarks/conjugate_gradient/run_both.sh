#!/bin/bash

# MTX="bcsstk01"
MTX=("bcsstk01.mtx" "bcsstk06.mtx")
# MTX=("2cubes_sphere.mtx" "2cubes_sphere_b.mtx" "Andrews.mtx" "BenElechi1.mtx" "BenElechi1_b.mtx" "Emilia_923.mtx" "G2_circuit.mtx" "G3_circuit.mtx" "Pres_Poisson.mtx" "Pres_Poisson_b.mtx" \
# 	"Queen_4147.mtx" "Serena.mtx" "StocF-1465.mtx" "bcsstk01.mtx" "bcsstk06.mtx" "bcsstk13.mtx" "bcsstk15.mtx" "bcsstk16.mtx" "bcsstk19.mtx" "bcsstk24.mtx" "bcsstk28.mtx" "bcsstm25.mtx" \
# 	"bcsstm27.mtx" "bone010.mtx" "bone010_B.mtx" "bone010_C.mtx" "bone010_M.mtx" "inline_1.mtx" "nasa2910.mtx" "nasa2910_b.mtx" "obstclae.mtx" "offshore.mtx" "offshore_b.mtx" "plat1919.mtx" \
# 	 "thermomech_dM.mtx" "tmt_sym.mtx" "x104.mtx")
TOL=1e-7
CG_ITER=2000
IR_ITER=10
STEP_CHECK="10"

INDICES=(0 1)
# INDICES=(0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36)

for i in ${INDICES[@]}
do
	echo "ref start ${MTX[$i]}"
	date
	./mpir/sparsesolve ./data/${MTX[$i]} $IR_ITER $TOL $CG_ITER $TOL $STEP_CHECK > ./results/ref_${MTX[$i]}.out
	date
	echo "ref end ${MTX[$i]}"

	echo "manseg start ${MTX[i]}"
	date
	./mpir_class_manseg/sparsesolve ./data/${MTX[$i]} $IR_ITER $TOL $CG_ITER $TOL $STEP_CHECK > ./results/manseg_${MTX[$i]}.out
	date
	echo "manseg end ${MTX[$i]}"
done