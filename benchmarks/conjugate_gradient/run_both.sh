#!/bin/bash

MTX="bcsstk01 bcsstk06 bcsstk13"
# MTX="bcsstk01 plat1919 nasa2910 bcsstm25 Pres_Poisson obstclae Andrews BenElechi1 tmt_sym G3_Circuit"
TOL=1e-7
CG_ITER=2000
IR_ITER=10
STEP_CHECK=10

for mat in $MTX
do
    echo "ref start ${mat}"
    date
    ./mpir/sparsesolve ./data/${mat}.mtx $IR_ITER $TOL $CG_ITER $TOL $STEP_CHECK > ./results/ref_${mat}.out
    date
    echo "ref end ${mat}"
    echo ""

    echo "manseg start ${mat}"
    date
    ./mpir_class_manseg/sparsesolve ./data/${mat}.mtx $IR_ITER $TOL $CG_ITER $TOL $STEP_CHECK > ./results/manseg_${mat}.out
    date
    echo "manseg end ${mat}"
    echo ""
done