#!/bin/bash

# comment this out if running on kelvin
cd benchmarks

## this is actual script from here

echo "starting ligra"
date
cd ligra-partition
./pagerank_bench.sh > pagerank_results.txt 2>&1
cd ..
date
echo "finished ligra"

echo "starting jacobi"
date
cd jacobi_stencil
./test_all.sh > jacobi_results.txt 2>&1
cd ..
date
echo "finished jacobi"

echo "starting cg"
date
cd conjugate_gradient
./test.sh > cg_results.txt 2>&1
cd ..
date
echo "finished cg"
