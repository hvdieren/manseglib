#!/bin/bash

## default small graphs ##
# std

export OMP_NUM_THREADS=8

echo "std version"
# COO
# ./pagerank_parallel default COO "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/smallGraph.coo"
# ./pagerank_parallel default COO "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/orkut_undir.coo"
# ./pagerank_parallel default COO "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/LiveJournal_dir.coo"
# ./pagerank_parallel default COO "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/rMatGraph_J_5_100.coo"

# CSR
# ./pagerank_parallel default CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/smallGraph.csr"
# ./pagerank_parallel default CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/orkut_undir.csr"
# ./pagerank_parallel default CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/LiveJournal_dir.csr"
# ./pagerank_parallel default CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/rMatGraph_J_5_100.csr"

# CSC
# ./pagerank_parallel default CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/smallGraph.csc"
# ./pagerank_parallel default CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/orkut_undir.csc"
./pagerank_parallel default CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/LiveJournal_dir.csc"
# ./pagerank_parallel default CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/rMatGraph_J_5_100.csc"

echo ""
echo "msa version"
# msa version
# COO
# ./msa_pagerank default COO "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/smallGraph.coo"
# ./msa_pagerank default COO "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/orkut_undir.coo"
# ./msa_pagerank default COO "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/LiveJournal_dir.coo"
# ./msa_pagerank default COO "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/rMatGraph_J_5_100.coo"

# CSR
# ./msa_pagerank default CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/smallGraph.csr"
# ./msa_pagerank default CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/orkut_undir.csr"
# ./msa_pagerank default CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/LiveJournal_dir.csr"
# ./msa_pagerank default CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/rMatGraph_J_5_100.csr"

# CSC
# ./msa_pagerank default CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/smallGraph.csc"
# ./msa_pagerank default CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/orkut_undir.csc"
./msa_pagerank_parallel default CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/LiveJournal_dir.csc"
# ./msa_pagerank default CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/rMatGraph_J_5_100.csc" 




## SNAP ##
# echo "std version"
# std
# COO
# ./pagerank SNAP COO "/mnt/d/Jorta/Documents/Uni/5thYear/CSC4006_MEngProj/Project_Stuff/graphs/Yahoo.txt"
# ./pagerank SNAP COO "/mnt/d/Jorta/Documents/Uni/5thYear/CSC4006_MEngProj/Project_Stuff/graphs/smallGraph.txt"

# CSR
# ./pagerank SNAP CSR "/mnt/d/Jorta/Documents/Uni/5thYear/CSC4006_MEngProj/Project_Stuff/graphs/Yahoo.txt"
# ./pagerank SNAP CSR "/mnt/d/Jorta/Documents/Uni/5thYear/CSC4006_MEngProj/Project_Stuff/graphs/smallGraph.txt"

# CSC
# ./pagerank SNAP CSC "/mnt/d/Jorta/Documents/Uni/5thYear/CSC4006_MEngProj/Project_Stuff/graphs/Yahoo.txt"
# ./pagerank SNAP CSC "/mnt/d/Jorta/Documents/Uni/5thYear/CSC4006_MEngProj/Project_Stuff/graphs/smallGraph.txt"

# echo ""
# echo "msa version"
# msa
# COO
# ./msa_pagerank SNAP COO "/mnt/d/Jorta/Documents/Uni/5thYear/CSC4006_MEngProj/Project_Stuff/graphs/Yahoo.txt"
