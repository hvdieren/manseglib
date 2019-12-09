#!/bin/bash

rm pagerank pagerank.o
rm msa_pagerank msa_pagerank.o

make pagerank msa_pagerank

# std
echo "std version"
# COO
# ./pagerank COO "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/smallGraph.coo"
# ./pagerank COO "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/orkut_undir.coo"
./pagerank COO "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/LiveJournal_dir.coo" >results/std_livejournal.coo.txt
# ./pagerank COO "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/rMatGraph_J_5_100.coo"

# CSR
# ./pagerank CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/smallGraph.csr"
# ./pagerank CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/orkut_undir.csr"
./pagerank CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/LiveJournal_dir.csr" >results/std_livejournal.csr.txt
# ./pagerank CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/rMatGraph_J_5_100.csr"

# CSC
# ./pagerank CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/smallGraph.csc"
# ./pagerank CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/orkut_undir.csc"
# ./pagerank CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/LiveJournal_dir.csc"
# ./pagerank CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/rMatGraph_J_5_100.csc"

echo ""
echo "msa version"
# msa version
# COO
# ./pagerank COO "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/smallGraph.coo"
# ./pagerank COO "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/orkut_undir.coo"
./pagerank COO "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/LiveJournal_dir.coo" >results/msa_livejournal.coo.txt
# ./pagerank COO "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/rMatGraph_J_5_100.coo"

# CSR
# ./pagerank CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/smallGraph.csr"
# ./pagerank CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/orkut_undir.csr"
./pagerank CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/LiveJournal_dir.csr" >results/msa_livejournal.csr.txt
# ./pagerank CSR "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/rMatGraph_J_5_100.csr"

# CSC
# ./pagerank CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/smallGraph.csc"
# ./pagerank CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/orkut_undir.csc"
# ./msa_pagerank CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/LiveJournal_dir.csc"
# ./pagerank CSC "/mnt/d/Jorta/Documents/Uni/4th Year/CSC3021_Concurrent_Programming/Assignments/graphs/rMatGraph_J_5_100.csc"
