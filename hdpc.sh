#!/bin/bash

# note: --all reserves all cores on the machine, and is important for performance measurements, for the sake of testing, 
# note: -n 1 ensures only one instance of script running
# note: using --cpus-per-task=1 is fine for dev and testing (remove for actual benchmarking)

## todo: include modules/add for gcc


# note: should probably try using echo ./<program> >> ./x.txt
# or something similar (see ligragraph or whatever .sh files for inspiration)

## use command srun -p xeon -w jacob -t 30 -n 1 --all ./hpdc.sh

#SRUN --ntasks=1
#SRUN --cpus-per-task=1

#SRUN -o pagerank.out
#SRUN -e pagerank.err

rm pagerank pagerank.o msa_pagerank msa_pagerank.o
# rm msa_pagerank msa_pagerank.o

make pagerank msa_pagerank

## default small graphs ##
# std
echo "std version"
# COO
./pagerank default COO "../graphs/smol/LiveJournal_dir.coo"

# CSR
./pagerank default CSR "../graphs/smol/LiveJournal_dir.csr"

# CSC
./pagerank default CSC "../graphs/smol/LiveJournal_dir.csc"

echo ""
echo "msa version"
# msa version
# COO
./msa_pagerank default COO "../graphs/smol/LiveJournal_dir.coo"

# CSR
./msa_pagerank default CSR "../graphs/smol/LiveJournal_dir.csr"

# CSC
./msa_pagerank default CSC "../graphs/smol/LiveJournal_dir.csc"
