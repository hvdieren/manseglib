CC=gcc
FLAGS= -O3 -mavx2
OMPFLAGS= -O3 -mavx2 -fopenmp

CCX=g++
CCXFLAGS= -std=c++11 -O3 -mavx2
CCXOMPFLAGS= -std=c++11 -O3 -mavx2 -fopenmp

ALL = jacobi_omp jacobi_mod_omp

all: $(ALL)

## jacobi update with omp *only*
jacobi_omp: jacobi_omp.c
	$(CC) $(OMPFLAGS) $< -o $@

## jacobi using manseg library with omp
.PHONY: jacobi_mod_omp.o
jacobi_mod_omp.o: jacobi_mod_omp.cpp ../../../../mantissaSegmentation_dev.hpp
	$(CCX) $(CCXOMPFLAGS) -c jacobi_mod_omp.cpp

jacobi_mod_omp: jacobi_mod_omp.o
	$(CCX) -fopenmp -o jacobi_mod_omp jacobi_mod_omp.o

clean:
	rm -f $(ALL) *.o
