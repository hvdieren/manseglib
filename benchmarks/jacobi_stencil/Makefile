CC=gcc
FLAGS= -O3 -mavx2
OMPFLAGS= -O3 -mavx2 -fopenmp

CCX=g++
CCXFLAGS= -std=c++11 -O3 -mavx2
CCXOMPFLAGS= -std=c++11 -O3 -mavx2 -fopenmp

ALL = J_jacobi_up jacobi_omp jacobi_mod jacobi_mod_omp

all: $(ALL)

## std jacobi
J_jacobi: J_jacobi.c
	$(CC) $(FLAGS) $< -o $@


## jacobi update (reference)
J_jacobi_up: J_jacobi_up.c
	$(CC) $(FLAGS) $< -o $@


## jacobi update with omp *only*
jacobi_omp: jacobi_omp.c
	$(CC) $(OMPFLAGS) $< -o $@


## jacobi with manseg library *only*
.PHONY: jacobi_mod.o
jacobi_mod.o: jacobi_mod.cpp ../../mantissaSegmentation.hpp
	$(CCX) $(CCXFLAGS) -c jacobi_mod.cpp

jacobi_mod: jacobi_mod.o
	$(CCX) -o jacobi_mod jacobi_mod.o


## jacobi using manseg library with omp
.PHONY: jacobi_mod_omp.o
jacobi_mod_omp.o: jacobi_mod_omp.cpp ../../mantissaSegmentation.hpp
	$(CCX) $(CCXOMPFLAGS) -c jacobi_mod_omp.cpp

jacobi_mod_omp: jacobi_mod_omp.o
	$(CCX) -fopenmp -o jacobi_mod_omp jacobi_mod_omp.o


# ## jacobi with manseg library, using sse2 vector ext.
# .PHONY: jacobi_mod_f.o
# jacobi_mod_f.o: jacobi_mod_f.cpp ../../mantissaSegmentation_f.hpp
# 	$(CCX) $(CCXFLAGS) -msse2 -c jacobi_mod_f.cpp

# jacobi_mod_f: jacobi_mod_f.o
# 	$(CCX) -o jacobi_mod_f jacobi_mod_f.o


# ## jacobi with manseg library, using vector ext, with omp
# .PHONY: jacobi_mod_f_omp.o
# jacobi_mod_f_omp.o: jacobi_mod_f_omp.cpp ../../mantissaSegmentation_f.hpp
# 	$(CCX) $(CCXOMPFLAGS) -msse2 -c jacobi_mod_f_omp.cpp

# jacobi_mod_f_omp: jacobi_mod_f_omp.o
# 	$(CCX) -fopenmp -o jacobi_mod_f_omp jacobi_mod_f_omp.o

clean:
	rm -f $(ALL) *.o
