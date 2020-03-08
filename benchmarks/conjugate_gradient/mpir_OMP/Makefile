CCX = g++
CCXFLAGS = -std=c++14 -O3 -fopenmp

.SECONDARY:

.PHONY: build

build: sparsesolve

sparsesolve: sparsesolve.o mmio.o matrix.o cg.o ir.o
	$(CCX) -fopenmp $^ -o $@

clean:
	rm -rf *.o sparsesolve

test: sparsesolve
	./sparsesolve ../data/bcsstk01.mtx 1 1e-7 10000 1e-7 100

%.o: %.cpp cg.h vector.h matrix.h
	$(CCX) $(CCXFLAGS) -c $< -o $@
