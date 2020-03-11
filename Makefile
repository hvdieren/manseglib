
CCX=g++
CCXFLAGS=-O3 -std=c++11 -march=native

all: 
	make main

main: main.o matrix.o
	${CCX} -o main main.o matrix.o

.PHONY: main.o
main.o: main.cpp mantissaSegmentation.hpp matrix.h vector.h
	${CCX} ${CCXFLAGS} -c main.cpp 

# mantissaSegmentation.o: mantissaSegmentation.cpp
# 	${CCX} ${CCXFLAGS} -c mantissaSegmentation.cpp

matrix.o: matrix.cpp mantissaSegmentation.hpp
	${CCX} ${CCXFLAGS} -c matrix.cpp

tests: tests.o
	g++ -o tests tests.o

.PHONY: clean
clean:
	rm main  *.o