
CCX=g++
CCXFLAGS=-O3 -std=c++11 -march=native

all: 
	make main

main: main.o
	${CCX} -o main main.o

.PHONY: main.o
main.o: main.cpp mantissaSegmentation.hpp
	${CCX} ${CCXFLAGS} -c main.cpp 

tests: tests.o
	g++ -o tests tests.o

.PHONY: clean
clean:
	rm main  *.o