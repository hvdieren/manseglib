
CCX=g++
CCXFLAGS=-O3 -std=c++17 -mavx2

all: 
	make tests

tests: tests.o
	g++ -o tests tests.o

.PHONY: tests.o
tests.o: tests.cpp ../mantissaSegmentation_dev.hpp
	${CCX} ${CCXFLAGS} -c tests.cpp 

.PHONY: clean
clean:
	rm tests *.o