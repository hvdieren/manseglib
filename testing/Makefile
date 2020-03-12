
CCX=g++
CCXFLAGS=-O3 -std=c++14 -msse

all: 
	make tests

main: main.o
	${CCX} -o main main.o

tests: tests.o
	g++ -o tests tests.o

.PHONY: tests.o
tests.o: tests.cpp ../mantissaSegmentation_dev.hpp
	${CCX} ${CCXFLAGS} -c tests.cpp 

.PHONY: clean
clean:
	rm tests *.o