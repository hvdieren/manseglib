
CCX=g++
CCXFLAGS=-O3 -std=c++17 -march=native  -msse2
CCXFLAGS_DEBUG=-O0 -std=c++17 -march=native

all: 
	make main tests pagerank pagerank_check msa_pagerank

main: main.o
	${CCX} -o main main.o

.PHONY: main.o
main.o: main.cpp mantissaSegmentation.hpp
	${CCX} ${CCXFLAGS_DEBUG} -c main.cpp 

tests: tests.o
	g++ -o tests tests.o

.PHONY: tests.o
tests.o: tests.cpp mantissaSegmentation.hpp mantissaSegmentation_f.hpp
	${CCX} ${CCXFLAGS} -c tests.cpp 

pagerank: pagerank.o
	${CCX} -o pagerank pagerank.o

.PHONY: pagerank.o
pagerank.o: pagerank.cpp quicksort.h
	${CCX} ${CCXFLAGS} -c pagerank.cpp

pagerank_check: pagerank_check.o
	${CCX} -o pagerank_check pagerank_check.o

.PHONY: pagerank_check.o
pagerank_check.o: pagerank_check.cpp quicksort.h
	${CCX} ${CCXFLAGS} -c pagerank_check.cpp

msa_pagerank: msa_pagerank.o
	${CCX} -o msa_pagerank msa_pagerank.o

msa_pagerank.o: msa_pagerank.cpp mantissaSegmentation.hpp mantissaSegmentation_f.hpp quicksort.h
	${CCX} ${CCXFLAGS} -c msa_pagerank.cpp

.PHONY: clean
clean:
	rm main tests pagerank pagerank_check msa_pagerank *.o