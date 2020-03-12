
CCX=g++
CCXFLAGS=-O3 -std=c++11

all: 
	make pagerank pagerank_check msa_pagerank

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

msa_pagerank.o: msa_pagerank.cpp ../../../mantissaSegmentation_dev.hpp quicksort.h
	${CCX} ${CCXFLAGS} -c msa_pagerank.cpp

.PHONY: clean
clean:
	rm pagerank pagerank_check msa_pagerank *.o