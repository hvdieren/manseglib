main: main.o
	g++ -o main main.o

all: 
	make main tests pagerank pagerank_check msa_pagerank

.PHONY: main.o
main.o: main.cpp mantissaSegmentation.hpp
	g++ -c main.cpp 

tests: tests.o
	g++ -o tests tests.o

.PHONY: tests.o
tests.o: tests.cpp mantissaSegmentation.hpp
	g++ -O4 -c tests.cpp 

pagerank: pagerank.o
	g++ -o pagerank pagerank.o

.PHONY: pagerank.o
pagerank.o: pagerank.cpp quicksort.h
	g++ -O4 -c pagerank.cpp

pagerank_check: pagerank_check.o
	g++ -o pagerank_check pagerank_check.o

.PHONY: pagerank_check.o
pagerank_check.o: pagerank_check.cpp quicksort.h
	g++ -O4 -c pagerank_check.cpp

#pagerank_parallel: pagerank_parallel.o
#	g++ -fopenmp -o pagerank_parallel pagerank_parallel.o

#pagerank_parallel.o: pagerank_parallel.cpp quicksort.h
#	g++ -O4 -fopenmp -c pagerank_parallel.cpp

msa_pagerank: msa_pagerank.o
	g++ -o msa_pagerank msa_pagerank.o

msa_pagerank.o: msa_pagerank.cpp mantissaSegmentation.hpp
	g++ -O4 -c msa_pagerank.cpp

.PHONY: clean
clean:
	rm main tests pagerank pagerank_check msa_pagerank *.o