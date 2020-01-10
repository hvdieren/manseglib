main: main.o
	g++ -o main main.o

main.o: main.cpp mantissaSegmentation.hpp
	g++ -O2 -c main.cpp 

pagerank: pagerank.o
	g++ -o pagerank pagerank.o

pagerank.o: pagerank.cpp
	g++ -O4 -c pagerank.cpp

msa_pagerank: msa_pagerank.o
	g++ -o msa_pagerank msa_pagerank.o

msa_pagerank.o: msa_pagerank.cpp mantissaSegmentation.hpp
	g++ -O4 -c msa_pagerank.cpp

.PHONY: clean
clean:
	rm main main.o pagerank pagerank.o msa_pagerank msa_pagerank.o