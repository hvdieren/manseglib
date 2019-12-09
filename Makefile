main: main.o
	g++ -o main main.o

pagerank: pagerank.o
	g++ -O4 -o pagerank pagerank.o

msa_pagerank: msa_pagerank.o mantissaSegmentation.hpp
	g++ -O4 -o msa_pagerank msa_pagerank.o

main.o: main.cpp mantissaSegmentation.hpp
	g++ -O4 -c main.cpp 

.PHONY: clean
clean:
	rm main main.o pagerank pagerank.o msa_pagerank msa_pagerank.o