main: main.o
	g++ -o main main.o

main.o: main.cpp mantissaSegmentation.hpp
	g++ -O4 -c main.cpp 

.PHONY: clean
clean:
	rm main main.o