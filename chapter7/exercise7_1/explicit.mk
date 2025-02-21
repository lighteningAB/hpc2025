default: fibonacci_main

fibonacci_main.o: main.cpp fibonacci.h
	g++ -c main.cpp -o fibonacci_main.o

fibonacci.o: fibonacci.cpp fibonacci.h
	g++ -c fibonacci.cpp -o fibonacci.o

fibonacci_main: fibonacci_main.o fibonacci.o
	g++ fibonacci_main.o fibonacci.o -o fibonacci_main