CC=g++
CXXFLAGS=-std=c++14 -Wall -O2
HDRS=fibonacci.h 

default: fibonacci_main

# Use implicit rules to compile .cpp files into .o files
%.o: %.cpp $(HDRS)
	$(CC) $(CXXFLAGS) -c $< -o $@

# Linking the final executable
fibonacci_main: main.o fibonacci.o
	$(CC) $(CXXFLAGS) $^ -o $@

clean:
	rm -f fibonacci_main *.o