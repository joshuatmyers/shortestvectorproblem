CXX = g++
CXXFLAGS = -std=c++17

SRC = main.cpp
OBJ = $(SRC:.cpp=.o)
EXEC = runme

all: main.cpp main.o
	$(CXX) -o runme main.cpp

%.o: %.cpp main.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)

test:
	./$(EXEC) [15.0 23.0 11.0] [46.0 15.0 3.0] [32.0 1.0 1.0]
	
clean:
	rm -f main.o runme result.txt