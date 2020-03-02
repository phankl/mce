#Generic c++ Makefile

CXX = mpicxx
CXXFLAGS = -O3 -g -fopenmp -std=c++11

MAIN = mce

SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)

$(MAIN): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

.PHONY: clean
clean:
	rm -f $(OBJ) $(MAIN)
